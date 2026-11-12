tf_screening_validate_intervals <- function(target_interval, fallback_target_interval) {
  valid_level <- function(value) {
    is.numeric(value) && length(value) == 1L &&
      is.finite(value) && value > 0 && value < 1
  }
  if (!valid_level(target_interval)) {
    stop("`target_interval` must be one finite number strictly between zero and one.",
         call. = FALSE)
  }
  if (!is.null(fallback_target_interval) &&
      (!valid_level(fallback_target_interval) || fallback_target_interval >= target_interval)) {
    stop("`fallback_target_interval` must be NULL or a positive interval level below `target_interval`.",
         call. = FALSE)
  }
}

tf_screening_target_draws <- function(model, gene, inference, sampling = NULL) {
  valid_count <- function(value, minimum = 1L) {
    is.numeric(value) && length(value) == 1L && is.finite(value) &&
      value >= minimum && value == trunc(value) && value <= .Machine$integer.max
  }
  source_count <- model$diagnostics$posterior_draws
  draws <- model$target_screening_draws
  if (is.null(draws) && identical(inference, "mcmc")) {
    legacy <- model$stage2_parameter_draws
    count <- legacy$draw_count
    complete <- is.list(legacy) &&
      identical(legacy$interface_version, "tf_em_stage2_parameter_draws_v2") &&
      valid_count(count, 2L) && valid_count(source_count, 2L) && count == source_count &&
      length(legacy$source_draw_index) == count &&
      is.numeric(legacy$source_draw_index) &&
      all(is.finite(legacy$source_draw_index)) &&
      all(legacy$source_draw_index == seq_len(count)) &&
      is.matrix(legacy$beta) && nrow(legacy$beta) == count
    if (isTRUE(complete)) {
      target_column <- which(toupper(colnames(legacy$beta)) ==
        toupper(model$target_tf$tf[[1]]))
      if (length(target_column) == 1L) {
        draws <- list(interface_version = "tf_target_screening_draws_v1",
          draw_count = count, overall = as.numeric(legacy$beta[, target_column]))
        if (isTRUE(model$condition_model)) {
          draws$control <- legacy$beta_target_control
          draws$disease <- legacy$beta_target_disease
        }
      }
    }
  }
  required <- if (isTRUE(model$condition_model)) c("overall", "control", "disease") else "overall"
  valid <- is.list(draws) &&
    identical(draws$interface_version, "tf_target_screening_draws_v1") &&
    valid_count(draws$draw_count, 2L)
  if (!is.null(source_count)) {
    valid <- valid && valid_count(source_count, 2L) && draws$draw_count == source_count
  }
  sampling <- if (!is.null(model$sampling)) model$sampling else sampling
  if (identical(inference, "mcmc") && !is.null(sampling)) {
    if (!is.list(sampling)) {
      valid <- FALSE
    } else if (!is.null(sampling$chains) || !is.null(sampling$iter_sampling)) {
      valid <- valid && valid_count(sampling$chains) && valid_count(sampling$iter_sampling) &&
        draws$draw_count == as.double(sampling$chains) * sampling$iter_sampling
    }
  }
  if (isTRUE(valid)) {
    valid <- all(vapply(required, function(field) {
      value <- draws[[field]]
      is.numeric(value) && length(value) == draws$draw_count && all(is.finite(value))
    }, logical(1)))
  }
  if (!isTRUE(valid)) {
    stop(sprintf(paste0(
      "Cannot relax the target interval for `%s` from this cached result. ",
      "Complete target posterior draws are required. Rerun with `force_refit = TRUE`."
    ), gene), call. = FALSE)
  }
  draws
}

tf_screening_recompute_target <- function(model, gene, interval, inference, sampling = NULL) {
  draws <- tf_screening_target_draws(model, gene, inference, sampling)
  overall <- tf_stage1_screening_quantile_interval(draws$overall, interval)
  target <- model$target_tf
  target$screening_interval_level <- interval
  target$screening_interval_lower <- overall[[1]]
  target$screening_interval_upper <- overall[[2]]
  target$screening_interval_excludes_zero <- overall[[1]] > 0 || overall[[2]] < 0
  if (isTRUE(model$condition_model)) {
    for (condition in c("control", "disease")) {
      bounds <- tf_stage1_screening_quantile_interval(draws[[condition]], interval)
      target[[paste0("beta_", condition, "_interval_lower")]] <- bounds[[1]]
      target[[paste0("beta_", condition, "_interval_upper")]] <- bounds[[2]]
      target[[paste0(condition, "_interval_excludes_zero")]] <-
        bounds[[1]] > 0 || bounds[[2]] < 0
    }
  }
  model$target_tf <- target
  model
}

tf_screening_filter_with_fallback <- function(
  results, target_interval, fallback_target_interval, filter_once
) {
  tf_screening_validate_intervals(target_interval, fallback_target_interval)
  initial <- filter_once(results, target_interval)
  initial_summary <- attr(initial, "filter_summary")
  initial_count <- sum(initial_summary$target_gene_retained)
  successful <- initial_summary$fit_status == "ok"
  valid_models <- vapply(seq_along(results), function(index) {
    if (!successful[[index]]) return(FALSE)
    target <- results[[index]]$target_tf
    fields <- if (isTRUE(results[[index]]$condition_model)) {
      c("beta_control_interval_lower", "beta_control_interval_upper",
        "beta_disease_interval_lower", "beta_disease_interval_upper")
    } else c("screening_interval_lower", "screening_interval_upper")
    values <- unlist(target[fields], use.names = FALSE)
    length(values) == length(fields) && is.numeric(values) && all(is.finite(values)) &&
      all(values[seq.int(1L, length(values), by = 2L)] <=
            values[seq.int(2L, length(values), by = 2L)])
  }, logical(1))
  attempted <- !is.null(fallback_target_interval) && initial_count == 0L && any(valid_models)
  filtered <- initial
  effective_interval <- target_interval
  if (attempted) {
    invalid_genes <- names(results)[successful & !valid_models]
    if (length(invalid_genes)) {
      stop(sprintf(
        "Cannot relax screening with invalid successful-model intervals for: %s.",
        paste(invalid_genes, collapse = ", ")
      ), call. = FALSE)
    }
    copied_models <- results
    for (index in which(valid_models)) {
      copied_models[[index]] <- tf_screening_recompute_target(
        copied_models[[index]], names(results)[[index]], fallback_target_interval,
        attr(results, "inference"), attr(results, "sampling"))
    }
    effective_interval <- fallback_target_interval
    filtered <- filter_once(copied_models, effective_interval)
  }
  final_summary <- attr(filtered, "filter_summary")
  final_count <- sum(final_summary$target_gene_retained)
  reason <- if (is.null(fallback_target_interval)) {
    "disabled"
  } else if (initial_count > 0L) {
    "initial_selection_nonempty"
  } else if (!any(valid_models)) {
    "no_valid_successful_models"
  } else if (final_count > 0L) {
    "fallback_recovered_targets"
  } else "fallback_still_empty"
  audit <- list(
    requested_target_interval = target_interval, effective_target_interval = effective_interval,
    fallback_target_interval = fallback_target_interval, fallback_attempted = attempted,
    fallback_applied = attempted, fallback_succeeded = attempted && final_count > 0L,
    fallback_reason = reason, initial_retained_count = initial_count,
    final_retained_count = final_count
  )
  initial_index <- match(final_summary$target_gene, initial_summary$target_gene)
  final_summary$initial_target_gene_retained <- initial_summary$target_gene_retained[initial_index]
  final_summary$final_target_gene_retained <- final_summary$target_gene_retained
  final_summary$initial_target_supported_in <- initial_summary$target_supported_in[initial_index]
  for (field in names(audit)) {
    value <- audit[[field]]
    final_summary[[field]] <- if (is.null(value)) NA_real_ else value
    attr(filtered, field) <- value
  }
  attr(filtered, "target_interval") <- effective_interval
  attr(filtered, "initial_filter_counts") <- attr(initial, "filter_counts")
  attr(filtered, "initial_filter_summary") <- initial_summary
  attr(filtered, "screening_audit") <- audit
  attr(filtered, "filter_summary") <- final_summary
  filtered
}
