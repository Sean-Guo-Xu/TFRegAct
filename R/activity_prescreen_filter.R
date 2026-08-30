#!/usr/bin/env Rscript

tf_prescreen_stop <- function(...) {
  stop(sprintf(...), call. = FALSE)
}

# Stage 1: retain a target-gene model when the target TF is supported in
# control or disease, and retain only confounders whose own interval excludes
# zero. Both intervals default to 90%.
filter_TF_prescreen_results <- function(
  prescreen_results,
  target_interval = 0.90,
  confounder_interval = 0.90,
  output_file = NULL
) {
  if (!is.list(prescreen_results) || !length(prescreen_results)) {
    tf_prescreen_stop("`prescreen_results` must be a non-empty result list.")
  }
  if (!identical(attr(prescreen_results, "inference"), "variational") ||
      !identical(attr(prescreen_results, "prior_family"), "laplace")) {
    tf_prescreen_stop(
      "Stage 1 prescreening requires variational inference with a Laplace prior."
    )
  }

  filtered <- list()
  summary_rows <- vector("list", length(prescreen_results))
  confounder_rows <- list()

  for (i in seq_along(prescreen_results)) {
    gene <- names(prescreen_results)[[i]]
    model <- prescreen_results[[i]]
    if (!is.list(model) || !identical(model$status, "ok") ||
        nrow(model$target_tf) != 1L) {
      # A single VI failure (commonly a near-all-zero response gene) must not
      # abort screening for every other downstream target.  Preserve its
      # reason in the audit table and exclude it from Stage 2.
      summary_rows[[i]] <- data.frame(
        target_gene = gene,
        target_tf = if (is.list(model) && !is.null(model$target_tf_name)) {
          as.character(model$target_tf_name[[1]])
        } else {
          NA_character_
        },
        fit_status = if (is.list(model) && !is.null(model$status)) {
          as.character(model$status[[1]])
        } else {
          "invalid_result"
        },
        fit_error = if (is.list(model) && !is.null(model$error)) {
          as.character(model$error[[1]])
        } else {
          "Missing or invalid prescreen result."
        },
        target_interval = target_interval,
        target_overall_interval_lower = NA_real_,
        target_overall_interval_upper = NA_real_,
        target_control_interval_lower = NA_real_,
        target_control_interval_upper = NA_real_,
        target_disease_interval_lower = NA_real_,
        target_disease_interval_upper = NA_real_,
        target_supported_in = "fit_failed",
        target_gene_retained = FALSE,
        confounder_interval = confounder_interval,
        confounders_before = 0L,
        confounders_retained = 0L,
        confounders_removed_by_interval = 0L,
        stringsAsFactors = FALSE
      )
      next
    }
    target_level_ok <- isTRUE(all.equal(
      model$target_tf$screening_interval_level[[1]],
      target_interval,
      tolerance = 1e-12
    ))
    if (!target_level_ok) {
      tf_prescreen_stop(
        "Target interval in `%s` is not the requested %.0f%% interval.",
        gene,
        100 * target_interval
      )
    }
    if (nrow(model$confounders)) {
      level_ok <- vapply(
        model$confounders$screening_interval_level,
        function(x) isTRUE(all.equal(
          x, confounder_interval, tolerance = 1e-12
        )),
        logical(1)
      )
      if (!all(level_ok)) {
        tf_prescreen_stop(
          "Confounder intervals in `%s` are not all %.0f%% intervals.",
          gene,
          100 * confounder_interval
        )
      }
    }

    has_condition <- isTRUE(model$condition_model)
    keep_overall <- isTRUE(
      model$target_tf$screening_interval_excludes_zero[[1]]
    )
    keep_control <- has_condition && isTRUE(
      model$target_tf$control_interval_excludes_zero[[1]]
    )
    keep_disease <- has_condition && isTRUE(
      model$target_tf$disease_interval_excludes_zero[[1]]
    )
    keep_gene <- if (has_condition) {
      keep_control || keep_disease
    } else {
      keep_overall
    }
    support <- if (!has_condition && keep_overall) {
      "overall"
    } else if (!has_condition) {
      "neither"
    } else if (keep_control && keep_disease) {
      "both"
    } else if (keep_control) {
      "control_only"
    } else if (keep_disease) {
      "disease_only"
    } else {
      "neither"
    }
    keep_confounder <- if (nrow(model$confounders)) {
      !is.na(model$confounders$screening_interval_excludes_zero) &
        model$confounders$screening_interval_excludes_zero
    } else {
      logical(0)
    }
    retained_confounders <- model$confounders[keep_confounder, , drop = FALSE]

    if (keep_gene) {
      filtered[[gene]] <- list(
        target_tf = model$target_tf,
        confounders = retained_confounders,
        control_level = model$control_level,
        disease_level = model$disease_level
      )
    }

    summary_rows[[i]] <- data.frame(
      target_gene = gene,
      target_tf = model$target_tf$tf[[1]],
      fit_status = "ok",
      fit_error = NA_character_,
      target_interval = target_interval,
      target_overall_interval_lower =
        model$target_tf$screening_interval_lower[[1]],
      target_overall_interval_upper =
        model$target_tf$screening_interval_upper[[1]],
      target_control_interval_lower =
        if (has_condition) model$target_tf$beta_control_interval_lower[[1]] else NA_real_,
      target_control_interval_upper =
        if (has_condition) model$target_tf$beta_control_interval_upper[[1]] else NA_real_,
      target_disease_interval_lower =
        if (has_condition) model$target_tf$beta_disease_interval_lower[[1]] else NA_real_,
      target_disease_interval_upper =
        if (has_condition) model$target_tf$beta_disease_interval_upper[[1]] else NA_real_,
      target_supported_in = support,
      target_gene_retained = keep_gene,
      confounder_interval = confounder_interval,
      confounders_before = nrow(model$confounders),
      confounders_retained = if (keep_gene) sum(keep_confounder) else 0L,
      confounders_removed_by_interval = if (keep_gene) {
        sum(!keep_confounder)
      } else {
        0L
      },
      stringsAsFactors = FALSE
    )
    if (nrow(model$confounders)) {
      confounder_rows[[gene]] <- data.frame(
        target_gene = gene,
        tf = model$confounders$tf,
        interval_level = model$confounders$screening_interval_level,
        interval_lower = model$confounders$screening_interval_lower,
        interval_upper = model$confounders$screening_interval_upper,
        retained = keep_gene & keep_confounder,
        gene_model_retained = keep_gene,
        stringsAsFactors = FALSE
      )
    }
  }

  filter_summary <- do.call(rbind, summary_rows)
  rownames(filter_summary) <- NULL
  confounder_summary <- if (length(confounder_rows)) {
    result <- do.call(rbind, confounder_rows)
    rownames(result) <- NULL
    result
  } else {
    data.frame()
  }
  filter_counts <- c(
    gene_models_before = nrow(filter_summary),
    gene_models_failed = sum(filter_summary$fit_status != "ok"),
    gene_models_retained = sum(filter_summary$target_gene_retained),
    gene_models_removed = sum(!filter_summary$target_gene_retained),
    confounders_before = sum(filter_summary$confounders_before),
    confounders_retained = sum(filter_summary$confounders_retained),
    confounders_removed_by_interval =
      sum(filter_summary$confounders_removed_by_interval),
    confounders_in_removed_gene_models = sum(
      filter_summary$confounders_before[!filter_summary$target_gene_retained]
    )
  )

  class(filtered) <- c("FilteredTFPrescreenList", "list")
  attr(filtered, "pipeline_stage") <- "prescreen"
  attr(filtered, "target_tf") <- attr(prescreen_results, "target_tf")
  attr(filtered, "model_version") <- attr(prescreen_results, "model_version")
  attr(filtered, "prior_family") <- "laplace"
  attr(filtered, "target_interval") <- target_interval
  attr(filtered, "confounder_interval") <- confounder_interval
  attr(filtered, "confounders_filtered") <- TRUE
  attr(filtered, "filter_counts") <- filter_counts
  attr(filtered, "filter_summary") <- filter_summary
  attr(filtered, "confounder_summary") <- confounder_summary

  if (!is.null(output_file)) {
    output_file <- normalizePath(output_file, winslash = "/", mustWork = FALSE)
    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
    attr(filtered, "output_file") <- output_file
    saveRDS(filtered, output_file)
  }
  filtered
}

# Convert the Stage 1 selection into the input for Stage 2 MCMC. The target TF
# is always retained. Every predictor-aligned field is subset in the same order.
build_TF_mcmc_input_from_prescreen <- function(
  prescreen_filtered_results,
  screening_input,
  direction_effect = 0.2,
  beta_sd_floor = 0.5,
  stage1_sd_multiplier = 1.5,
  alpha_prior_sd = 1,
  output_file = NULL
) {
  if (!is.list(prescreen_filtered_results) ||
      !length(prescreen_filtered_results) ||
      !identical(attr(prescreen_filtered_results, "pipeline_stage"), "prescreen") ||
      !identical(attr(prescreen_filtered_results, "confounders_filtered"), TRUE)) {
    tf_prescreen_stop("Invalid Stage 1 prescreen result object.")
  }
  if (!is.list(screening_input) || !length(screening_input)) {
    tf_prescreen_stop("`screening_input` must be a non-empty nested list.")
  }
  if (!is.numeric(alpha_prior_sd) || length(alpha_prior_sd) != 1L ||
      !is.finite(alpha_prior_sd) || alpha_prior_sd <= 0) {
    tf_prescreen_stop("`alpha_prior_sd` must be one finite positive number.")
  }

  target_genes <- names(prescreen_filtered_results)
  if (any(!target_genes %in% names(screening_input))) {
    tf_prescreen_stop("One or more retained genes are absent from the input.")
  }
  result <- setNames(vector("list", length(target_genes)), target_genes)

  for (gene in target_genes) {
    entry <- screening_input[[gene]]
    selection <- prescreen_filtered_results[[gene]]
    if (!identical(entry$status, "ready")) {
      tf_prescreen_stop("Input model `%s` is not ready.", gene)
    }
    keep_tf <- c(
      as.character(selection$target_tf$tf[[1]]),
      as.character(selection$confounders$tf)
    )
    keep_index <- match(toupper(keep_tf), toupper(entry$feature_names))
    if (anyNA(keep_index) || anyDuplicated(keep_index)) {
      tf_prescreen_stop("Predictor selection is inconsistent for `%s`.", gene)
    }
    # Preserve the original design-matrix column order.
    keep_index <- sort(keep_index)
    old_target_index <- as.integer(entry$target_tf_index)
    if (!old_target_index %in% keep_index) {
      tf_prescreen_stop("Target TF was removed from `%s`.", gene)
    }

    entry$feature_names <- entry$feature_names[keep_index]
    entry$predictor_center <- entry$predictor_center[keep_index]
    entry$predictor_scale <- entry$predictor_scale[keep_index]
    edge_index <- match(
      toupper(entry$feature_names),
      toupper(entry$predictor_edges$tf)
    )
    if (anyNA(edge_index)) {
      tf_prescreen_stop("Predictor edges are inconsistent for `%s`.", gene)
    }
    entry$predictor_edges <- entry$predictor_edges[edge_index, , drop = FALSE]
    entry$stan_data$X <- entry$stan_data$X[, keep_index, drop = FALSE]
    entry$stan_data$confidence <- entry$stan_data$confidence[keep_index]
    entry$stan_data$direction <- entry$stan_data$direction[keep_index]
    entry$stan_data$P <- length(keep_index)
    entry$target_tf_index <- match(old_target_index, keep_index)
    entry$stan_data$target_tf_index <- entry$target_tf_index

    stage1_table <- rbind(selection$target_tf, selection$confounders)
    stage1_index <- match(toupper(entry$feature_names), toupper(stage1_table$tf))
    if (anyNA(stage1_index) ||
        !all(c("beta_mean", "beta_sd") %in% names(stage1_table))) {
      tf_prescreen_stop(
        "Stage 1 posterior summaries are incomplete for `%s`.", gene
      )
    }
    stage1_summary <- stage1_table[stage1_index, , drop = FALSE]
    stage2_prior <- tf_stage2_prior_from_stage1(
      stage1_beta_mean = stage1_summary$beta_mean,
      stage1_beta_sd = stage1_summary$beta_sd,
      direction = entry$stan_data$direction,
      direction_effect = direction_effect,
      beta_sd_floor = beta_sd_floor,
      stage1_sd_multiplier = stage1_sd_multiplier
    )
    entry$stan_data$beta_prior_mean <- stage2_prior$beta_prior_mean
    entry$stan_data$beta_prior_sd <- stage2_prior$beta_prior_sd
    entry$stan_data$beta_init <- stage2_prior$stage1_beta_mean
    entry$stan_data$alpha_prior_sd <- as.numeric(alpha_prior_sd)
    entry$stage1_beta_summary <- data.frame(
      tf = entry$feature_names,
      stage1_beta_mean = stage2_prior$stage1_beta_mean,
      stage1_beta_sd = stage2_prior$stage1_beta_sd,
      direction = stage2_prior$direction,
      stringsAsFactors = FALSE
    )
    entry$stage2_prior <- data.frame(
      tf = entry$feature_names,
      beta_prior_mean = stage2_prior$beta_prior_mean,
      beta_prior_sd = stage2_prior$beta_prior_sd,
      stringsAsFactors = FALSE
    )
    entry$prescreen_original_predictor_count <- length(keep_index) +
      attr(prescreen_filtered_results, "filter_summary")$
        confounders_removed_by_interval[
          match(gene, attr(prescreen_filtered_results, "filter_summary")$target_gene)
        ]
    entry$prescreen_retained_confounder_count <- nrow(selection$confounders)
    entry$adjustment_tf_count <- nrow(selection$confounders)
    result[[gene]] <- entry
  }

  source_attributes <- attributes(screening_input)
  for (attribute_name in setdiff(names(source_attributes), c("names", "class"))) {
    attr(result, attribute_name) <- source_attributes[[attribute_name]]
  }
  class(result) <- unique(c("TFMCMCInputFromPrescreen", class(screening_input)))
  attr(result, "pipeline_stage") <- "mcmc_input"
  attr(result, "prescreen_prior_family") <- "laplace"
  attr(result, "prescreen_target_interval") <-
    attr(prescreen_filtered_results, "target_interval")
  attr(result, "prescreen_confounder_interval") <-
    attr(prescreen_filtered_results, "confounder_interval")
  attr(result, "prescreen_filter_counts") <-
    attr(prescreen_filtered_results, "filter_counts")
  attr(result, "model_version") <- tf_shared_stage2_model_version()
  attr(result, "stan_file") <- "TF_stage2_directional_nb_model.stan"
  attr(result, "stage2_prior_config") <- list(
    direction_effect = as.numeric(direction_effect),
    beta_sd_floor = as.numeric(beta_sd_floor),
    stage1_sd_multiplier = as.numeric(stage1_sd_multiplier),
    alpha_prior_sd = as.numeric(alpha_prior_sd)
  )

  if (!is.null(output_file)) {
    output_file <- normalizePath(output_file, winslash = "/", mustWork = FALSE)
    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
    attr(result, "output_file") <- output_file
    saveRDS(result, output_file)
  }
  result
}
