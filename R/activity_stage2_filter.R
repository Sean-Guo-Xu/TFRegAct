#!/usr/bin/env Rscript

# Stage 2: keep a target-gene model when the target-TF posterior interval
# excludes zero in control, disease, or both. Confounders entering this MCMC
# stage are never filtered again and propagate unchanged to EM.
filter_TF_mcmc_gene_results <- function(
  screening_results,
  target_interval = 0.90,
  confounder_interval = 0.90,
  output_file = NULL
) {
  if (!is.list(screening_results) || !length(screening_results)) {
    stop("`screening_results` must be a non-empty Stage 1 result list.",
         call. = FALSE)
  }
  if (!identical(attr(screening_results, "inference"), "mcmc") ||
      !identical(attr(screening_results, "prior_family"), "normal")) {
    stop(
      "Stage 2 filtering requires MCMC results from the Normal-prior model.",
      call. = FALSE
    )
  }

  required_columns <- c(
    "tf", "role", "screening_interval_level",
    "screening_interval_lower", "screening_interval_upper",
    "screening_interval_excludes_zero"
  )
  required_target_columns <- c(
    "beta_control_interval_lower", "beta_control_interval_upper",
    "beta_disease_interval_lower", "beta_disease_interval_upper",
    "control_interval_excludes_zero",
    "disease_interval_excludes_zero"
  )
  filtered <- list()
  summary_rows <- vector("list", length(screening_results))

  for (i in seq_along(screening_results)) {
    gene <- names(screening_results)[[i]]
    model <- screening_results[[i]]
    if (!is.list(model) || !identical(model$status, "ok")) {
      # Do not let one failed MCMC regression abort the complete TF activity
      # analysis.  It is excluded from EM and retained in the audit summary.
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
          "Missing or invalid MCMC result."
        },
        target_interval = target_interval,
        target_control_interval_lower = NA_real_,
        target_control_interval_upper = NA_real_,
        target_disease_interval_lower = NA_real_,
        target_disease_interval_upper = NA_real_,
        target_supported_in = "fit_failed",
        target_gene_retained = FALSE,
        confounder_interval = confounder_interval,
        confounders_before = 0L,
        confounders_supported = 0L,
        confounders_retained = 0L,
        confounders_removed_by_interval = 0L,
        stringsAsFactors = FALSE
      )
      next
    }
    if (nrow(model$target_tf) != 1L ||
        !all(required_columns %in% names(model$target_tf)) ||
        !all(required_target_columns %in% names(model$target_tf)) ||
        !all(required_columns %in% names(model$confounders))) {
      stop(sprintf("Stage 1 model `%s` has an invalid parameter table.", gene),
           call. = FALSE)
    }
    if (!is.list(model$stage2_parameter_draws) ||
        !identical(
          model$stage2_parameter_draws$interface_version,
          "tf_em_stage2_parameter_draws_v1"
        )) {
      stop(sprintf(
        "Stage 1 model `%s` does not contain the EM Stage 2 draw interface.",
        gene
      ), call. = FALSE)
    }

    target_level_ok <- isTRUE(all.equal(
      model$target_tf$screening_interval_level[[1]],
      target_interval,
      tolerance = 1e-12
    ))
    if (!target_level_ok) {
      stop(sprintf(
        "Target interval in `%s` is not the requested %.0f%% interval.",
        gene, 100 * target_interval
      ), call. = FALSE)
    }

    if (nrow(model$confounders)) {
      confounder_level_ok <- vapply(
        model$confounders$screening_interval_level,
        function(level) isTRUE(all.equal(
          level, confounder_interval, tolerance = 1e-12
        )),
        logical(1)
      )
      if (!all(confounder_level_ok)) {
        stop(sprintf(
          "One or more confounder intervals in `%s` are not %.0f%% intervals.",
          gene, 100 * confounder_interval
        ), call. = FALSE)
      }
    }

    target_keep_control <- isTRUE(
      model$target_tf$control_interval_excludes_zero[[1]]
    )
    target_keep_disease <- isTRUE(
      model$target_tf$disease_interval_excludes_zero[[1]]
    )
    target_keep <- target_keep_control || target_keep_disease
    target_support <- if (target_keep_control && target_keep_disease) {
      "both"
    } else if (target_keep_control) {
      "control_only"
    } else if (target_keep_disease) {
      "disease_only"
    } else {
      "neither"
    }
    confounder_supported <- if (nrow(model$confounders)) {
      !is.na(model$confounders$screening_interval_excludes_zero) &
        model$confounders$screening_interval_excludes_zero
    } else {
      logical(0)
    }

    if (target_keep) {
      filtered[[gene]] <- list(
        target_tf = model$target_tf,
        confounders = model$confounders,
        stage2_parameter_draws = model$stage2_parameter_draws,
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
      target_control_interval_lower =
        model$target_tf$beta_control_interval_lower[[1]],
      target_control_interval_upper =
        model$target_tf$beta_control_interval_upper[[1]],
      target_disease_interval_lower =
        model$target_tf$beta_disease_interval_lower[[1]],
      target_disease_interval_upper =
        model$target_tf$beta_disease_interval_upper[[1]],
      target_supported_in = target_support,
      target_gene_retained = target_keep,
      confounder_interval = confounder_interval,
      confounders_before = nrow(model$confounders),
      confounders_supported = sum(confounder_supported),
      confounders_retained = if (target_keep) nrow(model$confounders) else 0L,
      confounders_removed_by_interval = 0L,
      stringsAsFactors = FALSE
    )
  }

  filter_summary <- do.call(rbind, summary_rows)
  rownames(filter_summary) <- NULL
  filter_counts <- c(
    gene_models_before = nrow(filter_summary),
    gene_models_failed = sum(filter_summary$fit_status != "ok"),
    gene_models_retained = sum(filter_summary$target_gene_retained),
    gene_models_removed = sum(!filter_summary$target_gene_retained),
    confounders_before = sum(filter_summary$confounders_before),
    confounders_retained = sum(filter_summary$confounders_retained),
    confounders_removed_by_interval = 0L,
    confounders_in_removed_gene_models = sum(
      filter_summary$confounders_before[!filter_summary$target_gene_retained]
    )
  )

  class(filtered) <- c("FilteredTFStage1ScreeningList", "list")
  attr(filtered, "pipeline_stage") <- "mcmc"
  attr(filtered, "target_tf") <- attr(screening_results, "target_tf")
  attr(filtered, "model_version") <- attr(screening_results, "model_version")
  attr(filtered, "target_interval") <- target_interval
  attr(filtered, "confounder_interval") <- confounder_interval
  attr(filtered, "confounders_filtered") <- FALSE
  attr(filtered, "prescreen_confounders_filtered") <- identical(
    attr(screening_results, "input_pipeline_stage"),
    "mcmc_input"
  )
  attr(filtered, "prescreen_prior_family") <-
    attr(screening_results, "prescreen_prior_family")
  attr(filtered, "prescreen_target_interval") <-
    attr(screening_results, "prescreen_target_interval")
  attr(filtered, "prescreen_confounder_interval") <-
    attr(screening_results, "prescreen_confounder_interval")
  attr(filtered, "prescreen_filter_counts") <-
    attr(screening_results, "prescreen_filter_counts")
  attr(filtered, "filter_counts") <- filter_counts
  attr(filtered, "filter_summary") <- filter_summary

  if (is.null(output_file)) {
    input_file <- attr(screening_results, "output_file")
    output_dir <- if (!is.null(input_file) && nzchar(input_file[[1]])) {
      dirname(input_file[[1]])
    } else {
      file.path(getwd(), "adjustment_output")
    }
    target_tf <- attr(screening_results, "target_tf")
    if (is.null(target_tf) || !nzchar(target_tf[[1]])) target_tf <- "target_TF"
    output_file <- file.path(
      output_dir,
      sprintf(
        "%s_stage1_screening_stage2_ready_filtered_results.rds",
        target_tf[[1]]
      )
    )
  }
  output_file <- normalizePath(output_file, winslash = "/", mustWork = FALSE)
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  attr(filtered, "output_file") <- output_file
  saveRDS(filtered, output_file)
  filtered
}

# Backward-compatible name. Its behavior now follows the three-stage design:
# it accepts only Normal-prior MCMC results and never filters confounders.
filter_TF_stage1_screening_results <- function(
  screening_results,
  target_interval = 0.90,
  confounder_interval = 0.90,
  output_file = NULL
) {
  filter_TF_mcmc_gene_results(
    screening_results = screening_results,
    target_interval = target_interval,
    confounder_interval = confounder_interval,
    output_file = output_file
  )
}
