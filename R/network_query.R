#!/usr/bin/env Rscript

# Find every direct (distance = 1) downstream gene of one target TF.
#
# This script deliberately does not use expression, Seurat detection rate, or
# biological-importance filters. Confidence scores are kept on the original
# network scale (for example, values up to 12 are not capped or rescaled).

.query_TF_direct_target_cache <- new.env(parent = emptyenv())

query_TF_direct_target_stop <- function(..., call. = FALSE) {
  stop(sprintf(...), call. = call.)
}

query_TF_direct_target_key <- function(x) {
  toupper(trimws(as.character(x)))
}

query_TF_direct_target_script_dir <- local({
  source_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  source_dir <- if (is.null(source_file) || !nzchar(as.character(source_file[[1]]))) {
    getwd()
  } else {
    dirname(normalizePath(source_file[[1]], winslash = "/", mustWork = FALSE))
  }
  function() source_dir
})

query_TF_direct_target_load_bundle <- function(path, use_cache = TRUE) {
  normalized_path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  cache_key <- paste0("bundle::", normalized_path)

  if (isTRUE(use_cache) && exists(cache_key, envir = .query_TF_direct_target_cache, inherits = FALSE)) {
    return(get(cache_key, envir = .query_TF_direct_target_cache, inherits = FALSE))
  }

  loaded <- new.env(parent = emptyenv())
  load(normalized_path, envir = loaded)
  if (!exists("dagitty_network_bundle", envir = loaded, inherits = FALSE)) {
    query_TF_direct_target_stop(
      "RData file does not contain `dagitty_network_bundle`: %s",
      normalized_path
    )
  }
  bundle <- get("dagitty_network_bundle", envir = loaded, inherits = FALSE)
  if (!is.list(bundle) || !is.data.frame(bundle$edges)) {
    query_TF_direct_target_stop("`dagitty_network_bundle$edges` must be a data.frame.")
  }

  if (isTRUE(use_cache)) {
    assign(cache_key, bundle, envir = .query_TF_direct_target_cache)
  }
  bundle
}

query_TF_direct_target_read_file <- function(path, use_cache = TRUE) {
  normalized_path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  if (grepl("\\.RData$", normalized_path, ignore.case = TRUE)) {
    return(query_TF_direct_target_load_bundle(normalized_path, use_cache = use_cache)$edges)
  }

  cache_key <- paste0("edges::", normalized_path)
  if (isTRUE(use_cache) && exists(cache_key, envir = .query_TF_direct_target_cache, inherits = FALSE)) {
    return(get(cache_key, envir = .query_TF_direct_target_cache, inherits = FALSE))
  }

  extension <- tolower(tools::file_ext(normalized_path))
  edges <- if (identical(extension, "csv")) {
    utils::read.csv(normalized_path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    utils::read.delim(normalized_path, stringsAsFactors = FALSE, check.names = FALSE)
  }
  if (isTRUE(use_cache)) {
    assign(cache_key, edges, envir = .query_TF_direct_target_cache)
  }
  edges
}

query_TF_direct_target_extract_edges <- function(
  network = NULL,
  edge_file = .tfregact_default_network_file(),
  use_cache = TRUE
) {
  if (is.data.frame(network)) {
    return(network)
  }

  if (is.list(network) && !is.null(network$dagitty_bundle)) {
    if (!is.data.frame(network$dagitty_bundle$edges)) {
      query_TF_direct_target_stop("`network$dagitty_bundle$edges` must be a data.frame.")
    }
    return(network$dagitty_bundle$edges)
  }

  if (is.list(network) && is.data.frame(network$edges)) {
    return(network$edges)
  }

  if (is.list(network) && !is.null(network$files$dagitty_rdata)) {
    path <- as.character(network$files$dagitty_rdata[[1]])
    if (nzchar(path) && file.exists(path)) {
      return(query_TF_direct_target_read_file(path, use_cache = use_cache))
    }
  }

  if (is.list(network) && !is.null(network$files$csv)) {
    path <- as.character(network$files$csv[[1]])
    if (nzchar(path) && file.exists(path)) {
      return(query_TF_direct_target_read_file(path, use_cache = use_cache))
    }
  }

  if (is.character(network) && length(network) == 1L && nzchar(network) && file.exists(network)) {
    return(query_TF_direct_target_read_file(network, use_cache = use_cache))
  }

  if (!is.null(network)) {
    query_TF_direct_target_stop(
      paste0(
        "`network` must be an edge data.frame, a bundle/build_tf_network() result, ",
        "or a path to an edge file."
      )
    )
  }

  configured_network <- getOption("TFRegAct.network_path", Sys.getenv("TFREGACT_NETWORK_PATH", unset = ""))
  candidate_paths <- unique(c(
    configured_network,
    edge_file,
    file.path(query_TF_direct_target_script_dir(), edge_file),
    .tfregact_default_network_file(),
    file.path(query_TF_direct_target_script_dir(), "tf_union_output", "TF_Full_Map.RData")
  ))
  candidate_paths <- candidate_paths[!is.na(candidate_paths) & nzchar(candidate_paths)]
  existing <- candidate_paths[file.exists(candidate_paths)]
  if (length(existing) == 0L) {
    query_TF_direct_target_stop(
      "No TF network was found. Supply `network`/`edge_file`, or set options(TFRegAct.network_path = '<path>')."
    )
  }
  query_TF_direct_target_read_file(existing[[1]], use_cache = use_cache)
}

query_TF_direct_target_normalize_edges <- function(edges, target_tf = NULL) {
  required <- c("tf", "target", "confidence_score")
  missing_columns <- setdiff(required, colnames(edges))
  if (length(missing_columns) > 0L) {
    query_TF_direct_target_stop(
      "Network edges are missing required column(s): %s",
      paste(missing_columns, collapse = ", ")
    )
  }

  if (!is.null(target_tf)) {
    tf_values <- trimws(as.character(edges$tf))
    keep_tf <- query_TF_direct_target_key(tf_values) == query_TF_direct_target_key(target_tf)
    keep_tf[is.na(keep_tf)] <- FALSE
    edges <- edges[keep_tf, , drop = FALSE]
    if (nrow(edges) == 0L) {
      query_TF_direct_target_stop("Target TF `%s` was not found in the network.", target_tf)
    }
  }

  keep <- intersect(
    c("tf", "target", "effect", "confidence_score", "supporting_databases"),
    colnames(edges)
  )
  edges <- edges[, keep, drop = FALSE]
  if (!("effect" %in% colnames(edges))) {
    edges$effect <- "unknown"
  }
  if (!("supporting_databases" %in% colnames(edges))) {
    edges$supporting_databases <- ""
  }

  edges$tf <- trimws(as.character(edges$tf))
  edges$target <- trimws(as.character(edges$target))
  edges$effect <- trimws(as.character(edges$effect))
  edges$confidence_score <- suppressWarnings(as.numeric(edges$confidence_score))
  edges$supporting_databases <- trimws(as.character(edges$supporting_databases))
  valid <- !is.na(edges$tf) & !is.na(edges$target) & edges$tf != "" & edges$target != ""
  edges <- edges[valid, , drop = FALSE]
  edges[query_TF_direct_target_key(edges$tf) != query_TF_direct_target_key(edges$target), , drop = FALSE]
}

query_TF_direct_target_effect <- function(x) {
  effect <- tolower(trimws(as.character(x)))
  effect[is.na(effect) | effect == ""] <- "unknown"
  informative <- unique(effect[effect != "unknown"])
  if (length(informative) == 0L) {
    return("unknown")
  }
  if (length(informative) == 1L) {
    return(informative[[1]])
  }
  "mixed"
}

query_TF_direct_target_direction <- function(effect) {
  effect <- tolower(trimws(as.character(effect)))
  out <- integer(length(effect))
  out[effect %in% c("activation", "activate", "positive", "up")] <- 1L
  out[effect %in% c("repression", "repress", "negative", "down")] <- -1L
  out
}

query_TF_direct_target_databases <- function(x) {
  tokens <- unlist(strsplit(as.character(x), "\\s*;\\s*"), use.names = FALSE)
  tokens <- sort(unique(trimws(tokens[!is.na(tokens) & nzchar(trimws(tokens))])))
  paste(tokens, collapse = ";")
}

query_TF_direct_target_confidence <- function(x, target_gene) {
  x <- suppressWarnings(as.numeric(x))
  if (any(is.nan(x)) || any(is.infinite(x) & x < 0)) {
    query_TF_direct_target_stop(
      "Invalid confidence for direct edge to `%s`: NaN and negative Inf are not allowed.",
      target_gene
    )
  }
  finite <- x[is.finite(x)]
  if (length(finite) > 0L) {
    if (any(finite < 0)) {
      query_TF_direct_target_stop(
        "Invalid confidence for direct edge to `%s`: finite values must be non-negative.",
        target_gene
      )
    }
    return(max(finite))
  }
  if (any(is.infinite(x) & x > 0)) {
    return(Inf)
  }
  query_TF_direct_target_stop(
    "Direct edge to `%s` has no usable confidence value.",
    target_gene
  )
}

query_TF_direct_target_collapse <- function(direct_edges, resolved_tf) {
  gene_key <- query_TF_direct_target_key(direct_edges$target)

  if (!anyDuplicated(gene_key)) {
    confidence <- suppressWarnings(as.numeric(direct_edges$confidence_score))
    if (any(is.na(confidence) | is.nan(confidence))) {
      bad_gene <- direct_edges$target[which(is.na(confidence) | is.nan(confidence))[[1]]]
      query_TF_direct_target_stop(
        "Direct edge to `%s` has no usable confidence value.",
        bad_gene
      )
    }
    if (any(is.infinite(confidence) & confidence < 0) ||
        any(is.finite(confidence) & confidence < 0)) {
      bad_gene <- direct_edges$target[
        which((is.infinite(confidence) & confidence < 0) |
          (is.finite(confidence) & confidence < 0))[[1]]
      ]
      query_TF_direct_target_stop(
        "Invalid confidence for direct edge to `%s`: values must be non-negative.",
        bad_gene
      )
    }
    effect <- tolower(trimws(as.character(direct_edges$effect)))
    effect[is.na(effect) | effect == ""] <- "unknown"
    databases <- trimws(as.character(direct_edges$supporting_databases))
    databases[is.na(databases)] <- ""
    return(data.frame(
      target_tf = resolved_tf,
      target_gene = direct_edges$target,
      distance = 1L,
      effect = effect,
      direction = query_TF_direct_target_direction(effect),
      supporting_databases = databases,
      query_confidence = confidence,
      stringsAsFactors = FALSE
    ))
  }

  groups <- split(seq_len(nrow(direct_edges)), gene_key)

  rows <- lapply(groups, function(index) {
    block <- direct_edges[index, , drop = FALSE]
    gene_names <- block$target
    preferred_gene <- names(sort(table(gene_names), decreasing = TRUE))[[1]]
    collapsed_effect <- query_TF_direct_target_effect(block$effect)
    data.frame(
      target_tf = resolved_tf,
      target_gene = preferred_gene,
      distance = 1L,
      effect = collapsed_effect,
      direction = query_TF_direct_target_direction(collapsed_effect),
      supporting_databases = query_TF_direct_target_databases(block$supporting_databases),
      query_confidence = query_TF_direct_target_confidence(
        block$confidence_score,
        target_gene = preferred_gene
      ),
      stringsAsFactors = FALSE
    )
  })
  rownames_out <- NULL
  out <- do.call(rbind, rows)
  rownames(out) <- rownames_out
  out
}

query_TF_direct_target_override_map <- function(x, target_genes) {
  out <- rep(NA_real_, length(target_genes))
  names(out) <- target_genes
  if (is.null(x)) {
    return(out)
  }

  if (is.data.frame(x)) {
    gene_columns <- intersect(c("target_gene", "gene", "target"), colnames(x))
    value_columns <- intersect(
      c("manual_confidence", "confidence", "confidence_override"),
      colnames(x)
    )
    if (length(gene_columns) == 0L || length(value_columns) == 0L) {
      query_TF_direct_target_stop(
        paste0(
          "`target_confidence_override` data.frame must contain a target-gene column ",
          "and a confidence column."
        )
      )
    }
    supplied_genes <- trimws(as.character(x[[gene_columns[[1]]]]))
    supplied_values <- suppressWarnings(as.numeric(x[[value_columns[[1]]]]))
  } else if (is.numeric(x) || is.integer(x)) {
    supplied_values <- as.numeric(x)
    supplied_genes <- names(x)
    if (is.null(supplied_genes) || any(!nzchar(supplied_genes))) {
      if (length(supplied_values) == length(target_genes)) {
        supplied_genes <- target_genes
      } else if (length(supplied_values) == 1L && length(target_genes) == 1L) {
        supplied_genes <- target_genes
      } else {
        query_TF_direct_target_stop(
          paste0(
            "Unnamed `target_confidence_override` must contain one value per discovered ",
            "target gene; otherwise use a named vector."
          )
        )
      }
    }
  } else {
    query_TF_direct_target_stop(
      "`target_confidence_override` must be NULL, a named numeric vector, or a data.frame."
    )
  }

  supplied_keys <- query_TF_direct_target_key(supplied_genes)
  target_keys <- query_TF_direct_target_key(target_genes)
  if (any(is.na(supplied_keys)) || any(supplied_keys == "") || anyDuplicated(supplied_keys)) {
    query_TF_direct_target_stop(
      "`target_confidence_override` contains empty or duplicated target-gene names."
    )
  }
  match_index <- match(supplied_keys, target_keys)
  if (any(is.na(match_index))) {
    query_TF_direct_target_stop(
      "Unknown target gene(s) in `target_confidence_override`: %s",
      paste(supplied_genes[is.na(match_index)], collapse = ", ")
    )
  }
  if (length(supplied_values) != length(match_index)) {
    query_TF_direct_target_stop(
      "Gene and confidence columns in `target_confidence_override` have different lengths."
    )
  }
  out[match_index] <- supplied_values
  out
}

#' Query all directly regulated genes for one TF
#'
#' @param target_tf One TF symbol.
#' @param network Optional edge data.frame, network bundle/build result, or file path.
#' @param edge_file Default network bundle used when `network` is NULL.
#' @param confidence_threshold Minimum final confidence; defaults to 4.
#' @param target_confidence_override Optional named numeric vector or data.frame.
#'   A supplied value replaces both the finite query value and the automatic
#'   no-path default for that gene.
#' @param edit_target_confidence Open an interactive confidence-editing window.
#' @param return_all If FALSE, return only genes passing the confidence threshold;
#'   if TRUE, return all direct genes and retain the `included` audit column.
#' @param use_cache Cache file-backed network data during the current R session.
#'
#' @return A data.frame containing direct TF-gene edges and confidence auditing.
query_TF_direct_target_genes <- function(
  target_tf,
  network = NULL,
  edge_file = .tfregact_default_network_file(),
  confidence_threshold = 4,
  target_confidence_override = NULL,
  edit_target_confidence = FALSE,
  return_all = FALSE,
  use_cache = TRUE
) {
  if (length(target_tf) != 1L || is.na(target_tf) || !nzchar(trimws(as.character(target_tf)))) {
    query_TF_direct_target_stop("`target_tf` must be one non-empty TF symbol.")
  }
  target_tf <- trimws(as.character(target_tf))
  confidence_threshold <- suppressWarnings(as.numeric(confidence_threshold[[1]]))
  if (!is.finite(confidence_threshold) || confidence_threshold < 0) {
    query_TF_direct_target_stop("`confidence_threshold` must be finite and non-negative.")
  }

  direct_edges <- query_TF_direct_target_normalize_edges(
    edges = query_TF_direct_target_extract_edges(
      network = network,
      edge_file = edge_file,
      use_cache = use_cache
    ),
    target_tf = target_tf
  )
  if (nrow(direct_edges) == 0L) {
    empty <- data.frame(
      target_tf = character(0),
      target_gene = character(0),
      distance = integer(0),
      effect = character(0),
      direction = integer(0),
      supporting_databases = character(0),
      query_confidence = numeric(0),
      automatic_confidence = numeric(0),
      manual_confidence = numeric(0),
      final_confidence = numeric(0),
      confidence_source = character(0),
      included = logical(0),
      stringsAsFactors = FALSE
    )
    attr(empty, "query_summary") <- list(
      target_tf = target_tf,
      direct_genes_found = 0L,
      genes_returned = 0L,
      confidence_threshold = confidence_threshold,
      no_path_default = 2,
      expression_filter_used = FALSE
    )
    return(empty)
  }
  tf_names <- direct_edges$tf
  exact <- unique(tf_names[tf_names == target_tf])
  resolved_tf <- if (length(exact) > 0L) {
    exact[[1]]
  } else {
    names(sort(table(tf_names), decreasing = TRUE))[[1]]
  }

  audit <- query_TF_direct_target_collapse(direct_edges, resolved_tf = resolved_tf)
  no_path <- is.infinite(audit$query_confidence) & audit$query_confidence > 0
  audit$automatic_confidence <- audit$query_confidence
  audit$automatic_confidence[no_path] <- 2

  manual <- query_TF_direct_target_override_map(
    target_confidence_override,
    target_genes = audit$target_gene
  )
  editor_table <- data.frame(
    target_gene = audit$target_gene,
    query_confidence = audit$query_confidence,
    automatic_confidence = audit$automatic_confidence,
    manual_confidence = manual,
    stringsAsFactors = FALSE
  )

  if (isTRUE(edit_target_confidence)) {
    if (!interactive()) {
      query_TF_direct_target_stop(
        "`edit_target_confidence = TRUE` requires an interactive R session."
      )
    }
    edited <- utils::edit(editor_table)
    if (!is.data.frame(edited) || !("target_gene" %in% colnames(edited)) ||
        !("manual_confidence" %in% colnames(edited)) ||
        !identical(
          query_TF_direct_target_key(edited$target_gene),
          query_TF_direct_target_key(audit$target_gene)
        )) {
      query_TF_direct_target_stop(
        "The confidence editor must preserve every target gene and its order."
      )
    }
    manual <- suppressWarnings(as.numeric(edited$manual_confidence))
  }

  manual_used <- !is.na(manual)
  if (any(manual_used & (!is.finite(manual) | manual < 0))) {
    query_TF_direct_target_stop(
      "Manual confidence values must be finite and non-negative; leave unused entries as NA."
    )
  }

  audit$manual_confidence <- unname(manual)
  audit$final_confidence <- audit$automatic_confidence
  audit$final_confidence[manual_used] <- manual[manual_used]
  audit$confidence_source <- ifelse(no_path, "no_path_default", "query")
  audit$confidence_source[manual_used] <- "manual"
  audit$included <- audit$final_confidence >= confidence_threshold

  audit <- audit[order(-audit$final_confidence, audit$target_gene), , drop = FALSE]
  rownames(audit) <- NULL
  result <- if (isTRUE(return_all)) audit else audit[audit$included, , drop = FALSE]
  rownames(result) <- NULL

  attr(result, "query_summary") <- list(
    target_tf = resolved_tf,
    direct_genes_found = nrow(audit),
    genes_returned = nrow(result),
    confidence_threshold = confidence_threshold,
    no_path_default = 2,
    expression_filter_used = FALSE
  )
  result
}
