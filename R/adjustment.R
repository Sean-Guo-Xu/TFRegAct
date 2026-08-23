#!/usr/bin/env Rscript

find_adjustment_dagitty_get_arg <- function(args, flag, default = NULL) {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(hit) == 0) {
    return(default)
  }
  sub(paste0("^", flag, "="), "", hit[[1]])
}

find_adjustment_dagitty_install_if_missing <- function(pkg) {
  .tfregact_require(pkg)
}

find_adjustment_dagitty_load_bundle_file <- function(path) {
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  if (!exists("dagitty_network_bundle", envir = env, inherits = FALSE)) {
    stop(sprintf("Bundle file does not contain `dagitty_network_bundle`: %s", path))
  }
  get("dagitty_network_bundle", envir = env, inherits = FALSE)
}

find_adjustment_dagitty_extract_bundle <- function(network = NULL, edge_file = NULL) {
  if (is.null(edge_file) || !length(edge_file) || !nzchar(as.character(edge_file[[1]]))) {
    edge_file <- getOption("TFRegAct.network_path", Sys.getenv("TFREGACT_NETWORK_PATH", unset = ""))
  }
  if (!is.null(network)) {
    if (is.list(network) && !is.null(network$dagitty_bundle)) {
      return(network$dagitty_bundle)
    }

    if (is.list(network) && !is.null(network$files) && !is.null(network$files$dagitty_rdata)) {
      bundle_path <- as.character(network$files$dagitty_rdata[[1]])
      if (nzchar(bundle_path) && file.exists(bundle_path)) {
        return(find_adjustment_dagitty_load_bundle_file(bundle_path))
      }
    }

    if (is.character(network) && length(network) == 1 && grepl("\\.RData$", network, ignore.case = TRUE) && file.exists(network)) {
      return(find_adjustment_dagitty_load_bundle_file(network))
    }
  }

  if (!is.null(edge_file) && nzchar(edge_file) && grepl("\\.RData$", edge_file, ignore.case = TRUE) && file.exists(edge_file)) {
    return(find_adjustment_dagitty_load_bundle_file(edge_file))
  }

  default_bundle <- "tf_union_output/TF_Full_Map.RData"
  if (file.exists(default_bundle)) {
    return(find_adjustment_dagitty_load_bundle_file(default_bundle))
  }

  NULL
}

find_adjustment_dagitty_extract_edge_file <- function(network = NULL, edge_file = NULL) {
  if (is.null(edge_file) || !length(edge_file) || !nzchar(as.character(edge_file[[1]]))) {
    edge_file <- getOption("TFRegAct.network_path", Sys.getenv("TFREGACT_NETWORK_PATH", unset = ""))
  }
  if (!is.null(network)) {
    if (is.character(network) && length(network) == 1 && nzchar(network) && file.exists(network)) {
      if (grepl("\\.RData$", network, ignore.case = TRUE)) {
        return(NULL)
      }
      return(normalizePath(network, winslash = "/", mustWork = TRUE))
    }

    if (is.list(network) && !is.null(network$files) && !is.null(network$files$csv)) {
      csv_path <- as.character(network$files$csv[[1]])
      if (nzchar(csv_path) && file.exists(csv_path)) {
        return(normalizePath(csv_path, winslash = "/", mustWork = TRUE))
      }
    }

    stop("`network` must be a valid edge-file path or a build_tf_network() result with `files$csv`.")
  }

  if (!is.null(edge_file) && nzchar(edge_file) && file.exists(edge_file)) {
    if (grepl("\\.RData$", edge_file, ignore.case = TRUE)) {
      return(NULL)
    }
    return(normalizePath(edge_file, winslash = "/", mustWork = TRUE))
  }

  default_file <- "tf_union_output/tf_gene_merged_human_weighted_clean.csv"
  if (file.exists(default_file)) {
    return(normalizePath(default_file, winslash = "/", mustWork = TRUE))
  }

  if (exists("build_tf_network", mode = "function")) {
    build_res <- build_tf_network()
    if (is.list(build_res) && !is.null(build_res$files) && !is.null(build_res$files$csv)) {
      csv_path <- as.character(build_res$files$csv[[1]])
      if (nzchar(csv_path) && file.exists(csv_path)) {
        return(normalizePath(csv_path, winslash = "/", mustWork = TRUE))
      }
    }
  }

  stop("No usable TF network file was found. Supply `edge_file` explicitly or set options(TFRegAct.network_path = '<path>').")
}

find_adjustment_dagitty_read_edges <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "csv") {
    df <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    df <- utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  }

  need_cols <- c("tf", "target")
  if (!all(need_cols %in% colnames(df))) {
    stop("Edge file must contain at least columns: tf, target")
  }

  keep_cols <- intersect(c("tf", "target", "effect", "confidence_score", "supporting_databases"), colnames(df))
  df <- df[, keep_cols, drop = FALSE]
  if (!("effect" %in% colnames(df))) {
    df$effect <- "unknown"
  }
  if (!("confidence_score" %in% colnames(df))) {
    df$confidence_score <- NA_integer_
  }
  if (!("supporting_databases" %in% colnames(df))) {
    df$supporting_databases <- ""
  }

  df$tf <- trimws(as.character(df$tf))
  df$target <- trimws(as.character(df$target))
  df$effect <- trimws(as.character(df$effect))
  df$confidence_score <- suppressWarnings(as.integer(df$confidence_score))
  df$supporting_databases <- trimws(as.character(df$supporting_databases))

  if ("supporting_databases" %in% colnames(df)) {
    db_tokens <- strsplit(df$supporting_databases, "\\s*;\\s*")
    grouped_score <- vapply(db_tokens, function(x) {
      vals <- trimws(x)
      vals <- vals[vals != ""]
      has_trrust <- any(vals == "TRRUST")
      has_regnetwork <- any(vals == "RegNetwork")
      has_chip <- any(vals %in% c("ChEA", "ENCODE"))
      has_low <- any(vals %in% c("JASPAR"))
      4L * as.integer(has_trrust) +
        3L * as.integer(has_regnetwork) +
        2L * as.integer(has_chip) +
        1L * as.integer(has_low)
    }, integer(1))
    df$confidence_score <- as.integer(grouped_score)
  }

  df <- df[!is.na(df$tf) & !is.na(df$target) & df$tf != "" & df$target != "", , drop = FALSE]
  df <- df[df$tf != df$target, , drop = FALSE]
  unique(df)
}

find_adjustment_dagitty_sanitize_node <- function(x) {
  y <- gsub("[^A-Za-z0-9_]", "_", x)
  y <- gsub("_+", "_", y)
  y <- sub("^_+", "", y)
  ifelse(grepl("^[A-Za-z]", y), y, paste0("N_", y))
}

find_adjustment_dagitty_make_node_map <- function(nodes) {
  safe <- find_adjustment_dagitty_sanitize_node(nodes)
  out <- character(length(nodes))
  seen <- integer(0)
  names(seen) <- character(0)

  for (i in seq_along(nodes)) {
    base <- safe[[i]]
    if (!base %in% names(seen)) {
      seen[[base]] <- 1L
      out[[i]] <- base
    } else {
      seen[[base]] <- seen[[base]] + 1L
      out[[i]] <- sprintf("%s_%d", base, seen[[base]])
    }
  }

  data.frame(node = nodes, dagitty_id = out, stringsAsFactors = FALSE)
}

find_adjustment_dagitty_direct_distance_map <- function(graph, target_node) {
  suppressPackageStartupMessages(library(igraph))

  if (!target_node %in% igraph::V(graph)$name) {
    stop(sprintf("Target node '%s' not found in graph.", target_node))
  }

  rev_graph <- igraph::reverse_edges(graph)
  bfs_res <- igraph::bfs(
    rev_graph,
    root = target_node,
    mode = "out",
    unreachable = FALSE,
    dist = TRUE,
    order = FALSE
  )

  dist_vals <- as.numeric(bfs_res$dist)
  dist_vals[dist_vals < 0] <- Inf
  names(dist_vals) <- igraph::V(graph)$name
  dist_vals
}

find_adjustment_dagitty_build_local_dag <- function(edges, tf_query, gene_query, graph = NULL, reverse_graph = NULL) {
  suppressPackageStartupMessages(library(igraph))

  g <- if (is.null(graph)) {
    igraph::graph_from_data_frame(edges[, c("tf", "target"), drop = FALSE], directed = TRUE)
  } else {
    graph
  }
  rev_g <- if (is.null(reverse_graph)) {
    igraph::reverse_edges(g)
  } else {
    reverse_graph
  }

  anc_gene <- igraph::subcomponent(rev_g, gene_query, mode = "out")
  anc_tf <- igraph::subcomponent(rev_g, tf_query, mode = "out")

  # For backdoor adjustment we only need the ancestor graph of exposure/outcome.
  relevant_nodes <- unique(c(igraph::V(g)$name[anc_gene], igraph::V(g)$name[anc_tf], tf_query, gene_query))
  sub_g <- igraph::induced_subgraph(g, vids = relevant_nodes)
  sub_edges <- igraph::as_data_frame(sub_g, what = "edges")
  sub_edges <- merge(sub_edges, edges, by.x = c("from", "to"), by.y = c("tf", "target"), all.x = TRUE, sort = FALSE)
  colnames(sub_edges)[colnames(sub_edges) == "from"] <- "tf"
  colnames(sub_edges)[colnames(sub_edges) == "to"] <- "target"

  nodes <- sort(unique(c(sub_edges$tf, sub_edges$target, tf_query, gene_query)))
  dist_map_gene <- find_adjustment_dagitty_direct_distance_map(g, gene_query)
  dist_map_tf <- find_adjustment_dagitty_direct_distance_map(g, tf_query)
  d_to_gene <- as.numeric(dist_map_gene[nodes])
  d_to_tf <- as.numeric(dist_map_tf[nodes])

  node_meta <- data.frame(
    node = nodes,
    d_to_gene = d_to_gene,
    d_to_tf = d_to_tf,
    stringsAsFactors = FALSE
  )

  node_meta$class_rank <- 2L
  node_meta$class_rank[is.finite(node_meta$d_to_tf) & is.finite(node_meta$d_to_gene) & node_meta$node != tf_query & node_meta$node != gene_query] <- 0L
  node_meta$class_rank[node_meta$node == tf_query] <- 1L
  node_meta$class_rank[node_meta$node == gene_query] <- 3L

  node_meta$order_score <- with(
    node_meta,
    class_rank * 1e6 +
      ifelse(is.finite(d_to_tf), d_to_tf, 999) * 1e4 +
      ifelse(is.finite(d_to_gene), d_to_gene, 999)
  )

  order_map <- setNames(rank(node_meta$order_score, ties.method = "first"), node_meta$node)
  sub_edges$from_rank <- order_map[sub_edges$tf]
  sub_edges$to_rank <- order_map[sub_edges$target]
  dag_edges <- sub_edges[sub_edges$from_rank < sub_edges$to_rank, , drop = FALSE]
  dag_edges <- dag_edges[!duplicated(dag_edges[, c("tf", "target")]), , drop = FALSE]

  list(
    dag_edges = dag_edges,
    node_meta = node_meta
  )
}

find_adjustment_dagitty_choose_recommended_set <- function(adj_sets_df, dag_edges_out, tf_query, gene_query) {
  if (nrow(adj_sets_df) == 0) {
    return(list(
      recommended_set_id = NA_integer_,
      recommended_variables = character(0),
      recommended_table = data.frame(
        set_id = integer(0),
        variable = character(0),
        variable_confidence = integer(0),
        stringsAsFactors = FALSE
      ),
      set_scores = data.frame(
        set_id = integer(0),
        set_size = integer(0),
        total_confidence = integer(0),
        stringsAsFactors = FALSE
      )
    ))
  }

  vars <- sort(unique(adj_sets_df$variable))
  variable_conf <- vapply(vars, function(v) {
    vals <- dag_edges_out$confidence_score[dag_edges_out$tf == v | dag_edges_out$target == v]
    vals <- vals[is.finite(vals)]
    if (length(vals) == 0) {
      0L
    } else {
      as.integer(max(vals))
    }
  }, integer(1))
  var_conf_df <- data.frame(
    variable = vars,
    variable_confidence = as.integer(variable_conf),
    stringsAsFactors = FALSE
  )

  scored_df <- merge(adj_sets_df, var_conf_df, by = "variable", all.x = TRUE, sort = FALSE)
  scored_df$variable_confidence[is.na(scored_df$variable_confidence)] <- 0L

  split_sets <- split(scored_df, scored_df$set_id)
  set_scores <- do.call(
    rbind,
    lapply(split_sets, function(chunk) {
      data.frame(
        set_id = chunk$set_id[[1]],
        set_size = nrow(chunk),
        total_confidence = sum(chunk$variable_confidence),
        stringsAsFactors = FALSE
      )
    })
  )
  set_scores <- set_scores[order(set_scores$set_size, -set_scores$total_confidence, set_scores$set_id), , drop = FALSE]
  best_id <- set_scores$set_id[[1]]
  recommended_table <- scored_df[scored_df$set_id == best_id, c("set_id", "variable", "variable_confidence"), drop = FALSE]
  recommended_table <- recommended_table[order(-recommended_table$variable_confidence, recommended_table$variable), , drop = FALSE]
  rownames(recommended_table) <- NULL

  list(
    recommended_set_id = best_id,
    recommended_variables = recommended_table$variable,
    recommended_table = recommended_table,
    set_scores = set_scores
  )
}

find_adjustment_dagitty_prepare_target_tree <- function(graph, target_node) {
  if (!target_node %in% igraph::V(graph)$name) {
    stop(sprintf("Target node '%s' not found in graph.", target_node))
  }

  bfs_res <- igraph::bfs(
    graph,
    root = target_node,
    mode = "all",
    unreachable = TRUE,
    dist = TRUE,
    parent = TRUE,
    order = FALSE
  )

  vertex_names <- igraph::V(graph)$name
  dist_vals <- bfs_res$dist
  names(dist_vals) <- vertex_names

  father_idx <- bfs_res$parent
  father_names <- rep(NA_character_, length(father_idx))
  valid_father <- !is.na(father_idx) & father_idx > 0
  father_names[valid_father] <- vertex_names[father_idx[valid_father]]
  names(father_names) <- vertex_names

  list(
    graph = graph,
    target = target_node,
    distances = dist_vals,
    fathers = father_names
  )
}

find_adjustment_dagitty_best_path_metric <- function(tree_info, from_node) {
  graph <- tree_info$graph
  target_node <- tree_info$target

  if (!from_node %in% names(tree_info$distances)) {
    return(list(distance = Inf, avg_confidence = 0))
  }

  if (identical(from_node, target_node)) {
    return(list(distance = 0, avg_confidence = 0))
  }

  dist_val <- tree_info$distances[[from_node]]
  if (!is.finite(dist_val)) {
    return(list(distance = Inf, avg_confidence = 0))
  }

  total_conf <- 0
  step_count <- 0
  current <- from_node

  while (!identical(current, target_node)) {
    parent <- tree_info$fathers[[current]]
    if (is.na(parent) || !nzchar(parent)) {
      return(list(distance = Inf, avg_confidence = 0))
    }

    edge_id <- igraph::get_edge_ids(graph, vp = c(current, parent), directed = FALSE)
    if (length(edge_id) == 0 || edge_id[[1]] <= 0) {
      return(list(distance = Inf, avg_confidence = 0))
    }

    conf_val <- igraph::E(graph)$confidence_score[[edge_id[[1]]]]
    if (is.finite(conf_val)) {
      total_conf <- total_conf + conf_val
    }

    step_count <- step_count + 1
    current <- parent
  }

  avg_conf <- if (step_count == 0) 0 else total_conf / step_count
  list(distance = step_count, avg_confidence = avg_conf)
}

find_adjustment_dagitty_build_distance_graph <- function(edges) {
  suppressPackageStartupMessages(library(data.table))
  edge_dt <- data.table::as.data.table(edges[, c("tf", "target", "confidence_score"), drop = FALSE])
  edge_dt[, node1 := ifelse(tf <= target, tf, target)]
  edge_dt[, node2 := ifelse(tf <= target, target, tf)]
  collapsed <- edge_dt[
    ,
    .(confidence_score = as.integer(max(confidence_score, na.rm = TRUE))),
    by = .(node1, node2)
  ]
  collapsed[!is.finite(confidence_score), confidence_score := 0L]
  suppressPackageStartupMessages(library(igraph))
  graph <- igraph::graph_from_data_frame(
    collapsed[, .(from = node1, to = node2, confidence_score)],
    directed = FALSE
  )
  igraph::E(graph)$confidence_score <- collapsed$confidence_score
  graph
}

find_adjustment_dagitty_direct_effect_to_gene <- function(edges, gene_query, regulators) {
  regulator_vec <- unique(as.character(regulators))
  out <- rep("unknown", length(regulator_vec))
  names(out) <- regulator_vec

  edge_sub <- edges[edges$target == gene_query & edges$tf %in% regulator_vec, c("tf", "effect"), drop = FALSE]
  if (nrow(edge_sub) == 0) {
    return(out)
  }

  split_effects <- split(trimws(as.character(edge_sub$effect)), edge_sub$tf)
  for (reg in names(split_effects)) {
    vals <- split_effects[[reg]]
    vals <- vals[!is.na(vals) & vals != ""]
    vals <- vals[vals %in% c("activation", "repression", "unknown", "mixed")]
    vals_known <- vals[vals %in% c("activation", "repression")]

    out[[reg]] <- if (length(vals_known) == 0) {
      "unknown"
    } else if (all(vals_known == "activation")) {
      "activation"
    } else if (all(vals_known == "repression")) {
      "repression"
    } else {
      "unknown"
    }
  }

  out
}

find_adjustment_dagitty_effect_to_direction <- function(effect) {
  effect <- tolower(trimws(as.character(effect)))
  out <- integer(length(effect))
  out[effect %in% c("activation", "activate", "positive", "up")] <- 1L
  out[effect %in% c("repression", "repress", "negative", "down")] <- -1L
  out
}

find_adjustment_dagitty_direct_edge_metrics <- function(
  edges,
  gene_query,
  regulators
) {
  regulators <- unique(trimws(as.character(regulators)))
  out <- data.frame(
    tf = regulators,
    distance_to_gene_A = rep(Inf, length(regulators)),
    effect_on_gene_A = rep("unknown", length(regulators)),
    direction_to_gene_A = integer(length(regulators)),
    avg_confidence_to_gene_A = numeric(length(regulators)),
    stringsAsFactors = FALSE
  )
  if (length(regulators) == 0L) {
    return(out)
  }

  direct <- edges[
    edges$target == gene_query & edges$tf %in% regulators,
    c("tf", "effect", "confidence_score"),
    drop = FALSE
  ]
  if (nrow(direct) == 0L) {
    return(out)
  }

  grouped <- split(seq_len(nrow(direct)), direct$tf)
  for (regulator in names(grouped)) {
    block <- direct[grouped[[regulator]], , drop = FALSE]
    directions <- unique(find_adjustment_dagitty_effect_to_direction(block$effect))
    known_directions <- unique(directions[directions != 0L])
    direction <- if (length(known_directions) == 1L) {
      known_directions[[1]]
    } else {
      0L
    }
    effect <- if (direction == 1L) {
      "activation"
    } else if (direction == -1L) {
      "repression"
    } else if (length(known_directions) > 1L) {
      "mixed"
    } else {
      "unknown"
    }
    confidence <- suppressWarnings(as.numeric(block$confidence_score))
    confidence <- confidence[is.finite(confidence)]
    row_index <- match(regulator, out$tf)
    out$distance_to_gene_A[[row_index]] <- 1
    out$effect_on_gene_A[[row_index]] <- effect
    out$direction_to_gene_A[[row_index]] <- as.integer(direction)
    out$avg_confidence_to_gene_A[[row_index]] <- if (length(confidence)) {
      max(confidence)
    } else {
      0
    }
  }
  out
}

find_adjustment_dagitty_build_recommended_metrics <- function(
  edges,
  tf_query,
  gene_query,
  recommended_variables,
  beta = 2,
  direct_edges = edges
) {
  beta <- suppressWarnings(as.numeric(beta[[1]]))
  if (!is.finite(beta) || beta < 0) {
    stop("`beta` must be a non-negative number.")
  }

  graph <- find_adjustment_dagitty_build_distance_graph(edges)
  tf_tree <- find_adjustment_dagitty_prepare_target_tree(graph, tf_query)

  vars <- unique(c(recommended_variables, tf_query, gene_query))
  direct_metrics <- find_adjustment_dagitty_direct_edge_metrics(
    edges = direct_edges,
    gene_query = gene_query,
    regulators = vars
  )
  out <- lapply(vars, function(v) {
    direct_index <- match(v, direct_metrics$tf)
    direct_distance <- direct_metrics$distance_to_gene_A[[direct_index]]
    direct_confidence <- direct_metrics$avg_confidence_to_gene_A[[direct_index]]
    to_tf <- find_adjustment_dagitty_best_path_metric(tf_tree, v)
    avg_to_gene <- if (is.finite(direct_distance)) direct_confidence else 0
    avg_to_tf <- if (is.finite(to_tf$distance)) to_tf$avg_confidence else 0
    overall_avg <- if (is.finite(direct_distance) && is.finite(to_tf$distance)) {
      (beta * avg_to_gene + avg_to_tf) / (1 + beta)
    } else if (is.finite(direct_distance)) {
      avg_to_gene
    } else if (is.finite(to_tf$distance)) {
      avg_to_tf
    } else {
      0
    }
    data.frame(
      tf = v,
      distance_to_gene_A = direct_distance,
      effect_on_gene_A = direct_metrics$effect_on_gene_A[[direct_index]],
      direction_to_gene_A = direct_metrics$direction_to_gene_A[[direct_index]],
      avg_confidence_to_gene_A = avg_to_gene,
      distance_to_tf_B = if (is.finite(to_tf$distance)) to_tf$distance else Inf,
      avg_confidence_to_tf_B = avg_to_tf,
      overall_avg_confidence = overall_avg,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  out <- out[order(out$tf), , drop = FALSE]
  rownames(out) <- NULL
  out
}

find_adjustment_dagitty_recommended_metrics <- function(
  tf,
  gene,
  network = NULL,
  edge_file = "tf_union_output/tf_gene_merged_human_weighted_clean.csv",
  recommended_set_file = NULL,
  outdir = "adjustment_output",
  beta = 2,
  overall_confidence_threshold = 3,
  build_if_missing = TRUE,
  max_adjustment_sets = Inf,
  write_full_outputs = FALSE,
  write_files = TRUE
) {
  tf_query <- trimws(as.character(tf))
  gene_query <- trimws(as.character(gene))
  if (!nzchar(tf_query) || !nzchar(gene_query)) {
    stop("`tf` and `gene` must be non-empty strings.")
  }
  overall_confidence_threshold <- suppressWarnings(as.numeric(overall_confidence_threshold[[1]]))
  if (!is.finite(overall_confidence_threshold)) {
    stop("`overall_confidence_threshold` must be a finite number.")
  }

  local_dag_edge_file <- file.path(outdir, sprintf("dagitty_local_dag_edges_%s_%s.csv", tf_query, gene_query))
  bundle <- find_adjustment_dagitty_extract_bundle(network = network, edge_file = edge_file)
  use_default_edge <- is.null(edge_file) || identical(edge_file, "tf_union_output/tf_gene_merged_human_weighted_clean.csv")
  resolved_edge_file <- if (!is.null(bundle)) {
    if (!is.null(bundle$metadata$source_csv) && nzchar(bundle$metadata$source_csv)) {
      bundle$metadata$source_csv
    } else {
      edge_file
    }
  } else if (use_default_edge && file.exists(local_dag_edge_file)) {
    normalizePath(local_dag_edge_file, winslash = "/", mustWork = TRUE)
  } else {
    find_adjustment_dagitty_extract_edge_file(network = network, edge_file = edge_file)
  }
  edges <- if (!is.null(bundle)) bundle$edges else find_adjustment_dagitty_read_edges(resolved_edge_file)

  if (is.null(recommended_set_file) || !nzchar(recommended_set_file)) {
    recommended_set_file <- file.path(outdir, sprintf("dagitty_recommended_adjustment_set_%s_%s.csv", tf_query, gene_query))
  }
  if (!file.exists(recommended_set_file)) {
    if (!isTRUE(build_if_missing)) {
      stop(sprintf("Recommended adjustment set file not found: %s", recommended_set_file))
    }

    generated <- find_adjustment_dagitty_run(
      tf = tf_query,
      gene = gene_query,
      network = network,
      edge_file = edge_file,
      outdir = outdir,
      beta = beta,
      overall_confidence_threshold = overall_confidence_threshold,
      max_adjustment_sets = max_adjustment_sets,
      write_full_outputs = write_full_outputs,
      write_files = write_files
    )

    return(list(
      tf = tf_query,
      gene = gene_query,
      metrics = generated$recommended_adjustment_metrics,
      filtered_metrics = generated$filtered_recommended_adjustment_metrics,
      overall_confidence_threshold = overall_confidence_threshold,
      files = generated$files,
      generated_recommended_set = TRUE
    ))
  }

  rec_df <- utils::read.csv(recommended_set_file, stringsAsFactors = FALSE, check.names = FALSE)
  if (!("variable" %in% colnames(rec_df))) {
    stop("Recommended adjustment set file must contain a 'variable' column.")
  }

  recommended_variables <- unique(trimws(as.character(rec_df$variable)))
  recommended_variables <- recommended_variables[!is.na(recommended_variables) & recommended_variables != ""]
  metrics <- find_adjustment_dagitty_build_recommended_metrics(
    edges = edges,
    tf_query = tf_query,
    gene_query = gene_query,
    recommended_variables = recommended_variables,
    beta = beta,
    direct_edges = edges
  )
  keep_forced <- metrics$tf %in% c(tf_query, gene_query)
  keep_direct <- metrics$distance_to_gene_A == 1
  keep_conf <- keep_direct &
    metrics$overall_avg_confidence >= overall_confidence_threshold
  filtered_metrics <- metrics[keep_forced | keep_conf, , drop = FALSE]
  rownames(filtered_metrics) <- NULL

  outdir <- normalizePath(outdir, winslash = "/", mustWork = FALSE)
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  }
  metrics_file <- file.path(outdir, sprintf("dagitty_recommended_adjustment_metrics_%s_%s.csv", tf_query, gene_query))
  filtered_metrics_file <- file.path(
    outdir,
    sprintf("dagitty_recommended_adjustment_metrics_filtered_%s_%s.csv", tf_query, gene_query)
  )
  if (write_files) {
    utils::write.csv(filtered_metrics, file = filtered_metrics_file, row.names = FALSE)
    if (isTRUE(write_full_outputs)) {
      utils::write.csv(metrics, file = metrics_file, row.names = FALSE)
    }
  }

  list(
    tf = tf_query,
    gene = gene_query,
    metrics = metrics,
    filtered_metrics = filtered_metrics,
    overall_confidence_threshold = overall_confidence_threshold,
    files = list(
      recommended_set = normalizePath(recommended_set_file, winslash = "/", mustWork = TRUE),
      recommended_adjustment_metrics = normalizePath(metrics_file, winslash = "/", mustWork = FALSE),
      recommended_adjustment_metrics_filtered = normalizePath(filtered_metrics_file, winslash = "/", mustWork = FALSE),
      edge_file_used = resolved_edge_file
    )
  )
}

find_adjustment_dagitty_run <- function(
  tf,
  gene,
  network = NULL,
  edge_file = "tf_union_output/tf_gene_merged_human_weighted_clean.csv",
  outdir = "adjustment_output",
  beta = 2,
  overall_confidence_threshold = 3,
  max_adjustment_sets = Inf,
  include_candidate_confounders = TRUE,
  write_full_outputs = FALSE,
  write_files = TRUE
) {
  find_adjustment_dagitty_install_if_missing("igraph")
  find_adjustment_dagitty_install_if_missing("dagitty")
  suppressPackageStartupMessages(library(igraph))
  suppressPackageStartupMessages(library(dagitty))

  tf_query <- trimws(as.character(tf))
  gene_query <- trimws(as.character(gene))
  if (!nzchar(tf_query) || !nzchar(gene_query)) {
    stop("`tf` and `gene` must be non-empty strings.")
  }
  if (length(max_adjustment_sets) == 0 || is.null(max_adjustment_sets) || !is.finite(max_adjustment_sets[[1]])) {
    max_adjustment_sets <- Inf
  } else {
    max_adjustment_sets <- max(1L, as.integer(max_adjustment_sets[[1]]))
  }

  bundle <- find_adjustment_dagitty_extract_bundle(network = network, edge_file = edge_file)
  edge_file <- if (is.null(bundle)) {
    find_adjustment_dagitty_extract_edge_file(network = network, edge_file = edge_file)
  } else if (!is.null(bundle$metadata$source_csv) && nzchar(bundle$metadata$source_csv)) {
    bundle$metadata$source_csv
  } else {
    ""
  }
  outdir <- normalizePath(outdir, winslash = "/", mustWork = FALSE)
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  }

  edges <- if (!is.null(bundle)) bundle$edges else find_adjustment_dagitty_read_edges(edge_file)
  all_nodes <- sort(unique(c(edges$tf, edges$target)))

  if (!tf_query %in% all_nodes) {
    stop(sprintf("TF '%s' not found in edge file.", tf_query))
  }
  if (!gene_query %in% all_nodes) {
    stop(sprintf("Gene '%s' not found in edge file.", gene_query))
  }

  local_res <- find_adjustment_dagitty_build_local_dag(
    edges,
    tf_query,
    gene_query,
    graph = if (!is.null(bundle)) bundle$graph else NULL,
    reverse_graph = if (!is.null(bundle)) bundle$reverse_graph else NULL
  )
  dag_edges <- local_res$dag_edges
  node_meta <- local_res$node_meta

  local_nodes <- sort(unique(c(dag_edges$tf, dag_edges$target, tf_query, gene_query)))
  node_map <- find_adjustment_dagitty_make_node_map(local_nodes)
  id_map <- setNames(node_map$dagitty_id, node_map$node)
  rev_map <- setNames(node_map$node, node_map$dagitty_id)

  dagitty_lines <- c("dag {")
  for (id in node_map$dagitty_id) {
    tags <- character(0)
    if (id == id_map[[tf_query]]) {
      tags <- c(tags, "exposure")
    }
    if (id == id_map[[gene_query]]) {
      tags <- c(tags, "outcome")
    }
    if (length(tags) > 0) {
      dagitty_lines <- c(dagitty_lines, sprintf("  %s [%s]", id, paste(tags, collapse = ",")))
    }
  }
  if (nrow(dag_edges) > 0) {
    edge_strings <- sprintf("  %s -> %s", id_map[dag_edges$tf], id_map[dag_edges$target])
    dagitty_lines <- c(dagitty_lines, edge_strings)
  }
  dagitty_lines <- c(dagitty_lines, "}")
  dagitty_string <- paste(dagitty_lines, collapse = "\n")
  g_dag <- dagitty::dagitty(dagitty_string)

  adj_sets_ids <- dagitty::adjustmentSets(
    g_dag,
    exposure = id_map[[tf_query]],
    outcome = id_map[[gene_query]],
    effect = "total",
    max.results = max_adjustment_sets
  )
  if (length(adj_sets_ids) == 0) {
    adj_sets_df <- data.frame(set_id = integer(0), variable = character(0), stringsAsFactors = FALSE)
  } else {
    adj_sets_df <- do.call(
      rbind,
      lapply(seq_along(adj_sets_ids), function(i) {
        vals <- as.character(adj_sets_ids[[i]])
        data.frame(
          set_id = rep(i, length(vals)),
          variable = unname(rev_map[vals]),
          stringsAsFactors = FALSE
        )
      })
    )
  }

  candidate_df <- if (isTRUE(include_candidate_confounders)) {
    candidate_confounders_ids <- dagitty::adjustmentSets(
      g_dag,
      exposure = id_map[[tf_query]],
      outcome = id_map[[gene_query]],
      type = "canonical"
    )
    if (length(candidate_confounders_ids) == 0) {
      data.frame(variable = character(0), stringsAsFactors = FALSE)
    } else {
      data.frame(
        variable = unname(rev_map[as.character(candidate_confounders_ids[[1]])]),
        stringsAsFactors = FALSE
      )
    }
  } else {
    data.frame(variable = character(0), stringsAsFactors = FALSE)
  }

  dag_edges_out <- dag_edges[, c("tf", "target", "effect", "confidence_score", "supporting_databases"), drop = FALSE]
  dag_edges_out$direction <- find_adjustment_dagitty_effect_to_direction(
    dag_edges_out$effect
  )
  node_meta_out <- merge(node_meta, node_map, by = "node", all.x = TRUE, sort = FALSE)
  recommended_res <- find_adjustment_dagitty_choose_recommended_set(adj_sets_df, dag_edges_out, tf_query, gene_query)
  recommended_metrics <- find_adjustment_dagitty_build_recommended_metrics(
    edges = dag_edges_out,
    tf_query = tf_query,
    gene_query = gene_query,
    recommended_variables = recommended_res$recommended_variables,
    beta = beta,
    direct_edges = edges
  )
  keep_forced <- recommended_metrics$tf %in% c(tf_query, gene_query)
  keep_direct <- recommended_metrics$distance_to_gene_A == 1
  keep_conf <- keep_direct &
    recommended_metrics$overall_avg_confidence >= overall_confidence_threshold
  filtered_recommended_metrics <- recommended_metrics[keep_forced | keep_conf, , drop = FALSE]
  rownames(filtered_recommended_metrics) <- NULL
  node_meta_out$role <- "other"
  node_meta_out$role[node_meta_out$node == tf_query] <- "exposure"
  node_meta_out$role[node_meta_out$node == gene_query] <- "outcome"
  node_meta_out$role[node_meta_out$node %in% candidate_df$variable] <- "candidate_confounder"
  node_meta_out$role[node_meta_out$node %in% recommended_res$recommended_variables] <- "recommended_adjustment"
  node_meta_out$role[node_meta_out$node == tf_query] <- "exposure"
  node_meta_out$role[node_meta_out$node == gene_query] <- "outcome"

  file_map <- list(
    dag_edges = file.path(outdir, sprintf("dagitty_local_dag_edges_%s_%s.csv", tf_query, gene_query)),
    dag_nodes = file.path(outdir, sprintf("dagitty_local_dag_nodes_%s_%s.csv", tf_query, gene_query)),
    candidate_confounders = file.path(outdir, sprintf("dagitty_candidate_confounders_%s_%s.csv", tf_query, gene_query)),
    minimal_adjustment_sets = file.path(outdir, sprintf("dagitty_minimal_adjustment_sets_%s_%s.csv", tf_query, gene_query)),
    recommended_adjustment_set = file.path(outdir, sprintf("dagitty_recommended_adjustment_set_%s_%s.csv", tf_query, gene_query)),
    recommended_adjustment_metrics = file.path(outdir, sprintf("dagitty_recommended_adjustment_metrics_%s_%s.csv", tf_query, gene_query)),
    recommended_adjustment_metrics_filtered = file.path(
      outdir,
      sprintf("dagitty_recommended_adjustment_metrics_filtered_%s_%s.csv", tf_query, gene_query)
    ),
    adjustment_set_scores = file.path(outdir, sprintf("dagitty_adjustment_set_scores_%s_%s.csv", tf_query, gene_query)),
    dagitty_graph = file.path(outdir, sprintf("dagitty_graph_%s_%s.txt", tf_query, gene_query)),
    summary = file.path(outdir, sprintf("dagitty_summary_%s_%s.txt", tf_query, gene_query))
  )

  if (write_files) {
    utils::write.csv(filtered_recommended_metrics, file = file_map$recommended_adjustment_metrics_filtered, row.names = FALSE)
    if (isTRUE(write_full_outputs)) {
      utils::write.csv(dag_edges_out, file = file_map$dag_edges, row.names = FALSE)
      utils::write.csv(node_meta_out, file = file_map$dag_nodes, row.names = FALSE)
      utils::write.csv(candidate_df, file = file_map$candidate_confounders, row.names = FALSE)
      utils::write.csv(adj_sets_df, file = file_map$minimal_adjustment_sets, row.names = FALSE)
      utils::write.csv(recommended_res$recommended_table, file = file_map$recommended_adjustment_set, row.names = FALSE)
      utils::write.csv(recommended_metrics, file = file_map$recommended_adjustment_metrics, row.names = FALSE)
      utils::write.csv(recommended_res$set_scores, file = file_map$adjustment_set_scores, row.names = FALSE)
      writeLines(dagitty_string, con = file_map$dagitty_graph)
    }
  }

  summary_lines <- c(
    sprintf("TF: %s", tf_query),
    sprintf("Gene: %s", gene_query),
    sprintf("Edge file: %s", edge_file),
    sprintf("Local DAG nodes: %d", length(local_nodes)),
    sprintf("Local DAG edges: %d", nrow(dag_edges_out)),
    sprintf("Candidate confounders: %d", nrow(candidate_df)),
    sprintf("Minimal adjustment sets: %d", if (nrow(adj_sets_df) == 0) 0 else length(unique(adj_sets_df$set_id))),
    sprintf("Max adjustment sets: %s", if (is.infinite(max_adjustment_sets)) "Inf" else as.character(max_adjustment_sets)),
    sprintf("Beta: %s", format(beta, trim = TRUE)),
    sprintf("Overall confidence threshold: %s", format(overall_confidence_threshold, trim = TRUE)),
    sprintf(
      "Recommended adjustment set: %s",
      if (length(recommended_res$recommended_variables) == 0) "none" else paste(recommended_res$recommended_variables, collapse = ";")
    ),
    sprintf("DAG edge file: %s", file_map$dag_edges),
    sprintf("Adjustment set file: %s", file_map$minimal_adjustment_sets),
    sprintf("Recommended set file: %s", file_map$recommended_adjustment_set),
    sprintf("Recommended metrics file: %s", file_map$recommended_adjustment_metrics),
    sprintf("Filtered recommended metrics file: %s", file_map$recommended_adjustment_metrics_filtered)
  )

  if (write_files && isTRUE(write_full_outputs)) {
    writeLines(summary_lines, con = file_map$summary)
  }

  result <- list(
    tf = tf_query,
    gene = gene_query,
    network = network,
    edge_file = edge_file,
    dag_edges = dag_edges_out,
    dag_nodes = node_meta_out,
    candidate_confounders = candidate_df,
    minimal_adjustment_sets = adj_sets_df,
    recommended_adjustment_set = recommended_res$recommended_table,
    recommended_adjustment_metrics = recommended_metrics,
    filtered_recommended_adjustment_metrics = filtered_recommended_metrics,
    adjustment_set_scores = recommended_res$set_scores,
    dagitty_graph = dagitty_string,
    files = file_map,
    summary = summary_lines
  )
  class(result) <- c("find_adjustment_dagitty_result", class(result))
  result
}

find_adjustment_dagitty <- function(
  tf,
  gene,
  network = NULL,
  edge_file = "tf_union_output/tf_gene_merged_human_weighted_clean.csv",
  outdir = "adjustment_output",
  beta = 2,
  overall_confidence_threshold = 3,
  max_adjustment_sets = Inf,
  include_candidate_confounders = TRUE,
  write_full_outputs = FALSE,
  write_files = TRUE
) {
  find_adjustment_dagitty_run(
    tf = tf,
    gene = gene,
    network = network,
    edge_file = edge_file,
    outdir = outdir,
    beta = beta,
    overall_confidence_threshold = overall_confidence_threshold,
    max_adjustment_sets = max_adjustment_sets,
    include_candidate_confounders = include_candidate_confounders,
    write_full_outputs = write_full_outputs,
    write_files = write_files
  )
}

print.find_adjustment_dagitty_result <- function(x, ...) {
  cat(paste(x$summary, collapse = "\n"), "\n")
  invisible(x)
}

find_adjustment_dagitty_main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  required_args <- c("--tf", "--gene")
  missing_required <- required_args[!vapply(required_args, function(x) any(grepl(paste0("^", x, "="), args)), logical(1))]
  if (length(missing_required) > 0) {
    stop(sprintf("Missing required argument(s): %s", paste(missing_required, collapse = ", ")))
  }

  result <- find_adjustment_dagitty(
    tf = find_adjustment_dagitty_get_arg(args, "--tf"),
    gene = find_adjustment_dagitty_get_arg(args, "--gene"),
    edge_file = find_adjustment_dagitty_get_arg(args, "--edge_file", "tf_union_output/tf_gene_merged_human_weighted_clean.csv"),
    outdir = find_adjustment_dagitty_get_arg(args, "--outdir", "adjustment_output"),
    beta = as.numeric(find_adjustment_dagitty_get_arg(args, "--beta", "2")),
    overall_confidence_threshold = as.numeric(find_adjustment_dagitty_get_arg(args, "--overall_confidence_threshold", "3")),
    max_adjustment_sets = suppressWarnings(as.numeric(find_adjustment_dagitty_get_arg(args, "--max_adjustment_sets", "Inf"))),
    write_full_outputs = tolower(find_adjustment_dagitty_get_arg(args, "--write_full_outputs", "false")) == "true",
    write_files = TRUE
  )

  print(result)
}

if (sys.nframe() == 0) {
  find_adjustment_dagitty_main()
}
