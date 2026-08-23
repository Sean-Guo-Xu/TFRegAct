#!/usr/bin/env Rscript

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Package 'data.table' is required for build_tf_network.R.")
}

write_csv <- function(df, path) {
  data.table::fwrite(df, file = path, sep = ",", quote = TRUE, na = "")
}

write_tsv <- function(df, path) {
  data.table::fwrite(df, file = path, sep = "\t", quote = FALSE, na = "")
}

clean_name <- function(x) {
  tolower(gsub("[^a-z0-9]+", "", x))
}

match_col <- function(df, candidates) {
  cn <- colnames(df)
  cn_clean <- clean_name(cn)
  cand_clean <- clean_name(candidates)
  idx <- match(cand_clean, cn_clean)
  idx <- idx[!is.na(idx)]
  if (length(idx) == 0) {
    return(NULL)
  }
  cn[[idx[[1]]]]
}

read_table_auto <- function(path, header = TRUE) {
  ext <- tolower(tools::file_ext(path))
  con <- switch(
    ext,
    gz = gzfile(path, open = "rt"),
    path
  )

  on.exit({
    if (inherits(con, "connection")) {
      close(con)
    }
  }, add = TRUE)

  utils::read.delim(
    con,
    header = header,
    sep = "\t",
    stringsAsFactors = FALSE,
    quote = "",
    check.names = FALSE
  )
}

download_file_then_replace <- function(url, destfile, quiet = FALSE) {
  dir.create(dirname(destfile), recursive = TRUE, showWarnings = FALSE)
  temp_dest <- tempfile(
    pattern = paste0(basename(destfile), ".download-"),
    tmpdir = dirname(destfile)
  )
  on.exit(unlink(temp_dest, force = TRUE), add = TRUE)

  utils::download.file(
    url = url,
    destfile = temp_dest,
    mode = "wb",
    quiet = quiet,
    method = "libcurl"
  )

  if (!file.exists(temp_dest) || file.info(temp_dest)$size <= 0) {
    stop(sprintf("Downloaded file is empty: %s", url))
  }
  if (!file.copy(temp_dest, destfile, overwrite = TRUE)) {
    stop(sprintf("Could not replace cached file: %s", destfile))
  }

  invisible(destfile)
}

download_with_cache <- function(urls, destfiles, use_cache = TRUE, quiet = FALSE) {
  if (length(urls) == 1 && length(destfiles) == 1) {
    urls <- as.list(urls)
    destfiles <- as.list(destfiles)
  }

  last_error <- NULL

  for (i in seq_along(urls)) {
    url <- urls[[i]]
    destfile <- destfiles[[i]]

    if (use_cache && file.exists(destfile) && file.info(destfile)$size > 0) {
      if (!quiet) {
        message(sprintf("Using cached file: %s", destfile))
      }
      return(list(path = destfile, url = url, from_cache = TRUE))
    }

    dir.create(dirname(destfile), recursive = TRUE, showWarnings = FALSE)

    try_ok <- tryCatch(
      {
        download_file_then_replace(url, destfile, quiet = quiet)
        TRUE
      },
      error = function(e) {
        last_error <<- conditionMessage(e)
        FALSE
      },
      warning = function(w) {
        last_error <<- conditionMessage(w)
        FALSE
      }
    )

    if (isTRUE(try_ok) && file.exists(destfile) && file.info(destfile)$size > 0) {
      return(list(path = destfile, url = url, from_cache = FALSE))
    }
  }

  stop(sprintf("Download failed for all candidate URLs. Last error: %s", last_error))
}

effect_from_text <- function(x) {
  x <- tolower(trimws(as.character(x)))
  out <- rep("unknown", length(x))
  out[grepl("activ|stimulat|up", x)] <- "activation"
  out[grepl("repress|inhib|down|suppress", x)] <- "repression"
  out
}

collapse_unique_sorted <- function(x) {
  x <- unique(as.character(x))
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) {
    return("")
  }
  paste(sort(x), collapse = ";")
}

resolve_effect <- function(x) {
  x <- unique(trimws(tolower(as.character(x))))
  x <- x[x != "" & !is.na(x)]

  if (length(x) == 0) {
    return("unknown")
  }
  if ("mixed" %in% x) {
    return("mixed")
  }

  known <- setdiff(x, "unknown")
  if (length(known) == 0) {
    return("unknown")
  }
  if (all(known == "activation")) {
    return("activation")
  }
  if (all(known == "repression")) {
    return("repression")
  }
  "mixed"
}

normalize_edges <- function(df, source_db, source_class, species, tf_col, target_col,
                            effect_col = NULL, evidence_col = NULL, score_col = NULL,
                            tf_id_col = NULL, target_id_col = NULL) {
  out <- data.frame(
    tf = trimws(as.character(df[[tf_col]])),
    target = trimws(as.character(df[[target_col]])),
    effect = "unknown",
    source_db = source_db,
    source_class = source_class,
    species = species,
    evidence = NA_character_,
    score = NA_real_,
    tf_id = NA_character_,
    target_id = NA_character_,
    stringsAsFactors = FALSE
  )

  if (!is.null(effect_col) && effect_col %in% colnames(df)) {
    out$effect <- effect_from_text(df[[effect_col]])
  }

  if (!is.null(evidence_col) && evidence_col %in% colnames(df)) {
    out$evidence <- as.character(df[[evidence_col]])
  }

  if (!is.null(score_col) && score_col %in% colnames(df)) {
    out$score <- suppressWarnings(as.numeric(df[[score_col]]))
  }

  if (!is.null(tf_id_col) && tf_id_col %in% colnames(df)) {
    out$tf_id <- as.character(df[[tf_id_col]])
  }

  if (!is.null(target_id_col) && target_id_col %in% colnames(df)) {
    out$target_id <- as.character(df[[target_id_col]])
  }

  out <- out[!is.na(out$tf) & !is.na(out$target), , drop = FALSE]
  out <- out[out$tf != "" & out$target != "", , drop = FALSE]

  tf_clean <- tolower(out$tf)
  target_clean <- tolower(out$target)
  invalid_tokens <- c("na", "n/a", "null", "none", "genesym", "geneid", "source", "target", "attribute")
  out <- out[!(tf_clean %in% invalid_tokens | target_clean %in% invalid_tokens), , drop = FALSE]

  unique(out)
}

parse_trrust <- function(species, cache_dir, use_cache = TRUE, quiet = FALSE) {
  urls <- switch(
    species,
    human = c(
      "https://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv",
      "https://raw.githubusercontent.com/bioinfonerd/Transcription-Factor-Databases/master/Ttrust_v2/trrust_rawdata.human.tsv"
    ),
    mouse = c(
      "https://www.grnpedia.org/trrust/data/trrust_rawdata.mouse.tsv",
      "https://raw.githubusercontent.com/bioinfonerd/Transcription-Factor-Databases/master/Ttrust_v2/trrust_rawdata.mouse.tsv.gz"
    )
  )

  destfiles <- file.path(cache_dir, c(
    sprintf("trrust_%s_primary.tsv", species),
    sprintf("trrust_%s_fallback%s", species, if (species == "mouse") ".tsv.gz" else ".tsv")
  ))

  dl <- download_with_cache(urls, destfiles, use_cache = use_cache, quiet = quiet)
  df <- if (grepl("\\.gz$", dl$path, ignore.case = TRUE)) {
    read_table_auto(dl$path, header = FALSE)
  } else {
    utils::read.delim(
      dl$path,
      header = FALSE,
      sep = "\t",
      stringsAsFactors = FALSE,
      quote = ""
    )
  }

  if (ncol(df) < 3) {
    stop("TRRUST format is unexpected.")
  }

  colnames(df)[1:3] <- c("tf", "target", "effect_raw")
  if (ncol(df) >= 4) {
    colnames(df)[4] <- "reference"
  }

  data <- normalize_edges(
    df = df,
    source_db = "TRRUST",
    source_class = "curated_regulatory_db",
    species = species,
    tf_col = "tf",
    target_col = "target",
    effect_col = "effect_raw",
    evidence_col = "reference"
  )

  list(data = data, raw_path = dl$path, url = dl$url, from_cache = dl$from_cache)
}

find_7zip_executable <- function() {
  candidates <- unname(Sys.which(c("7z", "7zz")))
  if (.Platform$OS.type == "windows") {
    candidates <- c(
      candidates,
      file.path(Sys.getenv("ProgramFiles"), "7-Zip", "7z.exe"),
      file.path(Sys.getenv("ProgramFiles(x86)"), "7-Zip", "7z.exe")
    )
  }
  candidates <- unique(candidates[nzchar(candidates) & file.exists(candidates)])
  if (!length(candidates)) {
    stop(
      "RegNetwork 2025 is distributed as a .7z archive. Install 7-Zip (7z or 7zz) before downloading it."
    )
  }
  normalizePath(candidates[[1]], winslash = "/", mustWork = TRUE)
}

extract_7z_file_with_cache <- function(archive, expected_name, extract_dir,
                                       use_cache = TRUE) {
  extracted_file <- file.path(extract_dir, expected_name)
  if (use_cache && file.exists(extracted_file) && file.info(extracted_file)$size > 0) {
    return(normalizePath(extracted_file, winslash = "/", mustWork = TRUE))
  }

  dir.create(dirname(extract_dir), recursive = TRUE, showWarnings = FALSE)
  stage_dir <- tempfile("regnetwork_7z_extract_", tmpdir = dirname(extract_dir))
  dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(stage_dir, recursive = TRUE, force = TRUE), add = TRUE)

  seven_zip <- find_7zip_executable()
  archive <- normalizePath(archive, winslash = "/", mustWork = TRUE)
  stage_dir <- normalizePath(stage_dir, winslash = "/", mustWork = TRUE)
  output <- system2(
    seven_zip,
    args = c(
      "x",
      "-y",
      shQuote(archive),
      paste0("-o", shQuote(stage_dir)),
      shQuote(expected_name)
    ),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) {
    stop(sprintf("Could not extract RegNetwork 2025 archive: %s", paste(output, collapse = "\n")))
  }

  staged_file <- file.path(stage_dir, expected_name)
  if (!file.exists(staged_file) || file.info(staged_file)$size <= 0) {
    stop(sprintf("RegNetwork 2025 archive does not contain `%s`.", expected_name))
  }

  dir.create(extract_dir, recursive = TRUE, showWarnings = FALSE)
  if (!file.copy(staged_file, extracted_file, overwrite = TRUE)) {
    stop(sprintf("Could not replace extracted RegNetwork file: %s", extracted_file))
  }
  normalizePath(extracted_file, winslash = "/", mustWork = TRUE)
}

parse_regnetwork <- function(species, cache_dir, use_cache = TRUE, quiet = FALSE) {
  archive_name <- sprintf("%s_core_TF_Target.7z", species)
  source_url <- sprintf(
    "https://www.zpliulab.cn/RegNetwork/static/data/%s",
    archive_name
  )
  source_archive <- file.path(cache_dir, paste0("regnetwork_2025_", archive_name))
  source_dir <- file.path(cache_dir, sprintf("regnetwork_2025_%s_core_extract", species))
  source_name <- sprintf("%s_core_TF_Target.txt", species)

  source_dl <- download_with_cache(
    source_url,
    source_archive,
    use_cache = use_cache,
    quiet = quiet
  )
  source_file <- extract_7z_file_with_cache(
    archive = source_dl$path,
    expected_name = source_name,
    extract_dir = source_dir,
    use_cache = use_cache
  )

  source_df <- data.table::fread(
    source_file,
    sep = "\t",
    header = FALSE,
    data.table = FALSE,
    quote = "",
    showProgress = !quiet
  )
  if (ncol(source_df) != 6L) {
    stop(sprintf(
      "RegNetwork 2025 TF-Target format is unexpected: expected 6 columns, found %d.",
      ncol(source_df)
    ))
  }
  colnames(source_df) <- c(
    "tf", "tf_id", "target", "target_id", "regulator_type", "target_type"
  )

  source_df <- source_df[
    toupper(trimws(source_df$regulator_type)) == "TF" &
      toupper(trimws(source_df$target_type)) %in% c("TF", "GENE"),
    ,
    drop = FALSE
  ]
  source_df$tf_id <- sub("^NCBI:", "", source_df$tf_id, ignore.case = TRUE)
  source_df$target_id <- sub("^NCBI:", "", source_df$target_id, ignore.case = TRUE)
  source_df$evidence_raw <- "RegNetwork 2025 core TF-Target"

  data <- normalize_edges(
    df = source_df,
    source_db = "RegNetwork",
    source_class = "curated_regulatory_db",
    species = species,
    tf_col = "tf",
    target_col = "target",
    evidence_col = "evidence_raw",
    tf_id_col = "tf_id",
    target_id_col = "target_id"
  )

  list(
    data = data,
    raw_path = source_dl$path,
    url = source_dl$url,
    from_cache = source_dl$from_cache,
    dataset = "RegNetwork 2025 core TF-Target"
  )
}

parse_harmonizome <- function(species, source_db, source_class, url, cache_name,
                              cache_dir, use_cache = TRUE, quiet = FALSE) {
  if (species != "human") {
    stop(sprintf("%s currently only has the configured human download source.", source_db))
  }

  raw_path <- file.path(cache_dir, cache_name)
  dl <- download_with_cache(url, raw_path, use_cache = use_cache, quiet = quiet)
  df <- read_table_auto(dl$path, header = TRUE)

  source_col <- match_col(df, c("source", "gene", "genesymbol", "gene_symbol"))
  attribute_col <- match_col(df, c("target", "attribute", "tf", "regulator"))
  score_col <- match_col(df, c("weight", "score"))

  if (is.null(source_col) || is.null(attribute_col)) {
    stop(sprintf("%s format is unexpected.", source_db))
  }

  data <- normalize_edges(
    df = df,
    source_db = source_db,
    source_class = source_class,
    species = species,
    tf_col = attribute_col,
    target_col = source_col,
    score_col = score_col
  )

  list(data = data, raw_path = dl$path, url = dl$url, from_cache = dl$from_cache)
}

download_tf_sources <- function(species = "human",
                                databases = c("TRRUST", "RegNetwork", "ChEA", "ENCODE", "JASPAR"),
                                outdir = "tf_source_output",
                                use_cache = TRUE,
                                write_files = TRUE,
                                quiet = FALSE) {
  species <- tolower(species)
  if (!species %in% c("human", "mouse")) {
    stop("species must be 'human' or 'mouse'.")
  }

  db_order <- c("TRRUST", "RegNetwork", "ChEA", "ENCODE", "JASPAR")
  databases <- unique(databases)
  bad_dbs <- setdiff(databases, db_order)
  if (length(bad_dbs) > 0) {
    stop(sprintf("Unsupported databases: %s", paste(bad_dbs, collapse = ", ")))
  }
  databases <- db_order[db_order %in% databases]

  outdir <- normalizePath(outdir, winslash = "/", mustWork = FALSE)
  raw_dir <- file.path(outdir, "raw_cache")
  std_dir <- file.path(outdir, "standardized")
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(std_dir, recursive = TRUE, showWarnings = FALSE)
  options(timeout = max(300, getOption("timeout")))

  parsers <- list(
    TRRUST = function() parse_trrust(species, raw_dir, use_cache = use_cache, quiet = quiet),
    RegNetwork = function() parse_regnetwork(species, raw_dir, use_cache = use_cache, quiet = quiet),
    ChEA = function() parse_harmonizome(
      species = species,
      source_db = "ChEA",
      source_class = "chip_based_binding",
      url = "https://maayanlab.cloud/static/hdfs/harmonizome/data/cheappi/gene_attribute_edges.txt.gz",
      cache_name = sprintf("chea_%s_gene_attribute_edges.txt.gz", species),
      cache_dir = raw_dir,
      use_cache = use_cache,
      quiet = quiet
    ),
    ENCODE = function() parse_harmonizome(
      species = species,
      source_db = "ENCODE",
      source_class = "chip_based_binding",
      url = "https://maayanlab.cloud/static/hdfs/harmonizome/data/encodetfppi/gene_attribute_edges.txt.gz",
      cache_name = sprintf("encode_%s_gene_attribute_edges.txt.gz", species),
      cache_dir = raw_dir,
      use_cache = use_cache,
      quiet = quiet
    ),
    JASPAR = function() parse_harmonizome(
      species = species,
      source_db = "JASPAR",
      source_class = "motif_prediction",
      url = "https://maayanlab.cloud/static/hdfs/harmonizome/data/jasparpwm/gene_attribute_edges.txt.gz",
      cache_name = sprintf("jaspar_%s_gene_attribute_edges.txt.gz", species),
      cache_dir = raw_dir,
      use_cache = use_cache,
      quiet = quiet
    )
  )

  results <- vector("list", length(databases))
  names(results) <- databases

  for (db in databases) {
    tsv_path <- file.path(std_dir, sprintf("%s_%s_standardized.tsv", tolower(db), species))
    csv_path <- file.path(std_dir, sprintf("%s_%s_standardized.csv", tolower(db), species))
    cache_summary_path <- file.path(std_dir, sprintf("%s_%s_source_info.tsv", tolower(db), species))

    if (use_cache && file.exists(tsv_path)) {
      if (!quiet) {
        message(sprintf("Using cached standardized table: %s", tsv_path))
      }

      data <- utils::read.delim(
        tsv_path,
        sep = "\t",
        stringsAsFactors = FALSE,
        check.names = FALSE,
        quote = ""
      )

      raw_path <- NA_character_
      source_url <- NA_character_
      if (file.exists(cache_summary_path)) {
        info_df <- utils::read.delim(
          cache_summary_path,
          sep = "\t",
          stringsAsFactors = FALSE,
          check.names = FALSE,
          quote = ""
        )
        if (nrow(info_df) > 0) {
          raw_path <- info_df$raw_path[[1]]
          source_url <- info_df$source_url[[1]]
        }
      }

      results[[db]] <- list(
        data = data,
        raw_path = raw_path,
        source_url = source_url,
        from_cache = TRUE,
        tsv_path = tsv_path,
        csv_path = csv_path
      )
      next
    }

    if (!quiet) {
      message(sprintf("Preparing %s ...", db))
    }

    res <- parsers[[db]]()
    data <- res$data

    if (write_files) {
      write_tsv(data, tsv_path)
      write_csv(data, csv_path)
      write_tsv(
        data.frame(
          database = db,
          raw_path = paste(as.character(res$raw_path), collapse = ";"),
          source_url = paste(as.character(res$url), collapse = ";"),
          stringsAsFactors = FALSE
        ),
        cache_summary_path
      )
    }

    results[[db]] <- list(
      data = data,
      raw_path = res$raw_path,
      source_url = res$url,
      from_cache = res$from_cache,
      tsv_path = tsv_path,
      csv_path = csv_path
    )
  }

  summary_df <- data.frame(
    database = names(results),
    n_pairs = vapply(results, function(x) nrow(x$data), integer(1)),
    tsv_path = vapply(results, function(x) x$tsv_path, character(1)),
    csv_path = vapply(results, function(x) x$csv_path, character(1)),
    raw_path = vapply(results, function(x) paste(as.character(x$raw_path), collapse = ";"), character(1)),
    source_url = vapply(results, function(x) paste(as.character(x$source_url), collapse = ";"), character(1)),
    stringsAsFactors = FALSE
  )

  if (write_files) {
    write_tsv(summary_df, file.path(outdir, sprintf("download_summary_%s.tsv", species)))
    write_csv(summary_df, file.path(outdir, sprintf("download_summary_%s.csv", species)))
  }

  list(
    summary = summary_df,
    data = lapply(results, function(x) x$data),
    files = lapply(results, function(x) x[c("raw_path", "tsv_path", "csv_path")]),
    outdir = outdir
  )
}

read_source_table <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "csv") {
    data.table::fread(path, data.table = TRUE, showProgress = FALSE)
  } else {
    data.table::fread(path, sep = "\t", data.table = TRUE, showProgress = FALSE)
  }
}

merge_tf_sources <- function(source_dir = "tf_source_output/standardized",
                             species = "human",
                             databases = c("TRRUST", "RegNetwork", "ChEA", "ENCODE", "JASPAR"),
                             outdir = .tfregact_network_cache_dir(),
                             output_prefix = NULL,
                             write_files = TRUE) {
  species <- tolower(species)
  db_order <- c("TRRUST", "RegNetwork", "ChEA", "ENCODE", "JASPAR")
  databases <- unique(databases)
  bad_dbs <- setdiff(databases, db_order)
  if (length(bad_dbs) > 0) {
    stop(sprintf("Unsupported databases: %s", paste(bad_dbs, collapse = ", ")))
  }
  databases <- db_order[db_order %in% databases]

  source_dir <- normalizePath(source_dir, winslash = "/", mustWork = TRUE)
  outdir <- normalizePath(outdir, winslash = "/", mustWork = FALSE)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  weight_map <- c(
    TRRUST = 4L,
    RegNetwork = 3L,
    ChEA = 2L,
    ENCODE = 2L,
    JASPAR = 1L
  )

  source_files <- setNames(
    file.path(source_dir, sprintf("%s_%s_standardized.csv", tolower(databases), species)),
    databases
  )

  missing_files <- source_files[!file.exists(source_files)]
  if (length(missing_files) > 0) {
    stop(sprintf("Missing standardized files: %s", paste(unname(missing_files), collapse = ", ")))
  }

  tables <- lapply(names(source_files), function(db) {
    df <- read_source_table(source_files[[db]])
    keep_cols <- intersect(c("tf", "target", "effect", "source_db", "source_class"), colnames(df))
    df <- df[, ..keep_cols]
    if (!"source_db" %in% colnames(df)) {
      df[, source_db := db]
    }
    df
  })
  names(tables) <- names(source_files)

  all_edges <- data.table::rbindlist(tables, use.names = TRUE, fill = TRUE)
  all_edges <- unique(all_edges[, .(tf, target, effect, source_db, source_class)])

  merged_df <- all_edges[
    ,
    .(
      effect = resolve_effect(effect),
      # ChEA and ENCODE are correlated ChIP-based binding resources: record
      # one evidence-class score when either is present, rather than adding
      # two partially redundant scores.
      confidence_score = as.integer(
        4L * any(source_db == "TRRUST") +
          3L * any(source_db == "RegNetwork") +
          2L * any(source_db %in% c("ChEA", "ENCODE")) +
          1L * any(source_db == "JASPAR")
      ),
      supporting_databases = collapse_unique_sorted(source_db),
      source_class = collapse_unique_sorted(source_class),
      source_count = data.table::uniqueN(source_db),
      has_trrust = any(source_db == "TRRUST"),
      has_regnetwork = any(source_db == "RegNetwork"),
      has_chea = any(source_db == "ChEA"),
      has_encode = any(source_db == "ENCODE"),
      has_jaspar = any(source_db == "JASPAR")
    ),
    by = .(tf, target)
  ]
  data.table::setorder(merged_df, tf, target)

  summary_df <- data.frame(
    database = names(source_files),
    file = unname(source_files),
    weight = as.integer(weight_map[names(source_files)]),
    n_pairs = vapply(tables, nrow, integer(1)),
    stringsAsFactors = FALSE
  )

  if (is.null(output_prefix) || output_prefix == "") {
    output_prefix <- sprintf("tf_gene_merged_%s_weighted", species)
  }

  csv_path <- file.path(outdir, sprintf("%s.csv", output_prefix))
  tsv_path <- file.path(outdir, sprintf("%s.tsv", output_prefix))
  summary_path <- file.path(outdir, sprintf("%s_summary.tsv", output_prefix))

  if (write_files) {
    write_csv(merged_df, csv_path)
    write_tsv(merged_df, tsv_path)
    write_tsv(summary_df, summary_path)
  }

  list(
    merged = merged_df,
    summary = summary_df,
    files = list(csv = csv_path, tsv = tsv_path, summary = summary_path)
  )
}

download_gene_info_with_cache <- function(species, cache_dir, use_cache = TRUE) {
  species <- tolower(species)
  url <- switch(
    species,
    human = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz",
    mouse = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Mus_musculus.gene_info.gz",
    stop("species must be 'human' or 'mouse'.")
  )
  dest <- file.path(cache_dir, basename(url))
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  if (!use_cache || !file.exists(dest) || file.info(dest)$size <= 0) {
    download_file_then_replace(url, dest, quiet = FALSE)
  }
  dest
}

load_gene_symbol_maps <- function(species, cache_dir, use_cache = TRUE) {
  gene_info_path <- download_gene_info_with_cache(species, cache_dir, use_cache = use_cache)
  gi <- data.table::fread(
    cmd = sprintf("gzip -dc %s", shQuote(gene_info_path)),
    sep = "\t",
    header = TRUE,
    quote = "",
    showProgress = FALSE
  )

  gi <- gi[type_of_gene == "protein-coding" & Symbol != "-" & Symbol != ""]
  id_map <- gi[, .(gene_id = as.character(GeneID), symbol = as.character(Symbol))]
  symbol_map <- gi[, .(alias = as.character(Symbol), symbol = as.character(Symbol))]

  synonym_rows <- gi[Synonyms != "-" & Synonyms != ""]
  if (nrow(synonym_rows) > 0) {
    syn_list <- strsplit(as.character(synonym_rows$Synonyms), "\\|", perl = TRUE)
    syn_dt <- data.table::rbindlist(
      lapply(seq_along(syn_list), function(i) {
        vals <- trimws(syn_list[[i]])
        vals <- vals[vals != "" & vals != "-"]
        if (length(vals) == 0) {
          return(NULL)
        }
        data.table::data.table(
          alias = vals,
          symbol = rep(as.character(synonym_rows$Symbol[[i]]), length(vals))
        )
      }),
      use.names = TRUE,
      fill = TRUE
    )
    if (!is.null(syn_dt) && nrow(syn_dt) > 0) {
      syn_dt <- unique(syn_dt)
      syn_dt <- syn_dt[, if (data.table::uniqueN(symbol) == 1) .(symbol = symbol[[1]]) else .SD[0], by = .(alias)]
      symbol_map <- unique(data.table::rbindlist(list(symbol_map, syn_dt), use.names = TRUE, fill = TRUE))
    }
  }

  list(id_map = unique(id_map), symbol_map = unique(symbol_map))
}

map_gene_values <- function(values, maps) {
  vals <- trimws(as.character(values))
  out <- rep(NA_character_, length(vals))

  is_numeric_id <- grepl("^[0-9]+$", vals)
  if (any(is_numeric_id)) {
    idx <- match(vals[is_numeric_id], maps$id_map$gene_id)
    out[is_numeric_id] <- maps$id_map$symbol[idx]
  }

  non_numeric <- !is_numeric_id
  if (any(non_numeric)) {
    idx <- match(vals[non_numeric], maps$symbol_map$alias)
    out[non_numeric] <- maps$symbol_map$symbol[idx]
  }

  out
}

clean_tf_network <- function(
  input_file = "tf_union_output/tf_gene_merged_human_weighted.csv",
  species = "human",
  outdir = .tfregact_network_cache_dir(),
  output_prefix = NULL,
  write_files = TRUE,
  use_cache = TRUE
) {
  input_file <- normalizePath(input_file, winslash = "/", mustWork = TRUE)
  outdir <- normalizePath(outdir, winslash = "/", mustWork = FALSE)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  cache_dir <- file.path(outdir, "annotation_cache")
  maps <- load_gene_symbol_maps(species, cache_dir, use_cache = use_cache)

  dt <- data.table::fread(input_file, data.table = TRUE, showProgress = FALSE)
  n_input_edges <- nrow(dt)
  required_cols <- c(
    "tf", "target", "effect", "confidence_score", "supporting_databases",
    "source_class", "source_count", "has_trrust", "has_regnetwork",
    "has_chea", "has_encode", "has_jaspar"
  )
  missing_cols <- setdiff(required_cols, colnames(dt))
  if (length(missing_cols) > 0) {
    stop(sprintf("Input file is missing required columns: %s", paste(missing_cols, collapse = ", ")))
  }

  dt[, tf_original := tf]
  dt[, target_original := target]
  dt[, tf := map_gene_values(tf_original, maps)]
  dt[, target := map_gene_values(target_original, maps)]

  dt <- dt[!is.na(tf) & !is.na(target) & tf != "" & target != ""]
  dt <- dt[tf != target]

  logical_cols <- c("has_trrust", "has_regnetwork", "has_chea", "has_encode", "has_jaspar")
  for (col in logical_cols) {
    dt[[col]] <- as.logical(dt[[col]])
    dt[[col]][is.na(dt[[col]])] <- FALSE
  }

  cleaned <- dt[
    ,
    .(
      effect = resolve_effect(effect),
      has_trrust = any(has_trrust),
      has_regnetwork = any(has_regnetwork),
      has_chea = any(has_chea),
      has_encode = any(has_encode),
      has_jaspar = any(has_jaspar),
      tf_original = collapse_unique_sorted(tf_original),
      target_original = collapse_unique_sorted(target_original)
    ),
    by = .(tf, target)
  ]

  # Maximum score is 10: TRRUST (4), RegNetwork (3), ChIP evidence from
  # ChEA and/or ENCODE (2, counted once), and JASPAR motif evidence (1).
  cleaned[, confidence_score := 4L * as.integer(has_trrust) +
                                3L * as.integer(has_regnetwork) +
                                2L * as.integer(has_chea | has_encode) +
                                1L * as.integer(has_jaspar)]

  cleaned[, supporting_databases := vapply(seq_len(.N), function(i) {
    dbs <- c(
      if (cleaned$has_trrust[[i]]) "TRRUST" else NULL,
      if (cleaned$has_regnetwork[[i]]) "RegNetwork" else NULL,
      if (cleaned$has_chea[[i]]) "ChEA" else NULL,
      if (cleaned$has_encode[[i]]) "ENCODE" else NULL,
      if (cleaned$has_jaspar[[i]]) "JASPAR" else NULL
    )
    collapse_unique_sorted(dbs)
  }, character(1))]

  cleaned[, source_class := vapply(seq_len(.N), function(i) {
    classes <- c(
      if (cleaned$has_trrust[[i]] || cleaned$has_regnetwork[[i]]) "curated_regulatory_db" else NULL,
      if (cleaned$has_chea[[i]] || cleaned$has_encode[[i]]) "chip_based_binding" else NULL,
      if (cleaned$has_jaspar[[i]]) "motif_prediction" else NULL
    )
    collapse_unique_sorted(classes)
  }, character(1))]

  cleaned[, source_count := as.integer(has_trrust) + as.integer(has_regnetwork) +
                            as.integer(has_chea) + as.integer(has_encode) + as.integer(has_jaspar)]

  data.table::setcolorder(cleaned, c(
    "tf", "target", "effect", "confidence_score", "supporting_databases",
    "source_class", "source_count", "has_trrust", "has_regnetwork",
    "has_chea", "has_encode", "has_jaspar", "tf_original", "target_original"
  ))
  data.table::setorder(cleaned, tf, target)

  summary_df <- data.table::data.table(
    input_file = input_file,
    species = species,
    n_input_edges = n_input_edges,
    n_clean_edges = nrow(cleaned),
    numeric_id_examples_mapped = sum(grepl("^[0-9]+$", cleaned$tf_original)) +
      sum(grepl("^[0-9]+$", cleaned$target_original))
  )

  if (is.null(output_prefix) || output_prefix == "") {
    output_prefix <- sprintf("tf_gene_merged_%s_weighted_clean", tolower(species))
  }

  csv_path <- file.path(outdir, sprintf("%s.csv", output_prefix))
  tsv_path <- file.path(outdir, sprintf("%s.tsv", output_prefix))
  summary_path <- file.path(outdir, sprintf("%s_summary.tsv", output_prefix))

  if (write_files) {
    write_csv(cleaned, csv_path)
    write_tsv(cleaned, tsv_path)
    write_tsv(summary_df, summary_path)
  }

  list(
    cleaned = cleaned,
    summary = summary_df,
    files = list(csv = csv_path, tsv = tsv_path, summary = summary_path)
  )
}

build_tf_network_prepare_dagitty_bundle <- function(cleaned_df, csv_path = "") {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required to build the DAGitty network bundle.")
  }

  keep_cols <- intersect(
    c("tf", "target", "effect", "confidence_score", "supporting_databases"),
    colnames(cleaned_df)
  )
  edges <- as.data.frame(cleaned_df[, ..keep_cols], stringsAsFactors = FALSE)
  graph <- igraph::graph_from_data_frame(edges[, c("tf", "target"), drop = FALSE], directed = TRUE)
  reverse_graph <- igraph::reverse_edges(graph)

  list(
    edges = edges,
    graph = graph,
    reverse_graph = reverse_graph,
    metadata = list(
      source_csv = csv_path,
      created_at = as.character(Sys.time())
    )
  )
}

build_tf_network_save_dagitty_bundle <- function(bundle, path) {
  dagitty_network_bundle <- bundle
  save(dagitty_network_bundle, file = path, compress = "xz")
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

build_tf_network_bundle_from_csv <- function(csv_path, bundle_path) {
  dt <- data.table::fread(csv_path, data.table = TRUE, showProgress = FALSE)
  bundle <- build_tf_network_prepare_dagitty_bundle(cleaned_df = dt, csv_path = csv_path)
  build_tf_network_save_dagitty_bundle(bundle, bundle_path)
  bundle
}

build_tf_network_default_prefix <- function(species, databases) {
  all_dbs <- c("TRRUST", "RegNetwork", "ChEA", "ENCODE", "JASPAR")
  if (identical(sort(databases), sort(all_dbs))) {
    return(sprintf("tf_gene_merged_%s_weighted_clean", tolower(species)))
  }
  sprintf(
    "tf_gene_network_%s_%s",
    tolower(species),
    paste(tolower(databases), collapse = "_")
  )
}

build_tf_network <- function(
  databases = c("TRRUST", "RegNetwork", "ChEA", "ENCODE", "JASPAR"),
  species = "human",
  source_outdir = "tf_source_output",
  outdir = "tf_union_output",
  output_prefix = NULL,
  use_cache = TRUE,
  write_files = TRUE,
  force_rebuild = FALSE,
  update = FALSE,
  quiet = FALSE
) {
  if (length(update) != 1L || is.na(update) || !is.logical(update)) {
    stop("`update` must be TRUE or FALSE.")
  }

  # `update = TRUE` means a full refresh: do not reuse the final map,
  # standardized source tables, raw downloads, or the gene annotation cache.
  effective_use_cache <- isTRUE(use_cache) && !isTRUE(update)

  db_order <- c("TRRUST", "RegNetwork", "ChEA", "ENCODE", "JASPAR")
  databases <- unique(trimws(as.character(databases)))
  databases <- databases[databases != ""]
  bad_dbs <- setdiff(databases, db_order)
  if (length(bad_dbs) > 0) {
    stop(sprintf("Unsupported databases: %s", paste(bad_dbs, collapse = ", ")))
  }
  databases <- db_order[db_order %in% databases]
  species <- tolower(species)

  if (is.null(output_prefix) || !nzchar(output_prefix)) {
    output_prefix <- build_tf_network_default_prefix(species, databases)
  }

  outdir <- normalizePath(outdir, winslash = "/", mustWork = FALSE)
  source_outdir <- normalizePath(source_outdir, winslash = "/", mustWork = FALSE)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(source_outdir, recursive = TRUE, showWarnings = FALSE)

  final_csv <- file.path(outdir, sprintf("%s.csv", output_prefix))
  final_tsv <- file.path(outdir, sprintf("%s.tsv", output_prefix))
  final_summary <- file.path(outdir, sprintf("%s_summary.tsv", output_prefix))
  final_dagitty_rdata <- file.path(outdir, "TF_Full_Map.RData")

  if (effective_use_cache && !force_rebuild && write_files && file.exists(final_csv) && file.info(final_csv)$size > 0 &&
      file.exists(final_dagitty_rdata) && file.info(final_dagitty_rdata)$size > 0) {
    if (!quiet) {
      message(sprintf("Using cached final table: %s", final_csv))
    }
    return(list(
      cleaned = NULL,
      files = list(csv = final_csv, tsv = final_tsv, summary = final_summary, dagitty_rdata = final_dagitty_rdata),
      source = NULL,
      merged = NULL,
      dagitty_bundle = NULL,
      cached = TRUE
    ))
  }

  if (effective_use_cache && !force_rebuild && write_files && file.exists(final_csv) && file.info(final_csv)$size > 0 &&
      (!file.exists(final_dagitty_rdata) || file.info(final_dagitty_rdata)$size <= 0)) {
    if (!quiet) {
      message(sprintf("Using cached final table and creating DAGitty bundle: %s", final_csv))
    }
    dagitty_bundle <- build_tf_network_bundle_from_csv(final_csv, final_dagitty_rdata)
    return(list(
      cleaned = NULL,
      files = list(csv = final_csv, tsv = final_tsv, summary = final_summary, dagitty_rdata = final_dagitty_rdata),
      source = NULL,
      merged = NULL,
      dagitty_bundle = dagitty_bundle,
      cached = TRUE
    ))
  }

  if (!quiet) {
    message("Step 1/3: download and standardize selected TF databases")
  }
  source_res <- download_tf_sources(
    species = species,
    databases = databases,
    outdir = source_outdir,
    use_cache = effective_use_cache,
    write_files = write_files,
    quiet = quiet
  )

  merged_prefix <- if (identical(output_prefix, sprintf("tf_gene_merged_%s_weighted_clean", species))) {
    sprintf("tf_gene_merged_%s_weighted", species)
  } else {
    paste0(output_prefix, "_merged")
  }

  if (!quiet) {
    message("Step 2/3: merge standardized TF edges and compute confidence scores")
  }
  merge_res <- merge_tf_sources(
    source_dir = file.path(source_outdir, "standardized"),
    species = species,
    databases = databases,
    outdir = outdir,
    output_prefix = merged_prefix,
    write_files = write_files
  )

  if (!quiet) {
    message("Step 3/3: clean TF network and keep gene-symbol version")
  }
  clean_res <- clean_tf_network(
    input_file = merge_res$files$csv,
    species = species,
    outdir = outdir,
    output_prefix = output_prefix,
    write_files = write_files,
    use_cache = effective_use_cache
  )

  dagitty_bundle <- build_tf_network_prepare_dagitty_bundle(
    cleaned_df = clean_res$cleaned,
    csv_path = clean_res$files$csv
  )
  if (write_files) {
    build_tf_network_save_dagitty_bundle(dagitty_bundle, final_dagitty_rdata)
  }
  clean_res$files$dagitty_rdata <- final_dagitty_rdata

  list(
    cleaned = clean_res$cleaned,
    files = clean_res$files,
    source = source_res,
    merged = merge_res,
    dagitty_bundle = dagitty_bundle,
    cached = FALSE
  )
}

build_tf_network_parse_args <- function(args) {
  get_arg <- function(flag, default = NULL) {
    hit <- grep(paste0("^", flag, "="), args, value = TRUE)
    if (length(hit) == 0) {
      return(default)
    }
    sub(paste0("^", flag, "="), "", hit[[1]])
  }

  dbs <- get_arg("--databases", "TRRUST,RegNetwork,ChEA,ENCODE,JASPAR")
  dbs <- trimws(unlist(strsplit(dbs, ",", fixed = TRUE)))
  dbs <- dbs[dbs != ""]

  list(
    databases = dbs,
    species = get_arg("--species", "human"),
    source_outdir = get_arg("--source_outdir", "tf_source_output"),
    outdir = get_arg("--outdir", "tf_union_output"),
    output_prefix = get_arg("--output_prefix", ""),
    use_cache = tolower(get_arg("--use_cache", "true")) != "false",
    write_files = tolower(get_arg("--write_files", "true")) != "false",
    force_rebuild = tolower(get_arg("--force_rebuild", "false")) == "true",
    update = tolower(get_arg("--update", "false")) == "true"
  )
}

if (sys.nframe() == 0) {
  args <- build_tf_network_parse_args(commandArgs(trailingOnly = TRUE))
  res <- build_tf_network(
    databases = args$databases,
    species = args$species,
    source_outdir = args$source_outdir,
    outdir = args$outdir,
    output_prefix = args$output_prefix,
    use_cache = args$use_cache,
    write_files = args$write_files,
    force_rebuild = args$force_rebuild,
    update = args$update,
    quiet = FALSE
  )

  message("Done.")
  message(sprintf("Final CSV file: %s", res$files$csv))
}
