# =============================================================================
# utils.R  —  Internal helper functions (not exported)
# =============================================================================

# -----------------------------------------------------------------------------
# Utility operators
# -----------------------------------------------------------------------------

#' Null-coalescing operator
#' @keywords internal
#' @noRd
`%||%` <- function(a, b) if (!is.null(a)) a else b

# -----------------------------------------------------------------------------
# Parallel backend — cross-platform (Windows + Linux)
# -----------------------------------------------------------------------------

#' Run lapply in parallel using a cluster (Windows + Linux compatible)
#' Falls back to sequential on n_cores == 1 or single item.
#' @keywords internal
#' @noRd
.par_lapply <- function(x, FUN, n_cores, ...) {
  n <- length(x)
  if (n_cores <= 1L || n <= 1L) return(lapply(x, FUN, ...))
  n_workers <- min(n_cores, n)
  cl <- parallel::makeCluster(n_workers)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterExport(cl, ls(envir = parent.env(environment())),
                          envir = parent.env(environment()))
  parallel::parLapply(cl, x, FUN, ...)
}

#' Export objects to a cluster for use in parLapply workers
#' @keywords internal
#' @noRd
.cluster_export <- function(cl, vars, envir = parent.frame()) {
  parallel::clusterExport(cl, vars, envir = envir)
}

# -----------------------------------------------------------------------------
# File I/O helpers
# -----------------------------------------------------------------------------

#' @keywords internal
#' @noRd
.check_file <- function(path, arg) {
  if (!file.exists(path))
    stop(sprintf("File not found for '%s': %s", arg, path), call. = FALSE)
}

#' Open a read connection for plain or gzipped file (cross-platform)
#' @keywords internal
#' @noRd
.read_con <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt")
  else file(path, "rt")
}

#' Read TSV / TSV.gz via fread — cross-platform, no external commands
#' @keywords internal
#' @noRd
.read_tsv <- function(path, col.names = NULL, select = NULL,
                      colClasses = NULL, header = FALSE, nrows = Inf) {
  args <- list(
    file         = path,
    header       = FALSE,
    nrows        = nrows,
    showProgress = FALSE,
    data.table   = TRUE
  )
  if (!is.null(select))     args$select     <- select
  if (!is.null(col.names))  args$col.names  <- col.names
  if (!is.null(colClasses)) args$colClasses <- colClasses
  do.call(data.table::fread, args)
}

#' Read TSV in chunks via a persistent R connection — truly sequential,
#' works correctly for both plain and .gz files without re-seeking.
#' @keywords internal
#' @noRd
.read_tsv_chunked <- function(path, col.names = NULL, select = NULL,
                               colClasses = NULL, chunk_size = 1e6L,
                               FUN, header = FALSE) {
  con <- .read_con(path)
  on.exit(close(con), add = TRUE)

  results <- list()
  idx     <- 1L

  malformed_total <- 0L
  repeat {
    # Read chunk as text lines, then parse with fread(text=)
    lines <- readLines(con, n = as.integer(chunk_size), warn = FALSE)
    if (length(lines) == 0L) break
    n_chunk <- length(lines)

    # Pre-filter malformed lines BEFORE fread to avoid early-stop data loss.
    # fread silently truncates on field-count mismatch (only warns), so we
    # validate field count up front and drop bad lines explicitly.
    expected_cols <- if (!is.null(select)) max(select) else NULL
    if (!is.null(expected_cols)) {
      tab_counts <- lengths(regmatches(lines, gregexpr("\t", lines, fixed = TRUE)))
      good_mask  <- (tab_counts + 1L) >= expected_cols
      n_bad      <- sum(!good_mask)
      if (n_bad > 0L) {
        malformed_total <- malformed_total + n_bad
        lines <- lines[good_mask]
      }
      if (length(lines) == 0L) next
    }

    args <- list(
      text         = paste(lines, collapse = "\n"),
      header       = FALSE,
      fill         = TRUE,    # tolerate trailing-field variations
      showProgress = FALSE,
      data.table   = TRUE
    )
    if (!is.null(select))     args$select     <- select
    if (!is.null(col.names))  args$col.names  <- col.names
    if (!is.null(colClasses)) args$colClasses <- colClasses

    # Capture both errors and warnings so a single bad line does not
    # silently truncate an entire chunk (data.table::fread default).
    chunk <- tryCatch(
      withCallingHandlers(
        do.call(data.table::fread, args),
        warning = function(w) {
          if (grepl("Stopped early|Expected.*fields", conditionMessage(w))) {
            malformed_total <<- malformed_total + 1L
          }
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) NULL
    )
    rm(lines)

    if (!is.null(chunk) && nrow(chunk) > 0L) {
      res <- FUN(chunk)
      rm(chunk); gc(verbose = FALSE)
      if (!is.null(res) && nrow(res) > 0L) {
        results[[idx]] <- res
        idx <- idx + 1L
      }
    }

    if (n_chunk < chunk_size) break
  }

  if (malformed_total > 0L)
    warning(sprintf(
      paste0("%d malformed line(s) skipped while reading %s. ",
             "Check the file for inconsistent field counts ",
             "(extra/missing tabs)."),
      malformed_total, basename(path)
    ), call. = FALSE)

  if (length(results) == 0L) return(data.table::data.table())
  data.table::rbindlist(results, use.names = TRUE, fill = TRUE)
}

#' Open a write connection (plain or gzip)
#' @keywords internal
#' @noRd
.write_con <- function(path, compress) {
  if (compress) gzfile(path, "wt") else file(path, "wt")
}

#' Write data.table to TSV (optionally gzipped)
#' @keywords internal
#' @noRd
.write_tsv <- function(dt, path, compress = FALSE, col.names = TRUE) {
  con <- .write_con(path, compress)
  on.exit(close(con))
  utils::write.table(dt, con,
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = col.names)
  invisible(path)
}

#' File extension helper
#' @keywords internal
#' @noRd
.ext <- function(base, compress) {
  if (compress) paste0(base, ".gz") else base
}

#' Estimate file row count (cross-platform, uses R connection)
#' @keywords internal
#' @noRd
.estimate_rows <- function(path) {
  sz <- file.info(path)$size
  if (is.na(sz)) return(NA_integer_)
  gz <- grepl("\\.gz$", path, ignore.case = TRUE)
  con <- tryCatch(.read_con(path), error = function(e) NULL)
  if (is.null(con)) return(NA_integer_)
  on.exit(close(con), add = TRUE)
  sample_lines <- tryCatch(readLines(con, n = 200L, warn = FALSE),
                           error = function(e) NULL)
  if (is.null(sample_lines) || length(sample_lines) == 0L) return(NA_integer_)
  bytes_per_row <- nchar(paste(sample_lines, collapse = "\n")) /
                   length(sample_lines)
  if (gz) sz <- sz * 4L
  as.integer(sz / max(bytes_per_row, 1))
}

#' Human-readable file size string
#' @keywords internal
#' @noRd
.size_str <- function(path) {
  sz <- file.info(path)$size
  if (is.na(sz)) return("unknown size")
  if (sz >= 1e9) sprintf("%.1f GB", sz / 1e9)
  else if (sz >= 1e6) sprintf("%.0f MB", sz / 1e6)
  else sprintf("%.0f KB", sz / 1e3)
}

#' Human-readable row count string
#' @keywords internal
#' @noRd
.nrow_str <- function(n) {
  if (is.na(n)) return("? rows")
  if (n >= 1e9) sprintf("%.1fB", n / 1e9)
  else if (n >= 1e6) sprintf("%.0fM", n / 1e6)
  else if (n >= 1e3) sprintf("%.0fK", n / 1e3)
  else as.character(n)
}

# -----------------------------------------------------------------------------
# Large file join — R data.table, cross-platform, chunked
# -----------------------------------------------------------------------------

#' Join two TSV files on column 1 using R data.table chunked join.
#' Loads file_a into memory, streams file_b in chunks.
#' Cross-platform: no shell commands.
#' @keywords internal
#' @noRd
.shell_sort_join <- function(file_a, file_b, outfile,
                              compress_out = TRUE, verbose = TRUE) {
  if (verbose) message("  Loading counts file for join key ...")
  dt_a <- data.table::fread(file_a, header = FALSE,
                              showProgress = FALSE, data.table = TRUE)
  data.table::setnames(dt_a, c("read_id", "transcript_id"))
  data.table::setkey(dt_a, read_id)
  valid_reads <- dt_a$read_id

  if (verbose) message("  Streaming annotation file and joining ...")
  con_out <- if (compress_out) gzfile(outfile, "wt") else file(outfile, "wt")
  on.exit(close(con_out), add = TRUE)

  con_in     <- .read_con(file_b)
  on.exit(close(con_in), add = TRUE)
  chunk_size <- 2e6L
  first      <- TRUE

  repeat {
    lines <- readLines(con_in, n = as.integer(chunk_size), warn = FALSE)
    if (length(lines) == 0L) break
    n_chunk <- length(lines)

    chunk <- tryCatch(
      data.table::fread(text = paste(lines, collapse = "\n"),
                        header = FALSE, showProgress = FALSE, data.table = TRUE),
      error = function(e) data.table::data.table()
    )
    rm(lines)
    if (nrow(chunk) == 0L) { if (n_chunk < chunk_size) break; next }

    data.table::setnames(chunk, 1L, "read_id")
    chunk_f <- chunk[read_id %in% valid_reads]
    if (nrow(chunk_f) > 0L) {
      merged <- dt_a[chunk_f, on = "read_id", nomatch = NULL]
      utils::write.table(merged, con_out, sep = "\t", quote = FALSE,
                         row.names = FALSE, col.names = FALSE, append = !first)
      first <- FALSE
      rm(merged)
    }
    rm(chunk, chunk_f); gc(verbose = FALSE)
    if (n_chunk < chunk_size) break
  }
  invisible(outfile)
}

# -----------------------------------------------------------------------------
# Integer mapping helpers
# -----------------------------------------------------------------------------

#' @keywords internal
#' @noRd
.build_map <- function(x) {
  u <- unique(x)
  stats::setNames(seq_along(u), u)
}

#' @keywords internal
#' @noRd
.apply_map <- function(x, map) unname(map[x])

# -----------------------------------------------------------------------------
# GTF parser
# -----------------------------------------------------------------------------

#' Parse GTF: extract transcript_id, gene_id, is_novel
#' @keywords internal
#' @noRd
.parse_gtf <- function(gtf_file) {
  con <- .read_con(gtf_file)
  on.exit(close(con))
  lines    <- readLines(con)
  tx_lines <- lines[grepl("\ttranscript\t", lines)]
  if (length(tx_lines) == 0L)
    tx_lines <- lines[grepl("\texon\t", lines)]

  .extract_attr <- function(field, text) {
    pattern <- paste0(field, ' "([^"]+)"')
    m   <- regmatches(text, regexpr(pattern, text, perl = TRUE))
    out <- rep(NA_character_, length(text))
    hit <- nchar(m) > 0
    out[hit] <- sub(paste0('.*', field, ' "([^"]+)".*'), "\\1", m[hit])
    out
  }

  meta <- data.table::data.table(
    transcript_id = .extract_attr("transcript_id", tx_lines),
    gene_id       = .extract_attr("gene_id",       tx_lines),
    is_novel      = grepl("novel|NOVEL", tx_lines, ignore.case = FALSE)
  )
  meta <- unique(meta, by = "transcript_id")
  meta[!is.na(transcript_id)]
}

# -----------------------------------------------------------------------------
# Core EM algorithm
# -----------------------------------------------------------------------------

#' EM for one group (one sample or one cell)
#' @param ec_sub   data.table: ec_id | t_indices (list col) | count
#' @param n_tx     integer: total transcripts in universe
#' @param max_iter integer
#' @param tol      numeric convergence tolerance
#' @return list(counts, n_iter, converged, n_ec, n_unique_ec)
#' @keywords internal
#' @noRd
.em_core <- function(ec_sub, n_tx, max_iter, tol) {
  n_ec        <- nrow(ec_sub)
  ec_cnt_v    <- ec_sub$count
  ec_tx_list  <- ec_sub$t_indices
  ec_sizes    <- lengths(ec_tx_list)
  n_unique_ec <- sum(ec_sizes == 1L)

  counts    <- rep(1 / n_tx, n_tx)
  converged <- FALSE
  n_iter    <- 0L

  for (iter in seq_len(max_iter)) {
    n_iter <- iter
    alloc  <- numeric(n_tx)
    for (k in seq_len(n_ec)) {
      txs   <- ec_tx_list[[k]]
      w     <- counts[txs]
      w_sum <- sum(w)
      if (w_sum < .Machine$double.eps)
        w <- rep(1 / length(txs), length(txs))
      else
        w <- w / w_sum
      alloc[txs] <- alloc[txs] + w * ec_cnt_v[k]
    }
    total <- sum(alloc)
    if (total < .Machine$double.eps) break
    delta  <- sum(abs(alloc - counts)) / total
    counts <- alloc
    if (delta < tol) { converged <- TRUE; break }
  }

  list(counts = counts, n_iter = n_iter, converged = converged,
       n_ec = n_ec, n_unique_ec = n_unique_ec)
}

# -----------------------------------------------------------------------------
# outdir / temp dir management
# -----------------------------------------------------------------------------

#' @keywords internal
#' @noRd
.setup_dirs <- function(outdir) {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  tmp_dir <- file.path(outdir, "temp")
  if (!dir.exists(tmp_dir)) dir.create(tmp_dir)
  tmp_dir
}

#' @keywords internal
#' @noRd
.cleanup_temp <- function(outdir, keep_temp, verbose) {
  tmp_dir <- file.path(outdir, "temp")
  if (!keep_temp && dir.exists(tmp_dir)) {
    unlink(tmp_dir, recursive = TRUE)
    if (verbose) message("  Temp directory removed: ", tmp_dir)
  } else if (keep_temp && dir.exists(tmp_dir)) {
    if (verbose) message("  Temp directory kept: ", tmp_dir)
  }
}

# -----------------------------------------------------------------------------
# Regex extraction helpers
# -----------------------------------------------------------------------------

#' Extract capture group names from a regex pattern (e.g. (?P<name>...))
#' @keywords internal
#' @noRd
.get_capture_names <- function(pattern) {
  m      <- gregexpr("\\(\\?P<([^>]+)>", pattern, perl = TRUE)
  starts <- m[[1L]]
  if (starts[1L] == -1L) return(character(0))
  sapply(seq_along(starts), function(i) {
    sub("\\(\\?P<([^>]+)>.*", "\\1",
        substr(pattern, starts[i],
               starts[i] + attr(m[[1L]], "match.length")[i]))
  })
}

#' Extract named capture groups from a character vector using PCRE regex
#' Returns named list of character vectors, one per group name.
#' @keywords internal
#' @noRd
.extract_named_groups <- function(x, pattern, groups) {
  matched <- regexpr(pattern, x, perl = TRUE)
  lapply(stats::setNames(groups, groups), function(grp) {
    vals    <- rep(NA_character_, length(x))
    hit     <- matched > 0
    if (!any(hit)) return(vals)
    starts  <- attr(matched, "capture.start")
    lengths <- attr(matched, "capture.length")
    names_c <- attr(matched, "capture.names")
    idx     <- which(names_c == grp)
    if (length(idx) == 0L) return(vals)
    vals[hit] <- substr(x[hit],
                        starts[hit, idx],
                        starts[hit, idx] + lengths[hit, idx] - 1L)
    vals
  })
}

# -----------------------------------------------------------------------------
# gzip helper (no R.utils dependency)
# -----------------------------------------------------------------------------

#' @keywords internal
#' @noRd
.gzip_file <- function(src, dest) {
  buf <- readBin(src, "raw", file.info(src)$size)
  con <- gzfile(dest, "wb")
  writeBin(buf, con)
  close(con)
  file.remove(src)
  invisible(dest)
}
