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

#' Estimate file row count
#'
#' For gzipped files, uses a fast file-size heuristic (no decompression)
#' to avoid the ~30s overhead of opening a gzip R connection.
#' For plain files, samples a few lines for a more accurate estimate.
#' @keywords internal
#' @noRd
.estimate_rows <- function(path) {
  sz <- file.info(path)$size
  if (is.na(sz)) return(NA_integer_)
  gz <- grepl("\\.gz$", path, ignore.case = TRUE)

  if (gz) {
    # Fast heuristic: assume ~4x compression, ~60 bytes/row (typical TSV)
    # Avoids opening gzip connection which is very slow on large files
    as.integer(sz * 4 / 60)
  } else {
    con <- tryCatch(file(path, "rt"), error = function(e) NULL)
    if (is.null(con)) return(NA_integer_)
    on.exit(close(con), add = TRUE)
    sample_lines <- tryCatch(readLines(con, n = 200L, warn = FALSE),
                             error = function(e) NULL)
    if (is.null(sample_lines) || length(sample_lines) == 0L) return(NA_integer_)
    bytes_per_row <- nchar(paste(sample_lines, collapse = "\n")) /
                     length(sample_lines)
    as.integer(sz / max(bytes_per_row, 1))
  }
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
#'
#' Uses readLines + grepl to filter to transcript lines, then vectorised
#' sub() directly on the lines for attribute extraction (no intermediate
#' fread/paste step).
#' @keywords internal
#' @noRd
.parse_gtf <- function(gtf_file) {
  con <- .read_con(gtf_file)
  on.exit(close(con))
  lines <- readLines(con)

  # Filter: drop comments, keep only transcript (or exon) features
  lines <- lines[!startsWith(lines, "#")]
  tx_lines <- lines[grepl("\ttranscript\t", lines, fixed = FALSE)]
  if (length(tx_lines) == 0L)
    tx_lines <- lines[grepl("\texon\t", lines, fixed = FALSE)]
  rm(lines)

  if (length(tx_lines) == 0L)
    return(data.table::data.table(
      transcript_id = character(), gene_id = character(),
      is_novel = logical()
    ))

  # Vectorised attribute extraction — apply sub() directly to full lines
  tx_id <- sub('.*transcript_id "([^"]+)".*', "\\1", tx_lines, perl = TRUE)
  g_id  <- sub('.*gene_id "([^"]+)".*',       "\\1", tx_lines, perl = TRUE)
  rm(tx_lines)

  # sub() returns the original string when no match — mark as NA
  long_mask <- nchar(tx_id) > 200L
  tx_id[long_mask] <- NA_character_
  long_mask <- nchar(g_id) > 200L
  g_id[long_mask]  <- NA_character_

  meta <- data.table::data.table(
    transcript_id = tx_id,
    gene_id       = g_id,
    is_novel      = grepl("novel", tx_id, ignore.case = TRUE)
  )
  meta <- unique(meta, by = "transcript_id")
  meta[!is.na(transcript_id)]
}

# -----------------------------------------------------------------------------
# Core EM algorithm
# -----------------------------------------------------------------------------

#' EM for one group (one sample or one cell)
#'
#' Vectorised implementation: unique ECs (ec_size == 1) are pre-computed once
#' (no iteration needed). Multi-mapping ECs are flattened into parallel arrays
#' and processed with rowsum() — no per-EC R for-loop.
#'
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

  # --- Pre-compute unique-EC allocations (constant across iterations) --------
  uniq_mask  <- ec_sizes == 1L
  uniq_alloc <- numeric(n_tx)
  if (any(uniq_mask)) {
    uniq_t <- unlist(ec_tx_list[uniq_mask], use.names = FALSE)
    uniq_c <- ec_cnt_v[uniq_mask]
    agg    <- rowsum(uniq_c, uniq_t, reorder = FALSE)
    uniq_alloc[as.integer(rownames(agg))] <- agg[, 1L]
  }

  # --- Multi-mapping ECs: flatten for vectorised ops -------------------------
  multi_idx <- which(!uniq_mask)
  n_multi   <- length(multi_idx)

  if (n_multi == 0L) {
    return(list(counts = uniq_alloc, n_iter = 0L, converged = TRUE,
                n_ec = n_ec, n_unique_ec = n_unique_ec))
  }

  multi_list  <- ec_tx_list[multi_idx]
  multi_cnt   <- ec_cnt_v[multi_idx]
  multi_sizes <- ec_sizes[multi_idx]

  flat_t      <- unlist(multi_list, use.names = FALSE)      # transcript indices
  flat_k      <- rep.int(seq_len(n_multi), multi_sizes)     # which multi-EC
  flat_c      <- rep.int(multi_cnt, multi_sizes)            # count per entry
  flat_inv_sz <- 1 / rep.int(multi_sizes, multi_sizes)      # uniform fallback

  # --- EM iterations (only on multi-mapping ECs) -----------------------------
  counts    <- rep(1 / n_tx, n_tx)
  converged <- FALSE
  n_iter    <- 0L

  for (iter in seq_len(max_iter)) {
    n_iter <- iter

    # E-step: vectorised weight computation
    flat_w    <- counts[flat_t]
    ec_sums   <- rowsum(flat_w, flat_k, reorder = FALSE)[, 1L]
    ec_sums_e <- ec_sums[flat_k]

    # Normalise (uniform fallback for zero-sum ECs)
    flat_wn   <- flat_w / ec_sums_e
    zero_mask <- ec_sums_e < .Machine$double.eps
    if (any(zero_mask)) flat_wn[zero_mask] <- flat_inv_sz[zero_mask]

    # M-step: scatter-add fractional allocations via rowsum
    alloc_raw   <- rowsum(flat_wn * flat_c, flat_t, reorder = FALSE)
    multi_alloc <- numeric(n_tx)
    multi_alloc[as.integer(rownames(alloc_raw))] <- alloc_raw[, 1L]

    alloc <- uniq_alloc + multi_alloc
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
