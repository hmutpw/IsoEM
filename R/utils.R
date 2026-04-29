# =============================================================================
# utils.R  —  Internal helper functions (not exported)
# =============================================================================

# -----------------------------------------------------------------------------
# File I/O helpers
# -----------------------------------------------------------------------------

#' @keywords internal
#' @noRd
.check_file <- function(path, arg) {
  if (!file.exists(path))
    stop(sprintf("File not found for '%s': %s", arg, path), call. = FALSE)
}

#' Read TSV / TSV.gz via fread; uses zcat pipe for .gz to avoid 2GB limit
#' @keywords internal
#' @noRd
.read_tsv <- function(path, col.names = NULL, select = NULL,
                      colClasses = NULL, header = FALSE, nrows = Inf) {
  gz  <- grepl("\\.gz$", path, ignore.case = TRUE)
  args <- list(header=FALSE, nrows=nrows, showProgress=FALSE, data.table=TRUE)
  if (gz) args$cmd  <- paste("zcat", shQuote(path))
  else    args$file <- path
  if (!is.null(select))     args$select     <- select
  if (!is.null(col.names))  args$col.names  <- col.names
  if (!is.null(colClasses)) args$colClasses <- colClasses
  do.call(data.table::fread, args)
}


.read_tsv_chunked <- function(path, col.names = NULL, select = NULL,
                               colClasses = NULL, chunk_size = 1e6L,
                               FUN, header = FALSE) {
  gz  <- grepl("\\.gz$", path, ignore.case = TRUE)
  cmd <- if (gz) paste("zcat", shQuote(path)) else NULL

  # Build fread call dynamically — avoids NULL col.names/select causing errors
  .fread_chunk <- function(skip_n) {
    args <- list(
      header       = FALSE,
      nrows        = as.integer(chunk_size),
      skip         = skip_n,
      showProgress = FALSE,
      data.table   = TRUE
    )
    if (!is.null(cmd))        args$cmd        <- cmd
    else              args$file <- path
    if (!is.null(select))     args$select     <- select
    if (!is.null(col.names))  args$col.names  <- col.names
    if (!is.null(colClasses)) args$colClasses <- colClasses
    tryCatch(do.call(data.table::fread, args), error = function(e) NULL)
  }

  results <- list()
  skip    <- 0L
  idx     <- 1L

  repeat {
    chunk <- .fread_chunk(skip)
    if (is.null(chunk) || nrow(chunk) == 0L) break
    n_chunk <- nrow(chunk)

    res <- FUN(chunk)
    rm(chunk); gc(verbose = FALSE)

    if (!is.null(res) && nrow(res) > 0L) {
      results[[idx]] <- res
      idx <- idx + 1L
    }

    if (n_chunk < chunk_size) break
    skip <- skip + chunk_size
  }

  if (length(results) == 0L) return(data.table::data.table())
  data.table::rbindlist(results, use.names = TRUE, fill = TRUE)
}


#' Open a write connection (plain or gzip)
#' @keywords internal
#' @noRd
.write_con <- function(path, compress) {
  if (compress) gzfile(path, "wt") else file(path, "wt")
}

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

#' Estimate file row count from size (rough heuristic for progress reporting)
#' @keywords internal
#' @noRd
.estimate_rows <- function(path) {
  sz <- file.info(path)$size
  if (is.na(sz)) return(NA_integer_)
  # sample first 1000 lines to estimate bytes/row
  gz  <- grepl("\\.gz$", path, ignore.case = TRUE)
  cmd <- if (gz) sprintf("zcat %s | head -1000", shQuote(path))
         else    sprintf("head -1000 %s", shQuote(path))
  sample_lines <- tryCatch(system(cmd, intern = TRUE), error = function(e) NULL)
  if (is.null(sample_lines) || length(sample_lines) == 0L) return(NA_integer_)
  bytes_per_row <- nchar(paste(sample_lines, collapse = "\n")) / length(sample_lines)
  if (gz) sz <- sz * 4L  # rough decompression ratio
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
# Shell sort+join for large file merging
# -----------------------------------------------------------------------------

#' Join two TSV files on column 1 using R data.table (chunked for memory efficiency)
#' Loads file_a fully (counts, smaller), streams file_b in chunks (anno, larger)
#' @keywords internal
#' @noRd
.shell_sort_join <- function(file_a, file_b, outfile,
                              compress_out = TRUE, verbose = TRUE) {
  if (verbose) message("  Loading counts file for join key ...")

  # Load file_a fully (counts file: read_id | transcript_id)
  gz_a <- grepl("\\.gz$", file_a, ignore.case = TRUE)
  dt_a <- if (gz_a) {
    data.table::fread(cmd=paste("zcat", shQuote(file_a)),
      header=FALSE, showProgress=FALSE, data.table=TRUE)
  } else {
    data.table::fread(file_a, header=FALSE, showProgress=FALSE, data.table=TRUE)
  }
  data.table::setnames(dt_a, c("read_id", "transcript_id"))
  data.table::setkey(dt_a, read_id)
  valid_reads <- dt_a$read_id

  if (verbose) message("  Streaming annotation file and joining ...")

  # Stream file_b in chunks, join, write incrementally
  gz_b    <- grepl("\\.gz$", file_b, ignore.case = TRUE)
  con_out <- if (compress_out) gzfile(outfile, "wt") else file(outfile, "wt")
  on.exit(close(con_out), add = TRUE)

  chunk_size <- 2e6L
  skip       <- 0L
  first      <- TRUE

  repeat {
    n_chunk <- 0L
    chunk <- tryCatch({
      if (gz_b)
        data.table::fread(cmd=paste("zcat", shQuote(file_b)),
          header=FALSE, nrows=chunk_size, skip=skip,
          showProgress=FALSE, data.table=TRUE)
      else
        data.table::fread(file_b,
          header=FALSE, nrows=chunk_size, skip=skip,
          showProgress=FALSE, data.table=TRUE)
    }, error = function(e) data.table::data.table())

    if (nrow(chunk) == 0L) break
    n_chunk <- nrow(chunk)

    data.table::setnames(chunk, 1L, "read_id")
    chunk_f <- chunk[read_id %in% valid_reads]
    if (nrow(chunk_f) > 0L) {
      merged <- dt_a[chunk_f, on = "read_id", nomatch = NULL]
      utils::write.table(merged, con_out,
        sep="	", quote=FALSE, row.names=FALSE, col.names=FALSE,
        append=!first)
      first <- FALSE
      rm(merged)
    }
    rm(chunk, chunk_f); gc(verbose=FALSE)
    if (n_chunk < chunk_size) break
    skip <- skip + chunk_size
  }

    invisible(outfile)
}


# -----------------------------------------------------------------------------
# Integer mapping helpers
# -----------------------------------------------------------------------------

#' Build integer index map from a character vector
#' Returns named integer vector: name -> index
#' @keywords internal
#' @noRd
.build_map <- function(x) {
  u   <- unique(x)
  idx <- seq_along(u)
  stats::setNames(idx, u)
}

#' Apply integer map to a character vector
#' @keywords internal
#' @noRd
.apply_map <- function(x, map) {
  unname(map[x])
}

# -----------------------------------------------------------------------------
# GTF parser
# -----------------------------------------------------------------------------

#' Parse GTF: extract transcript_id, gene_id, is_novel
#' @keywords internal
#' @noRd
.parse_gtf <- function(gtf_file) {
  con <- if (grepl("\\.gz$", gtf_file)) gzfile(gtf_file, "rt")
         else file(gtf_file, "rt")
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
  n_ec       <- nrow(ec_sub)
  ec_cnt_v   <- ec_sub$count
  ec_tx_list <- ec_sub$t_indices   # list column

  ec_sizes   <- lengths(ec_tx_list)
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

#' Setup outdir and temp subdir, return temp path
#' @keywords internal
#' @noRd
.setup_dirs <- function(outdir) {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  tmp_dir <- file.path(outdir, "temp")
  if (!dir.exists(tmp_dir)) dir.create(tmp_dir)
  tmp_dir
}

#' Remove temp dir if keep_temp = FALSE
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
# Regex extraction helper
# -----------------------------------------------------------------------------

#' Extract named capture groups from a character vector using a regex
#' Returns a data.table with one column per capture group name
#' @keywords internal
#' @noRd
.regex_extract <- function(x, pattern) {
  m <- regmatches(x, regexpr(pattern, x, perl = TRUE))
  # extract named capture groups
  cap_names <- .get_capture_names(pattern)
  if (length(cap_names) == 0L)
    stop("pattern must contain named capture groups, e.g. (?P<barcode>[ACGT]+)",
         call. = FALSE)
  result <- lapply(cap_names, function(nm) {
    sub(paste0(".*(?P<", nm, ">([^)]+)).*"), "\\1", m, perl = TRUE)
  })
  names(result) <- cap_names
  as.data.table(result)
}

#' Extract capture group names from a regex pattern
#' @keywords internal
#' @noRd
.get_capture_names <- function(pattern) {
  m <- gregexpr("\\(\\?P<([^>]+)>", pattern, perl = TRUE)
  starts <- m[[1L]]
  if (starts[1L] == -1L) return(character(0))
  sapply(seq_along(starts), function(i) {
    sub("\\(\\?P<([^>]+)>.*", "\\1",
        substr(pattern, starts[i], starts[i] + attr(m[[1L]], "match.length")[i]))
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
