# =============================================================================
# build_ec.R  --  EC construction for bulk and single-cell modes
# =============================================================================

# -----------------------------------------------------------------------------
# build_ec()  --  Bulk
# -----------------------------------------------------------------------------

#' Build equivalence classes for bulk RNA-seq quantification
#'
#' Reads IsoQuant counts, merges group (sample) information, and constructs
#' an \code{IsoEMEC} object ready for \code{\link{run_em}}.
#'
#' @param input An \code{IsoEMInput} object from \code{\link{prepare_isoem}},
#'   with \code{mode} set to \code{"bulk_single"} or \code{"bulk_multi"}.
#' @param keep_temp Logical. Keep the \code{temp/} directory after EC
#'   construction (default FALSE).
#' @param chunk_size Integer. Rows per chunk when reading large files
#'   (default 5e6).
#' @param n_cores Integer. Parallel workers for EC construction across
#'   samples (default 1). Uses \code{parallel::makeCluster} --
#'   compatible with Windows and Linux.
#' @param verbose Logical (default TRUE).
#'
#' @return An \code{IsoEMEC} object containing:
#' \itemize{
#'   \item \code{ec_table}  -- data.table with columns \code{group_id} (sample),
#'     \code{ec_id}, \code{t_indices} (list of integer transcript indices),
#'     \code{count}
#'   \item \code{tx_map}    -- character vector mapping integer index to transcript ID
#'   \item \code{group_ids} -- all sample IDs
#'   \item \code{n_tx}      -- number of transcripts in the universe
#' }
#' The \code{ec_table} is carried through \code{\link{run_em}} into the
#' \code{IsoEMResult} and written to \code{ec_table.tsv} by
#' \code{\link{write_isoem}}.
#'
#' @seealso \code{\link{prepare_isoem}}, \code{\link{run_em}},
#'   \code{\link{write_isoem}}
#' @export
#'
#' @examples
#' counts_f <- system.file("extdata", "toy_counts.tsv.gz",  package = "IsoEM")
#' gtf_f    <- system.file("extdata", "toy.gtf.gz",         package = "IsoEM")
#' input    <- prepare_isoem(counts_f, gtf_f,
#'   mode = "bulk_single", sample_id = "toy")
#' ec <- build_ec(input)
#' print(ec)
#' # saveRDS(ec, "ec.rds")  # save for parameter tuning
build_ec <- function(input,
                     keep_temp  = FALSE,
                     chunk_size = 5e6L,
                     n_cores    = 1L,
                     verbose    = TRUE) {

  stopifnot(inherits(input, "IsoEMInput"))
  if (!input$validated)
    stop("IsoEMInput is not validated. Re-run prepare_isoem().", call. = FALSE)
  if (input$mode == "sc")
    stop("Use build_sc_ec() for sc mode.", call. = FALSE)

  .setup_dirs(input$outdir)
  on.exit(.cleanup_temp(input$outdir, keep_temp, verbose), add = TRUE)

  if (verbose) message("=== build_ec (bulk mode) ===")

  if (verbose) message("[1/3] Reading and integerising counts ...")
  int_data <- .integerise_counts(input$counts_file, verbose)

  if (verbose) message("[2/3] Preparing group information ...")
  dt_group <- if (!is.null(input$sample_id)) {
    .group_from_single(int_data$dt, input$sample_id)
  } else if (!is.null(input$anno_file)) {
    .group_from_anno(input$anno_file, int_data$read_map,
                     mode = "bulk_multi", chunk_size = chunk_size,
                     verbose = verbose)
  } else {
    .group_from_regex(input$counts_file, int_data$read_map,
                      pattern = input$pattern, bc_pattern = NULL,
                      umi_pattern = NULL, chunk_size = chunk_size,
                      mode = "bulk", verbose = verbose)
  }

  if (verbose) message("[3/3] Building equivalence classes ...")
  ec <- .build_ec_from_ints(
    dt_counts = int_data$dt, dt_group = dt_group,
    tx_map    = int_data$tx_map, mode = "bulk", unit = "read",
    gtf_file  = input$gtf_file, n_cores = n_cores, verbose = verbose
  )
  rm(int_data, dt_group); gc(verbose = FALSE)
  if (verbose) message("EC construction complete.")
  ec
}


# -----------------------------------------------------------------------------
# build_sc_ec()  --  Single-cell
# -----------------------------------------------------------------------------

#' Build equivalence classes for single-cell RNA-seq quantification
#'
#' Reads IsoQuant counts, merges barcode/UMI information, applies optional
#' barcode whitelist filtering, and constructs an \code{IsoEMEC} object
#' ready for \code{\link{run_em}}.
#'
#' @param input An \code{IsoEMInput} object from \code{\link{prepare_isoem}},
#'   with \code{mode = "sc"}.
#' @param keep_temp Logical. Keep \code{temp/} directory (default FALSE).
#' @param chunk_size Integer. Rows per chunk (default 5e6).
#' @param n_cores Integer. Parallel workers (default 1).
#'   Uses \code{parallel::makeCluster} -- compatible with Windows and Linux.
#' @param verbose Logical (default TRUE).
#'
#' @return An \code{IsoEMEC} object with \code{mode = "sc"} containing:
#' \itemize{
#'   \item \code{ec_table}  -- data.table with columns \code{group_id} (barcode),
#'     \code{ec_id}, \code{t_indices} (list of integer transcript indices),
#'     \code{count}
#'   \item \code{tx_map}    -- character vector mapping integer index to transcript ID
#'   \item \code{group_ids} -- all cell barcodes
#'   \item \code{n_tx}      -- number of transcripts in the universe
#' }
#' The \code{ec_table} is carried through \code{\link{run_em}} into the
#' \code{IsoEMSCResult} and written to \code{ec_table.tsv} by
#' \code{\link{write_sc_isoem}}.
#'
#' @seealso \code{\link{prepare_isoem}}, \code{\link{run_em}},
#'   \code{\link{write_sc_isoem}}
#' @export
#'
#' @examples
#' counts_f <- system.file("extdata", "toy_counts_multi.tsv", package = "IsoEM")
#' gtf_f    <- system.file("extdata", "toy.gtf.gz",            package = "IsoEM")
#' bc_f     <- system.file("extdata", "toy_bc_umi.tsv",        package = "IsoEM")
#' input    <- prepare_isoem(counts_f, gtf_f,
#'   mode = "sc", anno_file = bc_f, unit = "umi")
#' ec <- build_sc_ec(input)
#' print(ec)
#' # Save for reuse (avoid rebuilding on large datasets)
#' # saveRDS(ec, "ec.rds")
build_sc_ec <- function(input,
                         keep_temp  = FALSE,
                         chunk_size = 5e6L,
                         n_cores    = 1L,
                         verbose    = TRUE) {

  stopifnot(inherits(input, "IsoEMInput"))
  if (!input$validated)
    stop("IsoEMInput is not validated. Re-run prepare_isoem().", call. = FALSE)
  if (input$mode != "sc")
    stop("build_sc_ec() requires mode = 'sc'. Use build_ec() for bulk.",
         call. = FALSE)

  .setup_dirs(input$outdir)
  on.exit(.cleanup_temp(input$outdir, keep_temp, verbose), add = TRUE)

  if (verbose) message("=== build_sc_ec (single-cell mode) ===")
  if (verbose) message("  Unit: ", input$unit)

  if (verbose) message("[1/4] Reading and integerising counts ...")
  int_data <- .integerise_counts(input$counts_file, verbose)

  if (verbose) message("[2/4] Preparing barcode/UMI information ...")
  dt_group <- if (!is.null(input$anno_file)) {
    .group_from_anno(input$anno_file, int_data$read_map,
                     mode = "sc", chunk_size = chunk_size, verbose = verbose)
  } else {
    .group_from_regex(input$counts_file, int_data$read_map,
                      pattern = input$pattern, bc_pattern = input$bc_pattern,
                      umi_pattern = input$umi_pattern, chunk_size = chunk_size,
                      mode = "sc", verbose = verbose)
  }

  # barcode whitelist
  if (!is.null(input$barcodes_use)) {
    if (verbose) message("[3/4] Applying barcode whitelist ...")
    barcodes_vec <- if (length(input$barcodes_use) == 1L &&
                        file.exists(input$barcodes_use)) {
      bc_lines <- readLines(input$barcodes_use)
      bc_lines[nchar(trimws(bc_lines)) > 0]
    } else {
      input$barcodes_use
    }
    n_before <- data.table::uniqueN(dt_group$group_id)
    dt_group <- dt_group[group_id %in% barcodes_vec]
    n_after  <- data.table::uniqueN(dt_group$group_id)
    if (verbose)
      message(sprintf("  Whitelist: %d / %d barcodes retained.",
                      n_after, n_before))
  } else {
    if (verbose) message("[3/4] No whitelist -- using all barcodes.")
  }

  if (verbose) message("[4/4] Building equivalence classes ...")
  ec <- .build_ec_from_ints(
    dt_counts = int_data$dt, dt_group = dt_group,
    tx_map    = int_data$tx_map, mode = "sc", unit = input$unit,
    gtf_file  = input$gtf_file, n_cores = n_cores, verbose = verbose
  )
  rm(int_data, dt_group); gc(verbose = FALSE)
  if (verbose) message("EC construction complete.")
  ec
}


# =============================================================================
# Internal workers
# =============================================================================

#' Read counts + build integer maps for read_id and transcript_id
#'
#' Uses match() for integer encoding — faster than named-vector lookup
#' and avoids allocating a large named vector for read_map.
#' @keywords internal
#' @noRd
.integerise_counts <- function(counts_file, verbose) {
  dt <- .read_tsv(counts_file,
    col.names = c("read_id", "transcript_id"),
    select    = c(1L, 2L)
  )
  dt <- dt[transcript_id != "*"]
  if (nrow(dt) == 0L)
    stop("No valid records in counts_file after filtering '*'.", call. = FALSE)
  if (verbose)
    message(sprintf("  counts: %d records, %d unique reads, %d unique transcripts.",
                    nrow(dt), data.table::uniqueN(dt$read_id),
                    data.table::uniqueN(dt$transcript_id)))

  all_reads <- unique(dt$read_id)
  all_tx    <- sort(unique(dt$transcript_id))

  # match() builds an internal hash once — faster than named-vector [ lookup
  dt[, r_idx := match(read_id, all_reads)]
  dt[, t_idx := match(transcript_id, all_tx)]
  dt[, c("read_id", "transcript_id") := NULL]
  data.table::setkey(dt, r_idx)

  # read_map: named vector kept for downstream group-from-anno/regex lookups
  read_map <- stats::setNames(seq_along(all_reads), all_reads)
  tx_map   <- stats::setNames(seq_along(all_tx),    all_tx)

  list(dt = dt, read_map = read_map, tx_map = tx_map)
}

#' Prepare group table from an external annotation file
#'
#' Uses fread() directly for fast parallel gzip decompression (5-10x faster
#' than readLines-based chunked reading). Falls back to chunked reading if
#' direct read fails (e.g. malformed files, extreme memory pressure).
#'
#' sc:         anno_file cols = read_id | barcode | umi
#' bulk_multi: anno_file cols = read_id | sample_id
#' Returns data.table(r_idx, group_id [, umi])
#' @keywords internal
#' @noRd
.group_from_anno <- function(anno_file, read_map, mode, chunk_size, verbose) {
  is_sc    <- (mode == "sc")
  n_select <- if (is_sc) 3L else 2L

  # --- Fast path: direct fread (parallel gzip decompression) -----------------
  if (verbose) message("  Reading anno_file ...")
  dt <- tryCatch(
    .read_tsv(anno_file, select = seq_len(n_select)),
    error = function(e) NULL
  )

  if (is.null(dt) || nrow(dt) == 0L) {
    # Fallback to chunked reading for problematic files
    if (verbose)
      message("  Direct read failed, falling back to chunked reading ...")
    return(.group_from_anno_chunked(anno_file, read_map, mode,
                                     chunk_size, verbose))
  }

  if (is_sc)
    data.table::setnames(dt, seq_len(3L), c("read_id", "group_id", "umi"))
  else
    data.table::setnames(dt, seq_len(2L), c("read_id", "group_id"))

  # Filter to reads present in counts — use keyed data.table join (fast)
  read_dt <- data.table::data.table(
    read_id = names(read_map),
    r_idx   = unname(read_map)
  )
  data.table::setkey(read_dt, read_id)
  data.table::setkey(dt, read_id)
  dt <- read_dt[dt, on = "read_id", nomatch = NULL]
  dt[, read_id := NULL]
  rm(read_dt); gc(verbose = FALSE)

  if (verbose) message(sprintf("  anno_file: %d matching records.", nrow(dt)))
  dt
}

#' Fallback chunked reader for annotation files
#' @keywords internal
#' @noRd
.group_from_anno_chunked <- function(anno_file, read_map, mode,
                                      chunk_size, verbose) {
  valid_reads <- names(read_map)
  is_sc       <- (mode == "sc")
  n_select    <- if (is_sc) 3L else 2L

  if (verbose)
    message(sprintf("  Reading anno_file in chunks (select %d cols) ...",
                    n_select))

  .read_tsv_chunked(
    path       = anno_file,
    select     = seq_len(n_select),
    chunk_size = chunk_size,
    FUN = function(chunk) {
      if (is_sc)
        data.table::setnames(chunk, seq_len(3L), c("read_id", "group_id", "umi"))
      else
        data.table::setnames(chunk, seq_len(2L), c("read_id", "group_id"))
      chunk <- chunk[read_id %in% valid_reads]
      if (nrow(chunk) == 0L) return(NULL)
      chunk[, r_idx   := read_map[read_id]]
      chunk[, read_id := NULL]
      chunk
    }
  )
}

#' Prepare group table for bulk_single (trivial: all reads -> same sample_id)
#' @keywords internal
#' @noRd
.group_from_single <- function(dt_counts, sample_id) {
  data.table::data.table(
    r_idx    = unique(dt_counts$r_idx),
    group_id = sample_id
  )
}

#' Prepare group table via regex extraction from read_id
#' @keywords internal
#' @noRd
.group_from_regex <- function(counts_file, read_map, pattern,
                               bc_pattern, umi_pattern, chunk_size,
                               mode, verbose) {
  valid_reads <- names(read_map)
  is_sc       <- (mode == "sc")

  .read_tsv_chunked(
    path       = counts_file,
    select     = 1L,
    chunk_size = chunk_size,
    FUN = function(chunk) {
      data.table::setnames(chunk, 1L, "read_id")
      # dedup first -- counts has multiple rows per read (multi-mapping)
      chunk <- unique(chunk, by = "read_id")
      chunk <- chunk[read_id %in% valid_reads]
      if (nrow(chunk) == 0L) return(NULL)

      if (is_sc) {
        if (!is.null(pattern)) {
          extracted <- .extract_named_groups(chunk$read_id, pattern,
                                             c("barcode", "umi"))
          data.table::data.table(
            r_idx    = read_map[chunk$read_id],
            group_id = extracted[["barcode"]],
            umi      = extracted[["umi"]]
          )
        } else {
          # separate bc_pattern + umi_pattern
          bc_vals  <- .extract_named_groups(chunk$read_id,
                        paste0("(", bc_pattern,  ")"), "1")[["1"]]
          umi_vals <- .extract_named_groups(chunk$read_id,
                        paste0("(", umi_pattern, ")"), "1")[["1"]]
          # fallback: use gsub if no named groups
          if (all(is.na(bc_vals)))
            bc_vals  <- gsub(paste0(".*", bc_pattern,  ".*"), "\\1",
                             chunk$read_id, perl = TRUE)
          if (all(is.na(umi_vals)))
            umi_vals <- gsub(paste0(".*", umi_pattern, ".*"), "\\1",
                             chunk$read_id, perl = TRUE)
          data.table::data.table(
            r_idx    = read_map[chunk$read_id],
            group_id = bc_vals,
            umi      = umi_vals
          )
        }
      } else {
        # bulk: extract sample_id
        extracted <- .extract_named_groups(chunk$read_id, pattern, "sample_id")
        data.table::data.table(
          r_idx    = read_map[chunk$read_id],
          group_id = extracted[["sample_id"]]
        )
      }
    }
  )
}

#' Core EC construction from integer tables
#'
#' Optimised for memory and speed:
#' - Uses integer obs_id instead of string obs_key (avoids large string alloc)
#' - Fast path for single-mapping observations (no paste, just as.character)
#' - Strategic rm() + gc() to free intermediates
#' @keywords internal
#' @noRd
.build_ec_from_ints <- function(dt_counts, dt_group, tx_map,
                                 mode, unit, gtf_file = NULL,
                                 n_cores, verbose) {
  n_tx <- length(tx_map)

  # join counts with group info on r_idx (both already keyed)
  data.table::setkey(dt_group, r_idx)
  dt_full <- dt_counts[dt_group, on = "r_idx", nomatch = NULL]
  rm(dt_counts, dt_group); gc(verbose = FALSE)

  if (verbose)
    message(sprintf("  Joined: %d records, %d groups.",
                    nrow(dt_full),
                    data.table::uniqueN(dt_full$group_id)))

  # UMI dedup: key = (group_id, umi, t_idx), preserves multi-tx UMIs
  if (mode == "sc" && unit == "umi") {
    n_before <- nrow(dt_full)
    dt_full  <- unique(dt_full, by = c("group_id", "umi", "t_idx"))
    if (verbose)
      message(sprintf("  UMI dedup: %d -> %d records.",
                      n_before, nrow(dt_full)))
    # Integer obs_id (avoids large paste(group_id, umi) string allocation)
    dt_full[, obs_id := .GRP, by = .(group_id, umi)]
  } else {
    dt_full[, obs_id := r_idx]
  }

  all_groups <- unique(dt_full$group_id)
  if (verbose)
    message(sprintf("  Building ECs for %d groups ...", length(all_groups)))

  # ---- Vectorised EC construction ------------------------------------------
  # A: sort t_idx within obs so paste order is canonical
  data.table::setorder(dt_full, group_id, obs_id, t_idx)

  # B: count transcripts per observation (for fast path below)
  obs_n <- dt_full[, .N, by = .(group_id, obs_id)]

  # C: fast path — single-mapping observations (majority)
  #    EC key is just as.character(t_idx), no paste needed.
  single <- obs_n[N == 1L, .(group_id, obs_id)]
  single_ec <- dt_full[single, on = .(group_id, obs_id), nomatch = NULL
    ][, .(group_id, obs_id, ec_key = as.character(t_idx))]

  # D: multi-mapping observations — need paste(collapse)
  multi <- obs_n[N > 1L, .(group_id, obs_id)]
  if (nrow(multi) > 0L) {
    multi_ec <- dt_full[multi, on = .(group_id, obs_id), nomatch = NULL
      ][, .(ec_key = paste(t_idx, collapse = "|")), by = .(group_id, obs_id)]
    obs_ec <- data.table::rbindlist(list(single_ec, multi_ec), use.names = TRUE)
    rm(multi_ec)
  } else {
    obs_ec <- single_ec
  }
  rm(dt_full, single_ec, single, multi, obs_n); gc(verbose = FALSE)

  # E: assign EC id per (group_id, ec_key)
  obs_ec[, ec_id := .GRP, by = .(group_id, ec_key)]

  # F: EC count = number of observations per EC
  ec_counts <- obs_ec[, .(count = .N), by = .(group_id, ec_id, ec_key)]
  rm(obs_ec); gc(verbose = FALSE)

  # G: t_indices — parse ec_key string
  unique_ecs <- unique(ec_counts[, .(group_id, ec_id, ec_key)])
  unique_ecs[, t_indices := lapply(
    strsplit(ec_key, "|", fixed = TRUE), as.integer
  )]

  # H: assemble final EC table
  ec_table <- ec_counts[
    unique_ecs[, .(group_id, ec_id, t_indices)],
    on = c("group_id", "ec_id")
  ]
  data.table::setkey(ec_table, group_id)

  rm(ec_counts, unique_ecs); gc(verbose = FALSE)

  if (verbose)
    message(sprintf("  EC table: %d rows, %d unique ECs total.",
                    nrow(ec_table), data.table::uniqueN(ec_table$ec_id)))

  new_isoem_ec(
    ec_table     = ec_table,
    tx_map       = names(tx_map),
    group_ids    = all_groups,
    n_tx         = n_tx,
    mode         = mode,
    unit         = unit,
    gtf_file     = gtf_file,
    input_params = list(mode = mode, unit = unit, n_tx = n_tx)
  )
}
