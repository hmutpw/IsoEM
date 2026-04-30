# =============================================================================
# build_ec.R  —  EC construction for bulk and single-cell modes
# =============================================================================

# -----------------------------------------------------------------------------
# build_ec()  —  Bulk
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
#'   samples (default 1). Uses \code{parallel::makeCluster} —
#'   compatible with Windows and Linux.
#' @param verbose Logical (default TRUE).
#'
#' @return An \code{IsoEMEC} object.
#' @seealso \code{\link{prepare_isoem}}, \code{\link{run_em}}
#' @export
#'
#' @examples
#' counts_f <- system.file("extdata", "toy_counts.tsv",  package = "IsoEM")
#' gtf_f    <- system.file("extdata", "toy.gtf",         package = "IsoEM")
#' input    <- prepare_isoem(counts_f, gtf_f,
#'   mode = "bulk_single", sample_id = "toy")
#' ec <- build_ec(input)
#' print(ec)
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
# build_sc_ec()  —  Single-cell
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
#'   Uses \code{parallel::makeCluster} — compatible with Windows and Linux.
#' @param verbose Logical (default TRUE).
#'
#' @return An \code{IsoEMEC} object with \code{mode = "sc"}.
#' @seealso \code{\link{prepare_isoem}}, \code{\link{run_em}}
#' @export
#'
#' @examples
#' counts_f <- system.file("extdata", "toy_counts.tsv",  package = "IsoEM")
#' gtf_f    <- system.file("extdata", "toy.gtf",         package = "IsoEM")
#' bc_f     <- system.file("extdata", "toy_bc_umi.tsv",  package = "IsoEM")
#' input    <- prepare_isoem(counts_f, gtf_f,
#'   mode = "sc", anno_file = bc_f, unit = "umi")
#' ec <- build_sc_ec(input)
#' print(ec)
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
    if (verbose) message("[3/4] No whitelist — using all barcodes.")
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
  read_map  <- stats::setNames(seq_along(all_reads), all_reads)
  tx_map    <- stats::setNames(seq_along(all_tx),    all_tx)

  dt[, r_idx := read_map[read_id]]
  dt[, t_idx := tx_map[transcript_id]]
  dt[, c("read_id", "transcript_id") := NULL]
  data.table::setkey(dt, r_idx)

  list(dt = dt, read_map = read_map, tx_map = tx_map)
}

#' Prepare group table from an external annotation file (in-memory chunked join)
#' sc:         anno_file cols = read_id | barcode | umi
#' bulk_multi: anno_file cols = read_id | sample_id
#' Returns data.table(r_idx, group_id [, umi])
#' @keywords internal
#' @noRd
.group_from_anno <- function(anno_file, read_map, mode, chunk_size, verbose) {
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

#' Prepare group table for bulk_single (trivial: all reads → same sample_id)
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
      # dedup first — counts has multiple rows per read (multi-mapping)
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
#' Fully vectorised data.table operations; parallel across groups.
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
    # obs_key = unique observation unit per cell
    dt_full[, obs_key := paste(group_id, umi, sep = "__")]
  } else {
    dt_full[, obs_key := as.character(r_idx)]
  }

  all_groups <- unique(dt_full$group_id)
  if (verbose)
    message(sprintf("  Building ECs for %d groups ...", length(all_groups)))

  # ---- Vectorised EC construction (no per-group loop) ----------------------
  # A: sort t_idx within obs so paste order is canonical
  data.table::setorder(dt_full, group_id, obs_key, t_idx)

  # B: EC key per (group_id, obs_key) — one vectorised groupby
  obs_ec <- dt_full[,
    .(ec_key = paste(t_idx, collapse = "|")),
    by = .(group_id, obs_key)
  ]

  # C: global ec_id per (group_id, ec_key)
  obs_ec[, ec_id := .GRP, by = .(group_id, ec_key)]

  # D: EC count = number of observations per EC
  ec_counts <- obs_ec[, .(count = .N), by = .(group_id, ec_id, ec_key)]

  # E: t_indices — parse ec_key string (avoids re-joining dt_full)
  unique_ecs <- unique(ec_counts[, .(group_id, ec_id, ec_key)])
  unique_ecs[, t_indices := lapply(
    strsplit(ec_key, "|", fixed = TRUE), as.integer
  )]

  # F: assemble final EC table
  ec_table <- ec_counts[
    unique_ecs[, .(group_id, ec_id, t_indices)],
    on = c("group_id", "ec_id")
  ]
  data.table::setkey(ec_table, group_id)  # pre-key for fast run_em lookups

  rm(dt_full, obs_ec, ec_counts, unique_ecs); gc(verbose = FALSE)

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
