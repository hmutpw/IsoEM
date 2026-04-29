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
#'   samples (default 1).
#' @param verbose Logical (default TRUE).
#'
#' @return An \code{IsoEMEC} object.
#'
#' @seealso \code{\link{prepare_isoem}}, \code{\link{run_em}}
#'
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

  tmp_dir <- .setup_dirs(input$outdir)
  on.exit(.cleanup_temp(input$outdir, keep_temp, verbose), add = TRUE)

  if (verbose) message("=== build_ec (bulk mode) ===")

  # ---- Step 1: Read + integerise counts ------------------------------------
  if (verbose) message("[1/3] Reading and integerising counts ...")
  int_data <- .integerise_counts(input$counts_file, chunk_size, verbose)
  # int_data: list(dt = data.table(r_idx|t_idx), read_map, tx_map)

  # ---- Step 2: Prepare group info ------------------------------------------
  if (verbose) message("[2/3] Preparing group information ...")

  dt_group <- if (!is.null(input$sample_id)) {
    # bulk_single: assign fixed sample_id to all reads
    .group_from_single(int_data$dt, input$sample_id)

  } else if (!is.null(input$anno_file)) {
    # bulk_multi: shell sort+join anno_file with counts
    merged_f <- file.path(tmp_dir, "merged_counts.tsv.gz")
    .group_from_join(
      counts_file = input$counts_file,
      anno_file   = input$anno_file,
      read_map    = int_data$read_map,
      outfile     = merged_f,
      mode        = "bulk_multi",
      chunk_size  = chunk_size,
      verbose     = verbose
    )

  } else {
    # regex extraction
    .group_from_regex(
      counts_file    = input$counts_file,
      read_map       = int_data$read_map,
      pattern        = input$pattern,
      chunk_size     = chunk_size,
      mode           = "bulk",
      verbose        = verbose
    )
  }
  # dt_group: data.table(r_idx, group_id)

  # ---- Step 3: Build EC ----------------------------------------------------
  if (verbose) message("[3/3] Building equivalence classes ...")
  ec <- .build_ec_from_ints(
    dt_counts  = int_data$dt,
    dt_group   = dt_group,
    tx_map     = int_data$tx_map,
    mode       = "bulk",
    unit       = "read",
    n_cores    = n_cores,
    verbose    = verbose
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
#' @param verbose Logical (default TRUE).
#'
#' @return An \code{IsoEMEC} object with \code{mode = "sc"}.
#'
#' @seealso \code{\link{prepare_isoem}}, \code{\link{run_em}}
#'
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

  tmp_dir <- .setup_dirs(input$outdir)
  on.exit(.cleanup_temp(input$outdir, keep_temp, verbose), add = TRUE)

  if (verbose) message("=== build_sc_ec (single-cell mode) ===")
  if (verbose) message("  Unit: ", input$unit)

  # ---- Step 1: Read + integerise counts ------------------------------------
  if (verbose) message("[1/4] Reading and integerising counts ...")
  int_data <- .integerise_counts(input$counts_file, chunk_size, verbose)

  # ---- Step 2: Prepare barcode + UMI info ----------------------------------
  if (verbose) message("[2/4] Preparing barcode/UMI information ...")

  dt_group <- if (!is.null(input$anno_file)) {
    merged_f <- file.path(tmp_dir, "merged_sc_counts.tsv.gz")
    .group_from_join(
      counts_file = input$counts_file,
      anno_file   = input$anno_file,
      read_map    = int_data$read_map,
      outfile     = merged_f,
      mode        = "sc",
      chunk_size  = chunk_size,
      verbose     = verbose
    )
  } else {
    .group_from_regex(
      counts_file = input$counts_file,
      read_map    = int_data$read_map,
      pattern     = input$pattern,
      bc_pattern  = input$bc_pattern,
      umi_pattern = input$umi_pattern,
      chunk_size  = chunk_size,
      mode        = "sc",
      verbose     = verbose
    )
  }
  # dt_group: data.table(r_idx, group_id [=barcode], umi)

  # ---- Step 3: Barcode whitelist filter ------------------------------------
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

  # ---- Step 4: Build EC ----------------------------------------------------
  if (verbose) message("[4/4] Building equivalence classes ...")
  ec <- .build_ec_from_ints(
    dt_counts  = int_data$dt,
    dt_group   = dt_group,
    tx_map     = int_data$tx_map,
    mode       = "sc",
    unit       = input$unit,
    n_cores    = n_cores,
    verbose    = verbose
  )

  rm(int_data, dt_group); gc(verbose = FALSE)
  if (verbose) message("EC construction complete.")
  ec
}


# =============================================================================
# Internal workers
# =============================================================================

#' Read counts file and build integer maps for read_id and transcript_id
#' @keywords internal
#' @noRd
.integerise_counts <- function(counts_file, chunk_size, verbose) {
  # First pass: build string -> int maps via chunk reading
  read_map <- character(0)   # will grow
  tx_map   <- character(0)

  # Read full file to build global maps
  # (counts file is typically smaller than bc_umi file)
  dt <- .read_tsv(counts_file,
    col.names  = c("read_id", "transcript_id"),
    select     = c(1L, 2L),
    colClasses = c("character", "character")
  )
  dt <- dt[transcript_id != "*"]

  if (nrow(dt) == 0L)
    stop("No valid records in counts_file after filtering '*'.", call. = FALSE)

  if (verbose)
    message(sprintf("  counts: %d records, %d unique reads, %d unique transcripts.",
                    nrow(dt),
                    data.table::uniqueN(dt$read_id),
                    data.table::uniqueN(dt$transcript_id)))

  # build maps
  all_reads <- unique(dt$read_id)
  all_tx    <- sort(unique(dt$transcript_id))
  read_map  <- stats::setNames(seq_along(all_reads), all_reads)
  tx_map    <- stats::setNames(seq_along(all_tx),    all_tx)

  # convert to integer table
  dt[, r_idx := read_map[read_id]]
  dt[, t_idx := tx_map[transcript_id]]
  dt[, c("read_id", "transcript_id") := NULL]

  list(
    dt       = dt,          # data.table(r_idx, t_idx)
    read_map = read_map,    # named int: read_id -> r_idx
    tx_map   = tx_map       # named int: transcript_id -> t_idx
  )
}

#' Prepare group table via shell sort+join
#' Returns data.table(r_idx, group_id) for bulk
#'                  or data.table(r_idx, group_id, umi) for sc
#' @keywords internal
#' @noRd
.group_from_join <- function(counts_file, anno_file, read_map,
                              outfile, mode, chunk_size, verbose) {
  # shell sort+join: counts col1 (read_id) joined with anno col1 (read_id)
  .shell_sort_join(counts_file, anno_file, outfile,
                   compress_out = TRUE, verbose = verbose)

  # read result in chunks, immediately map read_id -> r_idx
  col_names  <- if (mode == "sc") c("read_id", "transcript_id", "barcode", "umi")
                else               c("read_id", "transcript_id", "sample_id")
  select_idx <- if (mode == "sc") c(1L, 3L, 4L) else c(1L, 3L)
  col_out    <- if (mode == "sc") c("read_id", "group_id", "umi")
                else               c("read_id", "group_id")

  valid_reads <- names(read_map)

  result <- .read_tsv_chunked(
    path       = outfile,
    col.names  = NULL,        # do NOT pass col.names when using select
    select     = select_idx,
    colClasses = NULL,        # do NOT pass colClasses with select (fread applies to all cols first)
    chunk_size = chunk_size,
    FUN = function(chunk) {
      # rename selected columns to col_out names
      data.table::setnames(chunk, seq_along(col_out), col_out)
      # filter to valid read_ids and convert to integer
      chunk <- chunk[read_id %in% valid_reads]
      chunk[, r_idx := read_map[read_id]]
      chunk[, read_id := NULL]
      chunk
    }
  )
  result
}

#' Prepare group table for bulk_single (no join needed)
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
                               bc_pattern = NULL, umi_pattern = NULL,
                               chunk_size, mode, verbose) {
  valid_reads <- names(read_map)

  result <- .read_tsv_chunked(
    path       = counts_file,
    col.names  = NULL,   # select=1L picks first col; rename in FUN
    select     = 1L,
    colClasses = NULL,
    chunk_size = chunk_size,
    FUN = function(chunk) {
      data.table::setnames(chunk, 1L, "read_id")
      # deduplicate read_ids first (counts has multiple rows per read due to multi-mapping)
      chunk <- unique(chunk, by = "read_id")
      # only process reads that have a transcript assignment
      chunk <- chunk[read_id %in% valid_reads]
      if (nrow(chunk) == 0L) return(NULL)

      if (mode == "sc") {
        if (!is.null(pattern)) {
          # single pattern with named groups
          bc  <- regmatches(chunk$read_id,
            regexpr(sub(".*\\(\\?P<barcode>", "(?P<barcode>", pattern, perl=TRUE),
                    chunk$read_id, perl = TRUE))
          umi <- regmatches(chunk$read_id,
            regexpr(sub(".*\\(\\?P<umi>", "(?P<umi>", pattern, perl=TRUE),
                    chunk$read_id, perl = TRUE))
          # extract values from matched strings
          bc_val  <- sub(".*(?P<barcode>).*", "\\1",
            regmatches(chunk$read_id,
              regexpr(pattern, chunk$read_id, perl=TRUE)), perl=TRUE)
          # simpler: use dedicated helper
          extracted <- .extract_named_groups(chunk$read_id, pattern,
                                              c("barcode", "umi"))
          data.table::data.table(
            r_idx    = read_map[chunk$read_id],
            group_id = extracted[["barcode"]],
            umi      = extracted[["umi"]]
          )
        } else {
          # separate patterns
          bc_val  <- gsub(paste0(".*", bc_pattern,  ".*"), "\\1",
                          chunk$read_id, perl = TRUE)
          umi_val <- gsub(paste0(".*", umi_pattern, ".*"), "\\1",
                          chunk$read_id, perl = TRUE)
          data.table::data.table(
            r_idx    = read_map[chunk$read_id],
            group_id = bc_val,
            umi      = umi_val
          )
        }
      } else {
        # bulk: extract sample_id
        sid_val <- .extract_named_groups(chunk$read_id, pattern, "sample_id")
        data.table::data.table(
          r_idx    = read_map[chunk$read_id],
          group_id = sid_val[["sample_id"]]
        )
      }
    }
  )
  result
}

#' Extract named capture groups from read_id vector using a regex
#' @keywords internal
#' @noRd
.extract_named_groups <- function(x, pattern, groups) {
  result <- lapply(groups, function(grp) {
    pat <- paste0("(?<=(?P<", grp, ">))[^_]+")
    # use PCRE named backreference approach
    m <- regmatches(x, regexpr(
      paste0("(?P<", grp, ">[^_]+)"),
      x, perl = TRUE
    ))
    # fallback: just run the full pattern and extract group
    vals <- rep(NA_character_, length(x))
    matched <- regexpr(pattern, x, perl = TRUE)
    hit <- matched > 0
    if (any(hit)) {
      starts <- attr(matched, "capture.start")
      lengths <- attr(matched, "capture.length")
      names_cap <- attr(matched, "capture.names")
      grp_idx <- which(names_cap == grp)
      if (length(grp_idx) > 0L && any(hit)) {
        vals[hit] <- substr(x[hit],
                            starts[hit, grp_idx],
                            starts[hit, grp_idx] + lengths[hit, grp_idx] - 1L)
      }
    }
    vals
  })
  names(result) <- groups
  result
}

#' Build IsoEMEC from integer dt_counts + dt_group
#' @keywords internal
#' @noRd
.build_ec_from_ints <- function(dt_counts, dt_group, tx_map,
                                 mode, unit, n_cores, verbose) {
  n_tx <- length(tx_map)

  # join counts with group info on r_idx
  data.table::setkey(dt_counts, r_idx)
  data.table::setkey(dt_group,  r_idx)
  dt_full <- dt_counts[dt_group, on = "r_idx", nomatch = NULL]
  rm(dt_counts, dt_group); gc(verbose = FALSE)

  if (verbose)
    message(sprintf("  Joined: %d records, %d groups.",
                    nrow(dt_full),
                    data.table::uniqueN(dt_full$group_id)))

  # UMI dedup for sc/umi mode
  if (mode == "sc" && unit == "umi") {
    n_before <- nrow(dt_full)
    # dedup key: (group_id, umi, t_idx) — preserves multi-transcript UMIs
    dt_full <- unique(dt_full, by = c("group_id", "umi", "t_idx"))
    if (verbose)
      message(sprintf("  UMI dedup: %d -> %d records.", n_before, nrow(dt_full)))
    # replace r_idx with umi_key for EC construction
    dt_full[, obs_key := paste(group_id, umi, sep = "__")]
  } else {
    dt_full[, obs_key := as.character(r_idx)]
  }

  # ---- EC construction per group ------------------------------------------
  all_groups <- unique(dt_full$group_id)

  if (verbose)
    message(sprintf("  Building ECs for %d groups ...", length(all_groups)))

  # split by group
  dt_split <- split(dt_full[, .(obs_key, t_idx)], dt_full$group_id)
  rm(dt_full); gc(verbose = FALSE)

  # batch parallel: divide groups into n_cores batches
  build_one_group <- function(gdt) {
    # gdt: data.table(obs_key, t_idx)
    data.table::setorder(gdt, obs_key, t_idx)
    # build EC key per obs_key
    obs_ec <- gdt[, .(
      ec_key = paste(sort(t_idx), collapse = "|")
    ), by = obs_key]
    obs_ec[, ec_id := .GRP, by = ec_key]
    # EC counts
    ec_counts <- obs_ec[, .(count = .N), by = .(ec_id, ec_key)]
    # EC transcript membership as list column
    ec_tx <- gdt[obs_ec[, .(obs_key, ec_id)], on = "obs_key"][,
      .(t_indices = list(unique(t_idx))), by = ec_id
    ]
    ec_counts[ec_tx, on = "ec_id"]
  }

  if (n_cores > 1L && .Platform$OS.type != "windows") {
    batches <- split(names(dt_split),
                     cut(seq_along(dt_split), n_cores, labels = FALSE))
    ec_list_batched <- parallel::mclapply(batches, function(batch) {
      lapply(stats::setNames(batch, batch), function(grp) {
        res <- build_one_group(dt_split[[grp]])
        res[, group_id := grp]
        res
      })
    }, mc.cores = n_cores)
    ec_list <- unlist(ec_list_batched, recursive = FALSE)
  } else {
    ec_list <- lapply(stats::setNames(all_groups, all_groups), function(grp) {
      res <- build_one_group(dt_split[[grp]])
      res[, group_id := grp]
      res
    })
  }

  ec_table <- data.table::rbindlist(ec_list, use.names = TRUE)
  rm(ec_list, dt_split); gc(verbose = FALSE)

  # reassign global ec_id (per-group ec_ids may collide)
  ec_table[, ec_id := .GRP, by = .(group_id, ec_key)]
  ec_table[, ec_key := NULL]

  if (verbose)
    message(sprintf("  EC table: %d rows, %d unique ECs total.",
                    nrow(ec_table),
                    data.table::uniqueN(ec_table$ec_id)))

  new_isoem_ec(
    ec_table     = ec_table,
    tx_map       = names(tx_map),   # index -> transcript_id (1-based)
    group_ids    = all_groups,
    n_tx         = n_tx,
    mode         = mode,
    unit         = unit,
    input_params = list(mode = mode, unit = unit, n_tx = n_tx)
  )
}
