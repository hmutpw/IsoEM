# =============================================================================
# prepare_isoem.R  —  Input validation + IsoEMInput object creation
# =============================================================================

#' Prepare and validate IsoEM input files
#'
#' Creates a lightweight \code{IsoEMInput} object that stores file paths and
#' validated metadata without loading any data into memory. This object is
#' passed to \code{\link{build_ec}} or \code{\link{build_sc_ec}} to begin
#' quantification.
#'
#' @param counts_file Character. Path to IsoQuant
#'   \code{transcript_model_reads.tsv[.gz]}. Two columns, no header:
#'   \code{read_id | transcript_id}.
#' @param gtf_file Character. Path to IsoQuant
#'   \code{transcript_models.gtf[.gz]}.
#' @param mode Character. One of \code{"bulk_single"}, \code{"bulk_multi"},
#'   or \code{"sc"}.
#' @param anno_file Character or NULL. Path to annotation file.
#'   \itemize{
#'     \item \code{bulk_multi}: \code{read_id | sample_id} (2 cols)
#'     \item \code{sc}: \code{read_id | barcode | umi} (3 cols)
#'   }
#'   Cannot be set together with \code{sample_id} or \code{pattern}.
#' @param sample_id Character or NULL. For \code{bulk_single} only.
#'   Fixed sample ID assigned to all reads. Cannot be set with
#'   \code{anno_file} or \code{pattern}.
#' @param pattern Character or NULL. Single regex with named capture groups
#'   to extract group info from read_id. Example:
#'   \code{".*_bc=(?P<barcode>[ACGT]+).*_sq=(?P<umi>[ACGT]+)"}.
#'   Cannot be set with \code{anno_file} or \code{sample_id}.
#' @param bc_pattern Character or NULL. Separate barcode regex (used with
#'   \code{umi_pattern} when \code{pattern} not provided).
#' @param umi_pattern Character or NULL. Separate UMI regex.
#' @param barcodes_use Character vector or file path or NULL.
#'   Cell barcode whitelist for \code{sc} mode.
#' @param unit Character. \code{"umi"} (default for sc) or \code{"read"}.
#'   Controls whether EC construction uses UMI-level or read-level counts.
#' @param outdir Character. Output directory. A \code{temp/} subdirectory
#'   will be created here for intermediate files. Default: \code{"IsoEM_out"}.
#' @param verbose Logical (default TRUE).
#'
#' @return An \code{IsoEMInput} object.
#'
#' @seealso \code{\link{build_ec}}, \code{\link{build_sc_ec}}
#'
#' @export
#'
#' @examples
#' counts_f <- system.file("extdata", "toy_counts.tsv.gz",  package = "IsoEM")
#' gtf_f    <- system.file("extdata", "toy.gtf.gz",         package = "IsoEM")
#'
#' # bulk single-sample
#' input <- prepare_isoem(counts_f, gtf_f,
#'   mode = "bulk_single", sample_id = "toy_sample")
#' print(input)
#'
#' # single-cell with anno_file
#' counts_sc_f <- system.file("extdata", "toy_counts_multi.tsv", package = "IsoEM")
#' bc_f        <- system.file("extdata", "toy_bc_umi.tsv",       package = "IsoEM")
#' input_sc    <- prepare_isoem(counts_sc_f, gtf_f,
#'   mode = "sc", anno_file = bc_f, unit = "umi")
#' print(input_sc)
prepare_isoem <- function(
    counts_file,
    gtf_file,
    mode           = c("bulk_single", "bulk_multi", "sc"),
    anno_file      = NULL,
    sample_id      = NULL,
    pattern        = NULL,
    bc_pattern     = NULL,
    umi_pattern    = NULL,
    barcodes_use   = NULL,
    unit           = c("umi", "read"),
    outdir         = "IsoEM_out",
    verbose        = TRUE
) {
  mode <- match.arg(mode)
  unit <- match.arg(unit)

  if (verbose) message("=== IsoEM input preparation ===")
  if (verbose) message("  Mode: ", mode)

  # ---- 1. Validate required files ------------------------------------------
  .check_file(counts_file, "counts_file")
  .check_file(gtf_file,    "gtf_file")

  # ---- 2. Validate group info sources (mutually exclusive) -----------------
  n_sources <- sum(
    !is.null(anno_file),
    !is.null(sample_id),
    !is.null(pattern) || (!is.null(bc_pattern) || !is.null(umi_pattern))
  )

  if (mode == "bulk_single") {
    if (n_sources == 0)
      stop("bulk_single mode requires 'sample_id' or 'pattern'.", call. = FALSE)
    if (!is.null(anno_file))
      stop("bulk_single mode does not use 'anno_file'. Use 'sample_id'.",
           call. = FALSE)
  } else {
    # bulk_multi / sc
    if (n_sources == 0)
      stop(sprintf("%s mode requires 'anno_file' or 'pattern'.", mode),
           call. = FALSE)
  }
  if (n_sources > 1)
    stop("Only one of 'anno_file', 'sample_id', or 'pattern' can be provided.",
         call. = FALSE)

  # ---- 3. Validate anno_file if provided ------------------------------------
  if (!is.null(anno_file)) {
    .check_file(anno_file, "anno_file")
    # peek at first few lines to verify column count
    peek <- data.table::fread(anno_file, header=FALSE, nrows=5L,
                              showProgress=FALSE, data.table=TRUE)
    expected_cols <- if (mode == "sc") 3L else 2L
    if (ncol(peek) < expected_cols)
      stop(sprintf(
        "anno_file for mode '%s' must have >= %d columns (read_id | %s).",
        mode, expected_cols,
        if (mode == "sc") "barcode | umi" else "sample_id"),
        call. = FALSE)
  }

  # ---- 4. Validate regex patterns ------------------------------------------
  if (!is.null(pattern)) {
    cap_names <- .get_capture_names(pattern)
    if (length(cap_names) == 0L)
      stop("'pattern' must contain named capture groups, e.g. (?P<barcode>...)",
           call. = FALSE)
    if (mode == "sc" && !all(c("barcode", "umi") %in% cap_names))
      stop("For sc mode, 'pattern' must capture groups named 'barcode' and 'umi'.",
           call. = FALSE)
    if (mode != "sc" && !"sample_id" %in% cap_names)
      stop("For bulk mode, 'pattern' must capture a group named 'sample_id'.",
           call. = FALSE)
    # test on a few real read_ids
    peek <- .read_tsv(counts_file,
      col.names  = NULL,    # avoid col.names+select length mismatch
      select     = 1L,
      colClasses = NULL,
      nrows      = 20L
    )
    data.table::setnames(peek, 1L, "read_id")
    n_match <- sum(grepl(pattern, peek$read_id, perl = TRUE))
    if (n_match == 0L)
      warning("'pattern' matched 0 / ", nrow(peek),
              " sampled read_ids. Check your regex.", call. = FALSE)
    else if (verbose)
      message(sprintf("  Pattern test: matched %d / %d sampled read_ids.",
                      n_match, nrow(peek)))
  }

  if (!is.null(bc_pattern) || !is.null(umi_pattern)) {
    if (mode != "sc")
      stop("bc_pattern / umi_pattern are only for sc mode.", call. = FALSE)
    if (is.null(bc_pattern) || is.null(umi_pattern))
      stop("Both 'bc_pattern' and 'umi_pattern' must be provided together.",
           call. = FALSE)
  }

  # ---- 5. Validate barcodes_use --------------------------------------------
  barcodes_n <- NULL
  if (!is.null(barcodes_use)) {
    if (mode != "sc")
      warning("'barcodes_use' is only used in sc mode.", call. = FALSE)
    if (length(barcodes_use) == 1L && file.exists(barcodes_use)) {
      # it's a file path — validate
      .check_file(barcodes_use, "barcodes_use")
      bc_lines <- readLines(barcodes_use)
      barcodes_n <- length(bc_lines[nchar(trimws(bc_lines)) > 0])
    } else {
      barcodes_n <- length(barcodes_use)
    }
    if (verbose)
      message(sprintf("  Whitelist: %d barcodes.", barcodes_n))
  }

  # ---- 6. Validate counts_file format --------------------------------------
  counts_peek <- .read_tsv(counts_file,
    col.names  = NULL,
    select     = c(1L, 2L),
    colClasses = NULL,
    nrows      = 5L
  )
  if (ncol(counts_peek) < 2L)
    stop("counts_file must have at least 2 columns (read_id | transcript_id).",
         call. = FALSE)
  if (verbose) message("  counts_file format: OK")

  # ---- 7. Collect file info ------------------------------------------------
  counts_nrow <- .estimate_rows(counts_file)
  anno_nrow   <- if (!is.null(anno_file)) .estimate_rows(anno_file) else NULL

  # rough memory estimate: ~50 bytes/row for integer table after processing
  mem_est <- if (!is.na(counts_nrow))
    round(counts_nrow * 50 / 1e9, 1) else NA

  file_info <- list(
    counts = list(
      size_str = .size_str(counts_file),
      nrow_str = .nrow_str(counts_nrow)
    ),
    anno = if (!is.null(anno_file)) list(
      size_str = .size_str(anno_file),
      nrow_str = .nrow_str(anno_nrow)
    ) else NULL,
    barcodes_n  = barcodes_n,
    mem_est_gb  = mem_est
  )

  # ---- 8. Setup outdir -----------------------------------------------------
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  if (verbose) message("  outdir: ", outdir)

  # ---- 9. Build object -----------------------------------------------------
  input <- new_isoem_input(
    mode         = mode,
    counts_file  = counts_file,
    gtf_file     = gtf_file,
    anno_file    = anno_file,
    sample_id    = sample_id,
    pattern      = pattern,
    bc_pattern   = bc_pattern,
    umi_pattern  = umi_pattern,
    barcodes_use = barcodes_use,
    unit         = unit,
    outdir       = outdir,
    file_info    = file_info,
    validated    = TRUE
  )

  if (verbose) message("  Validation complete.")
  input
}
