# =============================================================================
# run_isoem.R  —  One-step wrappers
# =============================================================================

#' One-step bulk EM quantification
#'
#' Convenience wrapper that calls \code{\link{prepare_isoem}},
#' \code{\link{build_ec}}, and \code{\link{run_em}} in sequence.
#' For a modular workflow that allows saving the intermediate \code{IsoEMEC}
#' object, call those functions separately.
#'
#' @param counts_file Character. Path to IsoQuant
#'   \code{transcript_model_reads.tsv[.gz]}.
#'   Each row: \code{read_id [TAB] transcript_id}.
#' @param gtf_file Character. Path to IsoQuant \code{transcript_models.gtf[.gz]}.
#'   Used for gene annotation in the output.
#' @param mode Character. One of \code{"bulk_single"} (one sample per run)
#'   or \code{"bulk_multi"} (multiple samples, requires \code{anno_file}).
#' @param anno_file Character or NULL. For \code{"bulk_multi"}: path to a
#'   two-column TSV (\code{read_id | sample_id}), no header. For
#'   \code{"bulk_single"} with \code{sample_id = NULL}: path to a file
#'   whose first column is \code{read_id} (sample extracted via \code{pattern}).
#' @param sample_id Character or NULL. Sample name for \code{"bulk_single"} mode.
#'   If NULL, \code{pattern} must be supplied to extract it from read IDs.
#' @param pattern Character or NULL. A PCRE regex with a named capture group
#'   \code{(?P<sample_id>...)} to extract sample IDs from read names.
#'   Only used when \code{anno_file} is NULL.
#' @param outdir Character. Output directory. A \code{temp/} subdirectory is
#'   used for intermediate files and cleaned up unless \code{keep_temp = TRUE}.
#'   Default: \code{"IsoEM_out"}.
#' @param keep_temp Logical. Keep the \code{temp/} directory after EC
#'   construction (default FALSE).
#' @param chunk_size Integer. Rows per chunk when reading large annotation or
#'   counts files (default 5e6). Reduce to lower peak memory.
#' @param max_iter Integer. Maximum EM iterations (default 1000).
#' @param tol Numeric. EM convergence tolerance (default 1e-6). EM stops when
#'   the L1 change in the count vector normalised by total counts is below
#'   this threshold.
#' @param n_cores Integer. Parallel workers for EC construction and EM
#'   (default 1). Uses \code{parallel::makeCluster} internally.
#' @param verbose Logical. Print progress messages (default TRUE).
#' @param write Logical. Automatically write results to \code{outdir}
#'   after EM completes (default TRUE). Calls \code{\link{write_isoem}}
#'   with the parameters below. Set to FALSE to return the result object
#'   without writing, e.g. for interactive inspection before saving.
#' @param compress Logical. Write \code{.gz} compressed output files
#'   (default FALSE). Passed to \code{\link{write_isoem}}.
#' @param write_ec_table Logical. Write \code{ec_table.tsv} (default TRUE).
#'   Passed to \code{\link{write_isoem}}.
#' @param write_sharing_table Logical. Write \code{sharing_table.tsv}
#'   (default TRUE). Passed to \code{\link{write_isoem}}.
#' @param min_sharing Numeric. Minimum sharing fraction for sharing table
#'   (default 0.5). Passed to \code{\link{write_isoem}}.
#' @param max_ec_size_sharing Integer. Max EC size for pairwise expansion
#'   (default 50). Passed to \code{\link{write_isoem}}.
#'
#' @return An \code{IsoEMResult} (single-sample) or \code{IsoEMDataset}
#'   (multi-sample), invisibly if \code{write = TRUE}.
#'
#' @seealso \code{\link{prepare_isoem}}, \code{\link{build_ec}},
#'   \code{\link{run_em}}, \code{\link{write_isoem}}
#'
#' @export
#'
#' @examples
#' counts_f <- system.file("extdata", "toy_counts.tsv.gz", package = "IsoEM")
#' gtf_f    <- system.file("extdata", "toy.gtf.gz",        package = "IsoEM")
#'
#' # Single sample — write automatically
#' result <- run_isoem(counts_f, gtf_f,
#'   mode = "bulk_single", sample_id = "toy_sample",
#'   outdir = tempdir())
#'
#' # Inspect before writing
#' result <- run_isoem(counts_f, gtf_f,
#'   mode = "bulk_single", sample_id = "toy_sample",
#'   write = FALSE)
#' print(result)
#' write_isoem(result, outdir = tempdir())
#'
#' # Multi-sample
#' multi_f <- system.file("extdata", "toy_counts_multi.tsv", package = "IsoEM")
#' sid_f   <- system.file("extdata", "toy_sample_ids.tsv",   package = "IsoEM")
#' result  <- run_isoem(multi_f, gtf_f,
#'   mode = "bulk_multi", anno_file = sid_f, outdir = tempdir())
#' print(result)
run_isoem <- function(
    counts_file,
    gtf_file,
    mode                = c("bulk_single", "bulk_multi"),
    anno_file           = NULL,
    sample_id           = NULL,
    pattern             = NULL,
    outdir              = "IsoEM_out",
    keep_temp           = FALSE,
    chunk_size          = 5e6L,
    max_iter            = 1000L,
    tol                 = 1e-6,
    n_cores             = 1L,
    verbose             = TRUE,
    write               = TRUE,
    compress            = FALSE,
    write_ec_table      = TRUE,
    write_sharing_table = TRUE,
    min_sharing         = 0.5,
    max_ec_size_sharing = 50L
) {
  mode <- match.arg(mode)

  input <- prepare_isoem(
    counts_file = counts_file,
    gtf_file    = gtf_file,
    mode        = mode,
    anno_file   = anno_file,
    sample_id   = sample_id,
    pattern     = pattern,
    outdir      = outdir,
    verbose     = verbose
  )

  ec <- build_ec(
    input      = input,
    keep_temp  = keep_temp,
    chunk_size = chunk_size,
    n_cores    = n_cores,
    verbose    = verbose
  )

  result <- run_em(
    ec       = ec,
    max_iter = max_iter,
    tol      = tol,
    n_cores  = n_cores,
    verbose  = verbose
  )

  if (write) {
    write_isoem(
      result,
      outdir              = outdir,
      compress            = compress,
      write_ec_table      = write_ec_table,
      write_sharing_table = write_sharing_table,
      min_sharing         = min_sharing,
      max_ec_size_sharing = as.integer(max_ec_size_sharing),
      verbose             = verbose
    )
    return(invisible(result))
  }
  result
}


#' One-step single-cell EM quantification
#'
#' Convenience wrapper that calls \code{\link{prepare_isoem}},
#' \code{\link{build_sc_ec}}, and \code{\link{run_em}} in sequence.
#' For a modular workflow that allows saving the intermediate \code{IsoEMEC}
#' object (recommended for large datasets), call those functions separately.
#'
#' @param counts_file Character. Path to IsoQuant
#'   \code{transcript_model_reads.tsv[.gz]}.
#'   Each row: \code{read_id [TAB] transcript_id}. Reads mapping to \code{*}
#'   are discarded automatically.
#' @param gtf_file Character. Path to IsoQuant \code{transcript_models.gtf[.gz]}.
#'   Used for gene and novelty annotation in the output.
#' @param anno_file Character or NULL. Path to a three-column TSV
#'   (\code{read_id | barcode | umi}), no header, plain or gzip-compressed.
#'   Mutually exclusive with \code{pattern} / \code{bc_pattern} +
#'   \code{umi_pattern}. One of \code{anno_file} or \code{pattern} is
#'   required.
#' @param pattern Character or NULL. A single PCRE regex with two named
#'   capture groups \code{(?P<barcode>...)} and \code{(?P<umi>...)} to
#'   extract barcode and UMI directly from read IDs. Used when
#'   \code{anno_file} is NULL.
#' @param bc_pattern Character or NULL. Separate barcode regex. Used
#'   together with \code{umi_pattern} as an alternative to \code{pattern}.
#' @param umi_pattern Character or NULL. Separate UMI regex. Used together
#'   with \code{bc_pattern}.
#' @param barcodes_use Character vector, file path, or NULL. Cell barcode
#'   whitelist. If a file path, the file must contain one barcode per line
#'   (plain or gzip). Reads whose barcodes are not in the whitelist are
#'   discarded before EC construction (default NULL = use all barcodes).
#' @param unit Character. Observational unit for EC construction:
#'   \code{"umi"} (default, recommended) or \code{"read"}.
#'   In UMI mode, each \code{(barcode, UMI)} pair is a single observation;
#'   PCR duplicates are collapsed. In read mode, each aligned read is an
#'   independent observation.
#' @param outdir Character. Output directory (default \code{"IsoEM_out"}).
#' @param keep_temp Logical. Keep intermediate temp files (default FALSE).
#' @param chunk_size Integer. Rows per chunk for reading large files
#'   (default 5e6). Reduce to lower peak memory usage.
#' @param max_iter Integer. Maximum EM iterations per cell (default 500).
#' @param tol Numeric. EM convergence tolerance (default 1e-6).
#' @param n_cores Integer. Parallel workers (default 1). Parallelisation is
#'   applied across cells in the EM step and across groups in EC construction.
#' @param verbose Logical. Print progress messages (default TRUE).
#' @param write Logical. Automatically write results to \code{outdir}
#'   after EM completes (default TRUE). Calls \code{\link{write_sc_isoem}}
#'   with the parameters below. Set to FALSE to return the result object
#'   without writing.
#' @param compress Logical. Write \code{.gz} compressed output files
#'   (default FALSE). Passed to \code{\link{write_sc_isoem}}.
#' @param write_ec_table Logical. Write \code{ec_table.tsv} (default TRUE).
#'   Passed to \code{\link{write_sc_isoem}}.
#' @param write_sharing_table Logical. Write \code{sharing_table.tsv}
#'   (default TRUE). Passed to \code{\link{write_sc_isoem}}.
#' @param min_sharing Numeric. Minimum sharing fraction for sharing table
#'   (default 0.5). Passed to \code{\link{write_sc_isoem}}.
#' @param max_ec_size_sharing Integer. Max EC size for pairwise expansion
#'   (default 50). Passed to \code{\link{write_sc_isoem}}.
#'
#' @return An \code{IsoEMSCResult} containing:
#' \itemize{
#'   \item \code{count_matrix} — sparse Matrix (transcript x cell)
#'   \item \code{cell_qc}      — data.table of per-cell QC metrics
#'   \item \code{gtf_meta}     — transcript annotation from GTF
#'   \item \code{ec_table}     — raw EC data.table (used by \code{write_sc_isoem})
#' }
#' Returned invisibly if \code{write = TRUE}.
#'
#' @seealso \code{\link{prepare_isoem}}, \code{\link{build_sc_ec}},
#'   \code{\link{run_em}}, \code{\link{write_sc_isoem}}
#'
#' @export
#'
#' @examples
#' counts_f <- system.file("extdata", "toy_counts_multi.tsv", package = "IsoEM")
#' gtf_f    <- system.file("extdata", "toy.gtf.gz",          package = "IsoEM")
#' bc_f     <- system.file("extdata", "toy_bc_umi.tsv",   package = "IsoEM")
#'
#' # anno_file mode — write automatically
#' result <- run_sc_isoem(counts_f, gtf_f,
#'   anno_file = bc_f, unit = "umi", outdir = tempdir())
#'
#' # Inspect before writing
#' result <- run_sc_isoem(counts_f, gtf_f,
#'   anno_file = bc_f, unit = "umi", write = FALSE)
#' print(result)
#' write_sc_isoem(result, outdir = tempdir())
#'
#' # regex mode
#' regex_f <- system.file("extdata", "toy_counts_regex.tsv", package = "IsoEM")
#' result  <- run_sc_isoem(regex_f, gtf_f,
#'   pattern = ".*_bc=(?P<barcode>[A-Z0-9]+)_sq=(?P<umi>[ACGT]+)",
#'   unit = "umi", outdir = tempdir())
run_sc_isoem <- function(
    counts_file,
    gtf_file,
    anno_file           = NULL,
    pattern             = NULL,
    bc_pattern          = NULL,
    umi_pattern         = NULL,
    barcodes_use        = NULL,
    unit                = c("umi", "read"),
    outdir              = "IsoEM_out",
    keep_temp           = FALSE,
    chunk_size          = 5e6L,
    max_iter            = 500L,
    tol                 = 1e-6,
    n_cores             = 1L,
    verbose             = TRUE,
    write               = TRUE,
    compress            = FALSE,
    write_ec_table      = TRUE,
    write_sharing_table = TRUE,
    min_sharing         = 0.5,
    max_ec_size_sharing = 50L
) {
  unit <- match.arg(unit)

  input <- prepare_isoem(
    counts_file  = counts_file,
    gtf_file     = gtf_file,
    mode         = "sc",
    anno_file    = anno_file,
    pattern      = pattern,
    bc_pattern   = bc_pattern,
    umi_pattern  = umi_pattern,
    barcodes_use = barcodes_use,
    unit         = unit,
    outdir       = outdir,
    verbose      = verbose
  )

  ec <- build_sc_ec(
    input      = input,
    keep_temp  = keep_temp,
    chunk_size = chunk_size,
    n_cores    = n_cores,
    verbose    = verbose
  )

  result <- run_em(
    ec       = ec,
    max_iter = max_iter,
    tol      = tol,
    n_cores  = n_cores,
    verbose  = verbose
  )

  if (write) {
    write_sc_isoem(
      result,
      outdir              = outdir,
      compress            = compress,
      write_ec_table      = write_ec_table,
      write_sharing_table = write_sharing_table,
      min_sharing         = min_sharing,
      max_ec_size_sharing = as.integer(max_ec_size_sharing),
      verbose             = verbose
    )
    return(invisible(result))
  }
  result
}
