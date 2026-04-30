# =============================================================================
# run_isoem.R  —  One-step wrappers
# =============================================================================

#' One-step bulk EM quantification
#'
#' Convenience wrapper that calls \code{\link{prepare_isoem}},
#' \code{\link{build_ec}}, and \code{\link{run_em}} in sequence.
#'
#' @inheritParams prepare_isoem
#' @inheritParams build_ec
#' @inheritParams run_em
#'
#' @return An \code{IsoEMResult} (single sample) or
#'   \code{IsoEMDataset} (multi-sample).
#'
#' @seealso \code{\link{prepare_isoem}}, \code{\link{build_ec}},
#'   \code{\link{run_em}}, \code{\link{write_isoem}}
#'
#' @export
#'
#' @examples
#' counts_f <- system.file("extdata", "toy_counts.tsv", package = "IsoEM")
#' gtf_f    <- system.file("extdata", "toy.gtf",        package = "IsoEM")
#' result   <- run_isoem(counts_f, gtf_f,
#'   mode = "bulk_single", sample_id = "toy_sample")
#' print(result)
run_isoem <- function(
    counts_file,
    gtf_file,
    mode       = c("bulk_single", "bulk_multi"),
    anno_file  = NULL,
    sample_id  = NULL,
    pattern    = NULL,
    outdir     = "IsoEM_out",
    keep_temp  = FALSE,
    chunk_size = 5e6L,
    max_iter   = 1000L,
    tol        = 1e-6,
    n_cores    = 1L,
    verbose    = TRUE
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

  run_em(
    ec       = ec,
    max_iter = max_iter,
    tol      = tol,
    n_cores  = n_cores,
    verbose  = verbose
  )
}


#' One-step single-cell EM quantification
#'
#' Convenience wrapper that calls \code{\link{prepare_isoem}},
#' \code{\link{build_sc_ec}}, and \code{\link{run_em}} in sequence.
#'
#' @inheritParams prepare_isoem
#' @inheritParams build_sc_ec
#' @inheritParams run_em
#'
#' @return An \code{IsoEMSCResult}.
#'
#' @seealso \code{\link{prepare_isoem}}, \code{\link{build_sc_ec}},
#'   \code{\link{run_em}}, \code{\link{write_sc_isoem}}
#'
#' @export
#'
#' @examples
#' counts_f <- system.file("extdata", "toy_counts.tsv",   package = "IsoEM")
#' gtf_f    <- system.file("extdata", "toy.gtf",          package = "IsoEM")
#' bc_f     <- system.file("extdata", "toy_bc_umi.tsv",   package = "IsoEM")
#' result   <- run_sc_isoem(counts_f, gtf_f,
#'   anno_file = bc_f, unit = "umi")
#' print(result)
run_sc_isoem <- function(
    counts_file,
    gtf_file,
    anno_file    = NULL,
    pattern      = NULL,
    bc_pattern   = NULL,
    umi_pattern  = NULL,
    barcodes_use = NULL,
    unit         = c("umi", "read"),
    outdir       = "IsoEM_out",
    keep_temp    = FALSE,
    chunk_size   = 5e6L,
    max_iter     = 500L,
    tol          = 1e-6,
    n_cores      = 1L,
    verbose      = TRUE
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

  run_em(
    ec       = ec,
    max_iter = max_iter,
    tol      = tol,
    n_cores  = n_cores,
    verbose  = verbose
  )
}
