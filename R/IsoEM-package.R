#' IsoEM: EM-Based Transcript Quantification for Long-Read RNA-seq
#'
#' IsoEM provides a modular pipeline for IsoQuant output quantification.
#' The workflow: \code{\link{prepare_isoem}} validates inputs;
#' \code{\link{build_ec}} / \code{\link{build_sc_ec}} construct equivalence
#' classes; \code{\link{run_em}} performs EM quantification;
#' \code{\link{write_isoem}} / \code{\link{write_sc_isoem}} write results.
#' One-step wrappers \code{\link{run_isoem}} and \code{\link{run_sc_isoem}}
#' are also provided.
#'
#' @docType package
#' @name IsoEM-package
#' @aliases IsoEM
#'
#' @import data.table
#' @importFrom Matrix sparseMatrix nnzero writeMM
#' @importFrom methods is
#' @importFrom parallel mclapply
#' @importFrom stats median setNames
#' @importFrom utils write.table
#' @importFrom tools file_path_sans_ext
"_PACKAGE"
