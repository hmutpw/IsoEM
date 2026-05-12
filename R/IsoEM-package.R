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
#' @importFrom parallel makeCluster stopCluster clusterExport parLapply
#' @importFrom stats median setNames
#' @importFrom utils write.table
#' @importFrom tools file_path_sans_ext
"_PACKAGE"

# Suppress R CMD check NOTEs for data.table non-standard evaluation
utils::globalVariables(c(
  ".", "N", "certainty", "count", "ec_id", "ec_key", "ec_size", "ec_type",
  "em_count", "feature_type", "gene_id", "group_id", "is_novel",
  "multimapping_rate", "obs_id", "r_idx", "read_id", "recommendation",
  "shared_reads", "sharing_fraction", "t1", "t2", "t_idx", "t_indices",
  "t_key", "total_count", "total_reads_tx1", "total_reads_tx2",
  "transcript_1", "transcript_2", "transcript_id", "transcripts",
  "umi", "unique_count"
))
