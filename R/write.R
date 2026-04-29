# =============================================================================
# write.R  —  Write IsoEM results to disk
# =============================================================================

# -----------------------------------------------------------------------------
# write_isoem()  —  bulk
# -----------------------------------------------------------------------------

#' Write bulk IsoEM results to disk
#'
#' Writes all output files for an \code{IsoEMResult} (single sample) or
#' \code{IsoEMDataset} (multi-sample) object.
#'
#' Output files:
#' \itemize{
#'   \item \code{counts.tsv[.gz]} — EM-estimated transcript counts
#'   \item \code{qc.tsv[.gz]}     — QC summary metrics
#' }
#'
#' @param result An \code{IsoEMResult} or \code{IsoEMDataset}.
#' @param outdir Character. Output directory (created if absent).
#' @param compress Logical. Write \code{.gz} files (default FALSE).
#' @param verbose Logical (default TRUE).
#'
#' @return Invisibly returns \code{outdir}.
#' @export
#'
#' @examples
#' counts_f <- system.file("extdata", "toy_counts.tsv", package = "IsoEM")
#' gtf_f    <- system.file("extdata", "toy.gtf",        package = "IsoEM")
#' result   <- run_isoem(counts_f, gtf_f, sample_id = "toy")
#' tmp      <- tempdir()
#' write_isoem(result, outdir = tmp)
write_isoem <- function(result, outdir, compress = FALSE, verbose = TRUE) {
  UseMethod("write_isoem")
}

#' @export
write_isoem.IsoEMResult <- function(result, outdir, compress = FALSE,
                                     verbose = TRUE) {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

  # counts
  cnt_path <- file.path(outdir, .ext("counts.tsv", compress))
  .write_tsv(result$counts, cnt_path, compress)

  # qc
  qc_dt <- data.table::as.data.table(
    as.list(result$qc)[c("total_reads", "unique_assignment_rate",
                         "n_transcripts_detected", "n_unique_ecs",
                         "n_multi_ecs", "n_iter", "converged")]
  )
  qc_dt[, sample_id := result$sample_id]
  qc_path <- file.path(outdir, .ext("qc.tsv", compress))
  .write_tsv(qc_dt, qc_path, compress)

  if (verbose) {
    message("Wrote IsoEM results -> ", outdir)
    message("  counts : ", basename(cnt_path))
    message("  qc     : ", basename(qc_path))
  }
  invisible(outdir)
}

#' @export
write_isoem.IsoEMDataset <- function(result, outdir, compress = FALSE,
                                      verbose = TRUE) {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  for (s in names(result$results)) {
    sdir <- file.path(outdir, s)
    write_isoem(result$results[[s]], sdir, compress, verbose = FALSE)
  }
  # combined count matrix
  cnt_mat <- as.data.frame(result$count_matrix)
  cnt_mat <- cbind(transcript_id = rownames(cnt_mat), cnt_mat)
  mat_path <- file.path(outdir, .ext("count_matrix.tsv", compress))
  .write_tsv(cnt_mat, mat_path, compress)

  if (verbose) {
    message("Wrote IsoEMDataset results -> ", outdir)
    message("  per-sample subdirs + count_matrix.tsv")
  }
  invisible(outdir)
}


# -----------------------------------------------------------------------------
# write_sc_isoem()  —  single-cell
# -----------------------------------------------------------------------------

#' Write single-cell IsoEM results to disk
#'
#' Writes all output files for an \code{IsoEMSCResult} object:
#' \itemize{
#'   \item \code{matrix/matrix.mtx[.gz]}    — sparse EM count matrix
#'   \item \code{matrix/features.tsv[.gz]}  — transcript annotations (rows)
#'   \item \code{matrix/barcodes.tsv[.gz]}  — cell barcodes (columns)
#'   \item \code{cell_qc.tsv[.gz]}          — per-cell QC metrics
#' }
#'
#' The \code{matrix/} subdirectory follows the 10x Genomics / STARsolo
#' convention and can be loaded directly with \code{Seurat::Read10X()}.
#'
#' @param result An \code{IsoEMSCResult} from \code{\link{run_sc_isoem}}.
#' @param outdir Character. Output directory (created if absent).
#' @param compress Logical. Write \code{.gz} compressed files (default FALSE).
#' @param verbose Logical (default TRUE).
#'
#' @return Invisibly returns \code{outdir}.
#' @export
#'
#' @examples
#' counts_f <- system.file("extdata", "toy_counts.tsv", package = "IsoEM")
#' gtf_f    <- system.file("extdata", "toy.gtf",        package = "IsoEM")
#' bc_f     <- system.file("extdata", "toy_bc_umi.tsv", package = "IsoEM")
#' result   <- run_sc_isoem(counts_f, gtf_f, bc_f)
#' write_sc_isoem(result, outdir = tempdir())
write_sc_isoem <- function(result, outdir, compress = FALSE, verbose = TRUE) {
  stopifnot(inherits(result, "IsoEMSCResult"))
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

  # ---- matrix/ subdirectory ------------------------------------------------
  mat_dir <- file.path(outdir, "matrix")
  if (!dir.exists(mat_dir)) dir.create(mat_dir)

  # matrix.mtx
  mtx_path <- file.path(mat_dir, .ext("matrix.mtx", compress))
  if (compress) {
    tmp <- tempfile(fileext = ".mtx")
    Matrix::writeMM(result$count_matrix, file = tmp)
    .gzip_file(tmp, mtx_path)
  } else {
    Matrix::writeMM(result$count_matrix, file = mtx_path)
  }

  # features.tsv: transcript_id, gene_id, is_novel
  feat <- data.table::data.table(
    transcript_id = rownames(result$count_matrix)
  )
  if (!is.null(result$gtf_meta) && nrow(result$gtf_meta) > 0L) {
    feat <- result$gtf_meta[feat, on = "transcript_id"]
    feat[is.na(gene_id),  gene_id  := "unknown"]
    feat[is.na(is_novel), is_novel := FALSE]
  }
  feat_path <- file.path(mat_dir, .ext("features.tsv", compress))
  .write_tsv(feat, feat_path, compress, col.names = FALSE)

  # barcodes.tsv
  bc_path <- file.path(mat_dir, .ext("barcodes.tsv", compress))
  bc_con  <- .write_con(bc_path, compress)
  writeLines(colnames(result$count_matrix), bc_con)
  close(bc_con)

  # ---- cell_qc.tsv ---------------------------------------------------------
  qc_path <- file.path(outdir, .ext("cell_qc.tsv", compress))
  .write_tsv(result$cell_qc, qc_path, compress, col.names = TRUE)

  if (verbose) {
    message("Wrote sc-IsoEM results -> ", outdir)
    message("  matrix/matrix.mtx   : ", basename(mtx_path))
    message("  matrix/features.tsv : ", basename(feat_path))
    message("  matrix/barcodes.tsv : ", basename(bc_path))
    message("  cell_qc.tsv         : ", basename(qc_path))
  }
  invisible(outdir)
}

# helper: plain gzip without R.utils dependency
.gzip_file <- function(src, dest) {
  buf <- readBin(src, "raw", file.info(src)$size)
  con <- gzfile(dest, "wb")
  writeBin(buf, con)
  close(con)
  file.remove(src)
  invisible(dest)
}
