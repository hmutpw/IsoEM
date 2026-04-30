# =============================================================================
# objects.R  —  S3 class constructors + print/summary methods
# =============================================================================

# -----------------------------------------------------------------------------
# IsoEMInput  —  lightweight input descriptor (paths + validation only)
# -----------------------------------------------------------------------------

#' Constructor for IsoEMInput
#' @keywords internal
#' @noRd
new_isoem_input <- function(mode, counts_file, gtf_file,
                             anno_file, sample_id,
                             pattern, bc_pattern, umi_pattern,
                             barcodes_use, unit,
                             outdir, file_info, validated) {
  structure(
    list(
      mode        = mode,
      counts_file = counts_file,
      gtf_file    = gtf_file,
      anno_file   = anno_file,
      sample_id   = sample_id,
      pattern     = pattern,
      bc_pattern  = bc_pattern,
      umi_pattern = umi_pattern,
      barcodes_use = barcodes_use,
      unit        = unit,
      outdir      = outdir,
      file_info   = file_info,
      validated   = validated
    ),
    class = "IsoEMInput"
  )
}

#' Print an IsoEMInput object
#' @param x An \code{IsoEMInput} object.
#' @param ... Ignored.
#' @export
print.IsoEMInput <- function(x, ...) {
  cat("IsoEMInput\n")
  cat(sprintf("  Mode         : %s\n", x$mode))
  cat(sprintf("  counts_file  : %s\n", basename(x$counts_file)))
  if (!is.null(x$file_info$counts))
    cat(sprintf("                 (%s, ~%s rows)\n",
                x$file_info$counts$size_str,
                x$file_info$counts$nrow_str))
  cat(sprintf("  gtf_file     : %s\n", basename(x$gtf_file)))
  if (!is.null(x$anno_file))
    cat(sprintf("  anno_file    : %s\n", basename(x$anno_file)))
  if (!is.null(x$file_info$anno))
    cat(sprintf("                 (%s, ~%s rows)\n",
                x$file_info$anno$size_str,
                x$file_info$anno$nrow_str))
  if (!is.null(x$sample_id))
    cat(sprintf("  sample_id    : %s\n", x$sample_id))
  if (!is.null(x$pattern))
    cat(sprintf("  pattern      : %s\n", x$pattern))
  if (!is.null(x$bc_pattern))
    cat(sprintf("  bc_pattern   : %s\n", x$bc_pattern))
  if (!is.null(x$umi_pattern))
    cat(sprintf("  umi_pattern  : %s\n", x$umi_pattern))
  if (!is.null(x$barcodes_use)) {
    n_bc <- if (is.character(x$barcodes_use) && length(x$barcodes_use) == 1 &&
                file.exists(x$barcodes_use))
      x$file_info$barcodes_n
    else length(x$barcodes_use)
    cat(sprintf("  barcodes_use : %d barcodes\n", n_bc))
  }
  if (x$mode %in% c("sc"))
    cat(sprintf("  unit         : %s\n", x$unit))
  cat(sprintf("  outdir       : %s\n", x$outdir))
  cat(sprintf("  Validated    : %s\n", x$validated))
  if (!is.null(x$file_info$mem_est_gb))
    cat(sprintf("  Est. memory  : ~%.1f GB peak\n", x$file_info$mem_est_gb))
  invisible(x)
}

#' @export
summary.IsoEMInput <- function(object, ...) print(object, ...)


# -----------------------------------------------------------------------------
# IsoEMEC  —  equivalence class table (intermediate, saveable)
# -----------------------------------------------------------------------------

#' Constructor for IsoEMEC
#' @keywords internal
#' @noRd
new_isoem_ec <- function(ec_table, tx_map, group_ids, n_tx,
                          mode, unit, gtf_file, input_params) {
  structure(
    list(
      ec_table     = ec_table,    # data.table: group_id|ec_id|t_indices|count
      tx_map       = tx_map,      # character vector: t_idx -> transcript_id
      group_ids    = group_ids,   # all group IDs
      n_tx         = n_tx,        # total transcript count
      mode         = mode,        # "bulk" / "sc"
      unit         = unit,        # "read" / "umi"
      gtf_file     = gtf_file,    # path to GTF for annotation in run_em()
      input_params = input_params # snapshot of key params for reproducibility
    ),
    class = "IsoEMEC"
  )
}

#' Print an IsoEMEC object
#' @param x An \code{IsoEMEC} object.
#' @param ... Ignored.
#' @export
print.IsoEMEC <- function(x, ...) {
  cat("IsoEMEC\n")
  cat(sprintf("  Mode          : %s\n", x$mode))
  cat(sprintf("  Unit          : %s\n", x$unit))
  cat(sprintf("  GTF           : %s\n", basename(x$gtf_file %||% "not set")))
  cat(sprintf("  Groups        : %d\n", length(x$group_ids)))
  cat(sprintf("  Transcripts   : %d\n", x$n_tx))
  cat(sprintf("  Total ECs     : %d\n", nrow(x$ec_table)))
  med_ec <- stats::median(
    x$ec_table[, .N, by = group_id]$N
  )
  cat(sprintf("  Median EC/grp : %.0f\n", med_ec))
  invisible(x)
}

#' @export
summary.IsoEMEC <- function(object, ...) print(object, ...)


# -----------------------------------------------------------------------------
# IsoEMResult  —  bulk single-sample result
# -----------------------------------------------------------------------------

#' Constructor for IsoEMResult
#' @keywords internal
#' @noRd
new_isoem_result <- function(sample_id, counts, qc, gtf_meta) {
  structure(
    list(sample_id = sample_id,
         counts    = counts,
         qc        = qc,
         gtf_meta  = gtf_meta),
    class = "IsoEMResult"
  )
}

#' Print an IsoEMResult object
#' @param x An \code{IsoEMResult} object.
#' @param ... Ignored.
#' @export
print.IsoEMResult <- function(x, ...) {
  cat("IsoEMResult\n")
  cat(sprintf("  Sample             : %s\n", x$sample_id))
  cat(sprintf("  Transcripts        : %d\n", nrow(x$counts)))
  cat(sprintf("  Novel transcripts  : %d\n",
              sum(x$counts$is_novel, na.rm = TRUE)))
  cat(sprintf("  Total counts       : %d\n", x$qc$total_reads))
  cat(sprintf("  Unique assign rate : %.1f%%\n",
              x$qc$unique_assignment_rate * 100))
  cat(sprintf("  Converged          : %s (%d iterations)\n",
              x$qc$converged, x$qc$n_iter))
  invisible(x)
}

#' @export
summary.IsoEMResult <- function(object, ...) print(object, ...)


# -----------------------------------------------------------------------------
# IsoEMDataset  —  bulk multi-sample result
# -----------------------------------------------------------------------------

#' Constructor for IsoEMDataset
#' @keywords internal
#' @noRd
new_isoem_dataset <- function(results, count_matrix, tpm_matrix) {
  structure(
    list(results      = results,
         count_matrix = count_matrix,
         tpm_matrix   = tpm_matrix,
         n_samples    = length(results)),
    class = "IsoEMDataset"
  )
}

#' Print an IsoEMDataset object
#' @param x An \code{IsoEMDataset} object.
#' @param ... Ignored.
#' @export
print.IsoEMDataset <- function(x, ...) {
  cat("IsoEMDataset\n")
  cat(sprintf("  Samples      : %d\n", x$n_samples))
  cat(sprintf("  Transcripts  : %d\n", nrow(x$count_matrix)))
  cat("  Sample IDs   :", paste(colnames(x$count_matrix), collapse = ", "), "\n")
  invisible(x)
}

#' @export
summary.IsoEMDataset <- function(object, ...) print(object, ...)


# -----------------------------------------------------------------------------
# IsoEMSCResult  —  single-cell result
# -----------------------------------------------------------------------------

#' Constructor for IsoEMSCResult
#' @keywords internal
#' @noRd
new_isoem_sc_result <- function(count_matrix, cell_qc, gtf_meta, params) {
  structure(
    list(count_matrix = count_matrix,
         cell_qc      = cell_qc,
         gtf_meta     = gtf_meta,
         params       = params),
    class = "IsoEMSCResult"
  )
}

#' Print an IsoEMSCResult object
#' @param x An \code{IsoEMSCResult} object.
#' @param ... Ignored.
#' @export
print.IsoEMSCResult <- function(x, ...) {
  mat <- x$count_matrix
  qc  <- x$cell_qc
  cat("IsoEMSCResult\n")
  cat(sprintf("  Cells              : %d\n",  ncol(mat)))
  cat(sprintf("  Transcripts        : %d\n",  nrow(mat)))
  cat(sprintf("  Non-zero entries   : %d\n",  Matrix::nnzero(mat)))
  cat(sprintf("  Median UMI / cell  : %.0f\n",
              stats::median(Matrix::colSums(mat))))
  cat(sprintf("  Cells converged    : %d / %d\n",
              sum(qc$em_converged, na.rm = TRUE), nrow(qc)))
  invisible(x)
}

#' @export
summary.IsoEMSCResult <- function(object, ...) print(object, ...)
