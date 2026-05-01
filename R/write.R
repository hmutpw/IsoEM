# =============================================================================
# write.R  --  Write IsoEM results to disk
# =============================================================================

# -----------------------------------------------------------------------------
# write_isoem()  --  bulk
# -----------------------------------------------------------------------------

#' Write bulk IsoEM results to disk
#'
#' Writes all output files for an \code{IsoEMResult} (single sample) or
#' \code{IsoEMDataset} (multi-sample) object.
#'
#' Output files (all optional outputs require EC data to be present,
#' i.e. the result must have been produced by \code{\link{run_em}} on an
#' \code{IsoEMEC} object from \code{\link{build_ec}}):
#' \itemize{
#'   \item \code{counts.tsv[.gz]}        -- EM-estimated transcript counts
#'   \item \code{qc.tsv[.gz]}            -- QC summary metrics
#'   \item \code{ec_table.tsv[.gz]}      -- Equivalence class table (if EC data present)
#'   \item \code{sharing_table.tsv[.gz]} -- Transcript sharing pairs (if EC data present)
#' }
#'
#' For \code{IsoEMDataset} (multi-sample), per-sample subdirectories are
#' created alongside a combined \code{count_matrix.tsv}.
#'
#' @param result   An \code{IsoEMResult} or \code{IsoEMDataset}.
#' @param outdir   Character. Output directory (created if absent).
#' @param compress Logical. Write \code{.gz} compressed files (default FALSE).
#' @param write_ec_table Logical. Write \code{ec_table.tsv} if EC data is
#'   available in the result object (default TRUE). The EC table contains one
#'   row per equivalence class with columns \code{sample_id}, \code{ec_id},
#'   \code{transcripts} (pipe-separated), \code{count}, \code{ec_size},
#'   \code{ec_type} (\code{"unique"} or \code{"multi"}).
#' @param write_sharing_table Logical. Write \code{sharing_table.tsv} if EC
#'   data is available (default TRUE). The sharing table identifies transcript
#'   pairs that share a large fraction of reads through multi-mapping ECs.
#'   See \code{min_sharing} and \code{max_ec_size_sharing}.
#' @param min_sharing Numeric in [0, 1]. Only transcript pairs with
#'   \code{sharing_fraction >= min_sharing} are reported in the sharing table
#'   (default 0.5). \code{sharing_fraction} is defined as
#'   \code{shared_reads / min(total_reads_tx1, total_reads_tx2)}.
#'   Lower values produce larger output files.
#' @param max_ec_size_sharing Integer. ECs with more than this many transcripts
#'   are excluded from pairwise sharing computation (default 50). This prevents
#'   combinatorial explosion at highly repetitive TE loci (e.g. L1 LINE
#'   families with 500+ near-identical copies).
#' @param verbose  Logical. Print progress messages (default TRUE).
#'
#' @return Invisibly returns \code{outdir}.
#'
#' @seealso \code{\link{run_isoem}}, \code{\link{run_em}},
#'   \code{\link{write_sc_isoem}}
#'
#' @export
#'
#' @examples
#' counts_f <- system.file("extdata", "toy_counts.tsv.gz", package = "IsoEM")
#' gtf_f    <- system.file("extdata", "toy.gtf.gz",        package = "IsoEM")
#' result   <- run_isoem(counts_f, gtf_f,
#'   mode = "bulk_single", sample_id = "toy")
#' tmp <- tempdir()
#' write_isoem(result, outdir = tmp)
#'
#' # Write with all outputs, compressed
#' write_isoem(result, outdir = tmp, compress = TRUE,
#'             write_ec_table      = TRUE,
#'             write_sharing_table = TRUE,
#'             min_sharing         = 0.5,
#'             max_ec_size_sharing = 50L)
write_isoem <- function(result, outdir, compress = FALSE,
                        write_ec_table       = TRUE,
                        write_sharing_table  = TRUE,
                        min_sharing          = 0.5,
                        max_ec_size_sharing  = 50L,
                        verbose              = TRUE) {
  UseMethod("write_isoem")
}

#' @export
write_isoem.IsoEMResult <- function(result, outdir, compress = FALSE,
                                     write_ec_table      = TRUE,
                                     write_sharing_table = TRUE,
                                     min_sharing         = 0.5,
                                     max_ec_size_sharing = 50L,
                                     verbose             = TRUE) {
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

  # ec_table
  ec_path <- NULL
  if (write_ec_table) {
    if (!is.null(result$ec_table) && nrow(result$ec_table) > 0L) {
      ec_path <- file.path(outdir, .ext("ec_table.tsv", compress))
      .write_tsv(.format_ec_bulk(result$ec_table, result$counts$transcript_id,
                                 result$sample_id),
                 ec_path, compress)
    } else if (verbose) {
      message("  EC data not available -- skipping ec_table.tsv")
    }
  }

  # sharing_table
  sh_path <- NULL
  if (write_sharing_table) {
    if (!is.null(result$ec_table) && nrow(result$ec_table) > 0L) {
      sh_dt   <- .compute_sharing(result$ec_table,
                                  result$counts$transcript_id,
                                  result$sample_id,
                                  min_sharing, as.integer(max_ec_size_sharing))
      sh_path <- file.path(outdir, .ext("sharing_table.tsv", compress))
      .write_tsv(sh_dt, sh_path, compress)
    } else if (verbose) {
      message("  EC data not available -- skipping sharing_table.tsv")
    }
  }

  if (verbose) {
    message("Wrote IsoEM results -> ", outdir)
    message("  counts  : ", basename(cnt_path))
    message("  qc      : ", basename(qc_path))
    if (!is.null(ec_path)) message("  ec_table: ", basename(ec_path))
    if (!is.null(sh_path)) message("  sharing : ", basename(sh_path))
  }
  invisible(outdir)
}

#' @export
write_isoem.IsoEMDataset <- function(result, outdir, compress = FALSE,
                                      write_ec_table      = TRUE,
                                      write_sharing_table = TRUE,
                                      min_sharing         = 0.5,
                                      max_ec_size_sharing = 50L,
                                      verbose             = TRUE) {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  for (s in names(result$results)) {
    sdir <- file.path(outdir, s)
    write_isoem(result$results[[s]], sdir, compress,
                write_ec_table      = write_ec_table,
                write_sharing_table = write_sharing_table,
                min_sharing         = min_sharing,
                max_ec_size_sharing = as.integer(max_ec_size_sharing),
                verbose             = FALSE)
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
# write_sc_isoem()  --  single-cell
# -----------------------------------------------------------------------------

#' Write single-cell IsoEM results to disk
#'
#' Writes all output files for an \code{IsoEMSCResult} object.
#'
#' Output files:
#' \itemize{
#'   \item \code{matrix/matrix.mtx[.gz]}    -- sparse EM count matrix (transcript x cell)
#'   \item \code{matrix/features.tsv[.gz]}  -- transcript annotations: transcript_id, gene_id, is_novel
#'   \item \code{matrix/barcodes.tsv[.gz]}  -- one cell barcode per line
#'   \item \code{cell_qc.tsv[.gz]}          -- per-cell QC metrics
#'   \item \code{ec_table.tsv[.gz]}         -- equivalence class table (if EC data present)
#'   \item \code{sharing_table.tsv[.gz]}    -- transcript sharing pairs pooled across cells (if EC data present)
#' }
#'
#' The \code{matrix/} subdirectory follows the 10x Genomics / STARsolo
#' convention and can be loaded directly with \code{Seurat::Read10X()}.
#'
#' The \code{ec_table.tsv} for single-cell data contains one row per
#' barcode-EC combination; the \code{sharing_table.tsv} aggregates counts
#' across all cells before computing sharing fractions.
#'
#' \code{cell_qc.tsv} columns:
#' \itemize{
#'   \item \code{barcode}                -- cell barcode
#'   \item \code{n_transcripts_detected} -- transcripts with EM count > 0.01
#'   \item \code{em_converged}           -- whether EM converged
#'   \item \code{em_n_iter}              -- EM iterations used
#'   \item \code{n_ec}                  -- total equivalence classes
#'   \item \code{n_unique_ec}            -- single-transcript ECs
#'   \item \code{total_reads}            -- total UMIs in this cell
#'   \item \code{unique_reads}           -- UMIs from unique-mapping ECs
#'   \item \code{unique_read_frac}       -- unique_reads / total_reads
#' }
#'
#' @param result   An \code{IsoEMSCResult} from \code{\link{run_sc_isoem}}
#'   or \code{\link{run_em}}.
#' @param outdir   Character. Output directory (created if absent).
#' @param compress Logical. Write \code{.gz} compressed files (default FALSE).
#' @param write_ec_table Logical. Write \code{ec_table.tsv} if EC data is
#'   available in the result object (default TRUE). Columns: \code{group_id}
#'   (barcode), \code{ec_id}, \code{transcripts} (pipe-separated),
#'   \code{count}, \code{ec_size}, \code{ec_type} (\code{"unique"} or
#'   \code{"multi"}).
#' @param write_sharing_table Logical. Write \code{sharing_table.tsv} if EC
#'   data is available (default TRUE). Sharing is computed by pooling UMI
#'   counts across all barcodes first. See \code{min_sharing} and
#'   \code{max_ec_size_sharing}.
#' @param min_sharing Numeric in [0, 1]. Only transcript pairs with
#'   \code{sharing_fraction >= min_sharing} are reported (default 0.5).
#'   \code{sharing_fraction = shared_reads / min(total_reads_tx1, total_reads_tx2)}.
#' @param max_ec_size_sharing Integer. ECs with more than this many transcripts
#'   are excluded from pairwise sharing computation (default 50). Prevents
#'   combinatorial explosion at repetitive TE loci.
#' @param verbose  Logical. Print progress messages (default TRUE).
#'
#' @return Invisibly returns \code{outdir}.
#'
#' @seealso \code{\link{run_sc_isoem}}, \code{\link{run_em}},
#'   \code{\link{write_isoem}}
#'
#' @export
#'
#' @examples
#' counts_f <- system.file("extdata", "toy_counts_multi.tsv", package = "IsoEM")
#' gtf_f    <- system.file("extdata", "toy.gtf.gz",            package = "IsoEM")
#' bc_f     <- system.file("extdata", "toy_bc_umi.tsv",        package = "IsoEM")
#' result   <- run_sc_isoem(counts_f, gtf_f, anno_file = bc_f, unit = "umi")
#' write_sc_isoem(result, outdir = tempdir())
#'
#' # Write with all outputs, gzip compressed
#' write_sc_isoem(result, outdir = tempdir(),
#'                compress            = TRUE,
#'                write_ec_table      = TRUE,
#'                write_sharing_table = TRUE,
#'                min_sharing         = 0.5,
#'                max_ec_size_sharing = 50L)
write_sc_isoem <- function(result, outdir, compress = FALSE,
                            write_ec_table      = TRUE,
                            write_sharing_table = TRUE,
                            min_sharing         = 0.5,
                            max_ec_size_sharing = 50L,
                            verbose             = TRUE) {
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

  # ec_table
  ec_path <- NULL
  if (write_ec_table) {
    if (!is.null(result$ec_table) && nrow(result$ec_table) > 0L) {
      tx_ids  <- rownames(result$count_matrix)
      ec_path <- file.path(outdir, .ext("ec_table.tsv", compress))
      .write_tsv(.format_ec_sc(result$ec_table, tx_ids),
                 ec_path, compress)
    } else if (verbose) {
      message("  EC data not available -- skipping ec_table.tsv")
    }
  }

  # sharing_table
  sh_path <- NULL
  if (write_sharing_table) {
    if (!is.null(result$ec_table) && nrow(result$ec_table) > 0L) {
      tx_ids  <- rownames(result$count_matrix)
      sh_dt   <- .compute_sharing_sc(result$ec_table, tx_ids,
                                      min_sharing,
                                      as.integer(max_ec_size_sharing))
      sh_path <- file.path(outdir, .ext("sharing_table.tsv", compress))
      .write_tsv(sh_dt, sh_path, compress)
    } else if (verbose) {
      message("  EC data not available -- skipping sharing_table.tsv")
    }
  }

  if (verbose) {
    message("Wrote sc-IsoEM results -> ", outdir)
    message("  matrix/matrix.mtx   : ", basename(mtx_path))
    message("  matrix/features.tsv : ", basename(feat_path))
    message("  matrix/barcodes.tsv : ", basename(bc_path))
    message("  cell_qc.tsv         : ", basename(qc_path))
    if (!is.null(ec_path)) message("  ec_table.tsv        : ", basename(ec_path))
    if (!is.null(sh_path)) message("  sharing_table.tsv   : ", basename(sh_path))
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


# =============================================================================
# Internal helpers: EC formatting and sharing computation
# =============================================================================

#' Format bulk EC table for output
#' Translates integer t_indices -> pipe-separated transcript names,
#' adds ec_size and ec_type columns.
#' @keywords internal
#' @noRd
.format_ec_bulk <- function(ec_table, tx_ids, sample_id) {
  dt <- data.table::copy(ec_table)
  dt[, transcripts := vapply(t_indices, function(idx) {
    paste(tx_ids[idx], collapse = "|")
  }, character(1L))]
  dt[, ec_size := lengths(t_indices)]
  dt[, ec_type := data.table::fifelse(ec_size == 1L, "unique", "multi")]
  dt[, t_indices := NULL]
  dt[, sample_id := sample_id]
  data.table::setcolorder(dt, c("sample_id", "ec_id", "transcripts",
                                 "count", "ec_size", "ec_type"))
  dt
}

#' Format SC EC table for output
#' Translates integer t_indices -> pipe-separated transcript names.
#' group_id column (barcode) is preserved as-is.
#' @keywords internal
#' @noRd
.format_ec_sc <- function(ec_table, tx_ids) {
  dt <- data.table::copy(ec_table)
  dt[, transcripts := vapply(t_indices, function(idx) {
    paste(tx_ids[idx], collapse = "|")
  }, character(1L))]
  dt[, ec_size := lengths(t_indices)]
  dt[, ec_type := data.table::fifelse(ec_size == 1L, "unique", "multi")]
  dt[, t_indices := NULL]
  data.table::setcolorder(dt, c("group_id", "ec_id", "transcripts",
                                 "count", "ec_size", "ec_type"))
  dt
}

#' Compute pairwise transcript sharing table (bulk)
#'
#' For each pair of transcripts that co-occur in multi-mapping ECs,
#' computes the fraction of their reads that are shared.
#' ECs larger than max_ec_size_sharing are skipped (TE repeat clusters).
#'
#' Output columns: sample_id, transcript_1, transcript_2,
#'   shared_reads, total_reads_tx1, total_reads_tx2, sharing_fraction,
#'   recommendation
#' @keywords internal
#' @noRd
.compute_sharing <- function(ec_table, tx_ids, sample_id,
                               min_sharing, max_ec_size_sharing) {
  # multi-mapping ECs only; skip oversized (repetitive TE) ECs
  multi <- ec_table[lengths(t_indices) > 1L &
                    lengths(t_indices) <= max_ec_size_sharing]
  if (nrow(multi) == 0L) {
    return(data.table::data.table(
      sample_id = character(), transcript_1 = character(),
      transcript_2 = character(), shared_reads = integer(),
      total_reads_tx1 = numeric(), total_reads_tx2 = numeric(),
      sharing_fraction = numeric(), recommendation = character()
    ))
  }

  # Per-transcript total reads: flatten (t_idx, count) pairs across all ECs.
  # Build the flat table from raw vectors (not via data.table j) to avoid
  # length-mismatch errors when unlist(t_indices) != nrow(ec_table).
  flat <- data.table::data.table(
    t_idx = unlist(ec_table$t_indices, use.names = FALSE),
    count = rep(ec_table$count, lengths(ec_table$t_indices))
  )
  tx_total <- flat[, .(total_reads = sum(count)), by = t_idx]

  # Pairwise expansion of multi-ECs
  pairs_list <- vector("list", nrow(multi))
  for (k in seq_len(nrow(multi))) {
    idx <- multi$t_indices[[k]]
    n   <- length(idx)
    if (n < 2L) next
    cm  <- utils::combn(n, 2L)
    pairs_list[[k]] <- data.table::data.table(
      t1           = idx[cm[1L, ]],
      t2           = idx[cm[2L, ]],
      shared_reads = multi$count[k]
    )
  }
  pairs <- data.table::rbindlist(pairs_list, fill = TRUE)
  pairs <- pairs[!is.na(t1)]
  if (nrow(pairs) == 0L) return(.empty_sharing_dt(sample_id))

  # Aggregate shared reads per pair
  agg <- pairs[, .(shared_reads = sum(shared_reads)), by = .(t1, t2)]

  # Join total reads
  agg <- tx_total[agg, on = .(t_idx = t1)]
  data.table::setnames(agg, c("t_idx", "total_reads"), c("t1", "total_reads_tx1"))
  agg <- tx_total[agg, on = .(t_idx = t2)]
  data.table::setnames(agg, c("t_idx", "total_reads"), c("t2", "total_reads_tx2"))

  # sharing_fraction = shared / min(total_tx1, total_tx2)
  agg[, sharing_fraction := shared_reads /
        pmin(total_reads_tx1, total_reads_tx2, na.rm = TRUE)]
  agg <- agg[sharing_fraction >= min_sharing]
  if (nrow(agg) == 0L) return(.empty_sharing_dt(sample_id))

  agg[, transcript_1 := tx_ids[t1]]
  agg[, transcript_2 := tx_ids[t2]]
  agg[, sample_id := sample_id]
  agg[, recommendation := data.table::fifelse(
    sharing_fraction >= 0.9,
    "consider_merging",
    "ambiguous_quantification"
  )]
  data.table::setcolorder(agg, c("sample_id", "transcript_1", "transcript_2",
                                  "shared_reads", "total_reads_tx1",
                                  "total_reads_tx2", "sharing_fraction",
                                  "recommendation"))
  agg[, c("t1", "t2") := NULL]
  data.table::setorder(agg, -sharing_fraction)
  agg
}

#' SC variant of sharing computation -- aggregates across all cells
#' @keywords internal
#' @noRd
.compute_sharing_sc <- function(ec_table, tx_ids,
                                  min_sharing, max_ec_size_sharing) {
  # Pool counts across barcodes for sharing analysis
  pooled <- ec_table[, .(
    t_indices = t_indices,
    count     = count
  )]
  # Deduplicate by EC structure (group_id irrelevant for sharing)
  pooled <- ec_table[,
    .(count = sum(count)),
    by = .(t_key = vapply(t_indices,
                          function(x) paste(sort(x), collapse="|"),
                          character(1L)))
  ]
  pooled[, t_indices := lapply(
    strsplit(t_key, "|", fixed = TRUE), as.integer
  )]

  # Re-use bulk sharing logic (sample_id = "all_cells")
  .compute_sharing(pooled, tx_ids, "all_cells",
                   min_sharing, max_ec_size_sharing)
}

#' Empty sharing data.table with correct schema
#' @keywords internal
#' @noRd
.empty_sharing_dt <- function(sample_id = character()) {
  data.table::data.table(
    sample_id        = character(),
    transcript_1     = character(),
    transcript_2     = character(),
    shared_reads     = integer(),
    total_reads_tx1  = numeric(),
    total_reads_tx2  = numeric(),
    sharing_fraction = numeric(),
    recommendation   = character()
  )
}
