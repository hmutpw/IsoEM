# =============================================================================
# run_em.R  —  EM quantification from IsoEMEC object
# =============================================================================

#' Run EM quantification from an IsoEMEC object
#'
#' Performs equivalence-class-compressed EM for all groups in an
#' \code{IsoEMEC} object. Works for both bulk and single-cell modes.
#'
#' @param ec An \code{IsoEMEC} object from \code{\link{build_ec}} or
#'   \code{\link{build_sc_ec}}.
#' @param gtf_file Character. Path to IsoQuant \code{transcript_models.gtf[.gz]}.
#'   Used to annotate transcripts with gene_id and is_novel.
#' @param max_iter Integer. Maximum EM iterations per group (default: 1000
#'   for bulk, 500 for sc).
#' @param tol Numeric. Convergence tolerance (default 1e-6).
#' @param n_cores Integer. Parallel workers (default 1).
#' @param verbose Logical (default TRUE).
#'
#' @return
#' \describe{
#'   \item{bulk single-sample}{\code{IsoEMResult}}
#'   \item{bulk multi-sample}{\code{IsoEMDataset}}
#'   \item{sc}{\code{IsoEMSCResult}}
#' }
#'
#' @seealso \code{\link{build_ec}}, \code{\link{build_sc_ec}},
#'   \code{\link{write_isoem}}, \code{\link{write_sc_isoem}}
#'
#' @export
#'
#' @examples
#' counts_f <- system.file("extdata", "toy_counts.tsv",  package = "IsoEM")
#' gtf_f    <- system.file("extdata", "toy.gtf",         package = "IsoEM")
#'
#' # bulk
#' input  <- prepare_isoem(counts_f, gtf_f,
#'   mode = "bulk_single", sample_id = "toy")
#' ec     <- build_ec(input)
#' result <- run_em(ec, gtf_file = gtf_f)
#' print(result)
#'
#' # single-cell
#' bc_f     <- system.file("extdata", "toy_bc_umi.tsv", package = "IsoEM")
#' input_sc <- prepare_isoem(counts_f, gtf_f,
#'   mode = "sc", anno_file = bc_f, unit = "umi")
#' ec_sc    <- build_sc_ec(input_sc)
#' result_sc <- run_em(ec_sc, gtf_file = gtf_f)
#' print(result_sc)
run_em <- function(ec,
                   gtf_file,
                   max_iter = NULL,
                   tol      = 1e-6,
                   n_cores  = 1L,
                   verbose  = TRUE) {

  stopifnot(inherits(ec, "IsoEMEC"))
  .check_file(gtf_file, "gtf_file")
  if (.Platform$OS.type == "windows") n_cores <- 1L

  # default max_iter by mode
  if (is.null(max_iter))
    max_iter <- if (ec$mode == "sc") 500L else 1000L

  if (verbose) {
    message(sprintf("=== run_em (%s mode) ===", ec$mode))
    message(sprintf("  Groups     : %d", length(ec$group_ids)))
    message(sprintf("  Transcripts: %d", ec$n_tx))
    message(sprintf("  max_iter   : %d  tol: %g  n_cores: %d",
                    max_iter, tol, n_cores))
  }

  # ---- Parse GTF -----------------------------------------------------------
  if (verbose) message("Parsing GTF ...")
  gtf_meta <- .parse_gtf(gtf_file)

  # ---- Dispatch EM per group -----------------------------------------------
  if (verbose) message("Running EM ...")
  t0 <- proc.time()

  em_results <- .run_em_dispatch(ec, max_iter, tol, n_cores, verbose)

  elapsed <- round((proc.time() - t0)[["elapsed"]], 1)
  n_conv  <- sum(sapply(em_results, `[[`, "converged"))
  if (verbose)
    message(sprintf("  EM complete: %d / %d converged (%.1f s)",
                    n_conv, length(em_results), elapsed))

  # ---- Assemble result object ----------------------------------------------
  if (verbose) message("Assembling result ...")
  .assemble_result(em_results, ec, gtf_meta, verbose)
}


# =============================================================================
# Internal: EM dispatch
# =============================================================================

#' @keywords internal
#' @noRd
.run_em_dispatch <- function(ec, max_iter, tol, n_cores, verbose) {
  n_tx      <- ec$n_tx
  group_ids <- ec$group_ids

  # worker: EM for one group
  .run_one <- function(grp) {
    ec_sub <- ec$ec_table[group_id == grp,
                          .(ec_id, t_indices, count)]
    if (nrow(ec_sub) == 0L) {
      return(list(group_id = grp, counts = rep(0, n_tx),
                  n_iter = 0L, converged = FALSE,
                  n_ec = 0L, n_unique_ec = 0L))
    }
    res <- .em_core(ec_sub, n_tx, max_iter, tol)
    res$group_id <- grp
    res
  }

  if (n_cores > 1L) {
    # batch groups into n_cores buckets to minimise scheduling overhead
    batches <- split(group_ids,
                     cut(seq_along(group_ids), n_cores, labels = FALSE))
    em_batched <- parallel::mclapply(batches, function(batch) {
      lapply(batch, .run_one)
    }, mc.cores = n_cores)
    em_results <- unlist(em_batched, recursive = FALSE)
  } else {
    em_results <- lapply(group_ids, .run_one)
  }

  names(em_results) <- group_ids
  em_results
}


# =============================================================================
# Internal: result assembly
# =============================================================================

#' @keywords internal
#' @noRd
.assemble_result <- function(em_results, ec, gtf_meta, verbose) {
  tx_ids   <- ec$tx_map          # character vector, index = t_idx
  n_tx     <- ec$n_tx
  group_ids <- ec$group_ids

  if (ec$mode == "sc") {
    # ---- SC: sparse matrix + cell QC --------------------------------------
    i_vec <- integer(0); j_vec <- integer(0); x_vec <- numeric(0)
    qc_list <- vector("list", length(group_ids))

    for (j in seq_along(group_ids)) {
      grp <- group_ids[j]
      res <- em_results[[grp]]
      nz  <- which(res$counts > 0.01)
      if (length(nz) > 0L) {
        i_vec <- c(i_vec, nz)
        j_vec <- c(j_vec, rep(j, length(nz)))
        x_vec <- c(x_vec, res$counts[nz])
      }
      qc_list[[j]] <- data.table::data.table(
        barcode                = grp,
        n_transcripts_detected = length(nz),
        em_converged           = res$converged,
        em_n_iter              = res$n_iter,
        n_ec                   = res$n_ec,
        n_unique_ec            = res$n_unique_ec
      )
    }

    mat <- Matrix::sparseMatrix(
      i    = i_vec, j = j_vec, x = x_vec,
      dims = c(n_tx, length(group_ids)),
      dimnames = list(tx_ids, group_ids)
    )
    cell_qc <- data.table::rbindlist(qc_list)

    return(new_isoem_sc_result(
      count_matrix = mat,
      cell_qc      = cell_qc,
      gtf_meta     = gtf_meta,
      params       = ec$input_params
    ))
  }

  # ---- Bulk -----------------------------------------------------------------
  .make_counts_dt <- function(grp) {
    res    <- em_results[[grp]]
    counts <- res$counts
    total  <- sum(counts)

    # unique assignment rate: reads in size-1 ECs
    ec_sub    <- ec$ec_table[group_id == grp]
    uniq_cnt  <- ec_sub[lengths(t_indices) == 1L,
                        sum(count), by = .(t_idx = sapply(t_indices, `[[`, 1L))]
    total_cnt <- ec_sub[,
      unlist(t_indices, use.names = FALSE), by = count
    ][, .(total = sum(count)), by = V1]

    dt <- data.table::data.table(
      transcript_id      = tx_ids,
      sample_id          = grp,
      em_count           = counts,
      tpm                = if (total > 0) counts / total * 1e6 else 0
    )
    # merge unique / total counts
    dt[, t_idx := seq_len(n_tx)]
    # uniq_cnt: t_idx | V1 (count) -> rename
    if (nrow(uniq_cnt) > 0L) {
      data.table::setnames(uniq_cnt, c("t_idx", "unique_count"))
      dt <- uniq_cnt[dt, on = "t_idx"]
    } else {
      dt[, unique_count := 0L]
    }
    dt[is.na(unique_count), unique_count := 0L]
    # total_cnt: V1 (t_idx) | total -> rename
    if (nrow(total_cnt) > 0L) {
      data.table::setnames(total_cnt, c("t_idx", "total_count"))
      dt <- total_cnt[dt, on = "t_idx"]
    } else {
      dt[, total_count := 0L]
    }
    dt[is.na(total_count), total_count := 0L]
    dt[, certainty         := ifelse(total_count > 0,
                                      unique_count / total_count, 0)]
    dt[, multimapping_rate := 1 - certainty]

    # annotate from GTF
    dt <- gtf_meta[dt, on = "transcript_id"]
    dt[is.na(gene_id),  gene_id  := "unknown"]
    dt[is.na(is_novel), is_novel := FALSE]
    data.table::setorder(dt, -em_count)
    dt[, t_idx := NULL]
    dt
  }

  if (length(group_ids) == 1L) {
    # single sample
    grp    <- group_ids[1L]
    res    <- em_results[[grp]]
    counts <- .make_counts_dt(grp)
    n_reads <- sum(ec$ec_table[group_id == grp, count])
    n_uniq  <- sum(ec$ec_table[group_id == grp &
                               lengths(t_indices) == 1L, count])
    qc <- list(
      total_reads            = n_reads,
      unique_assignment_rate = round(n_uniq / max(n_reads, 1), 4),
      n_transcripts_detected = sum(counts$em_count > 0.01),
      n_ec                   = res$n_ec,
      n_unique_ec            = res$n_unique_ec,
      n_iter                 = res$n_iter,
      converged              = res$converged
    )
    return(new_isoem_result(grp, counts, qc, gtf_meta))
  }

  # multi-sample
  results_list <- lapply(group_ids, function(grp) {
    counts <- .make_counts_dt(grp)
    res    <- em_results[[grp]]
    n_reads <- sum(ec$ec_table[group_id == grp, count])
    n_uniq  <- sum(ec$ec_table[group_id == grp &
                               lengths(t_indices) == 1L, count])
    qc <- list(
      total_reads            = n_reads,
      unique_assignment_rate = round(n_uniq / max(n_reads, 1), 4),
      n_transcripts_detected = sum(counts$em_count > 0.01),
      n_ec                   = res$n_ec,
      n_unique_ec            = res$n_unique_ec,
      n_iter                 = res$n_iter,
      converged              = res$converged
    )
    new_isoem_result(grp, counts, qc, gtf_meta)
  })
  names(results_list) <- group_ids

  # assemble count + TPM matrices
  all_tx <- tx_ids
  .mat <- function(col) {
    m <- matrix(0, nrow = length(all_tx), ncol = length(group_ids),
                dimnames = list(all_tx, group_ids))
    for (grp in group_ids) {
      cnt <- results_list[[grp]]$counts
      m[cnt$transcript_id, grp] <- cnt[[col]]
    }
    m
  }
  new_isoem_dataset(results_list, .mat("em_count"), .mat("tpm"))
}
