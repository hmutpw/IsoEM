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
#' @param gtf_file Character or NULL. Path to IsoQuant
#'   \code{transcript_models.gtf[.gz]} for transcript annotation.
#'   If NULL (default), the path stored in \code{ec$gtf_file} is used
#'   automatically (set by \code{\link{prepare_isoem}}).
#'   Only needed if the EC object was built without \code{prepare_isoem()}.
#' @param max_iter Integer. Maximum EM iterations per group (default: 1000
#'   for bulk, 500 for sc).
#' @param tol Numeric. Convergence tolerance (default 1e-6).
#' @param n_cores Integer. Parallel workers (default 1).
#'   Uses \code{parallel::makeCluster} — compatible with Windows and Linux.
#' @param verbose Logical (default TRUE).
#'
#' @return Depending on \code{ec$mode}:
#' \describe{
#'   \item{bulk single-sample}{\code{IsoEMResult} with fields \code{sample_id},
#'     \code{counts} (data.table), \code{qc} (list), \code{gtf_meta},
#'     \code{ec_table} (data.table or NULL)}
#'   \item{bulk multi-sample}{\code{IsoEMDataset} containing a list of
#'     \code{IsoEMResult} objects plus \code{count_matrix} and \code{tpm_matrix}}
#'   \item{sc}{\code{IsoEMSCResult} with fields \code{count_matrix} (sparse
#'     Matrix, transcript x cell), \code{cell_qc} (data.table), \code{gtf_meta},
#'     \code{ec_table} (data.table or NULL)}
#' }
#'
#' The \code{ec_table} field is populated automatically from the \code{IsoEMEC}
#' object and is used by \code{\link{write_isoem}} / \code{\link{write_sc_isoem}}
#' to produce \code{ec_table.tsv} and \code{sharing_table.tsv}.
#'
#' @seealso \code{\link{build_ec}}, \code{\link{build_sc_ec}},
#'   \code{\link{write_isoem}}, \code{\link{write_sc_isoem}}
#' @export
#'
#' @examples
#' counts_f <- system.file("extdata", "toy_counts.tsv.gz",  package = "IsoEM")
#' gtf_f    <- system.file("extdata", "toy.gtf.gz",         package = "IsoEM")
#'
#' # Bulk single-sample
#' input  <- prepare_isoem(counts_f, gtf_f,
#'   mode = "bulk_single", sample_id = "toy")
#' ec     <- build_ec(input)
#' result <- run_em(ec)
#' print(result)
#' write_isoem(result, outdir = tempdir())
#'
#' # Single-cell: build EC once, rerun EM with different params
#' counts_sc_f <- system.file("extdata", "toy_counts_multi.tsv", package = "IsoEM")
#' bc_f      <- system.file("extdata", "toy_bc_umi.tsv", package = "IsoEM")
#' input_sc  <- prepare_isoem(counts_sc_f, gtf_f,
#'   mode = "sc", anno_file = bc_f, unit = "umi")
#' ec_sc     <- build_sc_ec(input_sc)
#' # saveRDS(ec_sc, "ec.rds")  # save for parameter tuning
#' result_sc <- run_em(ec_sc, max_iter = 500, tol = 1e-6)
#' print(result_sc)
#' write_sc_isoem(result_sc, outdir = tempdir())
run_em <- function(ec,
                   gtf_file = NULL,
                   max_iter = NULL,
                   tol      = 1e-6,
                   n_cores  = 1L,
                   verbose  = TRUE) {

  stopifnot(inherits(ec, "IsoEMEC"))
  gtf_file <- gtf_file %||% ec$gtf_file
  if (is.null(gtf_file))
    stop("gtf_file not found. Provide via run_em(gtf_file=) or ",
         "ensure prepare_isoem() was called before build_ec().", call. = FALSE)
  .check_file(gtf_file, "gtf_file")

  if (is.null(max_iter))
    max_iter <- if (ec$mode == "sc") 500L else 1000L

  if (verbose) {
    message(sprintf("=== run_em (%s mode) ===", ec$mode))
    message(sprintf("  Groups     : %d", length(ec$group_ids)))
    message(sprintf("  Transcripts: %d", ec$n_tx))
    message(sprintf("  max_iter   : %d  tol: %g  n_cores: %d",
                    max_iter, tol, n_cores))
  }

  if (verbose) message("Parsing GTF ...")
  gtf_meta <- .parse_gtf(gtf_file)

  if (verbose) message("Running EM ...")
  t0 <- proc.time()
  em_results <- .run_em_parallel(ec, max_iter, tol, n_cores)
  elapsed    <- round((proc.time() - t0)[["elapsed"]], 1)
  n_conv     <- sum(vapply(em_results, `[[`, FALSE, "converged"))
  if (verbose)
    message(sprintf("  EM complete: %d / %d converged (%.1f s)",
                    n_conv, length(em_results), elapsed))

  if (verbose) message("Assembling result ...")
  .assemble_result(em_results, ec, gtf_meta)
}


# =============================================================================
# Internal: parallel EM dispatch (Windows + Linux)
# =============================================================================

#' Run EM in parallel using makeCluster (Windows + Linux compatible).
#' Uses a smart threshold: only starts a cluster when estimated work
#' exceeds cluster setup overhead (~3-5 s). Exports only the bare minimum
#' (function + data) without loading the IsoEM package on workers.
#' @keywords internal
#' @noRd
.run_em_parallel <- function(ec, max_iter, tol, n_cores) {
  group_ids <- ec$group_ids
  n_groups  <- length(group_ids)
  n_tx      <- ec$n_tx

  # pre-split ec_table by group (O(1) lookup per worker)
  ec_by_grp <- split(
    ec$ec_table[, .(ec_id, t_indices, count)],
    ec$ec_table$group_id
  )

  # bare .em_core function reference (no package needed on workers)
  em_fn <- .em_core

  .run_one <- function(grp) {
    ec_sub <- ec_by_grp[[grp]]
    if (is.null(ec_sub) || nrow(ec_sub) == 0L)
      return(list(group_id = grp, counts = rep(0, n_tx),
                  n_iter = 0L, converged = FALSE,
                  n_ec = 0L, n_unique_ec = 0L))
    res <- em_fn(ec_sub, n_tx, max_iter, tol)
    res$group_id <- grp
    res
  }

  # Smart parallel threshold:
  # cluster setup costs ~3-5 s regardless of data size.
  # Only go parallel when total EC work justifies it.
  # Heuristic: n_groups * avg_n_ec > 500,000 EC-operations
  avg_ec   <- nrow(ec$ec_table) / max(n_groups, 1L)
  ec_ops   <- as.numeric(n_groups) * avg_ec
  use_par  <- n_cores > 1L && n_groups > 1L && ec_ops > 5e5

  if (use_par) {
    n_workers <- min(n_cores, n_groups)
    cl <- parallel::makeCluster(n_workers)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    # Export only function + pre-split data + scalars
    # Do NOT call clusterEvalQ(library(IsoEM)) — that's the main overhead
    parallel::clusterExport(
      cl,
      varlist = c("ec_by_grp", "em_fn", "n_tx", "max_iter", "tol",
                  ".run_one"),
      envir   = environment()
    )
    batches <- .split_batches(group_ids, n_workers)
    results_batched <- parallel::parLapply(cl, batches, function(batch) {
      lapply(batch, .run_one)
    })
    em_results <- unlist(results_batched, recursive = FALSE)
  } else {
    em_results <- lapply(group_ids, .run_one)
  }

  stats::setNames(em_results, group_ids)
}

#' Split a vector into n balanced batches
#' @keywords internal
#' @noRd
.split_batches <- function(x, n) {
  n <- min(n, length(x))
  split(x, cut(seq_along(x), n, labels = FALSE))
}


# =============================================================================
# Internal: result assembly
# =============================================================================

#' @keywords internal
#' @noRd
.assemble_result <- function(em_results, ec, gtf_meta) {
  tx_ids    <- ec$tx_map
  n_tx      <- ec$n_tx
  group_ids <- ec$group_ids

  if (ec$mode == "sc") {
    return(.assemble_sc(em_results, group_ids, n_tx, tx_ids, gtf_meta,
                        ec$input_params, ec_table = ec$ec_table))
  }
  .assemble_bulk(em_results, ec, group_ids, n_tx, tx_ids, gtf_meta)
}

#' Assemble SC sparse matrix result — pre-allocated for speed
#' @keywords internal
#' @noRd
.assemble_sc <- function(em_results, group_ids, n_tx, tx_ids,
                          gtf_meta, params, ec_table = NULL) {
  n_cells <- length(group_ids)
  threshold <- 0.01

  # pre-collect triplets (i, j, x) using pre-allocated lists
  i_list <- vector("list", n_cells)
  x_list <- vector("list", n_cells)
  qc_list <- vector("list", n_cells)

  for (j in seq_len(n_cells)) {
    grp <- group_ids[j]
    res <- em_results[[grp]]
    nz  <- which(res$counts > threshold)
    i_list[[j]] <- nz
    x_list[[j]] <- res$counts[nz]
    qc_list[[j]] <- list(
      barcode                = grp,
      n_transcripts_detected = length(nz),
      em_converged           = res$converged,
      em_n_iter              = res$n_iter,
      n_ec                   = res$n_ec,
      n_unique_ec            = res$n_unique_ec
    )
  }

  # build sparse matrix from pre-collected triplets
  i_vec <- unlist(i_list, use.names = FALSE)
  j_vec <- rep.int(seq_len(n_cells), lengths(i_list))
  x_vec <- unlist(x_list, use.names = FALSE)

  mat <- Matrix::sparseMatrix(
    i = i_vec, j = j_vec, x = x_vec,
    dims = c(n_tx, n_cells),
    dimnames = list(tx_ids, group_ids)
  )
  cell_qc <- data.table::rbindlist(
    lapply(qc_list, data.table::as.data.table)
  )

  new_isoem_sc_result(mat, cell_qc, gtf_meta, params,
                      ec_table = ec_table)
}

#' Assemble bulk result
#' @keywords internal
#' @noRd
.assemble_bulk <- function(em_results, ec, group_ids, n_tx, tx_ids, gtf_meta) {

  # Pre-compute per-group EC summaries once (avoid repeated table scans)
  ec_summary <- ec$ec_table[, .(
    total_count = sum(count),
    uniq_count  = sum(count[lengths(t_indices) == 1L])
  ), by = group_id]

  .make_counts_dt <- function(grp) {
    res    <- em_results[[grp]]
    em_cnt <- res$counts
    total  <- sum(em_cnt)

    # unique / total count per transcript from EC table
    ec_grp <- ec$ec_table[group_id == grp]

    uniq_per_tx <- ec_grp[
      lengths(t_indices) == 1L,
      .(unique_count = sum(count)),
      by = .(t_idx = vapply(t_indices, `[`, 1L, 1L))
    ]

    # expand t_indices to get total count per transcript
    total_per_tx <- ec_grp[,
      .(t_idx = unlist(t_indices, use.names = FALSE),
        count  = rep(count, lengths(t_indices)))
    ][, .(total_count = sum(count)), by = t_idx]

    dt <- data.table::data.table(
      transcript_id = tx_ids,
      sample_id     = grp,
      em_count      = em_cnt,
      tpm           = if (total > 0) em_cnt / total * 1e6 else 0,
      t_idx         = seq_len(n_tx)
    )
    dt <- uniq_per_tx[dt, on = "t_idx"]
    dt[is.na(unique_count), unique_count := 0L]
    dt <- total_per_tx[dt, on = "t_idx"]
    dt[is.na(total_count), total_count := 0L]
    dt[, certainty         := ifelse(total_count > 0,
                                      unique_count / total_count, 0)]
    dt[, multimapping_rate := 1 - certainty]
    dt <- gtf_meta[dt, on = "transcript_id"]
    dt[is.na(gene_id),  gene_id  := "unknown"]
    dt[is.na(is_novel), is_novel := FALSE]
    data.table::setorder(dt, -em_count)
    dt[, t_idx := NULL]
    dt
  }

  if (length(group_ids) == 1L) {
    grp    <- group_ids[1L]
    res    <- em_results[[grp]]
    counts <- .make_counts_dt(grp)
    summ   <- ec_summary[group_id == grp]
    qc <- list(
      total_reads            = summ$total_count,
      unique_assignment_rate = round(summ$uniq_count /
                                     max(summ$total_count, 1), 4),
      n_transcripts_detected = sum(counts$em_count > 0.01),
      n_ec                   = res$n_ec,
      n_unique_ec            = res$n_unique_ec,
      n_iter                 = res$n_iter,
      converged              = res$converged
    )
    ec_grp <- if (!is.null(ec$ec_table))
      ec$ec_table[group_id == grp, .(ec_id, t_indices, count)]
    else NULL
    return(new_isoem_result(grp, counts, qc, gtf_meta, ec_table = ec_grp))
  }

  results_list <- lapply(group_ids, function(grp) {
    res    <- em_results[[grp]]
    counts <- .make_counts_dt(grp)
    summ   <- ec_summary[group_id == grp]
    qc <- list(
      total_reads            = summ$total_count,
      unique_assignment_rate = round(summ$uniq_count /
                                     max(summ$total_count, 1), 4),
      n_transcripts_detected = sum(counts$em_count > 0.01),
      n_ec                   = res$n_ec,
      n_unique_ec            = res$n_unique_ec,
      n_iter                 = res$n_iter,
      converged              = res$converged
    )
    ec_grp <- if (!is.null(ec$ec_table))
      ec$ec_table[group_id == grp, .(ec_id, t_indices, count)]
    else NULL
    new_isoem_result(grp, counts, qc, gtf_meta, ec_table = ec_grp)
  })
  names(results_list) <- group_ids

  all_tx <- tx_ids
  .mat <- function(col) {
    m <- matrix(0, nrow = length(all_tx), ncol = length(group_ids),
                dimnames = list(all_tx, group_ids))
    for (grp in group_ids)
      m[results_list[[grp]]$counts$transcript_id, grp] <-
        results_list[[grp]]$counts[[col]]
    m
  }
  new_isoem_dataset(results_list, .mat("em_count"), .mat("tpm"))
}
