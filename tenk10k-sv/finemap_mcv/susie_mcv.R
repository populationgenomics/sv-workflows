#!/usr/bin/env Rscript
# SuSiE (multiple-causal-variant) fine-mapping of SNP + TR + SV dosages, per cohort and/or pooled.
#
# Usage: Rscript susie_mcv.R <prep_dir> <out_prefix> [L=10] [max_iter=100] [pip_floor=0.01] [mode=both]
#   mode: per_cohort | pooled | both
#
# <prep_dir> is finemap_prep.py's output: <prep_dir>/<gene>/<cohort>/{X.f32.gz, samples.txt, variants.tsv.gz,
# pheno_<cell_type>.tsv.gz} and <prep_dir>/<gene>/pooled_svs.tsv.
#
# per_cohort: for every gene x cohort x cell type
#   1. the cell type's samples and the variants passing associaTR's filters in that sample set,
#   2. missing calls mean-imputed within those samples, constant columns dropped,
#   3. y and every dosage column residualised on [1, covariates] (the SV eQTL model's covariates),
#   4. susie(X, y, L, standardize = TRUE); susie_get_cs(fit, X) (95% coverage, purity min |r| >= 0.5).
#
# pooled: for every gene x cell type present in both cohorts, one model on BioHEART + TOB stacked:
#   - candidates are the variants the meta-analysis could pool: SNPs/TRs (same IDs in both cohorts, jointly
#     called) and federated SV pairs (pooled_svs.tsv), each passing associaTR's filters in BOTH cohorts'
#     sample sets for this cell type, as the meta-analysis requires a locus tested in both;
#   - a federated SV is one column: the BioHEART SV's dosage for BioHEART samples, the TOB SV's for TOB samples;
#   - imputation and residualisation are done WITHIN each cohort, on that cohort's own covariates, and the
#     residuals are then stacked. This removes cohort mean/batch differences and the constant offset between a
#     biallelic DUP (alt count) and a CNV (copy number = 2 + alt count) in dup_cnv pairs, while keeping the
#     per-copy scale. Dosages are deliberately NOT standardised per cohort, so a column means the same per-unit
#     dosage in both cohorts.
#   - IDs duplicated within either cohort's table are dropped (meta_runner_sv.py drop_duplicates(keep = FALSE)).
#   Pooled SV rows carry the pair_id (federate_id, unique) as variant_id plus bioheart_id / tob_id / pair_class.
#
# Writes three gzipped TSVs (cohort = bioheart / tob / pooled):
#   <out_prefix>.variants.tsv.gz  variants with PIP >= pip_floor or in a credible set, and every SV;
#                                 with no signal PIPs spread at ~L/p (~0.001 for 10k variants), hence the 0.01 floor
#   <out_prefix>.cs.tsv.gz        one row per credible set: size, coverage, purity, lead variant, best SV,
#                                 summed alpha by variant type, members
#   <out_prefix>.runs.tsv.gz      one row per gene x cohort x cell type: n, variants by type, converged, n CS
suppressPackageStartupMessages(library(susieR))

args <- commandArgs(trailingOnly = TRUE)
prep_dir <- args[1]
out_prefix <- args[2]
L <- if (length(args) >= 3) as.integer(args[3]) else 10L
max_iter <- if (length(args) >= 4) as.integer(args[4]) else 100L
pip_floor <- if (length(args) >= 5) as.numeric(args[5]) else 0.01
mode <- if (length(args) >= 6) args[6] else "both"
stopifnot(mode %in% c("per_cohort", "pooled", "both"))

VCOLS <- c("variant_id", "variant_type", "chrom", "pos", "end", "motif", "svtype", "overlapping_trs",
           "bioheart_id", "tob_id", "pair_class")

read_X <- function(path, n, p) {
  con <- gzfile(path, "rb")
  on.exit(close(con))
  v <- readBin(con, what = "numeric", size = 4, n = n * p, endian = "little")
  stopifnot(length(v) == n * p)
  matrix(v, nrow = n, ncol = p, byrow = TRUE)   # written row-major by numpy
}

residualise <- function(M, C) {
  Q <- qr.Q(qr(cbind(1, C)))
  M - Q %*% crossprod(Q, M)
}

impute_mean <- function(X) {
  mu <- colMeans(X, na.rm = TRUE)
  miss <- which(is.na(X), arr.ind = TRUE)
  if (nrow(miss)) X[miss] <- mu[miss[, 2]]
  X
}

load_cohort <- function(d) {
  vars <- read.delim(file.path(d, "variants.tsv.gz"), colClasses = c(variant_id = "character", motif = "character",
                     svtype = "character", overlapping_trs = "character"))
  vars$overlapping_trs[is.na(vars$overlapping_trs)] <- ""
  for (col in c("bioheart_id", "tob_id", "pair_class")) vars[[col]] <- NA_character_
  samples <- readLines(file.path(d, "samples.txt"))
  list(vars = vars, samples = samples,
       X = if (nrow(vars)) read_X(file.path(d, "X.f32.gz"), length(samples), nrow(vars)) else NULL)
}

read_pheno <- function(d, ct) {
  ph <- read.delim(file.path(d, paste0("pheno_", ct, ".tsv.gz")), colClasses = c(sample = "character"))
  list(sample = ph$sample, y = ph$y, C = as.matrix(ph[, grep("^c[0-9]+$", names(ph)), drop = FALSE]))
}

# one cohort's residualised block for a cell type: rows = that cell type's samples, columns = `cols` of X
cohort_block <- function(co, ph, cols) {
  rows <- match(ph$sample, co$samples)
  stopifnot(!anyNA(rows))
  X <- impute_mean(co$X[rows, cols, drop = FALSE])
  list(X = residualise(X, ph$C), y = as.vector(residualise(matrix(ph$y), ph$C)), raw_sd = apply(X, 2, sd))
}

var_rows <- list(); cs_rows <- list(); run_rows <- list()

# fit SuSiE on residualised Xr / yr (columns described by v) and append to the three output tables
fit_and_record <- function(Xr, yr, v, gene, cohort, ct, n, note = "") {
  run <- data.frame(gene = gene, cohort = cohort, cell_type = ct, n = n, n_SNP = sum(v$variant_type == "SNP"),
                    n_TR = sum(v$variant_type == "TR"), n_SV = sum(v$variant_type == "SV"), converged = NA,
                    n_cs = 0L, n_cs_with_sv = 0L, sigma2 = NA_real_, pve_top_cs = NA_real_, note = note)
  if (ncol(Xr) < 2) {
    run$note <- "fewer than 2 variants"; run_rows[[length(run_rows) + 1]] <<- run; return(invisible())
  }
  colnames(Xr) <- v$variant_id
  fit <- tryCatch(susie(Xr, yr, L = L, max_iter = max_iter, standardize = TRUE),
                  error = function(e) { run$note <<- conditionMessage(e); NULL })
  if (is.null(fit)) { run_rows[[length(run_rows) + 1]] <<- run; return(invisible()) }
  cs <- susie_get_cs(fit, X = Xr)
  run$converged <- isTRUE(fit$converged); run$sigma2 <- fit$sigma2

  cs_id <- rep(NA_integer_, ncol(Xr))
  if (length(cs$cs)) {
    for (k in seq_along(cs$cs)) {
      idx <- cs$cs[[k]]
      cs_id[idx] <- k
      l <- cs$cs_index[k]
      a <- fit$alpha[l, idx]
      vt <- v$variant_type[idx]
      sv_idx <- idx[vt == "SV"]
      best_sv <- if (length(sv_idx)) sv_idx[which.max(fit$alpha[l, sv_idx])] else NA
      lead <- idx[which.max(a)]
      cs_rows[[length(cs_rows) + 1]] <<- data.frame(
        gene = gene, cohort = cohort, cell_type = ct, cs_id = k, susie_l = l, cs_size = length(idx),
        coverage = cs$coverage[k], min_abs_corr = cs$purity$min.abs.corr[k],
        mean_abs_corr = cs$purity$mean.abs.corr[k],
        lead_variant = v$variant_id[lead], lead_type = v$variant_type[lead], lead_alpha = fit$alpha[l, lead],
        lead_pip = fit$pip[lead],
        alpha_SNP = sum(a[vt == "SNP"]), alpha_TR = sum(a[vt == "TR"]), alpha_SV = sum(a[vt == "SV"]),
        n_SNP = sum(vt == "SNP"), n_TR = sum(vt == "TR"), n_SV = sum(vt == "SV"),
        best_sv = if (is.na(best_sv)) NA else v$variant_id[best_sv],
        best_sv_alpha = if (is.na(best_sv)) NA else fit$alpha[l, best_sv],
        best_sv_pip = if (is.na(best_sv)) NA else fit$pip[best_sv],
        best_sv_overlapping_trs = if (is.na(best_sv)) NA else v$overlapping_trs[best_sv],
        members = paste(v$variant_id[idx], collapse = ","))
    }
    run$n_cs <- length(cs$cs)
    run$n_cs_with_sv <- sum(sapply(cs$cs, function(idx) any(v$variant_type[idx] == "SV")))
    # fit$mu is on the standardised-X scale; divide by the column scale to get per-unit-dosage effects
    run$pve_top_cs <- max(sapply(cs$cs_index, function(l)
      var(as.vector(Xr %*% (fit$alpha[l, ] * fit$mu[l, ] / fit$X_column_scale_factors))) / var(yr)))
  }
  keep <- fit$pip >= pip_floor | !is.na(cs_id) | v$variant_type == "SV"
  if (any(keep)) {
    var_rows[[length(var_rows) + 1]] <<- data.frame(
      gene = gene, cohort = cohort, cell_type = ct, v[keep, VCOLS], pip = fit$pip[keep], cs_id = cs_id[keep],
      post_mean = coef(fit)[-1][keep], post_sd = susie_get_posterior_sd(fit)[keep], row.names = NULL)
  }
  run_rows[[length(run_rows) + 1]] <<- run
  cat(sprintf("%s %s %s: n=%d p=%d (SV %d) CS=%d (with SV %d)\n", gene, cohort, ct, n, ncol(Xr), run$n_SV,
              run$n_cs, run$n_cs_with_sv))
}

run_per_cohort <- function(gene, cohort, co) {
  d <- file.path(prep_dir, gene, cohort)
  for (pf in sort(list.files(d, pattern = "^pheno_.*\\.tsv\\.gz$"))) {
    ct <- sub("^pheno_(.*)\\.tsv\\.gz$", "\\1", pf)
    ph <- read_pheno(d, ct)
    cols <- which(as.logical(co$vars[[paste0("keep_", ct)]]))
    if (length(cols) < 2) {
      fit_and_record(matrix(0, length(ph$y), 0), ph$y, co$vars[0, ], gene, cohort, ct, length(ph$y)); next
    }
    blk <- cohort_block(co, ph, cols)
    ok <- is.finite(blk$raw_sd) & blk$raw_sd > 0
    fit_and_record(blk$X[, ok, drop = FALSE], blk$y, co$vars[cols[ok], ], gene, cohort, ct, length(ph$y))
  }
}

# key -> column index in one cohort's variant table; IDs duplicated in that table are dropped (keep = FALSE)
unique_index <- function(keys) {
  dup <- keys %in% keys[duplicated(keys)]
  setNames(which(!dup), keys[!dup])
}

run_pooled <- function(gene, bh, tob) {
  ps_file <- file.path(prep_dir, gene, "pooled_svs.tsv")
  ps <- if (file.exists(ps_file)) read.delim(ps_file, colClasses = "character") else NULL
  cts_bh <- sub("^pheno_(.*)\\.tsv\\.gz$", "\\1", list.files(file.path(prep_dir, gene, "bioheart"), "^pheno_"))
  cts_tob <- sub("^pheno_(.*)\\.tsv\\.gz$", "\\1", list.files(file.path(prep_dir, gene, "tob"), "^pheno_"))
  for (ct in sort(intersect(cts_bh, cts_tob))) {
    kcol <- paste0("keep_", ct)
    pass <- function(co) as.logical(co$vars[[kcol]])
    # SNPs / TRs: shared IDs (type-prefixed so a SNP and TR can never collide)
    key <- function(co) paste(co$vars$variant_type, co$vars$variant_id, sep = "|")
    ib <- unique_index(key(bh)); it <- unique_index(key(tob))
    st <- intersect(names(ib), names(it))
    st <- st[!startsWith(st, "SV|")]
    st <- st[pass(bh)[ib[st]] & pass(tob)[it[st]]]
    cols_b <- ib[st]; cols_t <- it[st]
    v <- bh$vars[cols_b, ]
    # SVs: federated pairs, matched on upper-cased VARIDs, both sides passing filters in this cell type
    n_sv_pairs <- 0L
    if (!is.null(ps) && nrow(ps)) {
      svb <- unique_index(toupper(bh$vars$variant_id[bh$vars$variant_type == "SV"]))
      svb[] <- which(bh$vars$variant_type == "SV")[svb]
      svt <- unique_index(toupper(tob$vars$variant_id[tob$vars$variant_type == "SV"]))
      svt[] <- which(tob$vars$variant_type == "SV")[svt]
      kb <- toupper(ps$bioheart_id); kt <- toupper(ps$tob_id)
      ok <- kb %in% names(svb) & kt %in% names(svt)
      ok[ok] <- pass(bh)[svb[kb[ok]]] & pass(tob)[svt[kt[ok]]]
      n_sv_pairs <- sum(ok)
      if (n_sv_pairs) {
        sv_b <- svb[kb[ok]]; sv_t <- svt[kt[ok]]
        vs <- bh$vars[sv_b, ]
        vs$variant_id <- ps$pair_id[ok]   # federate_id, made unique (finemap_prep.load_pairs)
        vs$bioheart_id <- ps$bioheart_id[ok]; vs$tob_id <- ps$tob_id[ok]; vs$pair_class <- ps$pair_class[ok]
        # an overlapping TR in either cohort's SV call counts
        ot <- tob$vars$overlapping_trs[sv_t]
        vs$overlapping_trs <- ifelse(vs$overlapping_trs == "", ot, vs$overlapping_trs)
        cols_b <- c(cols_b, sv_b); cols_t <- c(cols_t, sv_t)
        v <- rbind(v, vs)
      }
    }
    pb <- read_pheno(file.path(prep_dir, gene, "bioheart"), ct)
    pt <- read_pheno(file.path(prep_dir, gene, "tob"), ct)
    stopifnot(!length(intersect(pb$sample, pt$sample)))
    n <- length(pb$y) + length(pt$y)
    if (length(cols_b) < 2) { fit_and_record(matrix(0, n, 0), numeric(n), v[0, ], gene, "pooled", ct, n); next }
    blk_b <- cohort_block(bh, pb, cols_b)
    blk_t <- cohort_block(tob, pt, cols_t)
    # a column constant in either cohort's cell-type subset is not testable there (associaTR would crash on it)
    ok <- is.finite(blk_b$raw_sd) & blk_b$raw_sd > 0 & is.finite(blk_t$raw_sd) & blk_t$raw_sd > 0
    fit_and_record(rbind(blk_b$X, blk_t$X)[, ok, drop = FALSE], c(blk_b$y, blk_t$y), v[ok, ], gene, "pooled", ct, n,
                   note = sprintf("%d federated SV pairs usable", sum(ok[v$variant_type == "SV"])))
  }
}

for (gene in sort(list.dirs(prep_dir, recursive = FALSE, full.names = FALSE))) {
  cohorts <- sort(list.dirs(file.path(prep_dir, gene), recursive = FALSE, full.names = FALSE))
  co <- setNames(lapply(cohorts, function(c) load_cohort(file.path(prep_dir, gene, c))), cohorts)
  co <- co[sapply(co, function(x) nrow(x$vars) > 0)]
  if (mode %in% c("per_cohort", "both")) for (c in names(co)) run_per_cohort(gene, c, co[[c]])
  if (mode %in% c("pooled", "both") && all(c("bioheart", "tob") %in% names(co))) run_pooled(gene, co$bioheart, co$tob)
}

write_out <- function(rows, suffix) {
  df <- if (length(rows)) do.call(rbind, rows) else data.frame()
  con <- gzfile(paste0(out_prefix, suffix), "w")
  write.table(df, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  close(con)
}
write_out(var_rows, ".variants.tsv.gz")
write_out(cs_rows, ".cs.tsv.gz")
write_out(run_rows, ".runs.tsv.gz")
