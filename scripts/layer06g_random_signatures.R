## =====================================================================
## Layer 06g - How does the observed signal compare with random gene sets?
##
## Rationale. Layers 06d to 06f established that the published glycolytic
## signatures do not replicate across cohorts and that the surviving
## CGGA-325 signal (Chen 2017, HR 1.35 after age and MGMT adjustment) is
## not explained by any measurable difference in cohort composition. What
## remains untested is the more basic question: is a hazard ratio of that
## size unusual at all for a ten-gene score in this disease?
##
## Transcriptomes in glioblastoma are highly correlated with proliferation
## and with the mesenchymal axis, so a substantial share of arbitrary gene
## sets carries prognostic association. If the observed estimate sits
## inside the distribution produced by random gene sets, the claim changes
## from "this signature failed to replicate" to "this signature is not
## distinguishable from an arbitrary set of genes of comparable
## expression". That is a conceptual result rather than a negative one.
##
## Design.
##
##  1. MATCHED SAMPLING. Genes are not interchangeable: highly expressed,
##     highly variable genes are more likely to associate with outcome
##     than lowly expressed ones. Sampling uniformly from the
##     transcriptome would build a permissive null and would overstate how
##     special the real signature is. Each random gene is therefore drawn
##     from the neighbourhood of a real signature gene in the joint
##     ranking of mean expression and standard deviation, so the random
##     sets match the originals in size and in expression profile.
##
##  2. WEIGHT STRUCTURE PRESERVED. The published coefficients are applied
##     unchanged to the matched random genes, so the number of genes, the
##     sign pattern and the magnitude of the weights are held constant and
##     only gene identity varies. The null therefore isolates the question
##     of which genes were chosen.
##
##  3. SAME MODEL AS THE PRIMARY ANALYSIS. Each random score is
##     standardised within cohort and entered into the age-adjusted Cox
##     model used in Layers 06d to 06f, in the IDH-wildtype WHO grade IV
##     stratum. Nothing else changes.
##
##  4. TWO QUANTITIES ARE REPORTED. First, where the observed hazard ratio
##     falls in the null distribution. Second, and more informative for a
##     reader, the proportion of random sets that reach p < 0.05, which
##     measures how easily a "significant prognostic signature" can be
##     produced in this setting. A third quantity, the increment in
##     Harrell C over the age-only model, is included because the
##     manuscript's practical claim rests on discrimination rather than on
##     significance.
##
##  5. EXCLUSIONS. A signature's own genes are excluded from its null.
##     Genes of the other published signature are retained, because
##     removing every plausibly glycolytic gene would define a
##     non-glycolytic null and answer a different question. This choice is
##     stated in the footnote so the reader can judge it.
##
## Interpretation fixed in advance:
##   - Observed HR below the 95th percentile of the null -> the signature
##     is not distinguishable from matched random gene sets.
##   - Observed HR at or above the 95th percentile -> the specific genes
##     carry information beyond expression-matched chance.
##   - The proportion of random sets with p < 0.05 is reported whatever
##     the result, because it stands on its own.
##
## Usage:
##   setwd("~/Library/CloudStorage/Dropbox/HOPE-GBM"); source("layer06g_random_signatures.R")
## Requires Layers 06d and 06e.
## Runtime is a few minutes: N_PERM x 2 cohorts x 2 signatures Cox fits.
## =====================================================================

set.seed(1024)
suppressPackageStartupMessages({
  library(survival); library(dplyr); library(ggplot2)
  library(flextable); library(officer); library(readr)
})

PROJ    <- getwd()
DIR_RES <- file.path(PROJ, "results")
DIR_OUT <- file.path(PROJ, "Sep_06")
DIR_TAB <- file.path(DIR_OUT, "Tables")
DIR_FIG <- file.path(DIR_OUT, "Figures")
FIG_DPI <- 600
N_PERM  <- 1000          # resolution of the percentile is 1/N_PERM
N_NEIGH <- 100           # candidate pool per signature gene, in rank space

gen  <- readRDS(file.path(DIR_RES, "generality_results.rds"))
disc <- readRDS(file.path(DIR_RES, "discordance_results.rds"))
SIGS <- gen$signatures
G4   <- disc$G4
CGGA_DIR <- Find(function(p) dir.exists(p) &&
                   length(list.files(p, pattern = "mRNAseq.*genes.*[.]txt$")) > 0,
                 c(path.expand("~/Documents/CGGA_local"),
                   path.expand("~/Library/CloudStorage/Dropbox/00_DATA/CGGA"),
                   file.path(PROJ, "data_raw", "CGGA")))
stopifnot(!is.null(CGGA_DIR))

log_msg <- function(...) {
  m <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ...)
  cat(m, "\n"); cat(m, "\n", file = file.path(PROJ, "logs", "hope_run.log"), append = TRUE)
}
zvec <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
ft_style <- function(ft) {
  ft %>% flextable::font(fontname = "Helvetica", part = "all") %>%
    flextable::fontsize(size = 10, part = "all") %>% flextable::border_remove() %>%
    flextable::hline_top(part = "header", border = officer::fp_border(width = 0.5)) %>%
    flextable::hline_bottom(part = "header", border = officer::fp_border(width = 0.5)) %>%
    flextable::hline_bottom(part = "body", border = officer::fp_border(width = 0.5)) %>%
    flextable::align(align = "center", part = "all") %>%
    flextable::align(j = 1, align = "left", part = "all") %>% flextable::autofit()
}
save_fig <- function(obj, stem, width, height, dpi = FIG_DPI) {
  png(paste0(stem, ".png"), width = width, height = height, units = "in", res = dpi)
  print(obj); dev.off()
  tiff(paste0(stem, ".tif"), width = width, height = height, units = "in", res = dpi)
  print(obj); dev.off()
}
read_cgga_expr <- function(path) {
  con <- file(path, open = "rb")
  ln  <- readLines(con, skipNul = TRUE, warn = FALSE, encoding = "latin1"); close(con)
  ln  <- gsub("\r$", "", ln, useBytes = TRUE); ln <- ln[nzchar(ln, keepNA = FALSE)]
  fl  <- strsplit(ln, "\t", fixed = TRUE, useBytes = TRUE); nf <- lengths(fl)
  nc  <- as.integer(names(sort(table(nf[-1]), decreasing = TRUE))[1])
  hdr <- fl[[1]]; ids <- if (length(hdr) == nc) hdr[-1] else hdr
  body <- fl[-1][nf[-1] == nc]
  gene <- enc2utf8(vapply(body, `[`, "", 1L))
  expr <- matrix(suppressWarnings(as.numeric(unlist(lapply(body, `[`, -1L), use.names = FALSE))),
                 nrow = length(body), ncol = length(ids), byrow = TRUE, dimnames = list(NULL, ids))
  ok <- !is.na(gene) & nzchar(gene) & rowSums(!is.na(expr)) > 0
  expr <- expr[ok, , drop = FALSE]; gene <- gene[ok]
  ord <- order(rowMeans(expr, na.rm = TRUE), decreasing = TRUE)
  expr <- expr[ord, , drop = FALSE]; gene <- gene[ord]
  dup <- duplicated(gene); expr <- expr[!dup, , drop = FALSE]; rownames(expr) <- gene[!dup]
  expr
}
num_of <- function(p) regmatches(basename(p), regexpr("325|693", basename(p)))

expr_f <- list.files(CGGA_DIR, pattern = "mRNAseq.*genes.*[.]txt$", full.names = TRUE)
EXPR <- setNames(lapply(expr_f, read_cgga_expr), paste0("CGGA-", vapply(expr_f, num_of, "")))

## ---- Observed estimates, recomputed here so the comparison is exact ----------
observed <- function(cohort, sig) {
  d  <- G4[G4$cohort == cohort, ]
  lg <- log2(EXPR[[cohort]] + 1)[, d$patient, drop = FALSE]
  w  <- SIGS[[sig]]; gs <- intersect(names(w), rownames(lg))
  Z  <- t(apply(lg[gs, , drop = FALSE], 1, zvec))
  d$score <- zvec(as.numeric(colSums(Z * w[gs])))
  f0 <- coxph(Surv(OS.time, OS.event) ~ age, data = d)
  f1 <- coxph(Surv(OS.time, OS.event) ~ score + age, data = d)
  sm <- summary(f1)
  list(d = d, lg = lg, w = w, genes_used = gs,
       HR = unname(sm$coefficients["score", "exp(coef)"]),
       p  = unname(sm$coefficients["score", "Pr(>|z|)"]),
       c0 = unname(summary(f0)$concordance[1]),
       c1 = unname(sm$concordance[1]))
}

## ---- Expression- and variance-matched random gene sets ----------------------
## For each real signature gene, candidates are the N_NEIGH genes closest in the
## joint ranking of mean expression and standard deviation, excluding the
## signature's own genes. One candidate is drawn per real gene, without
## replacement across the set, so each random set has the same size.
make_matcher <- function(lg, exclude) {
  mu <- rowMeans(lg, na.rm = TRUE); sdv <- apply(lg, 1, sd, na.rm = TRUE)
  keep <- is.finite(mu) & is.finite(sdv) & sdv > 0
  mu <- mu[keep]; sdv <- sdv[keep]
  r_mu <- rank(mu, ties.method = "average"); r_sd <- rank(sdv, ties.method = "average")
  pool <- setdiff(names(mu), exclude)
  function(target) {
    if (!(target %in% names(r_mu))) return(NULL)
    dist <- sqrt((r_mu[pool] - r_mu[[target]])^2 + (r_sd[pool] - r_sd[[target]])^2)
    names(sort(dist))[seq_len(min(N_NEIGH, length(dist)))]
  }
}

run_null <- function(cohort, sig) {
  ob <- observed(cohort, sig)
  d  <- ob$d; lg <- ob$lg; w <- ob$w[ob$genes_used]
  match_fn <- make_matcher(lg, exclude = ob$genes_used)
  cand <- lapply(ob$genes_used, match_fn)
  ok   <- !vapply(cand, is.null, TRUE)
  cand <- cand[ok]; w_use <- w[ok]
  if (!length(cand)) return(NULL)
  ## Pre-standardise the whole matrix once; the same z-scores are reused.
  Zall <- t(apply(lg, 1, zvec))
  base <- coxph(Surv(OS.time, OS.event) ~ age, data = d)
  c0   <- unname(summary(base)$concordance[1])
  out <- vapply(seq_len(N_PERM), function(i) {
    gs <- character(length(cand))
    for (j in seq_along(cand)) {
      pool_j <- setdiff(cand[[j]], gs)
      gs[j]  <- if (length(pool_j)) sample(pool_j, 1) else sample(cand[[j]], 1)
    }
    sc <- zvec(as.numeric(colSums(Zall[gs, , drop = FALSE] * w_use)))
    if (!all(is.finite(sc))) return(c(NA, NA, NA))
    dd <- d; dd$score <- sc
    f <- tryCatch(coxph(Surv(OS.time, OS.event) ~ score + age, data = dd),
                  error = function(e) NULL)
    if (is.null(f)) return(c(NA, NA, NA))
    sm <- summary(f)
    c(unname(sm$coefficients["score", "exp(coef)"]),
      unname(sm$coefficients["score", "Pr(>|z|)"]),
      unname(sm$concordance[1]) - c0)
  }, numeric(3))
  HRs <- out[1, ]; ps <- out[2, ]; dC <- out[3, ]
  valid <- is.finite(HRs) & is.finite(ps)
  HRs <- HRs[valid]; ps <- ps[valid]; dC <- dC[valid]
  ## Directional percentile: how many random sets reached a hazard ratio at
  ## least as extreme as the observed one, in the observed direction.
  pct <- if (ob$HR >= 1) mean(HRs >= ob$HR) else mean(HRs <= ob$HR)
  list(cohort = cohort, signature = sig, n_valid = length(HRs),
       obs_HR = ob$HR, obs_p = ob$p, obs_dC = ob$c1 - ob$c0,
       null_HR = HRs, null_p = ps, null_dC = dC,
       null_median = median(HRs), null_q95 = quantile(HRs, 0.95),
       pct_more_extreme = pct, prop_sig = mean(ps < 0.05),
       prop_dC_above = mean(dC >= (ob$c1 - ob$c0)))
}

COMB <- expand.grid(cohort = levels(G4$cohort), sig = names(SIGS),
                    stringsAsFactors = FALSE)
RES <- Map(function(k, s) {
  log_msg(sprintf("06g: running %d matched random sets for %s / %s", N_PERM, k, s))
  run_null(k, s)
}, COMB$cohort, COMB$sig)
RES <- Filter(Negate(is.null), RES)

SUM <- do.call(rbind, lapply(RES, function(r) data.frame(
  cohort = r$cohort, signature = r$signature, n_random = r$n_valid,
  obs_HR = r$obs_HR, obs_p = r$obs_p, obs_dC = r$obs_dC,
  null_median_HR = r$null_median, null_q95_HR = unname(r$null_q95),
  pct_more_extreme = r$pct_more_extreme, prop_random_sig = r$prop_sig,
  prop_random_dC_ge_obs = r$prop_dC_above, row.names = NULL)))
print(SUM)
for (i in seq_len(nrow(SUM))) {
  r <- SUM[i, ]
  log_msg(sprintf("06g %s / %s: observed HR %.2f | null median %.2f, 95th pct %.2f | %.1f%% of random sets as extreme | %.1f%% of random sets p<0.05",
                  r$cohort, r$signature, r$obs_HR, r$null_median_HR, r$null_q95_HR,
                  100 * r$pct_more_extreme, 100 * r$prop_random_sig))
}

## ---- Table 13 ------------------------------------------------------------------
t13 <- data.frame(
  Cohort = SUM$cohort, Signature = SUM$signature,
  `Observed HR (p)` = sprintf("%.2f (p=%s)", SUM$obs_HR,
                              ifelse(SUM$obs_p < 0.001, "<0.001", sprintf("%.3f", SUM$obs_p))),
  `Random sets, median HR` = sprintf("%.2f", SUM$null_median_HR),
  `Random sets, 95th percentile HR` = sprintf("%.2f", SUM$null_q95_HR),
  `Random sets at least as extreme` = sprintf("%.1f%%", 100 * SUM$pct_more_extreme),
  `Random sets with p < 0.05` = sprintf("%.1f%%", 100 * SUM$prop_random_sig),
  `Random sets matching observed C increment` = sprintf("%.1f%%", 100 * SUM$prop_random_dC_ge_obs),
  check.names = FALSE)
ft13 <- flextable::flextable(t13) %>% ft_style() %>%
  flextable::add_footer_lines(values = c(
    sprintf("For each signature, %d random gene sets were generated in each cohort. Every random gene was drawn from the neighbourhood of a real signature gene in the joint ranking of mean expression and standard deviation, so that random and published sets are comparable in size and expression profile; the published coefficients were applied unchanged, leaving gene identity as the only quantity that varies.", N_PERM),
    "Each random score was standardised within cohort and entered into the same age-adjusted Cox model used for the published signatures, restricted to IDH-wildtype WHO grade IV gliomas.",
    "The proportion of random sets reaching p < 0.05 indicates how readily a statistically significant prognostic score can be obtained in this setting, irrespective of the genes selected. The final column gives the proportion of random sets whose increment in Harrell C over the age-only model equalled or exceeded that of the published signature.",
    "A signature's own genes were excluded from its null; genes belonging to the other published signature were retained, since removing all plausibly glycolytic genes would define a different comparison.",
    "C, concordance; HR, hazard ratio."))
flextable::save_as_docx(ft13, path = file.path(DIR_TAB, "Table13_random_signature_null.docx"))

## ---- Figure 9: null distributions with the observed estimate ------------------
pd <- do.call(rbind, lapply(RES, function(r) data.frame(
  cohort = r$cohort, signature = r$signature, HR = r$null_HR)))
pv <- do.call(rbind, lapply(RES, function(r) data.frame(
  cohort = r$cohort, signature = r$signature, HR = r$obs_HR,
  q95 = unname(r$null_q95))))
g9 <- ggplot(pd, aes(x = HR)) +
  geom_histogram(bins = 45, fill = "grey80", colour = "grey40", linewidth = 0.2) +
  geom_vline(xintercept = 1, linetype = "dotted", colour = "grey50") +
  geom_vline(data = pv, aes(xintercept = q95), linetype = "dashed", colour = "grey30") +
  geom_vline(data = pv, aes(xintercept = HR), linewidth = 0.8) +
  facet_grid(signature ~ cohort, scales = "free_y") +
  scale_x_log10() +
  labs(x = "HR per 1 SD, age-adjusted, IDH-wildtype WHO grade IV",
       y = "Random gene sets") +
  theme_classic(base_size = 11) +
  theme(strip.background = element_blank())
save_fig(g9, file.path(DIR_FIG, "Figure9_random_signature_null"), width = 7.2, height = 5.4)

## ---- Persist ---------------------------------------------------------------------
saveRDS(list(summary = SUM, detail = RES, n_perm = N_PERM, n_neigh = N_NEIGH),
        file.path(DIR_RES, "random_signature_null.rds"))
write.csv(SUM, file.path(DIR_RES, "table13_random_signature_null.csv"), row.names = FALSE)

nums <- c("", sprintf("--- Layer 06g: %d expression-matched random gene sets per cohort and signature ---", N_PERM),
  sprintf("%s / %s: observed HR %.2f (p=%.3f); null median %.2f, 95th percentile %.2f; %.1f%% of random sets at least as extreme; %.1f%% of random sets reached p<0.05; observed C increment %.3f, matched or exceeded by %.1f%% of random sets",
          SUM$cohort, SUM$signature, SUM$obs_HR, SUM$obs_p, SUM$null_median_HR,
          SUM$null_q95_HR, 100 * SUM$pct_more_extreme, 100 * SUM$prop_random_sig,
          SUM$obs_dC, 100 * SUM$prop_random_dC_ge_obs))
cat(nums, sep = "\n")
write(nums, file = file.path(DIR_OUT, "manuscript_numbers.txt"), append = TRUE)
log_msg("Layer 06g complete: Table 13 and Figure 9 exported to Sep_06/.")
