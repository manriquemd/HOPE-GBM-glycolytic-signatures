## =====================================================================
## Layer 06h - Calibration check on the random-signature null (Layer 06g)
##
## The problem. Layer 06g produced a type I error rate that varies far more
## across cells than a correctly specified null should:
##
##   CGGA-325 / Chen 2017   8.7% of random sets reached p < 0.05
##   CGGA-325 / Zhang 2021  7.3%
##   CGGA-693 / Zhang 2021  5.6%
##   CGGA-693 / Chen 2017   0.3%
##
## Under a well-behaved null the expectation is approximately 5% in every
## cell. A value of 0.3% in one cell of the same design, alongside 8.7% in
## another, indicates that the sampling procedure itself behaves
## differently depending on cohort and signature, not that the biology
## differs. Table 13 cannot be released until this is resolved, because a
## reviewer familiar with permutation methods will ask why the nominal
## error rate is not reproduced.
##
## Candidate explanations, each of which this layer measures directly:
##
##  1. VARIANCE COLLAPSE IN THE RANDOM SCORE. Chen 2017 carries three
##     negative coefficients. If the matched neighbourhoods of genes with
##     opposite signs overlap, a random set can place correlated genes on
##     opposite sides of the weighted sum, cancelling signal and producing
##     scores with little variance. A score with compressed variance,
##     once standardised, is dominated by noise and rarely reaches
##     significance. This would make the 0.3% an artefact of construction.
##
##  2. NEIGHBOURHOOD OVERLAP. If the candidate pools of the ten signature
##     genes intersect heavily, the effective diversity of random sets is
##     far lower than the nominal one thousand, and repeated genes appear
##     across draws.
##
##  3. INTERNAL CORRELATION OF THE RANDOM SETS. Genes matched on mean and
##     standard deviation are not independent; matching on expression may
##     inadvertently select co-expressed modules, and the degree to which
##     this happens can differ between cohorts.
##
##  4. GENUINE COHORT DIFFERENCE. CGGA-693 may simply show weaker coupling
##     between transcriptome and outcome, in which case the low rate is a
##     finding rather than a defect. This is distinguishable from the
##     others: if the random scores in that cell have normal variance and
##     normal internal correlation, and the rate is still low, the
##     explanation is the cohort.
##
## What this layer does NOT do. It does not change the Layer 06g result.
## It measures the properties of the null that Layer 06g already generated
## and, where a defect is found, re-runs a corrected variant so the two can
## be compared. If the corrected null reproduces the original, the original
## stands.
##
## Decision rule fixed in advance:
##   - Median random-score SD in the anomalous cell below 60% of the other
##     cells -> variance collapse; the corrected null replaces Table 13.
##   - SD comparable and Jaccard overlap of candidate pools comparable ->
##     the low rate is a property of CGGA-693 and is reported as such.
##   - Corrected and original nulls agree -> original retained, this layer
##     is cited in the supplement as a calibration check.
##
## Usage:
##   setwd("~/Library/CloudStorage/Dropbox/HOPE-GBM"); source("layer06h_null_calibration.R")
## Requires Layer 06g (results/random_signature_null.rds).
## =====================================================================

set.seed(1024)
suppressPackageStartupMessages({
  library(survival); library(dplyr); library(ggplot2)
  library(flextable); library(officer)
})

PROJ    <- getwd()
DIR_RES <- file.path(PROJ, "results")
DIR_OUT <- file.path(PROJ, "Sep_06")
DIR_TAB <- file.path(DIR_OUT, "Tables")
DIR_FIG <- file.path(DIR_OUT, "Figures")
FIG_DPI <- 600
N_PERM  <- 1000
N_NEIGH <- 100

gen  <- readRDS(file.path(DIR_RES, "generality_results.rds"))
disc <- readRDS(file.path(DIR_RES, "discordance_results.rds"))
g06g <- readRDS(file.path(DIR_RES, "random_signature_null.rds"))
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

jaccard <- function(a, b) length(intersect(a, b)) / length(union(a, b))

## ---- Core: rebuild the null while recording its internal properties ----------
## mode = "matched"     reproduces Layer 06g exactly (sampling without
##                      replacement across signature genes)
## mode = "independent" draws each gene from its own pool without the
##                      cross-gene exclusion, which is the suspected source
##                      of variance collapse
## mode = "signflip"    keeps the real signature genes and permutes the sign
##                      pattern of the weights; isolates the contribution of
##                      gene identity from that of weight structure
probe_cell <- function(cohort, sig, mode = "matched") {
  d  <- G4[G4$cohort == cohort, ]
  lg <- log2(EXPR[[cohort]] + 1)[, d$patient, drop = FALSE]
  w0 <- SIGS[[sig]]; gs0 <- intersect(names(w0), rownames(lg)); w0 <- w0[gs0]
  mu <- rowMeans(lg, na.rm = TRUE); sdv <- apply(lg, 1, sd, na.rm = TRUE)
  keep <- is.finite(mu) & is.finite(sdv) & sdv > 0
  mu <- mu[keep]; sdv <- sdv[keep]
  r_mu <- rank(mu, ties.method = "average"); r_sd <- rank(sdv, ties.method = "average")
  pool_all <- setdiff(names(mu), gs0)
  pools <- lapply(gs0, function(g) {
    dd <- sqrt((r_mu[pool_all] - r_mu[[g]])^2 + (r_sd[pool_all] - r_sd[[g]])^2)
    names(sort(dd))[seq_len(min(N_NEIGH, length(dd)))]
  })
  names(pools) <- gs0
  ## pairwise overlap of candidate pools, a direct measure of explanation 2
  jm <- outer(seq_along(pools), seq_along(pools),
              Vectorize(function(i, j) if (i < j) jaccard(pools[[i]], pools[[j]]) else NA_real_))
  Zall <- t(apply(lg, 1, zvec))
  base <- coxph(Surv(OS.time, OS.event) ~ age, data = d)
  c0 <- unname(summary(base)$concordance[1])
  draw <- function() {
    if (mode == "signflip") {
      list(genes = gs0, w = w0 * sample(c(-1, 1), length(w0), replace = TRUE))
    } else {
      g <- character(length(pools))
      for (j in seq_along(pools)) {
        p_j <- if (mode == "matched") setdiff(pools[[j]], g) else pools[[j]]
        g[j] <- if (length(p_j)) sample(p_j, 1) else sample(pools[[j]], 1)
      }
      list(genes = g, w = w0)
    }
  }
  out <- vapply(seq_len(N_PERM), function(i) {
    dr <- draw()
    sub <- Zall[dr$genes, , drop = FALSE]
    raw <- as.numeric(colSums(sub * dr$w))
    sd_raw <- sd(raw, na.rm = TRUE)
    ## mean absolute pairwise correlation among the selected genes
    cm <- suppressWarnings(cor(t(sub), use = "pairwise.complete.obs"))
    mabs <- mean(abs(cm[upper.tri(cm)]), na.rm = TRUE)
    if (!is.finite(sd_raw) || sd_raw <= 0) return(c(NA, NA, NA, NA, NA))
    dd <- d; dd$score <- zvec(raw)
    f <- tryCatch(coxph(Surv(OS.time, OS.event) ~ score + age, data = dd),
                  error = function(e) NULL)
    if (is.null(f)) return(c(NA, NA, NA, NA, NA))
    sm <- summary(f)
    c(unname(sm$coefficients["score", "exp(coef)"]),
      unname(sm$coefficients["score", "Pr(>|z|)"]),
      unname(sm$concordance[1]) - c0, sd_raw, mabs)
  }, numeric(5))
  ok <- is.finite(out[1, ]) & is.finite(out[2, ])
  list(cohort = cohort, signature = sig, mode = mode, n_valid = sum(ok),
       HR = out[1, ok], p = out[2, ok], dC = out[3, ok],
       sd_raw = out[4, ok], mean_abs_cor = out[5, ok],
       pool_jaccard_median = median(jm, na.rm = TRUE),
       pool_jaccard_max = suppressWarnings(max(jm, na.rm = TRUE)),
       n_genes = length(gs0), n_neg_weights = sum(w0 < 0),
       sd_observed = sd(as.numeric(colSums(Zall[gs0, , drop = FALSE] * w0)), na.rm = TRUE))
}

CELLS <- expand.grid(cohort = levels(G4$cohort), sig = names(SIGS), stringsAsFactors = FALSE)

## ---- Step 1. Reproduce the original null and record its properties ----------
MATCH <- Map(function(k, s) {
  log_msg(sprintf("06h: probing matched null, %s / %s", k, s))
  probe_cell(k, s, "matched")
}, CELLS$cohort, CELLS$sig)

prop_tab <- do.call(rbind, lapply(MATCH, function(r) data.frame(
  cohort = r$cohort, signature = r$signature,
  n_genes = r$n_genes, n_negative_weights = r$n_neg_weights,
  sd_observed_score = r$sd_observed,
  sd_random_median = median(r$sd_raw), sd_random_q10 = quantile(r$sd_raw, 0.10),
  sd_ratio_vs_observed = median(r$sd_raw) / r$sd_observed,
  mean_abs_cor_random = median(r$mean_abs_cor),
  pool_jaccard_median = r$pool_jaccard_median, pool_jaccard_max = r$pool_jaccard_max,
  prop_sig = mean(r$p < 0.05), row.names = NULL)))
print(prop_tab)
for (i in seq_len(nrow(prop_tab))) {
  r <- prop_tab[i, ]
  log_msg(sprintf("06h PROPS %s / %s: random score SD median %.2f (observed %.2f, ratio %.2f); mean |r| among genes %.2f; pool Jaccard median %.3f max %.3f; %.1f%% significant",
                  r$cohort, r$signature, r$sd_random_median, r$sd_observed_score,
                  r$sd_ratio_vs_observed, r$mean_abs_cor_random, r$pool_jaccard_median,
                  r$pool_jaccard_max, 100 * r$prop_sig))
}

## ---- Step 2. Apply the pre-specified rule -----------------------------------
ref_sd <- median(prop_tab$sd_random_median)
prop_tab$sd_vs_cells <- prop_tab$sd_random_median / ref_sd
collapsed <- prop_tab$sd_vs_cells < 0.60
if (any(collapsed)) {
  for (i in which(collapsed))
    log_msg(sprintf("06h FLAG: variance collapse in %s / %s (SD %.0f%% of the cross-cell median)",
                    prop_tab$cohort[i], prop_tab$signature[i], 100 * prop_tab$sd_vs_cells[i]))
} else {
  log_msg("06h: no variance collapse detected; random score SD comparable across all cells")
}

## ---- Step 3. Corrected null, drawn independently per gene -------------------
INDEP <- Map(function(k, s) {
  log_msg(sprintf("06h: probing independent-draw null, %s / %s", k, s))
  probe_cell(k, s, "independent")
}, CELLS$cohort, CELLS$sig)

## ---- Step 4. Sign-flip null, isolating weight structure from gene identity ---
FLIP <- Map(function(k, s) {
  log_msg(sprintf("06h: probing sign-flip null, %s / %s", k, s))
  probe_cell(k, s, "signflip")
}, CELLS$cohort, CELLS$sig)

obs_hr <- setNames(g06g$summary$obs_HR,
                   paste(g06g$summary$cohort, g06g$summary$signature))
compare <- do.call(rbind, lapply(seq_along(MATCH), function(i) {
  key <- paste(MATCH[[i]]$cohort, MATCH[[i]]$signature)
  oh  <- unname(obs_hr[key])
  f <- function(r) {
    pct <- if (oh >= 1) mean(r$HR >= oh) else mean(r$HR <= oh)
    c(median_HR = median(r$HR), q95 = unname(quantile(r$HR, 0.95)),
      prop_sig = mean(r$p < 0.05), pct_extreme = pct)
  }
  m <- f(MATCH[[i]]); n <- f(INDEP[[i]]); s <- f(FLIP[[i]])
  data.frame(cohort = MATCH[[i]]$cohort, signature = MATCH[[i]]$signature,
             observed_HR = oh,
             matched_median = m["median_HR"], matched_q95 = m["q95"],
             matched_sig = m["prop_sig"], matched_pct = m["pct_extreme"],
             indep_median = n["median_HR"], indep_q95 = n["q95"],
             indep_sig = n["prop_sig"], indep_pct = n["pct_extreme"],
             flip_median = s["median_HR"], flip_q95 = s["q95"],
             flip_sig = s["prop_sig"], flip_pct = s["pct_extreme"], row.names = NULL)
}))
print(compare)
for (i in seq_len(nrow(compare))) {
  r <- compare[i, ]
  log_msg(sprintf("06h COMPARE %s / %s: observed HR %.2f | matched pct %.1f%% (sig %.1f%%) | independent pct %.1f%% (sig %.1f%%) | sign-flip pct %.1f%% (sig %.1f%%)",
                  r$cohort, r$signature, r$observed_HR, 100 * r$matched_pct,
                  100 * r$matched_sig, 100 * r$indep_pct, 100 * r$indep_sig,
                  100 * r$flip_pct, 100 * r$flip_sig))
}

## ---- Verdict -------------------------------------------------------------------
verdict <- vapply(seq_len(nrow(compare)), function(i) {
  d_pct <- abs(compare$matched_pct[i] - compare$indep_pct[i])
  d_sig <- abs(compare$matched_sig[i] - compare$indep_sig[i])
  if (collapsed[i] && (d_pct > 0.05 || d_sig > 0.03))
    "Variance collapse confirmed; independent-draw null replaces the matched null"
  else if (collapsed[i])
    "Variance reduced but conclusions unchanged; matched null retained with a caveat"
  else if (d_pct <= 0.05 && d_sig <= 0.03)
    "Null well behaved; matched and independent draws agree"
  else
    "Sampling scheme materially affects the null; report both and prefer the independent draw"
}, "")
compare$verdict <- verdict
for (i in seq_len(nrow(compare)))
  log_msg(sprintf("06h VERDICT %s / %s: %s", compare$cohort[i], compare$signature[i], verdict[i]))

## ---- Tables and figure -----------------------------------------------------------
tS1 <- data.frame(
  Cohort = prop_tab$cohort, Signature = prop_tab$signature,
  `Genes (negative weights)` = sprintf("%d (%d)", prop_tab$n_genes, prop_tab$n_negative_weights),
  `SD of published score` = sprintf("%.2f", prop_tab$sd_observed_score),
  `Median SD of random scores` = sprintf("%.2f", prop_tab$sd_random_median),
  `Ratio to published` = sprintf("%.2f", prop_tab$sd_ratio_vs_observed),
  `Mean absolute correlation within random sets` = sprintf("%.2f", prop_tab$mean_abs_cor_random),
  `Median pool overlap (Jaccard)` = sprintf("%.3f", prop_tab$pool_jaccard_median),
  `Random sets with p < 0.05` = sprintf("%.1f%%", 100 * prop_tab$prop_sig),
  check.names = FALSE)
ftS1 <- flextable::flextable(tS1) %>% ft_style() %>%
  flextable::add_footer_lines(values = c(
    "Internal properties of the random gene sets generated in Layer 06g, reported so that the behaviour of the null can be judged independently of its result. The standard deviation of the random scores is given before standardisation; a value much smaller than that of the published score indicates cancellation between genes carrying weights of opposite sign.",
    "Pool overlap is the median Jaccard index between the candidate neighbourhoods of the signature genes, and measures how far the nominal number of random sets overstates their effective diversity.",
    "Under a correctly specified null approximately 5% of random sets are expected to reach p < 0.05 in each cell.",
    "SD, standard deviation."))
flextable::save_as_docx(ftS1, path = file.path(DIR_TAB, "TableS1_null_calibration.docx"))

tS2 <- data.frame(
  Cohort = compare$cohort, Signature = compare$signature,
  `Observed HR` = sprintf("%.2f", compare$observed_HR),
  `Matched draw: percentile (p<0.05 rate)` = sprintf("%.1f%% (%.1f%%)",
    100 * compare$matched_pct, 100 * compare$matched_sig),
  `Independent draw: percentile (p<0.05 rate)` = sprintf("%.1f%% (%.1f%%)",
    100 * compare$indep_pct, 100 * compare$indep_sig),
  `Sign-flip: percentile (p<0.05 rate)` = sprintf("%.1f%% (%.1f%%)",
    100 * compare$flip_pct, 100 * compare$flip_sig),
  check.names = FALSE)
ftS2 <- flextable::flextable(tS2) %>% ft_style() %>%
  flextable::add_footer_lines(values = c(
    "Three null constructions compared in each cell. The matched draw is the procedure of Layer 06g, in which genes are sampled without replacement across the signature. The independent draw removes the cross-gene exclusion. The sign-flip null retains the published genes and permutes the signs of their coefficients, isolating the contribution of gene identity from that of weight structure.",
    "The percentile is the proportion of random sets reaching a hazard ratio at least as extreme as the observed one in the observed direction; the figure in parentheses is the proportion reaching p < 0.05.",
    "Agreement between the first two columns indicates that the conclusions of Layer 06g do not depend on the sampling scheme.",
    "HR, hazard ratio."))
flextable::save_as_docx(ftS2, path = file.path(DIR_TAB, "TableS2_null_constructions.docx"))

pd <- do.call(rbind, lapply(seq_along(MATCH), function(i) rbind(
  data.frame(cohort = MATCH[[i]]$cohort, signature = MATCH[[i]]$signature,
             construction = "Matched draw", HR = MATCH[[i]]$HR),
  data.frame(cohort = INDEP[[i]]$cohort, signature = INDEP[[i]]$signature,
             construction = "Independent draw", HR = INDEP[[i]]$HR),
  data.frame(cohort = FLIP[[i]]$cohort, signature = FLIP[[i]]$signature,
             construction = "Sign-flip", HR = FLIP[[i]]$HR))))
pv <- data.frame(cohort = compare$cohort, signature = compare$signature,
                 HR = compare$observed_HR)
gS <- ggplot(pd, aes(x = HR, colour = construction)) +
  geom_density(linewidth = 0.6) +
  geom_vline(xintercept = 1, linetype = "dotted", colour = "grey60") +
  geom_vline(data = pv, aes(xintercept = HR), linewidth = 0.7, colour = "black") +
  facet_grid(signature ~ cohort) + scale_x_log10() +
  scale_colour_grey(start = 0.1, end = 0.7) +
  labs(x = "HR per 1 SD, age-adjusted, IDH-wildtype WHO grade IV",
       y = "Density", colour = NULL) +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom", strip.background = element_blank())
save_fig(gS, file.path(DIR_FIG, "FigureS1_null_constructions"), width = 7.4, height = 5.4)

## ---- Persist -----------------------------------------------------------------------
saveRDS(list(properties = prop_tab, comparison = compare,
             matched = MATCH, independent = INDEP, signflip = FLIP),
        file.path(DIR_RES, "null_calibration.rds"))
write.csv(prop_tab, file.path(DIR_RES, "tableS1_null_calibration.csv"), row.names = FALSE)
write.csv(compare,  file.path(DIR_RES, "tableS2_null_constructions.csv"), row.names = FALSE)

nums <- c("", "--- Layer 06h: calibration of the random-signature null ---",
  sprintf("%s / %s: random score SD median %.2f vs published %.2f (ratio %.2f); mean |r| within random sets %.2f; pool Jaccard median %.3f; %.1f%% of matched random sets significant",
          prop_tab$cohort, prop_tab$signature, prop_tab$sd_random_median,
          prop_tab$sd_observed_score, prop_tab$sd_ratio_vs_observed,
          prop_tab$mean_abs_cor_random, prop_tab$pool_jaccard_median,
          100 * prop_tab$prop_sig),
  "",
  sprintf("%s / %s: observed HR %.2f | matched %.1f%% extreme (%.1f%% sig) | independent %.1f%% extreme (%.1f%% sig) | sign-flip %.1f%% extreme (%.1f%% sig) | %s",
          compare$cohort, compare$signature, compare$observed_HR,
          100 * compare$matched_pct, 100 * compare$matched_sig,
          100 * compare$indep_pct, 100 * compare$indep_sig,
          100 * compare$flip_pct, 100 * compare$flip_sig, compare$verdict))
cat(nums, sep = "\n")
write(nums, file = file.path(DIR_OUT, "manuscript_numbers.txt"), append = TRUE)
log_msg("Layer 06h complete: Tables S1-S2 and Figure S1 exported to Sep_06/.")
