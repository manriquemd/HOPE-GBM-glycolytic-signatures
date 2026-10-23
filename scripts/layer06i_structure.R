## =====================================================================
## Layer 06i - Is the dominant expression structure of CGGA-693 technical?
##
## The finding to be explained. Layer 06h showed that ten genes drawn at
## random from the expression neighbourhood of a signature gene are far
## more correlated with one another in CGGA-693 (mean |r| = 0.47 to 0.48)
## than in CGGA-325 (0.17 to 0.18). This collapses the diversity of any
## ten-gene score in CGGA-693: nearly every random set measures the same
## underlying axis, the null distribution becomes narrow (95th percentile
## HR 1.10 against 1.29), and almost no random set reaches significance
## (0.3% against 8.7%).
##
## Two readings are possible and they lead to different manuscripts.
##
##   Biological. CGGA-693 is dominated by a real transcriptional programme
##   shared across tumours. The discordance is then a genuine difference
##   between patient populations and belongs in Results.
##
##   Technical. The dominant axis reflects sequencing depth, normalisation
##   or batch. Layer 06e already recorded a difference in the library-size
##   proxy (median 0.7 against 0.9 million, p < 0.001), which makes this
##   reading plausible. The discordance would then be partly an artefact
##   and belongs in Limitations, stated before a reviewer states it.
##
## The two are distinguishable. A technical axis correlates with depth and
## with the per-sample detection rate; a biological axis correlates with
## known tumour phenotypes and not with depth.
##
## What is measured.
##
##  1. Variance explained by the leading principal components in each
##     cohort, computed on the same gene space so the comparison is fair.
##  2. Correlation of PC1 and PC2 with the library-size proxy and with the
##     proportion of genes detected per sample.
##  3. Correlation of PC1 with age, sex, MGMT status and survival, to see
##     whether the axis carries clinical meaning.
##  4. The same mean absolute correlation statistic recomputed after
##     removing PC1, which answers the operative question directly: if the
##     elevated correlation in CGGA-693 disappears once the first component
##     is projected out, that single axis is responsible for it.
##  5. A sensitivity analysis: the observed Chen and Zhang hazard ratios
##     recomputed on PC1-removed expression. If the published estimates are
##     unchanged, the dominant axis, whatever its origin, does not drive
##     the prognostic results.
##
## Gene space. Both cohorts are restricted to their common genes and to
## those expressed in both, so a difference in the annotation or in the
## number of quantified genes cannot by itself produce a difference in
## principal component structure.
##
## Interpretation fixed in advance:
##   - |r| of PC1 with depth or detection >= 0.40 -> technical contribution
##     acknowledged in Limitations.
##   - Mean |r| in CGGA-693 falls to the level of CGGA-325 after removing
##     PC1 -> the difference is attributable to a single dominant axis.
##   - Published hazard ratios unchanged after removing PC1 -> the main
##     results do not depend on it, whatever its nature.
##
## This is the last analytical layer. Everything after this is REMARK
## completion and writing.
##
## Usage:
##   setwd("~/Library/CloudStorage/Dropbox/HOPE-GBM"); source("layer06i_structure.R")
## Requires Layers 06e and 06g.
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
N_RAND  <- 300     # random ten-gene sets used to recompute mean |r|
N_PC    <- 10      # leading components reported

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
RAW <- setNames(lapply(expr_f, read_cgga_expr), paste0("CGGA-", vapply(expr_f, num_of, "")))

## ---- Common gene space, restricted to the stratum under study ---------------
COH <- levels(G4$cohort)
common <- Reduce(intersect, lapply(RAW, rownames))
LG <- lapply(COH, function(k) {
  d <- G4[G4$cohort == k, ]
  m <- log2(RAW[[k]][common, d$patient, drop = FALSE] + 1)
  m[is.finite(rowSums(m)), , drop = FALSE]
})
names(LG) <- COH
common2 <- Reduce(intersect, lapply(LG, rownames))
## keep genes expressed in both cohorts, to avoid a detection-driven artefact
expressed <- Reduce(intersect, lapply(LG, function(m) rownames(m)[rowMeans(m > 0) >= 0.5]))
gene_space <- intersect(common2, expressed)
LG <- lapply(LG, function(m) m[gene_space, , drop = FALSE])
log_msg(sprintf("06i: common gene space = %d genes expressed in at least half of samples in both cohorts",
                length(gene_space)))

## ---- Sample-level technical covariates --------------------------------------
tech <- do.call(rbind, lapply(COH, function(k) {
  d <- G4[G4$cohort == k, ]
  raw <- RAW[[k]][gene_space, d$patient, drop = FALSE]
  data.frame(cohort = k, patient = d$patient,
             depth = colSums(raw, na.rm = TRUE),
             detection = colMeans(raw > 0, na.rm = TRUE), row.names = NULL)
}))

## ---- Principal components ----------------------------------------------------
pc_of <- function(k) {
  m <- LG[[k]]
  v <- apply(m, 1, var, na.rm = TRUE)
  m <- m[is.finite(v) & v > 0, , drop = FALSE]
  p <- prcomp(t(m), center = TRUE, scale. = TRUE)
  ve <- (p$sdev^2) / sum(p$sdev^2)
  list(pca = p, var_explained = ve, mat = m)
}
PC <- setNames(lapply(COH, pc_of), COH)

var_tab <- do.call(rbind, lapply(COH, function(k) data.frame(
  cohort = k, component = paste0("PC", seq_len(N_PC)),
  variance_explained = 100 * PC[[k]]$var_explained[seq_len(N_PC)],
  cumulative = 100 * cumsum(PC[[k]]$var_explained)[seq_len(N_PC)], row.names = NULL)))
for (k in COH)
  log_msg(sprintf("06i %s: PC1 explains %.1f%% of variance, PC1-PC5 %.1f%%",
                  k, 100 * PC[[k]]$var_explained[1],
                  100 * sum(PC[[k]]$var_explained[1:5])))

## ---- What does PC1 correlate with? ------------------------------------------
corr_tab <- do.call(rbind, lapply(COH, function(k) {
  d  <- G4[G4$cohort == k, ]
  tk <- tech[tech$cohort == k, ]
  tk <- tk[match(d$patient, tk$patient), ]
  sc <- PC[[k]]$pca$x
  sc <- sc[match(d$patient, rownames(sc)), , drop = FALSE]
  cs <- function(x, y) suppressWarnings(cor(x, y, method = "spearman", use = "complete.obs"))
  do.call(rbind, lapply(1:2, function(j) data.frame(
    cohort = k, component = paste0("PC", j),
    r_depth      = cs(sc[, j], tk$depth),
    r_detection  = cs(sc[, j], tk$detection),
    r_age        = cs(sc[, j], d$age),
    r_survival   = cs(sc[, j], d$OS.time),
    r_mgmt       = cs(sc[, j], as.numeric(d$MGMT == "Unmethylated")),
    row.names = NULL)))
}))
print(corr_tab)
for (i in seq_len(nrow(corr_tab))) {
  r <- corr_tab[i, ]
  log_msg(sprintf("06i CORR %s %s: depth rho=%.2f | detection rho=%.2f | age rho=%.2f | OS rho=%.2f | MGMT rho=%.2f",
                  r$cohort, r$component, r$r_depth, r$r_detection, r$r_age, r$r_survival, r$r_mgmt))
}

## ---- Mean |r| among random ten-gene sets, before and after removing PC1 ------
## Residualising the expression matrix on PC1 and repeating the statistic that
## first revealed the difference is the direct test of the hypothesis.
remove_pc1 <- function(k) {
  m <- PC[[k]]$mat
  s <- scale(t(m))                       # samples x genes, standardised
  p <- PC[[k]]$pca
  sc1 <- p$x[, 1, drop = FALSE]
  fitted <- sc1 %*% t(p$rotation[, 1, drop = FALSE])
  t(s - fitted)                          # genes x samples, PC1 projected out
}
mean_abs_r <- function(mat, n_sets = N_RAND, k_genes = 10) {
  v <- apply(mat, 1, var, na.rm = TRUE)
  pool <- rownames(mat)[is.finite(v) & v > 0]
  vapply(seq_len(n_sets), function(i) {
    g <- sample(pool, k_genes)
    cm <- suppressWarnings(cor(t(mat[g, , drop = FALSE]), use = "pairwise.complete.obs"))
    mean(abs(cm[upper.tri(cm)]), na.rm = TRUE)
  }, numeric(1))
}
struct_tab <- do.call(rbind, lapply(COH, function(k) {
  before <- mean_abs_r(PC[[k]]$mat)
  after  <- mean_abs_r(remove_pc1(k))
  data.frame(cohort = k,
             mean_abs_r_before = median(before), mean_abs_r_after = median(after),
             reduction_pct = 100 * (median(before) - median(after)) / median(before),
             pc1_variance = 100 * PC[[k]]$var_explained[1], row.names = NULL)
}))
print(struct_tab)
for (i in seq_len(nrow(struct_tab))) {
  r <- struct_tab[i, ]
  log_msg(sprintf("06i STRUCT %s: mean |r| among random 10-gene sets %.2f -> %.2f after removing PC1 (%.0f%% reduction); PC1 = %.1f%% of variance",
                  r$cohort, r$mean_abs_r_before, r$mean_abs_r_after,
                  r$reduction_pct, r$pc1_variance))
}

## ---- Sensitivity: published signatures on PC1-removed expression -------------
sens <- do.call(rbind, lapply(COH, function(k) {
  d <- G4[G4$cohort == k, ]
  m_raw <- PC[[k]]$mat
  m_res <- remove_pc1(k)
  do.call(rbind, lapply(names(SIGS), function(s) {
    w  <- SIGS[[s]]; gs <- intersect(names(w), rownames(m_raw))
    if (length(gs) < 3) return(NULL)
    fit <- function(mat) {
      Z <- t(apply(mat[gs, d$patient, drop = FALSE], 1, zvec))
      dd <- d; dd$score <- zvec(as.numeric(colSums(Z * w[gs])))
      f <- tryCatch(coxph(Surv(OS.time, OS.event) ~ score + age, data = dd),
                    error = function(e) NULL)
      if (is.null(f)) return(c(NA, NA))
      sm <- summary(f)
      c(unname(sm$coefficients["score", "exp(coef)"]),
        unname(sm$coefficients["score", "Pr(>|z|)"]))
    }
    a <- fit(m_raw); b <- fit(m_res)
    data.frame(cohort = k, signature = s, genes_used = length(gs),
               HR_original = a[1], p_original = a[2],
               HR_pc1_removed = b[1], p_pc1_removed = b[2], row.names = NULL)
  }))
}))
print(sens)
for (i in seq_len(nrow(sens))) {
  r <- sens[i, ]
  log_msg(sprintf("06i SENS %s / %s: HR %.2f (p=%.3f) -> %.2f (p=%.3f) after removing PC1",
                  r$cohort, r$signature, r$HR_original, r$p_original,
                  r$HR_pc1_removed, r$p_pc1_removed))
}

## ---- Verdict against the pre-specified rule ---------------------------------
tech_flag <- corr_tab %>% filter(component == "PC1") %>%
  mutate(technical = pmax(abs(r_depth), abs(r_detection)) >= 0.40)
for (i in seq_len(nrow(tech_flag)))
  log_msg(sprintf("06i VERDICT %s PC1: %s (max |rho| with depth or detection = %.2f)",
                  tech_flag$cohort[i],
                  ifelse(tech_flag$technical[i],
                         "technical contribution to be acknowledged in Limitations",
                         "no appreciable association with sequencing depth or detection"),
                  pmax(abs(tech_flag$r_depth[i]), abs(tech_flag$r_detection[i]))[1]))

## ---- Tables and figures ------------------------------------------------------
tS3 <- data.frame(
  Cohort = struct_tab$cohort,
  `Variance explained by PC1` = sprintf("%.1f%%", struct_tab$pc1_variance),
  `Mean absolute correlation, random ten-gene sets` = sprintf("%.2f", struct_tab$mean_abs_r_before),
  `After removing PC1` = sprintf("%.2f", struct_tab$mean_abs_r_after),
  `Reduction` = sprintf("%.0f%%", struct_tab$reduction_pct),
  `PC1 correlation with sequencing depth` = sprintf("%.2f",
    corr_tab$r_depth[corr_tab$component == "PC1"]),
  `PC1 correlation with detection rate` = sprintf("%.2f",
    corr_tab$r_detection[corr_tab$component == "PC1"]),
  check.names = FALSE)
ftS3 <- flextable::flextable(tS3) %>% ft_style() %>%
  flextable::add_footer_lines(values = c(
    sprintf("Principal component analysis was performed separately in each cohort on the %d genes quantified and expressed in at least half of the samples of both series, so that differences in annotation or detection cannot by themselves generate differences in component structure.", length(gene_space)),
    sprintf("The mean absolute correlation was computed across %d random sets of ten genes and is reported before and after projecting out the first principal component; this is the statistic that first revealed the difference between the two series.", N_RAND),
    "Correlations with sequencing depth, approximated by the column sum of the expression matrix, and with the proportion of genes detected per sample are given as Spearman coefficients and are intended to separate a technical axis from a biological one.",
    "PC, principal component."))
flextable::save_as_docx(ftS3, path = file.path(DIR_TAB, "TableS3_expression_structure.docx"))

tS4 <- data.frame(
  Cohort = sens$cohort, Signature = sens$signature,
  `Genes used` = sens$genes_used,
  `HR, original expression` = sprintf("%.2f (p=%s)", sens$HR_original,
    ifelse(sens$p_original < 0.001, "<0.001", sprintf("%.3f", sens$p_original))),
  `HR, first component removed` = sprintf("%.2f (p=%s)", sens$HR_pc1_removed,
    ifelse(sens$p_pc1_removed < 0.001, "<0.001", sprintf("%.3f", sens$p_pc1_removed))),
  check.names = FALSE)
ftS4 <- flextable::flextable(tS4) %>% ft_style() %>%
  flextable::add_footer_lines(values = c(
    "Sensitivity analysis in which each published signature was recomputed after projecting out the first principal component of the expression matrix, with the same age-adjusted Cox model and the same stratum. Estimates that are unchanged indicate that the dominant axis of expression, whatever its origin, does not account for the prognostic results.",
    "Hazard ratios are expressed per 1 standard deviation of the signature score.",
    "HR, hazard ratio; PC, principal component."))
flextable::save_as_docx(ftS4, path = file.path(DIR_TAB, "TableS4_PC1_sensitivity.docx"))

gS2 <- ggplot(var_tab, aes(x = component, y = variance_explained, fill = cohort)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.68) +
  scale_fill_grey(start = 0.35, end = 0.75) +
  scale_x_discrete(limits = paste0("PC", seq_len(N_PC))) +
  labs(x = NULL, y = "Variance explained (%)", fill = NULL) +
  theme_classic(base_size = 11) + theme(legend.position = "bottom")
save_fig(gS2, file.path(DIR_FIG, "FigureS2_variance_explained"), width = 7, height = 4.2)

## ---- Persist -------------------------------------------------------------------
saveRDS(list(variance = var_tab, correlations = corr_tab, structure = struct_tab,
             sensitivity = sens, gene_space_n = length(gene_space)),
        file.path(DIR_RES, "expression_structure.rds"))
write.csv(var_tab,    file.path(DIR_RES, "tableS3a_variance_explained.csv"), row.names = FALSE)
write.csv(corr_tab,   file.path(DIR_RES, "tableS3b_pc_correlations.csv"), row.names = FALSE)
write.csv(struct_tab, file.path(DIR_RES, "tableS3_expression_structure.csv"), row.names = FALSE)
write.csv(sens,       file.path(DIR_RES, "tableS4_pc1_sensitivity.csv"), row.names = FALSE)

nums <- c("", "--- Layer 06i: expression structure and its technical correlates ---",
  sprintf("Common gene space: %d genes", length(gene_space)),
  sprintf("%s: PC1 explains %.1f%% of variance; mean |r| among random ten-gene sets %.2f -> %.2f after removing PC1 (%.0f%% reduction)",
          struct_tab$cohort, struct_tab$pc1_variance, struct_tab$mean_abs_r_before,
          struct_tab$mean_abs_r_after, struct_tab$reduction_pct),
  "",
  sprintf("%s %s: rho with depth %.2f, with detection %.2f, with age %.2f, with OS %.2f, with unmethylated MGMT %.2f",
          corr_tab$cohort, corr_tab$component, corr_tab$r_depth, corr_tab$r_detection,
          corr_tab$r_age, corr_tab$r_survival, corr_tab$r_mgmt),
  "",
  sprintf("%s / %s: HR %.2f (p=%.3f) with original expression, %.2f (p=%.3f) after removing PC1",
          sens$cohort, sens$signature, sens$HR_original, sens$p_original,
          sens$HR_pc1_removed, sens$p_pc1_removed))
cat(nums, sep = "\n")
write(nums, file = file.path(DIR_OUT, "manuscript_numbers.txt"), append = TRUE)
log_msg("Layer 06i complete: Tables S3-S4 and Figure S2 exported to Sep_06/. Analysis phase closed.")
