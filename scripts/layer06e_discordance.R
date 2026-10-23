## =====================================================================
## Layer 06e - Why do CGGA-325 and CGGA-693 disagree?
##
## The problem. In the decisive stratum (IDH-wildtype, WHO grade IV) the
## two published signatures behave differently across cohorts:
##
##   CGGA-325 (n=74,  deaths=64): Chen HR 1.38 (p=0.016), Zhang 1.25 (p=0.071)
##   CGGA-693 (n=105, deaths=92): Chen HR 1.01 (p=0.96),  Zhang 1.08 (p=0.44)
##   TCGA-GBM (n=191, deaths=153): Chen HR 0.99 (p=0.91)
##
## The positive cohort is the smallest, so lack of power in CGGA-693 does
## not explain it. This layer tests the composition and technical
## explanations that can be checked with the data already on disk, and
## quantifies what the study could and could not have detected.
##
## Tests performed, restricted to the IDH-wildtype grade IV stratum unless
## stated otherwise:
##
##   A. Composition. Age, sex, MGMT promoter status, chemotherapy and
##      radiotherapy where recorded, survival distribution, and
##      transcriptomic subtype if available from the pipeline objects.
##   B. Score distribution. Mean, SD and range of each signature score
##      before standardization, plus the correlation between the two
##      signatures within each cohort. A signature that is compressed in
##      one cohort has less room to show an effect.
##   C. Technical. Library-size proxy and the distribution of the raw
##      expression of the signature genes, to detect a batch or platform
##      difference between the two CGGA series.
##   D. Formal heterogeneity. Fixed-effect and random-effects meta-analysis
##      of the three IDH-wildtype grade IV estimates per signature, with
##      Cochran Q and I-squared. This replaces the informal claim that the
##      cohorts "disagree" with a measured statement.
##   E. Power. The hazard ratio detectable with 80% power at alpha 0.05
##      given the observed number of events in each cohort, so the null
##      results are described by what they exclude rather than by p > 0.05.
##
## No new hypothesis is generated here; this layer is diagnostic and its
## results belong in Results (heterogeneity, power) and in Limitations
## (whatever remains unexplained).
##
## Usage:
##   setwd("~/Library/CloudStorage/Dropbox/HOPE-GBM"); source("layer06e_discordance.R")
## Requires Layer 06d to have been run (results/generality_results.rds).
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
CGGA_DIR <- Find(function(p) dir.exists(p) &&
                   length(list.files(p, pattern = "mRNAseq.*genes.*[.]txt$")) > 0,
                 c(path.expand("~/Documents/CGGA_local"),
                   path.expand("~/Library/CloudStorage/Dropbox/00_DATA/CGGA"),
                   file.path(PROJ, "data_raw", "CGGA")))
stopifnot(!is.null(CGGA_DIR))
gen <- readRDS(file.path(DIR_RES, "generality_results.rds"))
SIGS <- gen$signatures
FIG_DPI <- 600

log_msg <- function(...) {
  m <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ...)
  cat(m, "\n"); cat(m, "\n", file = file.path(PROJ, "logs", "hope_run.log"), append = TRUE)
}
zrow <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
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
  tiff(paste0(stem, ".tif"), width = width, height = height, units = "in", res = dpi,
       compression = "lzw"); print(obj); dev.off()
}
read_cgga_expr <- function(path, label) {
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
pick_col <- function(nm, exact, pattern) {
  hit <- exact[exact %in% nm]
  if (length(hit)) return(hit[1])
  grep(pattern, nm, value = TRUE, ignore.case = TRUE)[1]
}
harmonize_mgmt <- function(x) {
  x <- gsub("[^a-z]", "", tolower(as.character(x)))
  factor(ifelse(x == "methylated", "Methylated",
         ifelse(x == "unmethylated", "Unmethylated", NA_character_)),
         levels = c("Methylated", "Unmethylated"))
}

expr_f <- list.files(CGGA_DIR, pattern = "mRNAseq.*genes.*[.]txt$", full.names = TRUE)
clin_f <- list.files(CGGA_DIR, pattern = "clinical.*[.]txt$",       full.names = TRUE)
num_of <- function(p) regmatches(basename(p), regexpr("325|693", basename(p)))
clin_f <- clin_f[match(vapply(expr_f, num_of, ""), vapply(clin_f, num_of, ""))]

## ---- Build the IDH-wildtype grade IV stratum with full clinical detail -------
build <- function(expr_path, clin_path) {
  label <- paste0("CGGA-", num_of(expr_path))
  expr  <- read_cgga_expr(expr_path, label)
  clin  <- as.data.frame(readr::read_tsv(clin_path, show_col_types = FALSE, progress = FALSE))
  nm <- colnames(clin)
  cols <- list(
    id = pick_col(nm, "CGGA_ID", "CGGA|sample|ID"),
    prs = pick_col(nm, "PRS_type", "PRS"), idh = pick_col(nm, "IDH_mutation_status", "IDH"),
    grade = pick_col(nm, "Grade", "Grade|WHO"), age = pick_col(nm, "Age", "Age"),
    sex = pick_col(nm, "Gender", "Gender|Sex"), mgmt = pick_col(nm, "MGMTp_methylation_status", "MGMT"),
    os = pick_col(nm, "OS", "^OS$|overall"),
    chemo = grep("chemo", nm, value = TRUE, ignore.case = TRUE)[1],
    radio = grep("radio", nm, value = TRUE, ignore.case = TRUE)[1],
    x1p19q = grep("1p19q|codel", nm, value = TRUE, ignore.case = TRUE)[1])
  ev <- grep("Censor", nm, value = TRUE, ignore.case = TRUE)[1]
  log_msg(sprintf("06e %s: clinical columns detected -> chemo=%s, radio=%s, 1p19q=%s",
                  label, cols$chemo, cols$radio, cols$x1p19q))
  keep <- if (!is.na(cols$prs)) clin[[cols$prs]] == "Primary" else rep(TRUE, nrow(clin))
  keep[is.na(keep)] <- FALSE
  clin <- clin[keep, ]
  ids <- intersect(clin[[cols$id]], colnames(expr))
  expr <- expr[, ids, drop = FALSE]; clin <- clin[match(ids, clin[[cols$id]]), ]
  lg <- log2(expr + 1)

  d <- data.frame(
    patient = ids, cohort = label,
    OS.time = as.numeric(clin[[cols$os]]), OS.event = as.integer(clin[[ev]]),
    age = suppressWarnings(as.numeric(clin[[cols$age]])),
    sex = if (!is.na(cols$sex)) as.character(clin[[cols$sex]]) else NA_character_,
    MGMT = harmonize_mgmt(if (!is.na(cols$mgmt)) clin[[cols$mgmt]] else NA),
    chemo = if (!is.na(cols$chemo)) as.character(clin[[cols$chemo]]) else NA_character_,
    radio = if (!is.na(cols$radio)) as.character(clin[[cols$radio]]) else NA_character_,
    IDH = factor(ifelse(grepl("wildtype|^wt", tolower(clin[[cols$idh]])), "Wildtype",
                 ifelse(grepl("mutant|^mut", tolower(clin[[cols$idh]])), "Mutant", NA)),
                 levels = c("Wildtype", "Mutant")),
    grade = factor(gsub("\\s", "", as.character(clin[[cols$grade]]))),
    libsize = colSums(expr[, ids, drop = FALSE], na.rm = TRUE),
    stringsAsFactors = FALSE)
  for (s in names(SIGS)) {
    w  <- SIGS[[s]]; gs <- intersect(names(w), rownames(lg))
    Z  <- t(apply(lg[gs, , drop = FALSE], 1, zrow))
    d[[make.names(s)]] <- as.numeric(colSums(Z * w[gs]))[match(ids, colnames(Z))]
  }
  ## raw (unstandardized) mean expression of signature genes, for the technical check
  d$sig_gene_mean <- colMeans(lg[intersect(unique(unlist(lapply(SIGS, names))), rownames(lg)),
                                 ids, drop = FALSE], na.rm = TRUE)
  d <- d[!is.na(d$OS.time) & d$OS.time > 0 & !is.na(d$IDH) & !is.na(d$age), ]
  d
}
ALL <- do.call(rbind, Map(build, expr_f, clin_f))
G4 <- ALL[ALL$IDH == "Wildtype" & grepl("IV|4", ALL$grade), ]
G4$cohort <- factor(G4$cohort, levels = c("CGGA-325", "CGGA-693"))
log_msg(sprintf("06e: IDH-wildtype grade IV -> CGGA-325 n=%d (deaths %d), CGGA-693 n=%d (deaths %d)",
                sum(G4$cohort == "CGGA-325"), sum(G4$OS.event[G4$cohort == "CGGA-325"]),
                sum(G4$cohort == "CGGA-693"), sum(G4$OS.event[G4$cohort == "CGGA-693"])))

## ---- A. Composition ----------------------------------------------------------
num_row <- function(v, name) {
  a <- v[G4$cohort == "CGGA-325"]; b <- v[G4$cohort == "CGGA-693"]
  p <- tryCatch(wilcox.test(a, b)$p.value, error = function(e) NA_real_)
  data.frame(Variable = name,
             `CGGA-325` = sprintf("%.1f (%.1f, %.1f)", median(a, na.rm = TRUE),
                                  quantile(a, .25, na.rm = TRUE), quantile(a, .75, na.rm = TRUE)),
             `CGGA-693` = sprintf("%.1f (%.1f, %.1f)", median(b, na.rm = TRUE),
                                  quantile(b, .25, na.rm = TRUE), quantile(b, .75, na.rm = TRUE)),
             `p-value` = ifelse(is.na(p), "NA", sprintf("%.3f", p)), check.names = FALSE)
}
cat_row <- function(v, name, level) {
  tb <- table(G4$cohort, v)
  if (!(level %in% colnames(tb))) return(NULL)
  p <- tryCatch(chisq.test(tb)$p.value, error = function(e) NA_real_)
  n <- rowSums(tb)
  data.frame(Variable = sprintf("%s: %s, n (%%)", name, level),
             `CGGA-325` = sprintf("%d (%.0f%%)", tb["CGGA-325", level], 100 * tb["CGGA-325", level] / n["CGGA-325"]),
             `CGGA-693` = sprintf("%d (%.0f%%)", tb["CGGA-693", level], 100 * tb["CGGA-693", level] / n["CGGA-693"]),
             `p-value` = ifelse(is.na(p), "NA", sprintf("%.3f", p)), check.names = FALSE)
}
sexf <- ifelse(grepl("^f", tolower(G4$sex)), "Female", ifelse(grepl("^m", tolower(G4$sex)), "Male", NA))
comp <- rbind(
  num_row(G4$age, "Age, years, median (Q1, Q3)"),
  num_row(G4$OS.time / 30.4375, "Observed follow-up, months, median (Q1, Q3)"),
  cat_row(sexf, "Sex", "Female"),
  cat_row(as.character(G4$MGMT), "MGMT promoter", "Methylated"),
  if (!all(is.na(G4$chemo))) cat_row(as.character(G4$chemo), "Chemotherapy", names(sort(table(G4$chemo), decreasing = TRUE))[1]),
  if (!all(is.na(G4$radio))) cat_row(as.character(G4$radio), "Radiotherapy", names(sort(table(G4$radio), decreasing = TRUE))[1]),
  num_row(G4$Chen.2017, "Chen 2017 score, raw units"),
  num_row(G4$Zhang.2021, "Zhang 2021 score, raw units"),
  num_row(G4$sig_gene_mean, "Mean log2 expression of signature genes"),
  num_row(G4$libsize / 1e6, "Library size proxy, millions"))
## survival comparison uses log-rank, not a rank test on months
lr <- survdiff(Surv(OS.time, OS.event) ~ cohort, data = G4)
lr_p <- 1 - pchisq(lr$chisq, df = 1)
comp <- rbind(comp, data.frame(
  Variable = "Overall survival (log-rank between cohorts)",
  `CGGA-325` = sprintf("%d deaths / %d", sum(G4$OS.event[G4$cohort == "CGGA-325"]), sum(G4$cohort == "CGGA-325")),
  `CGGA-693` = sprintf("%d deaths / %d", sum(G4$OS.event[G4$cohort == "CGGA-693"]), sum(G4$cohort == "CGGA-693")),
  `p-value` = sprintf("%.3f", lr_p), check.names = FALSE))
print(comp)

## ---- B. Score behaviour ------------------------------------------------------
spread <- do.call(rbind, lapply(levels(G4$cohort), function(k) {
  d <- G4[G4$cohort == k, ]
  data.frame(cohort = k,
             sd_chen = sd(d$Chen.2017), sd_zhang = sd(d$Zhang.2021),
             iqr_chen = IQR(d$Chen.2017), iqr_zhang = IQR(d$Zhang.2021),
             rho_between_signatures = cor(d$Chen.2017, d$Zhang.2021, method = "spearman"),
             row.names = NULL)
}))
print(spread)
for (i in seq_len(nrow(spread)))
  log_msg(sprintf("06e %s: SD Chen %.2f, SD Zhang %.2f, rho(Chen,Zhang)=%.2f",
                  spread$cohort[i], spread$sd_chen[i], spread$sd_zhang[i],
                  spread$rho_between_signatures[i]))

## ---- D. Heterogeneity across the three IDH-wt grade IV estimates -------------
## TCGA estimate for the Chen signature comes from Layer 06 (published weights).
tab <- gen$TAB
g4rows <- tab[tab$stratum == "IDH-wildtype, WHO grade IV", ]
meta_one <- function(sig) {
  r <- g4rows[g4rows$signature == sig, ]
  est <- data.frame(cohort = r$cohort, beta = log(r$HR),
                    se = (log(r$hi) - log(r$lo)) / (2 * 1.959964))
  if (sig == "Chen 2017" && file.exists(file.path(DIR_RES, "table4_ggess_values.csv"))) {
    t4 <- read.csv(file.path(DIR_RES, "table4_ggess_values.csv"))
    t4 <- t4[t4$construction == "published" & t4$cohort == "TCGA-GBM", ]
    if (nrow(t4) == 1)
      est <- rbind(est, data.frame(cohort = "TCGA-GBM", beta = log(t4$HR_base),
                                   se = (log(t4$hi_base) - log(t4$lo_base)) / (2 * 1.959964)))
  }
  w <- 1 / est$se^2
  b_fe <- sum(w * est$beta) / sum(w); se_fe <- sqrt(1 / sum(w))
  Q <- sum(w * (est$beta - b_fe)^2); df <- nrow(est) - 1
  Qp <- 1 - pchisq(Q, df); I2 <- max(0, 100 * (Q - df) / Q)
  tau2 <- max(0, (Q - df) / (sum(w) - sum(w^2) / sum(w)))
  wr <- 1 / (est$se^2 + tau2)
  b_re <- sum(wr * est$beta) / sum(wr); se_re <- sqrt(1 / sum(wr))
  list(signature = sig, k = nrow(est), est = est,
       fe = c(HR = exp(b_fe), lo = exp(b_fe - 1.96 * se_fe), hi = exp(b_fe + 1.96 * se_fe)),
       re = c(HR = exp(b_re), lo = exp(b_re - 1.96 * se_re), hi = exp(b_re + 1.96 * se_re)),
       Q = Q, df = df, Qp = Qp, I2 = I2, tau2 = tau2)
}
metas <- lapply(unique(g4rows$signature), meta_one)
for (m in metas)
  log_msg(sprintf("06e META %s (k=%d): FE HR %.2f (%.2f-%.2f) | RE HR %.2f (%.2f-%.2f) | Q=%.2f df=%d p=%.3f | I2=%.0f%%",
                  m$signature, m$k, m$fe[1], m$fe[2], m$fe[3], m$re[1], m$re[2], m$re[3],
                  m$Q, m$df, m$Qp, m$I2))

## ---- E. Detectable effect size -----------------------------------------------
## Schoenfeld: events needed = (z_{1-a/2} + z_{1-b})^2 / (p(1-p) * log(HR)^2),
## with p = 0.5 for a continuous covariate expressed per 1 SD (variance 1),
## which reduces to events = (z1+z2)^2 / log(HR)^2.
detectable_hr <- function(events, power = 0.80, alpha = 0.05) {
  z1 <- qnorm(1 - alpha / 2); z2 <- qnorm(power)
  exp((z1 + z2) / sqrt(events))
}
pw <- do.call(rbind, lapply(levels(G4$cohort), function(k) {
  e <- sum(G4$OS.event[G4$cohort == k])
  data.frame(cohort = k, events = e, HR_detectable_80 = detectable_hr(e), row.names = NULL)
}))
tcga_ev <- tryCatch({
  t2 <- read.csv(file.path(DIR_RES, "table2_figure1_values.csv"))
  t2$deaths[grepl("TCGA", t2$cohort)][1] }, error = function(e) NA_integer_)
if (!is.na(tcga_ev))
  pw <- rbind(pw, data.frame(cohort = "TCGA-GBM", events = tcga_ev,
                             HR_detectable_80 = detectable_hr(tcga_ev)))
print(pw)
for (i in seq_len(nrow(pw)))
  log_msg(sprintf("06e POWER %s: %d events -> detectable HR per SD >= %.2f (80%% power, alpha 0.05)",
                  pw$cohort[i], pw$events[i], pw$HR_detectable_80[i]))

## ---- Tables and figure --------------------------------------------------------
ft9 <- flextable::flextable(comp) %>% ft_style() %>%
  flextable::add_footer_lines(values = c(
    "Comparison restricted to IDH-wildtype, WHO grade IV gliomas, the stratum in which the two CGGA series gave discordant results. Continuous variables were compared with the Wilcoxon rank-sum test and categorical variables with the chi-square test; overall survival was compared with the log-rank test.",
    "The library-size proxy is the column sum of the RSEM expression matrix and is reported to detect gross technical differences between the two series, not to substitute for a formal batch analysis.",
    "Q1, first quartile; Q3, third quartile."))
flextable::save_as_docx(ft9, path = file.path(DIR_TAB, "Table9_CGGA_discordance_composition.docx"))

meta_tab <- do.call(rbind, lapply(metas, function(m) data.frame(
  Signature = m$signature, `Cohorts (k)` = m$k,
  `Fixed-effect HR (95% CI)` = sprintf("%.2f (%.2f, %.2f)", m$fe[1], m$fe[2], m$fe[3]),
  `Random-effects HR (95% CI)` = sprintf("%.2f (%.2f, %.2f)", m$re[1], m$re[2], m$re[3]),
  `Cochran Q (df), p` = sprintf("%.2f (%d), p=%.3f", m$Q, m$df, m$Qp),
  `I-squared` = sprintf("%.0f%%", m$I2), check.names = FALSE)))
ft10 <- flextable::flextable(meta_tab) %>% ft_style() %>%
  flextable::add_footer_lines(values = c(
    "Meta-analysis of the cohort-specific hazard ratios for each published signature in the IDH-wildtype, WHO grade IV stratum. Estimates were combined on the log hazard ratio scale with inverse-variance weights; between-cohort variance was estimated by the method of moments.",
    "With so few cohorts, the heterogeneity statistics are imprecise and are reported to describe the dispersion of the estimates rather than to support a formal test.",
    "CI, confidence interval; HR, hazard ratio."))
flextable::save_as_docx(ft10, path = file.path(DIR_TAB, "Table10_G4_heterogeneity.docx"))

fp <- do.call(rbind, lapply(metas, function(m) {
  e <- m$est
  rbind(data.frame(signature = m$signature, label = e$cohort, HR = exp(e$beta),
                   lo = exp(e$beta - 1.96 * e$se), hi = exp(e$beta + 1.96 * e$se), kind = "Cohort"),
        data.frame(signature = m$signature, label = "Random-effects summary",
                   HR = m$re[1], lo = m$re[2], hi = m$re[3], kind = "Summary"))
}))
fp$label <- factor(fp$label, levels = rev(unique(fp$label)))
g7 <- ggplot(fp, aes(x = HR, y = label, shape = kind)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.16) +
  geom_point(size = 2.6) + facet_wrap(~signature, ncol = 1, scales = "free_y") +
  scale_x_log10() + scale_shape_manual(values = c(Cohort = 16, Summary = 18)) +
  labs(x = "HR per 1 SD of signature score (95% CI), IDH-wildtype WHO grade IV",
       y = NULL, shape = NULL) +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom", axis.line.y = element_blank(), axis.ticks.y = element_blank())
save_fig(g7, file.path(DIR_FIG, "Figure7_G4_heterogeneity"), width = 7, height = 5.2)

## ---- Persist -------------------------------------------------------------------
saveRDS(list(comp = comp, spread = spread, metas = metas, power = pw, G4 = G4),
        file.path(DIR_RES, "discordance_results.rds"))
write.csv(comp,     file.path(DIR_RES, "table9_discordance_composition.csv"), row.names = FALSE)
write.csv(meta_tab, file.path(DIR_RES, "table10_g4_heterogeneity.csv"), row.names = FALSE)
write.csv(pw,       file.path(DIR_RES, "power_detectable_hr.csv"), row.names = FALSE)

nums <- c("", "--- Layer 06e: discordance between CGGA-325 and CGGA-693 (IDH-wt, WHO grade IV) ---",
  sprintf("%s | CGGA-325: %s | CGGA-693: %s | p=%s", comp$Variable, comp$`CGGA-325`,
          comp$`CGGA-693`, comp$`p-value`),
  "",
  sprintf("Score spread %s: SD Chen %.2f, SD Zhang %.2f, Spearman rho between signatures %.2f",
          spread$cohort, spread$sd_chen, spread$sd_zhang, spread$rho_between_signatures),
  "",
  unlist(lapply(metas, function(m) sprintf(
    "META %s (k=%d): fixed-effect HR %.2f (%.2f-%.2f); random-effects HR %.2f (%.2f-%.2f); Q=%.2f df=%d p=%.3f; I2=%.0f%%",
    m$signature, m$k, m$fe[1], m$fe[2], m$fe[3], m$re[1], m$re[2], m$re[3], m$Q, m$df, m$Qp, m$I2))),
  "",
  sprintf("POWER %s: %d events, detectable HR per SD >= %.2f (80%% power, two-sided alpha 0.05)",
          pw$cohort, pw$events, pw$HR_detectable_80))
cat(nums, sep = "\n")
write(nums, file = file.path(DIR_OUT, "manuscript_numbers.txt"), append = TRUE)
log_msg("Layer 06e complete: Tables 9-10 and Figure 7 exported to Sep_06/.")
