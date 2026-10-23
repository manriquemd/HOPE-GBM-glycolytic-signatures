## =====================================================================
## Layer 06d - Generality of the IDH/grade explanation, and the anomaly
##
## Two questions left open by Layers 06b and 06c:
##
##   Q1 (generality). Is the behaviour of the Chen 2017 GGESS specific to
##       that signature, or shared by any glycolysis/glucose-metabolism
##       signature trained on mixed-grade, mixed-IDH glioma cohorts?
##       Tested with a second, independent, published signature: the ten
##       glucose-metabolism genes of Zhang et al. (BioMed Res Int 2021;
##       2021:8872977), with their published Cox coefficients.
##
##   Q2 (the anomaly). In Layer 06c the Chen signature remained prognostic
##       WITHIN IDH-mutant gliomas (HR 2.85 in CGGA-325, 1.84 in CGGA-693),
##       which a purely IDH-driven mechanism does not predict. IDH-mutant
##       gliomas span WHO grades II to IV, so the candidate explanation is
##       grade. Tested by stratifying the IDH-mutant subgroup by grade and
##       by adding grade to the within-stratum model.
##
## Also computed, for every stratum and signature: likelihood-ratio test and
## C-index of the covariate-only model versus the same model with the
## signature added, on identical complete cases. This quantifies "adds
## nothing" instead of inferring it from a confidence interval.
##
## Strata analysed per cohort (CGGA-325 and CGGA-693, never pooled):
##   S1 all primary gliomas, grades II-IV, both IDH states
##   S2 IDH-wildtype only
##   S3 IDH-mutant only
##   S4 IDH-wildtype, WHO grade IV only  (the WHO 2021 GBM population)
##
## Usage:
##   setwd("~/Library/CloudStorage/Dropbox/HOPE-GBM"); source("layer06d_generality.R")
## Outputs to Sep_06/: Table7, Table8, Figure6, and appended numbers.
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
stopifnot(!is.null(CGGA_DIR), dir.exists(DIR_TAB))

## ---- Published signatures --------------------------------------------------
## Chen C, et al. J Genet Genomics 2017;44:519-530, Table 1 (univariate Cox b).
## PKM2 in the paper is HGNC PKM.
SIG_CHEN  <- c(ALDOB = -0.414, ENO1 = 0.218, GALM = 0.242, GAPDH = 0.194,
               HK2 = 0.211, HK3 = 0.319, LDHA = 0.269, LDHB = -0.296,
               PKLR = -0.181, PKM = 0.213)
## Zhang C, et al. BioMed Res Int 2021;2021:8872977, Table 1 / risk-score
## formula in Results 3.2 (multivariable Cox b).
SIG_ZHANG <- c(B4GALT7 = 2.0604, CHST12 = 1.4322, G6PC2 = -2.8374, GALE = 1.4081,
               IL13RA1 = 1.0801, LDHB = -3.2119, SPAG4 = 0.3957, STC1 = 0.4413,
               TGFBI = -1.4198, TPBG = 0.5223)
SIGS <- list(`Chen 2017`  = SIG_CHEN, `Zhang 2021` = SIG_ZHANG)
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
  ncol_data <- as.integer(names(sort(table(nf[-1]), decreasing = TRUE))[1])
  hdr  <- fl[[1]]
  ids  <- if (length(hdr) == ncol_data) hdr[-1] else hdr
  body <- fl[-1][nf[-1] == ncol_data]
  stopifnot(length(ids) > 10, length(body) > 1000)
  gene <- enc2utf8(vapply(body, `[`, "", 1L))
  expr <- matrix(suppressWarnings(as.numeric(unlist(lapply(body, `[`, -1L), use.names = FALSE))),
                 nrow = length(body), ncol = length(ids), byrow = TRUE, dimnames = list(NULL, ids))
  ok   <- !is.na(gene) & nzchar(gene) & rowSums(!is.na(expr)) > 0
  expr <- expr[ok, , drop = FALSE]; gene <- gene[ok]
  ord  <- order(rowMeans(expr, na.rm = TRUE), decreasing = TRUE)
  expr <- expr[ord, , drop = FALSE]; gene <- gene[ord]
  dup  <- duplicated(gene)
  expr <- expr[!dup, , drop = FALSE]; rownames(expr) <- gene[!dup]
  log_msg(sprintf("06d %s: %d genes x %d samples", label, nrow(expr), ncol(expr)))
  expr
}
pick_col <- function(nm, exact, pattern) {
  hit <- exact[exact %in% nm]
  if (length(hit)) return(hit[1])
  grep(pattern, nm, value = TRUE, ignore.case = TRUE)[1]
}

expr_f <- list.files(CGGA_DIR, pattern = "mRNAseq.*genes.*[.]txt$", full.names = TRUE)
clin_f <- list.files(CGGA_DIR, pattern = "clinical.*[.]txt$",       full.names = TRUE)
num_of <- function(p) regmatches(basename(p), regexpr("325|693", basename(p)))
clin_f <- clin_f[match(vapply(expr_f, num_of, ""), vapply(clin_f, num_of, ""))]
stopifnot(!any(is.na(clin_f)))

## ---- Build one clinical + score frame per cohort ----------------------------
build <- function(expr_path, clin_path) {
  label <- paste0("CGGA-", num_of(expr_path))
  expr  <- read_cgga_expr(expr_path, label)
  clin  <- as.data.frame(readr::read_tsv(clin_path, show_col_types = FALSE, progress = FALSE))
  nm <- colnames(clin)
  id_col <- pick_col(nm, "CGGA_ID", "CGGA|sample|ID")
  prs    <- pick_col(nm, "PRS_type", "PRS")
  idhc   <- pick_col(nm, "IDH_mutation_status", "IDH")
  grc    <- pick_col(nm, "Grade", "Grade|WHO")
  agec   <- pick_col(nm, "Age", "Age")
  osc    <- pick_col(nm, "OS", "^OS$|overall")
  evc    <- grep("Censor", nm, value = TRUE, ignore.case = TRUE)[1]
  keep <- if (!is.na(prs)) clin[[prs]] == "Primary" else rep(TRUE, nrow(clin))
  keep[is.na(keep)] <- FALSE
  clin <- clin[keep, ]
  ids  <- intersect(clin[[id_col]], colnames(expr))
  expr <- expr[, ids, drop = FALSE]; clin <- clin[match(ids, clin[[id_col]]), ]
  lg <- log2(expr + 1)

  d <- data.frame(
    patient  = ids,
    OS.time  = as.numeric(clin[[osc]]),
    OS.event = as.integer(clin[[evc]]),
    age      = suppressWarnings(as.numeric(clin[[agec]])),
    IDH      = factor(ifelse(grepl("wildtype|^wt", tolower(clin[[idhc]])), "Wildtype",
                      ifelse(grepl("mutant|^mut",  tolower(clin[[idhc]])), "Mutant", NA)),
                      levels = c("Wildtype", "Mutant")),
    grade    = if (!is.na(grc)) factor(gsub("\\s", "", as.character(clin[[grc]]))) else NA,
    stringsAsFactors = FALSE)
  for (s in names(SIGS)) {
    w  <- SIGS[[s]]
    gs <- intersect(names(w), rownames(lg))
    if (length(gs) < length(w))
      log_msg(sprintf("06d %s / %s: genes absent: %s", label, s,
                      paste(setdiff(names(w), gs), collapse = ", ")))
    Z <- t(apply(lg[gs, , drop = FALSE], 1, zrow))
    d[[make.names(s)]] <- as.numeric(colSums(Z * w[gs]))[match(ids, colnames(Z))]
  }
  d <- d[!is.na(d$OS.time) & d$OS.time > 0 & !is.na(d$IDH) & !is.na(d$age), ]
  log_msg(sprintf("06d %s: analysed n=%d (wt %d, mut %d), deaths=%d | grades: %s", label, nrow(d),
                  sum(d$IDH == "Wildtype"), sum(d$IDH == "Mutant"), sum(d$OS.event),
                  paste(names(table(d$grade)), table(d$grade), sep = "=", collapse = ", ")))
  list(label = label, d = d)
}
CO <- Map(build, expr_f, clin_f); names(CO) <- vapply(CO, `[[`, "", "label")

## ---- Strata -----------------------------------------------------------------
grade4 <- function(g) grepl("IV|4", as.character(g))
strata_of <- function(d) list(
  `All gliomas (II-IV), both IDH states` = d,
  `IDH-wildtype`                         = d[d$IDH == "Wildtype", ],
  `IDH-mutant`                           = d[d$IDH == "Mutant", ],
  `IDH-wildtype, WHO grade IV`           = d[d$IDH == "Wildtype" & grade4(d$grade), ])

## ---- Core model fitting -----------------------------------------------------
## Within each stratum: covariates alone vs covariates + signature.
## Covariates are age, plus grade where the stratum still contains >1 grade.
fit_block <- function(d, sig_col, stratum, cohort, signature) {
  if (nrow(d) < 25 || sum(d$OS.event) < 12) {
    log_msg(sprintf("06d %s / %s / %s: skipped (n=%d, deaths=%d)",
                    cohort, signature, stratum, nrow(d), sum(d$OS.event)))
    return(NULL)
  }
  d$score <- as.numeric(scale(d[[sig_col]]))
  use_grade <- !all(is.na(d$grade)) && nlevels(droplevels(factor(d$grade))) > 1
  if (use_grade) d$grade <- droplevels(factor(d$grade))
  cc <- d[complete.cases(d[, c("score", "age", if (use_grade) "grade")]), ]
  rhs0 <- paste("age", if (use_grade) "+ grade")
  f0 <- coxph(as.formula(paste("Surv(OS.time, OS.event) ~", rhs0)), data = cc)
  f1 <- coxph(as.formula(paste("Surv(OS.time, OS.event) ~ score +", rhs0)), data = cc)
  lrt <- anova(f0, f1, test = "LRT")
  s1 <- summary(f1); ci <- s1$conf.int["score", ]
  data.frame(
    cohort = cohort, signature = signature, stratum = stratum,
    covariates = if (use_grade) "age + WHO grade" else "age",
    n = f1$n, deaths = f1$nevent,
    HR = unname(ci[1]), lo = unname(ci[3]), hi = unname(ci[4]),
    p  = s1$coefficients["score", "Pr(>|z|)"],
    c_base = unname(summary(f0)$concordance[1]),
    c_full = unname(s1$concordance[1]),
    lrt_chisq = unname(lrt$Chisq[2]), lrt_p = unname(lrt$`Pr(>|Chi|)`[2]),
    row.names = NULL)
}

TAB <- do.call(rbind, unlist(recursive = FALSE, lapply(CO, function(co)
  unlist(recursive = FALSE, lapply(names(SIGS), function(s)
    lapply(names(strata_of(co$d)), function(st)
      fit_block(strata_of(co$d)[[st]], make.names(s), st, co$label, s)))))))
rownames(TAB) <- NULL
print(TAB)
write.csv(TAB, file.path(DIR_RES, "table7_generality_values.csv"), row.names = FALSE)
for (i in seq_len(nrow(TAB)))
  log_msg(sprintf("06d %s | %s | %s (%s): HR %.2f (%.2f-%.2f) p=%.4f | C %.3f->%.3f | LRT p=%.4f [n=%d, d=%d]",
                  TAB$cohort[i], TAB$signature[i], TAB$stratum[i], TAB$covariates[i],
                  TAB$HR[i], TAB$lo[i], TAB$hi[i], TAB$p[i], TAB$c_base[i], TAB$c_full[i],
                  TAB$lrt_p[i], TAB$n[i], TAB$deaths[i]))

## ---- Q2: the IDH-mutant anomaly, resolved by grade --------------------------
anom <- do.call(rbind, lapply(CO, function(co) {
  d <- co$d[co$d$IDH == "Mutant", ]
  do.call(rbind, lapply(names(SIGS), function(s) {
    d$score <- as.numeric(scale(d[[make.names(s)]]))
    if (nrow(d) < 25 || sum(d$OS.event) < 12) return(NULL)
    f_age <- coxph(Surv(OS.time, OS.event) ~ score + age, data = d)
    gr_ok <- !all(is.na(d$grade)) && nlevels(droplevels(factor(d$grade))) > 1
    f_gr  <- if (gr_ok) coxph(Surv(OS.time, OS.event) ~ score + age + droplevels(factor(grade)),
                              data = d) else NULL
    b0 <- unname(coef(f_age)["score"]); b1 <- if (is.null(f_gr)) NA else unname(coef(f_gr)["score"])
    ci0 <- summary(f_age)$conf.int["score", ]
    ci1 <- if (is.null(f_gr)) rep(NA_real_, 4) else summary(f_gr)$conf.int["score", ]
    data.frame(cohort = co$label, signature = s, n = f_age$n, deaths = f_age$nevent,
               HR_age = unname(ci0[1]), lo_age = unname(ci0[3]), hi_age = unname(ci0[4]),
               p_age = summary(f_age)$coefficients["score", "Pr(>|z|)"],
               HR_grade = unname(ci1[1]), lo_grade = unname(ci1[3]), hi_grade = unname(ci1[4]),
               p_grade = if (is.null(f_gr)) NA_real_ else summary(f_gr)$coefficients["score", "Pr(>|z|)"],
               attenuation_by_grade = if (is.na(b1) || abs(b0) < 1e-8) NA_real_ else 1 - b1 / b0,
               row.names = NULL)
  }))
}))
rownames(anom) <- NULL
print(anom)
write.csv(anom, file.path(DIR_RES, "table8_idhmutant_grade_values.csv"), row.names = FALSE)
for (i in seq_len(nrow(anom)))
  log_msg(sprintf("06d ANOMALY %s | %s: age-only HR %.2f (p=%.4f) -> +grade HR %.2f (p=%.4f); attenuation %s",
                  anom$cohort[i], anom$signature[i], anom$HR_age[i], anom$p_age[i],
                  anom$HR_grade[i], anom$p_grade[i],
                  if (is.na(anom$attenuation_by_grade[i])) "NA"
                  else sprintf("%.0f%%", 100 * anom$attenuation_by_grade[i])))

## ---- Table 7 -----------------------------------------------------------------
t7 <- data.frame(
  Cohort = TAB$cohort, Signature = TAB$signature, Stratum = TAB$stratum,
  `Adjusted for` = TAB$covariates,
  `n (deaths)` = sprintf("%d (%d)", TAB$n, TAB$deaths),
  `HR per 1 SD (95% CI)` = sprintf("%.2f (%.2f, %.2f)", TAB$HR, TAB$lo, TAB$hi),
  `p-value` = ifelse(TAB$p < 0.001, "<0.001", sprintf("%.3f", TAB$p)),
  `C-index without / with signature` = sprintf("%.2f / %.2f", TAB$c_base, TAB$c_full),
  `LRT p-value` = ifelse(TAB$lrt_p < 0.001, "<0.001", sprintf("%.3f", TAB$lrt_p)),
  check.names = FALSE)
ft7 <- flextable::flextable(t7) %>% ft_style() %>%
  flextable::merge_v(j = c("Cohort", "Signature")) %>%
  flextable::add_footer_lines(values = c(
    "Two independently published signatures were reconstructed with their own published Cox coefficients and applied to log2-transformed, cohort-standardized expression: the ten glycolytic genes of Chen et al. (2017) and the ten glucose-metabolism genes of Zhang et al. (2021). Each score was standardized within stratum, so the hazard ratio is expressed per 1 standard deviation.",
    "Cohorts were analysed separately and never pooled. Within each stratum the signature was compared against a covariate-only model on identical complete cases; the likelihood-ratio test and the pair of C-indices quantify the information added by the signature beyond age and, where more than one grade was present, WHO grade.",
    "CI, confidence interval; HR, hazard ratio; LRT, likelihood-ratio test; SD, standard deviation."))
flextable::save_as_docx(ft7, path = file.path(DIR_TAB, "Table7_signature_generality_by_stratum.docx"))

## ---- Table 8 -----------------------------------------------------------------
t8 <- data.frame(
  Cohort = anom$cohort, Signature = anom$signature,
  `n (deaths)` = sprintf("%d (%d)", anom$n, anom$deaths),
  `HR, age only (95% CI)` = sprintf("%.2f (%.2f, %.2f); p=%s", anom$HR_age, anom$lo_age, anom$hi_age,
                                    ifelse(anom$p_age < 0.001, "<0.001", sprintf("%.3f", anom$p_age))),
  `HR, age + WHO grade (95% CI)` = ifelse(is.na(anom$HR_grade), "not estimable",
    sprintf("%.2f (%.2f, %.2f); p=%s", anom$HR_grade, anom$lo_grade, anom$hi_grade,
            ifelse(anom$p_grade < 0.001, "<0.001", sprintf("%.3f", anom$p_grade)))),
  `Attenuation by grade` = ifelse(is.na(anom$attenuation_by_grade), "NA",
                                  sprintf("%.0f%%", 100 * anom$attenuation_by_grade)),
  check.names = FALSE)
ft8 <- flextable::flextable(t8) %>% ft_style() %>%
  flextable::add_footer_lines(values = c(
    "Analysis restricted to IDH-mutant gliomas, which span WHO grades II to IV. Attenuation is the reduction of the log hazard ratio of the signature after WHO grade is added to the model.",
    "This stratum tests whether the residual prognostic value of the signatures within a single IDH class reflects tumour grade rather than glycolysis.",
    "CI, confidence interval; HR, hazard ratio."))
flextable::save_as_docx(ft8, path = file.path(DIR_TAB, "Table8_IDHmutant_grade_anomaly.docx"))

## ---- Figure 6 -----------------------------------------------------------------
fp <- TAB
fp$stratum <- factor(fp$stratum, levels = rev(c(
  "All gliomas (II-IV), both IDH states", "IDH-wildtype", "IDH-mutant",
  "IDH-wildtype, WHO grade IV")))
g6 <- ggplot(fp, aes(x = HR, y = stratum, colour = signature)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.18,
                 position = position_dodge(width = 0.55)) +
  geom_point(size = 2.4, position = position_dodge(width = 0.55)) +
  facet_wrap(~cohort, ncol = 1) + scale_x_log10() +
  scale_colour_manual(values = c("#0072B2", "#D55E00")) +
  labs(x = "Hazard ratio per 1 SD of signature score (95% CI)", y = NULL,
       colour = "Published signature") +
  theme_classic(base_size = 11) + theme(legend.position = "bottom")
save_fig(g6, file.path(DIR_FIG, "Figure6_generality_by_stratum"), width = 7.5, height = 6.5)

## ---- Persist -------------------------------------------------------------------
saveRDS(list(TAB = TAB, anom = anom, cohorts = CO, signatures = SIGS),
        file.path(DIR_RES, "generality_results.rds"))
nums <- c("", "--- Layer 06d: generality across two published signatures, and the IDH-mutant anomaly ---",
  sprintf("%s | %s | %s (adj: %s): HR %.2f (95%% CI %.2f-%.2f), p=%.4f | C %.3f->%.3f | LRT p=%.4f [n=%d, deaths=%d]",
          TAB$cohort, TAB$signature, TAB$stratum, TAB$covariates, TAB$HR, TAB$lo, TAB$hi, TAB$p,
          TAB$c_base, TAB$c_full, TAB$lrt_p, TAB$n, TAB$deaths),
  "",
  sprintf("ANOMALY (IDH-mutant only) %s | %s: age-only HR %.2f (p=%.4f); +grade HR %.2f (p=%.4f); attenuation by grade %s [n=%d, deaths=%d]",
          anom$cohort, anom$signature, anom$HR_age, anom$p_age, anom$HR_grade, anom$p_grade,
          ifelse(is.na(anom$attenuation_by_grade), "NA",
                 sprintf("%.0f%%", 100 * anom$attenuation_by_grade)), anom$n, anom$deaths))
cat(nums, sep = "\n")
write(nums, file = file.path(DIR_OUT, "manuscript_numbers.txt"), append = TRUE)
log_msg("Layer 06d complete: Table 7, Table 8 and Figure 6 exported to Sep_06/.")
