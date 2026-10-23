## =====================================================================
## Layer 06c - IDH explanation for the Chen 2017 GGESS, tested in CGGA
##             gliomas of all grades (adequate IDH-mutant representation)
##
## Why CGGA rather than TCGA. Layer 06b confirmed both premises in TCGA
## (GGESS lower in IDH-mutant tumours, IDH-mutant tumours survive longer)
## but could not demonstrate the consequence: only 18 of 209 patients with
## complete data were IDH-mutant, and the signature was not significant
## even in the mixed cohort. Chen et al. trained on 512 cases. The CGGA
## mRNAseq cohorts contain gliomas of WHO grades II to IV with IDH status,
## where IDH-mutant tumours are abundant, which is the composition needed
## for the dilution argument to be testable.
##
## Design, run separately in CGGA-693 and CGGA-325 (never pooled):
##   Cohort   - all gliomas with valid OS and known IDH status; primary
##              tumours only, matching the filter used upstream
##   Model 1  - GGESS + age                    (mixed-grade, mixed-IDH)
##   Model 2  - GGESS + age + IDH status
##   Model 3  - GGESS + age + IDH + WHO grade  (grade is the obvious confounder)
##   Model 4  - GGESS + age, IDH-wildtype only
##   Model 5  - GGESS + age, IDH-mutant only
##
## Interpretation. A strong, significant HR in Model 1 that attenuates in
## Model 2 and disappears in Model 4 supports the claim that the published
## signature tracks IDH status. Model 3 is the honest control: much of what
## IDH captures is grade, and the manuscript should say so rather than
## overclaim a purely IDH-driven mechanism.
##
## Attenuation is computed only when the base model is significant.
## MGMT is not used here: it is incomplete across mixed-grade CGGA and its
## inclusion would drop the cases the design depends on.
##
## Usage:
##   setwd("~/Library/CloudStorage/Dropbox/HOPE-GBM"); source("layer06c_cgga_idh.R")
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
CGGA_CANDIDATES <- c(path.expand("~/Documents/CGGA_local"),
                     path.expand("~/Library/CloudStorage/Dropbox/00_DATA/CGGA"),
                     file.path(PROJ, "data_raw", "CGGA"))
CGGA_DIR <- Find(function(p) dir.exists(p) &&
                   length(list.files(p, pattern = "mRNAseq.*genes.*[.]txt$")) > 0, CGGA_CANDIDATES)
stopifnot(!is.null(CGGA_DIR), dir.exists(DIR_TAB))

CHEN2017_GGESS_COEF <- c(ALDOB = -0.414, ENO1 = 0.218, GALM = 0.242,
                         GAPDH = 0.194, HK2 = 0.211, HK3 = 0.319,
                         LDHA = 0.269, LDHB = -0.296, PKLR = -0.181,
                         PKM = 0.213)   # Chen 2017, Table 1; PKM2 = HGNC PKM
FIG_DPI <- 600

log_msg <- function(...) {
  m <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ...)
  cat(m, "\n"); cat(m, "\n", file = file.path(PROJ, "logs", "hope_run.log"), append = TRUE)
}
zrow <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
ft_style <- function(ft) {
  ft %>% flextable::font(fontname = "Helvetica", part = "all") %>%
    flextable::fontsize(size = 11, part = "all") %>% flextable::border_remove() %>%
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
       compression = "lzw")
  print(obj); dev.off()
}
## Same byte-level reader validated in the v3.1 pipeline
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
  log_msg(sprintf("06c %s: %d genes x %d samples", label, nrow(expr), ncol(expr)))
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

## ---- Per-cohort analysis ---------------------------------------------------
analyse <- function(expr_path, clin_path) {
  label <- paste0("CGGA-", num_of(expr_path))
  expr  <- read_cgga_expr(expr_path, label)
  clin  <- as.data.frame(readr::read_tsv(clin_path, show_col_types = FALSE, progress = FALSE))
  nm <- colnames(clin)
  id_col   <- pick_col(nm, "CGGA_ID", "CGGA|sample|ID")
  prs_col  <- pick_col(nm, "PRS_type", "PRS")
  idh_col  <- pick_col(nm, "IDH_mutation_status", "IDH")
  grade_col<- pick_col(nm, "Grade", "Grade|WHO")
  age_col  <- pick_col(nm, "Age", "Age")
  os_col   <- pick_col(nm, "OS", "^OS$|overall")
  ev_col   <- grep("Censor", nm, value = TRUE, ignore.case = TRUE)[1]

  n0 <- nrow(clin)
  keep <- rep(TRUE, n0)
  if (!is.na(prs_col)) keep <- keep & clin[[prs_col]] == "Primary"
  keep[is.na(keep)] <- FALSE
  clin <- clin[keep, ]
  ids  <- intersect(clin[[id_col]], colnames(expr))
  expr <- expr[, ids, drop = FALSE]; clin <- clin[match(ids, clin[[id_col]]), ]

  genes <- intersect(names(CHEN2017_GGESS_COEF), rownames(expr))
  if (length(genes) < length(CHEN2017_GGESS_COEF))
    log_msg(sprintf("06c %s: GGESS genes absent: %s", label,
                    paste(setdiff(names(CHEN2017_GGESS_COEF), genes), collapse = ", ")))
  Z   <- t(apply(log2(expr[genes, , drop = FALSE] + 1), 1, zrow))
  raw <- as.numeric(colSums(Z * CHEN2017_GGESS_COEF[genes])); names(raw) <- colnames(Z)

  idh_raw <- tolower(as.character(clin[[idh_col]]))
  d <- data.frame(
    patient  = ids,
    OS.time  = as.numeric(clin[[os_col]]),
    OS.event = as.integer(clin[[ev_col]]),
    age      = suppressWarnings(as.numeric(clin[[age_col]])),
    IDH      = factor(ifelse(grepl("wildtype|^wt", idh_raw), "Wildtype",
                      ifelse(grepl("mutant|^mut", idh_raw), "Mutant", NA_character_)),
                      levels = c("Wildtype", "Mutant")),
    grade    = if (!is.na(grade_col)) factor(gsub("\\s", "", as.character(clin[[grade_col]]))) else NA,
    ggess    = as.numeric(raw[ids]),
    stringsAsFactors = FALSE)
  d <- d[!is.na(d$OS.time) & d$OS.time > 0 & !is.na(d$IDH) & !is.na(d$age), ]
  d$ggess <- as.numeric(scale(d$ggess))
  log_msg(sprintf("06c %s: %d rows -> analysed n=%d (wildtype %d, mutant %d), deaths=%d",
                  label, n0, nrow(d), sum(d$IDH == "Wildtype"), sum(d$IDH == "Mutant"),
                  sum(d$OS.event)))
  if (!all(is.na(d$grade)))
    log_msg(sprintf("06c %s: grade distribution %s", label,
                    paste(names(table(d$grade)), table(d$grade), sep = "=", collapse = ", ")))

  row_of <- function(fit, lab) {
    s <- summary(fit); ci <- s$conf.int["ggess", ]
    data.frame(cohort = label, model = lab, n = fit$n, deaths = fit$nevent,
               HR = unname(ci[1]), lo = unname(ci[3]), hi = unname(ci[4]),
               p = s$coefficients["ggess", "Pr(>|z|)"], beta = unname(coef(fit)["ggess"]),
               row.names = NULL)
  }
  m1 <- coxph(Surv(OS.time, OS.event) ~ ggess + age, data = d)
  m2 <- coxph(Surv(OS.time, OS.event) ~ ggess + age + IDH, data = d)
  out <- rbind(row_of(m1, "Mixed cohort: age"),
               row_of(m2, "Mixed cohort: age + IDH"))
  has_grade <- !all(is.na(d$grade)) && nlevels(droplevels(d$grade)) > 1
  if (has_grade) {
    m3 <- coxph(Surv(OS.time, OS.event) ~ ggess + age + IDH + grade, data = d)
    out <- rbind(out, row_of(m3, "Mixed cohort: age + IDH + grade"))
  }
  d_wt  <- d[d$IDH == "Wildtype", ]; d_mut <- d[d$IDH == "Mutant", ]
  if (nrow(d_wt)  >= 20 && sum(d_wt$OS.event)  >= 10)
    out <- rbind(out, row_of(coxph(Surv(OS.time, OS.event) ~ ggess + age, data = d_wt),
                             "IDH-wildtype only: age"))
  if (nrow(d_mut) >= 20 && sum(d_mut$OS.event) >= 10)
    out <- rbind(out, row_of(coxph(Surv(OS.time, OS.event) ~ ggess + age, data = d_mut),
                             "IDH-mutant only: age"))

  base_p <- out$p[1]
  att_idh <- if (base_p < 0.05) 1 - out$beta[2] / out$beta[1] else NA_real_
  att_gr  <- if (base_p < 0.05 && has_grade) 1 - out$beta[3] / out$beta[1] else NA_real_
  wil <- wilcox.test(ggess ~ IDH, data = d)
  cox_idh <- coxph(Surv(OS.time, OS.event) ~ IDH + age, data = d)
  ci_idh  <- summary(cox_idh)$conf.int["IDHMutant", ]
  log_msg(sprintf("06c %s: GGESS median wt=%.2f mut=%.2f, Wilcoxon p=%.3g | IDH-mut survival HR %.2f (%.2f-%.2f)",
                  label, median(d$ggess[d$IDH == "Wildtype"]), median(d$ggess[d$IDH == "Mutant"]),
                  wil$p.value, ci_idh[1], ci_idh[3], ci_idh[4]))
  for (i in seq_len(nrow(out)))
    log_msg(sprintf("06c %s | %s -> HR %.2f (%.2f-%.2f), p=%.4f [n=%d, d=%d]", label,
                    out$model[i], out$HR[i], out$lo[i], out$hi[i], out$p[i], out$n[i], out$deaths[i]))
  log_msg(sprintf("06c %s: attenuation by IDH %s; by IDH + grade %s", label,
                  if (is.na(att_idh)) "NOT REPORTED (base null)" else sprintf("%.0f%%", 100 * att_idh),
                  if (is.na(att_gr))  "NOT REPORTED" else sprintf("%.0f%%", 100 * att_gr)))
  list(label = label, out = out, d = d, att_idh = att_idh, att_gr = att_gr,
       wil = wil, ci_idh = ci_idh)
}

res <- Map(analyse, expr_f, clin_f)
names(res) <- vapply(res, `[[`, "", "label")
tab <- do.call(rbind, lapply(res, `[[`, "out"))

## ---- Table 6 ---------------------------------------------------------------
t6 <- data.frame(
  Cohort = tab$cohort, `Cox model` = tab$model,
  `n (deaths)` = sprintf("%d (%d)", tab$n, tab$deaths),
  `HR per 1 SD of GGESS (95% CI)` = sprintf("%.2f (%.2f, %.2f)", tab$HR, tab$lo, tab$hi),
  `p-value` = ifelse(tab$p < 0.001, "<0.001", sprintf("%.3f", tab$p)),
  check.names = FALSE)
foot <- c(
  "CGGA mRNAseq cohorts, primary gliomas of WHO grades II to IV with known IDH status and available overall survival, analysed separately and never pooled. The signature was computed with the ten published weights of Chen et al. (2017) applied to log2-transformed, cohort-standardized expression, then standardized; the hazard ratio is expressed per 1 standard deviation.",
  paste(vapply(res, function(r) sprintf(
    "%s: GGESS median %.2f in IDH-wildtype and %.2f in IDH-mutant gliomas (Wilcoxon p = %s); IDH-mutant status was associated with longer survival (HR %.2f, 95%% CI %.2f to %.2f).",
    r$label, median(r$d$ggess[r$d$IDH == "Wildtype"]), median(r$d$ggess[r$d$IDH == "Mutant"]),
    format.pval(r$wil$p.value, digits = 2, eps = 1e-16),
    r$ci_idh[1], r$ci_idh[3], r$ci_idh[4]), ""), collapse = " "),
  paste(vapply(res, function(r) if (is.na(r$att_idh))
      sprintf("%s: attenuation is not reported because the model without IDH status was not significant.", r$label)
    else sprintf("%s: adding IDH status reduced the log hazard ratio of the signature by %.0f%%%s.",
                 r$label, 100 * r$att_idh,
                 if (is.na(r$att_gr)) "" else sprintf(", and adding IDH status together with WHO grade reduced it by %.0f%%", 100 * r$att_gr)), ""),
    collapse = " "),
  "WHO grade is retained as a control covariate because IDH status and grade are strongly associated; the grade-adjusted model indicates how much of the attenuation is attributable to IDH status beyond tumour grade.",
  "CI, confidence interval; GGESS, glycolytic gene expression signature score; HR, hazard ratio; SD, standard deviation.")
ft6 <- flextable::flextable(t6) %>% ft_style() %>%
  flextable::merge_v(j = "Cohort") %>% flextable::add_footer_lines(values = foot)
flextable::save_as_docx(ft6, path = file.path(DIR_TAB, "Table6_GGESS_IDH_CGGA.docx"))

## ---- Figure 5 ---------------------------------------------------------------
dl <- dplyr::bind_rows(lapply(res, function(r)
  data.frame(cohort = r$label, IDH = r$d$IDH, ggess = r$d$ggess)))
p_left <- ggplot(dl, aes(IDH, ggess)) +
  geom_boxplot(outlier.shape = NA, width = 0.55) +
  geom_jitter(width = 0.14, alpha = 0.25, size = 0.7) +
  facet_wrap(~cohort) +
  labs(x = "IDH status", y = "GGESS, published weights (z)") +
  theme_classic(base_size = 11)

fp <- tab
fp$lab <- paste(fp$cohort, fp$model, sep = " | ")
fp$lab <- factor(fp$lab, levels = rev(fp$lab))
p_right <- ggplot(fp, aes(x = HR, y = lab)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.18) +
  geom_point(size = 2.4) + scale_x_log10() +
  labs(x = "HR per 1 SD of GGESS (95% CI)", y = NULL) +
  theme_classic(base_size = 10) +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())

g5 <- if (requireNamespace("patchwork", quietly = TRUE))
  patchwork::wrap_plots(p_left, p_right, widths = c(1, 1.3)) else { print(p_left); p_right }
save_fig(g5, file.path(DIR_FIG, "Figure5_GGESS_IDH_CGGA"),
         width = 10, height = max(3.6, 0.32 * nrow(fp) + 2))

## ---- Persist ----------------------------------------------------------------
saveRDS(res, file.path(DIR_RES, "ggess_idh_cgga_results.rds"))
write.csv(tab, file.path(DIR_RES, "table6_ggess_idh_cgga_values.csv"), row.names = FALSE)

nums <- c("", "--- Layer 06c: IDH explanation tested in CGGA gliomas (all grades) ---",
  unlist(lapply(res, function(r) c(
    sprintf("%s: n=%d (wildtype %d, mutant %d), deaths=%d", r$label, nrow(r$d),
            sum(r$d$IDH == "Wildtype"), sum(r$d$IDH == "Mutant"), sum(r$d$OS.event)),
    sprintf("%s: GGESS median wildtype %.2f, mutant %.2f; Wilcoxon p=%.3g", r$label,
            median(r$d$ggess[r$d$IDH == "Wildtype"]), median(r$d$ggess[r$d$IDH == "Mutant"]),
            r$wil$p.value),
    sprintf("%s: IDH-mutant vs wildtype survival HR %.2f (95%% CI %.2f-%.2f)", r$label,
            r$ci_idh[1], r$ci_idh[3], r$ci_idh[4]),
    sprintf("%s | %s: HR %.2f (95%% CI %.2f-%.2f), p=%.4f [n=%d, deaths=%d]",
            r$out$cohort, r$out$model, r$out$HR, r$out$lo, r$out$hi, r$out$p, r$out$n, r$out$deaths),
    sprintf("%s: attenuation by IDH %s; by IDH + grade %s", r$label,
            if (is.na(r$att_idh)) "not reported (base null)" else sprintf("%.0f%%", 100 * r$att_idh),
            if (is.na(r$att_gr))  "not reported" else sprintf("%.0f%%", 100 * r$att_gr))))))
cat(nums, sep = "\n")
write(nums, file = file.path(DIR_OUT, "manuscript_numbers.txt"), append = TRUE)
log_msg("Layer 06c complete: Table 6 and Figure 5 exported to Sep_06/.")
