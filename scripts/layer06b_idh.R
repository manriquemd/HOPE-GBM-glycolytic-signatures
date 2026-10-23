## =====================================================================
## Layer 06b - Does the prognostic value of the Chen 2017 GGESS reflect
##             IDH1 status rather than glycolysis?
##
## Rationale. Chen et al. (J Genet Genomics 2017;44:519-530) trained the
## GGESS on 512 TCGA GBM cases that, under the pre-2021 definition,
## included IDH-mutant tumours, and they themselves report that GGESS is
## lower in IDH1-mutant GBM. IDH-mutant tumours therefore carry both a low
## score and a markedly better survival, which is sufficient to generate a
## prognostic association in a mixed cohort without any contribution from
## glycolysis. Layer 06 showed the published signature is null in
## IDH-wildtype cohorts (HR 0.99 in TCGA, 1.00 in CGGA-693). This layer
## tests the proposed explanation directly.
##
## Design (TCGA only; CGGA-325 and CGGA-693 were already filtered to
## IDH-wildtype upstream and are not re-derived here):
##   Cohort A - all GBM with valid OS, IDH-wildtype and IDH-mutant pooled
##              (approximates the cohort Chen trained on)
##   Model 1  - GGESS + age + MGMT                     (mixed cohort)
##   Model 2  - GGESS + age + MGMT + IDH status        (mixed cohort)
##   Model 3  - GGESS + age + MGMT, IDH-wildtype only
##   Model 4  - GGESS + age + MGMT, IDH-mutant only    (reported only if n allows)
##
## If the signal is present in Model 1, attenuates in Model 2 and is absent
## in Model 3, the prognostic value of the signature is attributable to IDH
## status. Attenuation is reported ONLY when the base model is significant;
## the ratio is uninterpretable when the base coefficient is near zero.
##
## Usage:
##   setwd("~/Library/CloudStorage/Dropbox/HOPE-GBM"); source("layer06b_idh.R")
## Requires: a completed run of HOPE-GBM_11_Sep06.Rmd (data_raw/, results/).
## =====================================================================

set.seed(1024)
suppressPackageStartupMessages({
  library(SummarizedExperiment); library(DESeq2); library(survival)
  library(dplyr); library(ggplot2); library(flextable); library(officer)
})

PROJ     <- getwd()
DIR_RAW  <- file.path(PROJ, "data_raw")
DIR_PROC <- file.path(PROJ, "data_processed")
DIR_RES  <- file.path(PROJ, "results")
DIR_OUT  <- file.path(PROJ, "Sep_06")
DIR_TAB  <- file.path(DIR_OUT, "Tables")
DIR_FIG  <- file.path(DIR_OUT, "Figures")
stopifnot(dir.exists(DIR_RAW), dir.exists(DIR_RES))

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
harmonize_mgmt <- function(x) {
  x <- gsub("[^a-z]", "", tolower(as.character(x)))
  factor(ifelse(x == "methylated", "Methylated",
         ifelse(x == "unmethylated", "Unmethylated", NA_character_)),
         levels = c("Methylated", "Unmethylated"))
}
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

## ---- 1. Rebuild the TCGA cohort WITHOUT the IDH filter -------------------
se      <- readRDS(file.path(DIR_RAW, "tcga_gbm_SE.rds"))
subtype <- readRDS(file.path(DIR_RAW, "tcga_gbm_subtype.rds"))

se <- se[, substr(colnames(se), 14, 15) == "01"]
pat <- substr(colnames(se), 1, 12)
se <- se[, !duplicated(pat)]
colnames(se) <- substr(colnames(se), 1, 12)
log_msg(sprintf("06b: %d primary, de-duplicated TCGA patients (no IDH filter).", ncol(se)))

cd <- as.data.frame(colData(se))
os_time <- ifelse(!is.na(cd$days_to_death), cd$days_to_death, cd$days_to_last_follow_up)
clin <- data.frame(
  patient  = colnames(se),
  OS.time  = as.numeric(os_time),
  OS.event = ifelse(tolower(cd$vital_status) == "dead", 1L, 0L),
  age      = as.numeric(cd$age_at_index),
  stringsAsFactors = FALSE)

idh_raw <- as.character(subtype[match(clin$patient, subtype$patient), "IDH.status"])
clin$IDH <- factor(ifelse(grepl("^wt|wild", tolower(idh_raw)), "Wildtype",
                   ifelse(grepl("mut", tolower(idh_raw)), "Mutant", NA_character_)),
                   levels = c("Wildtype", "Mutant"))
mcol <- grep("MGMT", colnames(subtype), value = TRUE, ignore.case = TRUE)[1]
clin$MGMT <- harmonize_mgmt(subtype[match(clin$patient, subtype$patient), mcol])

clin <- clin[!is.na(clin$OS.time) & clin$OS.time > 0, ]
se   <- se[, clin$patient]
log_msg(sprintf("06b: valid OS n=%d | IDH wildtype=%d, mutant=%d, unknown=%d",
                nrow(clin), sum(clin$IDH == "Wildtype", na.rm = TRUE),
                sum(clin$IDH == "Mutant", na.rm = TRUE), sum(is.na(clin$IDH))))

## ---- 2. VST on the mixed cohort ------------------------------------------
cts <- assay(se, "unstranded"); mode(cts) <- "integer"
dds <- DESeqDataSetFromMatrix(cts, colData = DataFrame(clin), design = ~ 1)
dds <- dds[rowSums(counts(dds)) >= 10, ]
vst_all <- assay(vst(dds, blind = TRUE))
gsym <- rowData(se)$gene_name[match(rownames(vst_all), rownames(se))]
ok <- !is.na(gsym) & gsym != ""
vst_all <- vst_all[ok, ]; gsym <- gsym[ok]
ord <- order(rowMeans(vst_all), decreasing = TRUE)
vst_all <- vst_all[ord, ]; gsym <- gsym[ord]
vst_all <- vst_all[!duplicated(gsym), ]; rownames(vst_all) <- gsym[!duplicated(gsym)]
log_msg(sprintf("06b: VST matrix %d genes x %d patients.", nrow(vst_all), ncol(vst_all)))

## ---- 3. GGESS with published weights -------------------------------------
genes <- intersect(names(CHEN2017_GGESS_COEF), rownames(vst_all))
if (length(genes) < length(CHEN2017_GGESS_COEF))
  log_msg("06b: GGESS genes absent from the matrix: ",
          paste(setdiff(names(CHEN2017_GGESS_COEF), genes), collapse = ", "))
Z <- t(apply(vst_all[genes, , drop = FALSE], 1, zrow))
raw <- as.numeric(colSums(Z * CHEN2017_GGESS_COEF[genes]))
names(raw) <- colnames(Z)
clin$ggess <- as.numeric(scale(raw[clin$patient]))

## ---- 4. Models ------------------------------------------------------------
row_of <- function(fit, label, n_note = NA_character_) {
  s <- summary(fit); ci <- s$conf.int["ggess", ]
  data.frame(model = label, n = fit$n, deaths = fit$nevent,
             HR = unname(ci[1]), lo = unname(ci[3]), hi = unname(ci[4]),
             p = s$coefficients["ggess", "Pr(>|z|)"],
             beta = unname(coef(fit)["ggess"]),
             note = n_note, row.names = NULL)
}
d_mix <- clin[complete.cases(clin[, c("ggess", "age", "MGMT", "IDH")]), ]
log_msg(sprintf("06b: complete-case mixed cohort n=%d (wt=%d, mut=%d)", nrow(d_mix),
                sum(d_mix$IDH == "Wildtype"), sum(d_mix$IDH == "Mutant")))

m1 <- coxph(Surv(OS.time, OS.event) ~ ggess + age + MGMT, data = d_mix)
m2 <- coxph(Surv(OS.time, OS.event) ~ ggess + age + MGMT + IDH, data = d_mix)
d_wt  <- d_mix[d_mix$IDH == "Wildtype", ]
m3 <- coxph(Surv(OS.time, OS.event) ~ ggess + age + MGMT, data = d_wt)
d_mut <- d_mix[d_mix$IDH == "Mutant", ]
m4 <- if (nrow(d_mut) >= 20 && sum(d_mut$OS.event) >= 10)
        coxph(Surv(OS.time, OS.event) ~ ggess + age, data = d_mut) else NULL

res <- rbind(
  row_of(m1, "Mixed cohort: age + MGMT"),
  row_of(m2, "Mixed cohort: age + MGMT + IDH status"),
  row_of(m3, "IDH-wildtype only: age + MGMT"))
if (!is.null(m4)) res <- rbind(res, row_of(m4, "IDH-mutant only: age"))
if (is.null(m4))
  log_msg(sprintf("06b: IDH-mutant subgroup not modelled (n=%d, deaths=%d): insufficient events.",
                  nrow(d_mut), sum(d_mut$OS.event)))

## Attenuation reported ONLY when the base model is significant.
att <- if (res$p[1] < 0.05) 1 - res$beta[2] / res$beta[1] else NA_real_
log_msg(sprintf("06b: base p=%.3f -> attenuation %s", res$p[1],
                if (is.na(att)) "NOT REPORTED (base model null)" else sprintf("%.0f%%", 100 * att)))

## GGESS by IDH status, and survival by IDH status (the two premises)
kw_idh <- wilcox.test(ggess ~ IDH, data = d_mix)
cox_idh <- coxph(Surv(OS.time, OS.event) ~ IDH + age, data = d_mix)
ci_idh <- summary(cox_idh)$conf.int["IDHMutant", ]
log_msg(sprintf("06b: GGESS by IDH, Wilcoxon p=%.3g | median wt=%.2f, mut=%.2f",
                kw_idh$p.value, median(d_mix$ggess[d_mix$IDH == "Wildtype"]),
                median(d_mix$ggess[d_mix$IDH == "Mutant"])))
log_msg(sprintf("06b: IDH-mutant vs wildtype survival: HR %.2f (%.2f-%.2f), p=%.3g",
                ci_idh[1], ci_idh[3], ci_idh[4],
                summary(cox_idh)$coefficients["IDHMutant", "Pr(>|z|)"]))
for (i in seq_len(nrow(res)))
  log_msg(sprintf("06b: %s -> HR %.2f (%.2f-%.2f), p=%.3f, n=%d, deaths=%d",
                  res$model[i], res$HR[i], res$lo[i], res$hi[i], res$p[i], res$n[i], res$deaths[i]))

## ---- 5. Table 5 ------------------------------------------------------------
t5 <- data.frame(
  `Cox model` = res$model,
  `n (deaths)` = sprintf("%d (%d)", res$n, res$deaths),
  `HR per 1 SD of GGESS (95% CI)` = sprintf("%.2f (%.2f, %.2f)", res$HR, res$lo, res$hi),
  `p-value` = sprintf("%.3f", res$p),
  check.names = FALSE)
foot <- c(
  sprintf("TCGA-GBM, primary tumours with available overall survival, IDH status and MGMT promoter status (n = %d; IDH-wildtype %d, IDH-mutant %d). The signature was computed with the ten published weights of Chen et al. (2017) applied to cohort-standardized expression, then standardized, so the hazard ratio is expressed per 1 standard deviation.",
          nrow(d_mix), sum(d_mix$IDH == "Wildtype"), sum(d_mix$IDH == "Mutant")),
  sprintf("GGESS was lower in IDH-mutant than in IDH-wildtype tumours (Wilcoxon p = %.3g), and IDH-mutant status was associated with longer survival (HR %.2f, 95%% CI %.2f to %.2f).",
          kw_idh$p.value, ci_idh[1], ci_idh[3], ci_idh[4]),
  if (is.na(att))
    "Attenuation of the coefficient is not reported because the model without IDH status did not show a significant association; the ratio is uninterpretable when the base coefficient approaches zero."
  else sprintf("Adding IDH status to the model reduced the log hazard ratio of the signature by %.0f%%.", 100 * att),
  "CI, confidence interval; GGESS, glycolytic gene expression signature score; HR, hazard ratio; SD, standard deviation.")
ft5 <- flextable::flextable(t5) %>% ft_style() %>% flextable::add_footer_lines(values = foot)
flextable::save_as_docx(ft5, path = file.path(DIR_TAB, "Table5_GGESS_IDH_explanation.docx"))

## ---- 6. Figure 4: two panels, the premises and the models -----------------
p_left <- ggplot(d_mix, aes(IDH, ggess)) +
  geom_boxplot(outlier.shape = NA, width = 0.55) +
  geom_jitter(width = 0.14, alpha = 0.3, size = 0.8) +
  labs(x = "IDH status", y = "GGESS, published weights (z)",
       subtitle = sprintf("Wilcoxon p = %.2g", kw_idh$p.value)) +
  theme_classic(base_size = 11)

fp <- res
fp$model <- factor(fp$model, levels = rev(fp$model))
p_right <- ggplot(fp, aes(x = HR, y = model)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.18) +
  geom_point(size = 2.6) +
  scale_x_log10() +
  labs(x = "HR per 1 SD of GGESS (95% CI)", y = NULL,
       subtitle = "Prognostic value of the published signature") +
  theme_classic(base_size = 11) +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())

g4 <- if (requireNamespace("patchwork", quietly = TRUE)) {
  patchwork::wrap_plots(p_left, p_right, widths = c(1, 1.6))
} else { print(p_left); p_right }
save_fig(g4, file.path(DIR_FIG, "Figure4_GGESS_IDH_explanation"), width = 9, height = 3.4)

## ---- 7. Persist and append to the numbers file ----------------------------
saveRDS(list(res = res, att = att, kw_idh = kw_idh, cox_idh = cox_idh,
             n_mix = nrow(d_mix), n_wt = sum(d_mix$IDH == "Wildtype"),
             n_mut = sum(d_mix$IDH == "Mutant"), genes = genes),
        file.path(DIR_RES, "ggess_idh_results.rds"))
write.csv(res, file.path(DIR_RES, "table5_ggess_idh_values.csv"), row.names = FALSE)

nums <- c("", "--- Layer 06b: IDH explanation for the Chen 2017 GGESS ---",
  sprintf("Mixed TCGA cohort: n=%d (IDH-wildtype %d, IDH-mutant %d), deaths=%d",
          nrow(d_mix), sum(d_mix$IDH == "Wildtype"), sum(d_mix$IDH == "Mutant"), sum(d_mix$OS.event)),
  sprintf("GGESS median: wildtype %.2f, mutant %.2f; Wilcoxon p=%.3g",
          median(d_mix$ggess[d_mix$IDH == "Wildtype"]), median(d_mix$ggess[d_mix$IDH == "Mutant"]),
          kw_idh$p.value),
  sprintf("IDH-mutant vs wildtype survival: HR %.2f (95%% CI %.2f-%.2f), p=%.3g",
          ci_idh[1], ci_idh[3], ci_idh[4], summary(cox_idh)$coefficients["IDHMutant", "Pr(>|z|)"]),
  sprintf("%s: HR %.2f (95%% CI %.2f-%.2f), p=%.3f [n=%d, deaths=%d]",
          res$model, res$HR, res$lo, res$hi, res$p, res$n, res$deaths),
  if (is.na(att)) "Attenuation not reported: base model null (ratio uninterpretable)."
  else sprintf("Attenuation of log HR after adding IDH status: %.0f%%", 100 * att))
cat(nums, sep = "\n")
write(nums, file = file.path(DIR_OUT, "manuscript_numbers.txt"), append = TRUE)
log_msg("Layer 06b complete: Table 5 and Figure 4 exported to Sep_06/.")
