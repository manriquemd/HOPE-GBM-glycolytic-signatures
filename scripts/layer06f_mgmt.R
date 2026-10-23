## =====================================================================
## Layer 06f - Is the CGGA-325 signal explained by MGMT promoter status?
##
## Rationale. Layer 06e left MGMT as the only composition difference that
## met the pre-specified relevance threshold: 32% methylated in CGGA-325
## versus 57% in CGGA-693 (25 percentage points, p=0.003). The direction
## is coherent with the discordance, because CGGA-325 is enriched for the
## unmethylated, poorer-prognosis phenotype. If the glycolytic score is
## associated with MGMT status, the HR of 1.38 observed in CGGA-325 may be
## capturing MGMT indirectly rather than carrying independent prognostic
## information.
##
## The Layer 06d models adjusted for age only. This layer adds MGMT and
## measures what happens.
##
## Design decisions that matter for interpretation:
##
##  1. LIKE-FOR-LIKE BASELINE. MGMT is missing in a fraction of patients.
##     Adding MGMT to a model silently drops those patients, so a naive
##     comparison confounds adjustment with a change of sample. Every
##     model in this layer is fitted on the MGMT-complete subset, and the
##     age-only model is refitted on that same subset. The attenuation
##     ratio is therefore attributable to adjustment, not to attrition.
##     The original full-sample estimate is reported alongside so the
##     effect of restriction is visible.
##
##  2. ATTENUATION IS ONLY REPORTED WHEN THE BASE MODEL IS NON-NULL.
##     A ratio computed from a base log-HR near zero is an artefact, as
##     Layer 06 demonstrated (-1126%, -10606%). The threshold used here is
##     a base HR at least 1.10 with p < 0.10; otherwise the cell reads
##     "not applicable (base model null)".
##
##  3. POWER IS RECOMPUTED ON THE RESTRICTED SUBSET. Dropping patients
##     with unknown MGMT reduces events, so the detectable effect size
##     changes and must be restated rather than carried over from 06e.
##
##  4. THE INTERACTION TEST IS EXPLORATORY AND LABELLED AS SUCH. A score
##     by MGMT interaction in fewer than a hundred patients is
##     underpowered by construction. It is computed to be reported, not
##     to be believed.
##
## Interpretation set in advance, to avoid a post hoc reading:
##   - Attenuation >= 30% with loss of significance -> MGMT confounding is
##     the working explanation for the discordance.
##   - Attenuation < 30% and signal persists -> MGMT is not the
##     explanation; the discordance goes to Limitations unresolved.
##   - Score not associated with MGMT at all -> confounding is ruled out
##     regardless of what the adjusted model shows.
##
## Usage:
##   setwd("~/Library/CloudStorage/Dropbox/HOPE-GBM"); source("layer06f_mgmt.R")
## Requires Layer 06e (results/discordance_results.rds).
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

disc <- readRDS(file.path(DIR_RES, "discordance_results.rds"))
G4   <- disc$G4
SIG_COLS <- c("Chen.2017", "Zhang.2021")
SIG_LAB  <- c(Chen.2017 = "Chen 2017", Zhang.2021 = "Zhang 2021")
stopifnot(all(SIG_COLS %in% names(G4)))

log_msg <- function(...) {
  m <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ...)
  cat(m, "\n"); cat(m, "\n", file = file.path(PROJ, "logs", "hope_run.log"), append = TRUE)
}
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
detectable_hr <- function(events, power = 0.80, alpha = 0.05)
  exp((qnorm(1 - alpha / 2) + qnorm(power)) / sqrt(events))

## Standardize each score within cohort so HR is per 1 SD of that cohort.
for (s in SIG_COLS)
  G4[[paste0(s, ".z")]] <- ave(G4[[s]], G4$cohort,
                               FUN = function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))

## ---- Step 1. Missingness of MGMT, reported before anything is modelled -------
miss <- do.call(rbind, lapply(levels(G4$cohort), function(k) {
  d <- G4[G4$cohort == k, ]
  data.frame(cohort = k, n_total = nrow(d), n_mgmt_known = sum(!is.na(d$MGMT)),
             pct_missing = 100 * mean(is.na(d$MGMT)),
             events_total = sum(d$OS.event),
             events_mgmt_known = sum(d$OS.event[!is.na(d$MGMT)]), row.names = NULL)
}))
miss$HR_detectable_restricted <- detectable_hr(miss$events_mgmt_known)
print(miss)
for (i in seq_len(nrow(miss)))
  log_msg(sprintf("06f %s: MGMT known in %d/%d (%.0f%% missing); events %d -> %d; detectable HR >= %.2f",
                  miss$cohort[i], miss$n_mgmt_known[i], miss$n_total[i], miss$pct_missing[i],
                  miss$events_total[i], miss$events_mgmt_known[i], miss$HR_detectable_restricted[i]))

## ---- Step 2. Is the score associated with MGMT at all? -----------------------
## If it is not, confounding by MGMT cannot be the mechanism, whatever the
## adjusted model happens to show.
assoc <- do.call(rbind, lapply(levels(G4$cohort), function(k) {
  d <- G4[G4$cohort == k & !is.na(G4$MGMT), ]
  do.call(rbind, lapply(SIG_COLS, function(s) {
    x <- d[[paste0(s, ".z")]]
    m <- d$MGMT
    w <- tryCatch(wilcox.test(x ~ m), error = function(e) NULL)
    data.frame(cohort = k, signature = SIG_LAB[s],
               median_meth  = median(x[m == "Methylated"], na.rm = TRUE),
               median_unmeth = median(x[m == "Unmethylated"], na.rm = TRUE),
               delta_SD = median(x[m == "Unmethylated"], na.rm = TRUE) -
                          median(x[m == "Methylated"], na.rm = TRUE),
               p_wilcox = if (is.null(w)) NA_real_ else w$p.value,
               rho_pointbiserial = suppressWarnings(
                 cor(x, as.numeric(m == "Unmethylated"), method = "spearman", use = "complete.obs")),
               row.names = NULL)
  }))
}))
print(assoc)
for (i in seq_len(nrow(assoc)))
  log_msg(sprintf("06f ASSOC %s / %s: median difference unmethylated minus methylated = %.2f SD, Wilcoxon p=%.3f",
                  assoc$cohort[i], assoc$signature[i], assoc$delta_SD[i], assoc$p_wilcox[i]))

## ---- Step 3. Cox models, all fitted on the MGMT-complete subset ---------------
fit_one <- function(k, s) {
  d_all <- G4[G4$cohort == k, ]
  d     <- d_all[!is.na(d_all$MGMT), ]
  zz    <- paste0(s, ".z")
  grab <- function(fit, term) {
    if (is.null(fit)) return(c(NA, NA, NA, NA))
    sm <- summary(fit)$coefficients
    ci <- summary(fit)$conf.int
    c(sm[term, "coef"], ci[term, "lower .95"], ci[term, "upper .95"], sm[term, "Pr(>|z|)"])
  }
  f_full  <- tryCatch(coxph(as.formula(sprintf("Surv(OS.time, OS.event) ~ %s + age", zz)), data = d_all),
                      error = function(e) NULL)
  f_base  <- tryCatch(coxph(as.formula(sprintf("Surv(OS.time, OS.event) ~ %s + age", zz)), data = d),
                      error = function(e) NULL)
  f_mgmt  <- tryCatch(coxph(as.formula(sprintf("Surv(OS.time, OS.event) ~ %s + age + MGMT", zz)), data = d),
                      error = function(e) NULL)
  f_int   <- tryCatch(coxph(as.formula(sprintf("Surv(OS.time, OS.event) ~ %s * MGMT + age", zz)), data = d),
                      error = function(e) NULL)
  a <- grab(f_full, zz); b <- grab(f_base, zz); c3 <- grab(f_mgmt, zz)
  p_int <- if (is.null(f_int)) NA_real_ else {
    cn <- grep(":", rownames(summary(f_int)$coefficients), value = TRUE)
    if (length(cn)) summary(f_int)$coefficients[cn[1], "Pr(>|z|)"] else NA_real_
  }
  ph <- if (is.null(f_mgmt)) NA_real_ else
    tryCatch(cox.zph(f_mgmt)$table["GLOBAL", "p"], error = function(e) NA_real_)
  ## Attenuation only when the restricted base model is genuinely non-null.
  base_ok <- !is.na(b[1]) && exp(b[1]) >= 1.10 && b[4] < 0.10
  att <- if (base_ok) 100 * (b[1] - c3[1]) / b[1] else NA_real_
  data.frame(
    cohort = k, signature = SIG_LAB[s],
    n_full = nrow(d_all), n_restricted = nrow(d), events_restricted = sum(d$OS.event),
    HR_full_sample   = exp(a[1]), p_full_sample   = a[4],
    HR_base_restr    = exp(b[1]), lo_base = b[2], hi_base = b[3], p_base_restr = b[4],
    HR_mgmt_adj      = exp(c3[1]), lo_adj = c3[2], hi_adj = c3[3], p_mgmt_adj = c3[4],
    attenuation_pct  = att, base_nonnull = base_ok,
    p_interaction    = p_int, ph_global_p = ph, row.names = NULL)
}
COX <- do.call(rbind, lapply(levels(G4$cohort), function(k)
  do.call(rbind, lapply(SIG_COLS, function(s) fit_one(k, s)))))
print(COX[, c("cohort", "signature", "n_restricted", "HR_full_sample", "HR_base_restr",
              "HR_mgmt_adj", "p_mgmt_adj", "attenuation_pct")])
for (i in seq_len(nrow(COX))) {
  r <- COX[i, ]
  log_msg(sprintf("06f COX %s / %s: full n=%d HR %.2f | restricted n=%d HR %.2f (p=%.3f) | +MGMT HR %.2f (p=%.3f) | attenuation %s",
                  r$cohort, r$signature, r$n_full, r$HR_full_sample, r$n_restricted,
                  r$HR_base_restr, r$p_base_restr, r$HR_mgmt_adj, r$p_mgmt_adj,
                  if (is.na(r$attenuation_pct)) "not applicable (base null)" else
                    sprintf("%.0f%%", r$attenuation_pct)))
}

## ---- Step 4. Verdict against the pre-specified rule ---------------------------
verdict <- function(r, a) {
  if (!is.na(a$p_wilcox) && a$p_wilcox >= 0.05 && abs(a$delta_SD) < 0.3)
    return("Score not associated with MGMT; confounding by MGMT ruled out")
  if (!r$base_nonnull)
    return("Base model null in restricted subset; attenuation not evaluable")
  if (!is.na(r$attenuation_pct) && r$attenuation_pct >= 30 && r$p_mgmt_adj >= 0.05)
    return("MGMT confounding is the working explanation")
  if (!is.na(r$attenuation_pct) && r$attenuation_pct >= 30)
    return("Partial attenuation; signal persists")
  "MGMT does not explain the signal; discordance unresolved"
}
COX$verdict <- vapply(seq_len(nrow(COX)), function(i) {
  a <- assoc[assoc$cohort == COX$cohort[i] & assoc$signature == COX$signature[i], ][1, ]
  verdict(COX[i, ], a)
}, "")
for (i in seq_len(nrow(COX)))
  log_msg(sprintf("06f VERDICT %s / %s: %s", COX$cohort[i], COX$signature[i], COX$verdict[i]))

## ---- Tables -------------------------------------------------------------------
fmt <- function(hr, lo, hi, p) ifelse(is.na(hr), "NA",
  sprintf("%.2f (%.2f, %.2f), p=%s", hr, lo, hi, ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))))
t11 <- data.frame(
  Cohort = COX$cohort, Signature = COX$signature,
  `n analysed (events)` = sprintf("%d (%d)", COX$n_restricted, COX$events_restricted),
  `Age-adjusted HR (95% CI)` = fmt(COX$HR_base_restr, COX$lo_base, COX$hi_base, COX$p_base_restr),
  `Age- and MGMT-adjusted HR (95% CI)` = fmt(COX$HR_mgmt_adj, COX$lo_adj, COX$hi_adj, COX$p_mgmt_adj),
  `Attenuation` = ifelse(is.na(COX$attenuation_pct), "Not applicable",
                         sprintf("%.0f%%", COX$attenuation_pct)),
  `Interaction p (exploratory)` = sprintf("%.3f", COX$p_interaction), check.names = FALSE)
ft11 <- flextable::flextable(t11) %>% ft_style() %>%
  flextable::add_footer_lines(values = c(
    "Analysis restricted to IDH-wildtype, WHO grade IV gliomas with known MGMT promoter status. Both models were fitted on this same subset so that the change in the hazard ratio reflects adjustment rather than a change of sample; the estimate obtained in the full stratum is given in the text for comparison.",
    "Hazard ratios are expressed per 1 standard deviation of the signature score within each cohort. Attenuation is reported only where the age-adjusted model in the restricted subset was itself non-null, because a ratio computed from a null base estimate is uninformative.",
    "The interaction between signature score and MGMT status is reported for completeness; the available number of events is insufficient to test it reliably and it should not be interpreted as evidence of effect modification.",
    "CI, confidence interval; HR, hazard ratio; MGMT, O6-methylguanine-DNA methyltransferase."))
flextable::save_as_docx(ft11, path = file.path(DIR_TAB, "Table11_MGMT_adjustment.docx"))

t12 <- data.frame(
  Cohort = assoc$cohort, Signature = assoc$signature,
  `Median score, methylated (SD units)`   = sprintf("%.2f", assoc$median_meth),
  `Median score, unmethylated (SD units)` = sprintf("%.2f", assoc$median_unmeth),
  `Difference (SD units)` = sprintf("%.2f", assoc$delta_SD),
  `p-value` = ifelse(assoc$p_wilcox < 0.001, "<0.001", sprintf("%.3f", assoc$p_wilcox)),
  check.names = FALSE)
ft12 <- flextable::flextable(t12) %>% ft_style() %>%
  flextable::add_footer_lines(values = c(
    "Distribution of each signature score by MGMT promoter methylation status within the IDH-wildtype, WHO grade IV stratum. Scores were standardised within cohort, so a difference of 1.00 corresponds to one standard deviation of that cohort.",
    "Groups were compared with the Wilcoxon rank-sum test. A score that does not differ by MGMT status cannot be confounded by it.",
    "MGMT, O6-methylguanine-DNA methyltransferase; SD, standard deviation."))
flextable::save_as_docx(ft12, path = file.path(DIR_TAB, "Table12_score_by_MGMT.docx"))

## ---- Figure: age-adjusted versus MGMT-adjusted, side by side ------------------
pf <- rbind(
  data.frame(cohort = COX$cohort, signature = COX$signature, model = "Age-adjusted",
             HR = COX$HR_base_restr, lo = COX$lo_base, hi = COX$hi_base),
  data.frame(cohort = COX$cohort, signature = COX$signature, model = "Age- and MGMT-adjusted",
             HR = COX$HR_mgmt_adj, lo = COX$lo_adj, hi = COX$hi_adj))
pf$model <- factor(pf$model, levels = c("Age-adjusted", "Age- and MGMT-adjusted"))
g8 <- ggplot(pf, aes(x = HR, y = cohort, shape = model, group = model)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.14,
                 position = position_dodge(width = 0.55)) +
  geom_point(size = 2.6, position = position_dodge(width = 0.55)) +
  facet_wrap(~signature, ncol = 1) + scale_x_log10() +
  scale_shape_manual(values = c("Age-adjusted" = 16, "Age- and MGMT-adjusted" = 17)) +
  labs(x = "HR per 1 SD of signature score (95% CI), MGMT-complete subset",
       y = NULL, shape = NULL) +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom", axis.line.y = element_blank(), axis.ticks.y = element_blank())
save_fig(g8, file.path(DIR_FIG, "Figure8_MGMT_adjustment"), width = 7, height = 5)

## ---- Persist -------------------------------------------------------------------
saveRDS(list(missingness = miss, association = assoc, cox = COX),
        file.path(DIR_RES, "mgmt_adjustment_results.rds"))
write.csv(miss,  file.path(DIR_RES, "table_mgmt_missingness.csv"), row.names = FALSE)
write.csv(assoc, file.path(DIR_RES, "table12_score_by_mgmt.csv"), row.names = FALSE)
write.csv(COX,   file.path(DIR_RES, "table11_mgmt_adjustment.csv"), row.names = FALSE)

nums <- c("", "--- Layer 06f: MGMT adjustment in the IDH-wildtype WHO grade IV stratum ---",
  sprintf("%s: MGMT known in %d of %d (%.0f%% missing); events %d -> %d; detectable HR per SD >= %.2f",
          miss$cohort, miss$n_mgmt_known, miss$n_total, miss$pct_missing,
          miss$events_total, miss$events_mgmt_known, miss$HR_detectable_restricted),
  "",
  sprintf("Score by MGMT, %s / %s: methylated %.2f vs unmethylated %.2f SD units (difference %.2f), Wilcoxon p=%.3f",
          assoc$cohort, assoc$signature, assoc$median_meth, assoc$median_unmeth,
          assoc$delta_SD, assoc$p_wilcox),
  "",
  sprintf("%s / %s: full-stratum HR %.2f | restricted age-adjusted HR %.2f (%.2f-%.2f) p=%.3f | + MGMT HR %.2f (%.2f-%.2f) p=%.3f | attenuation %s | PH global p=%.3f",
          COX$cohort, COX$signature, COX$HR_full_sample, COX$HR_base_restr, COX$lo_base,
          COX$hi_base, COX$p_base_restr, COX$HR_mgmt_adj, COX$lo_adj, COX$hi_adj,
          COX$p_mgmt_adj,
          ifelse(is.na(COX$attenuation_pct), "not applicable (base null)",
                 sprintf("%.0f%%", COX$attenuation_pct)), COX$ph_global_p),
  "",
  sprintf("VERDICT %s / %s: %s", COX$cohort, COX$signature, COX$verdict))
cat(nums, sep = "\n")
write(nums, file = file.path(DIR_OUT, "manuscript_numbers.txt"), append = TRUE)
log_msg("Layer 06f complete: Tables 11-12 and Figure 8 exported to Sep_06/.")
