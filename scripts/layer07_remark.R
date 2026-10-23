## =====================================================================
## Layer 07 - REMARK items 10, 11 and 13
##
## Three reporting requirements remain open in the REMARK checklist. None
## of them is an analysis in the ordinary sense; each is a quantity that a
## prognostic marker study is expected to report and whose absence is a
## routine reason for a methodological reviewer to ask for revision.
##
##   Item 10. The effect size that the study was able to detect, given the
##            number of events observed. A null result is uninterpretable
##            unless the reader is told what it excludes. This layer
##            expresses the null findings as the range of hazard ratios
##            that each cohort could and could not have detected.
##
##   Item 11. A flow diagram, per cohort, from the assembled series to the
##            analysed stratum, with the reason and the count at every
##            exclusion. The figure is generated from the same objects the
##            analysis used, so the numbers cannot drift from the Results.
##
##   Item 13. Median follow-up among survivors, estimated by the reverse
##            Kaplan-Meier method, in which the censoring indicator is
##            reversed so that the estimator describes observation time
##            rather than survival time. The naive median of observed
##            follow-up, which is what most manuscripts report, is
##            downward-biased because patients who die stop contributing.
##
## Flow counts are reconstructed from the source files using the same
## filters as the pipeline, so the diagram reports the study as executed
## rather than as remembered.
##
## Usage:
##   setwd("~/Library/CloudStorage/Dropbox/HOPE-GBM"); source("layer07_remark.R")
## Requires Layers 06e and 06g.
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

disc <- readRDS(file.path(DIR_RES, "discordance_results.rds"))
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
pick_col <- function(nm, exact, pattern) {
  hit <- exact[exact %in% nm]; if (length(hit)) return(hit[1])
  grep(pattern, nm, value = TRUE, ignore.case = TRUE)[1]
}
read_cgga_ids <- function(path) {
  con <- file(path, open = "rb")
  ln  <- readLines(con, n = 1, skipNul = TRUE, warn = FALSE, encoding = "latin1"); close(con)
  f <- strsplit(gsub("\r$", "", ln, useBytes = TRUE), "\t", fixed = TRUE)[[1]]
  f[-1]
}
num_of <- function(p) regmatches(basename(p), regexpr("325|693", basename(p)))

## =====================================================================
## REMARK item 11: flow per cohort, reconstructed from source
## =====================================================================
expr_f <- list.files(CGGA_DIR, pattern = "mRNAseq.*genes.*[.]txt$", full.names = TRUE)
clin_f <- list.files(CGGA_DIR, pattern = "clinical.*[.]txt$", full.names = TRUE)
clin_f <- clin_f[match(vapply(expr_f, num_of, ""), vapply(clin_f, num_of, ""))]

flow_one <- function(ep, cp) {
  label <- paste0("CGGA-", num_of(ep))
  ids_expr <- read_cgga_ids(ep)
  clin <- as.data.frame(readr::read_tsv(cp, show_col_types = FALSE, progress = FALSE))
  nm <- colnames(clin)
  cid   <- pick_col(nm, "CGGA_ID", "CGGA|sample|ID")
  cprs  <- pick_col(nm, "PRS_type", "PRS")
  cidh  <- pick_col(nm, "IDH_mutation_status", "IDH")
  cgr   <- pick_col(nm, "Grade", "Grade|WHO")
  cage  <- pick_col(nm, "Age", "Age")
  cos   <- pick_col(nm, "OS", "^OS$|overall")
  cev   <- grep("Censor", nm, value = TRUE, ignore.case = TRUE)[1]

  steps <- list(); n <- function(x) length(unique(x))
  s <- clin[[cid]]
  steps[["Patients in the clinical file"]] <- n(s)
  s1 <- intersect(s, ids_expr)
  steps[["With RNA sequencing data available"]] <- n(s1)
  k <- clin[clin[[cid]] %in% s1, ]
  prs_ok <- if (!is.na(cprs)) !is.na(k[[cprs]]) & k[[cprs]] == "Primary" else rep(TRUE, nrow(k))
  steps[["Primary tumour (recurrent and secondary excluded)"]] <- sum(prs_ok)
  k <- k[prs_ok, ]
  surv_ok <- !is.na(k[[cos]]) & suppressWarnings(as.numeric(k[[cos]])) > 0 & !is.na(k[[cev]])
  steps[["Overall survival recorded and greater than zero"]] <- sum(surv_ok)
  k <- k[surv_ok, ]
  age_ok <- !is.na(suppressWarnings(as.numeric(k[[cage]])))
  steps[["Age recorded"]] <- sum(age_ok)
  k <- k[age_ok, ]
  idh_known <- !is.na(k[[cidh]]) & grepl("mut|wild|wt", tolower(k[[cidh]]))
  steps[["IDH status known"]] <- sum(idh_known)
  k <- k[idh_known, ]
  idh_wt <- grepl("wildtype|^wt", tolower(k[[cidh]]))
  steps[["IDH-wildtype"]] <- sum(idh_wt)
  k <- k[idh_wt, ]
  g4 <- grepl("IV|4", gsub("\\s", "", as.character(k[[cgr]])))
  steps[["WHO grade IV"]] <- sum(g4)
  k <- k[g4, ]

  df <- data.frame(cohort = label, step = names(steps),
                   n = unlist(steps, use.names = FALSE), row.names = NULL)
  df$excluded <- c(NA, -diff(df$n))
  df$order <- seq_len(nrow(df))
  ## reconcile with the analysed object
  analysed <- sum(G4$cohort == label)
  df$analysed_in_pipeline <- analysed
  if (tail(df$n, 1) != analysed)
    log_msg(sprintf("07 WARNING %s: flow reconstruction ends at n=%d but the analysed object holds n=%d; reconcile before submission",
                    label, tail(df$n, 1), analysed))
  df
}
FLOW <- do.call(rbind, Map(flow_one, expr_f, clin_f))
print(FLOW[, c("cohort", "step", "n", "excluded")])
for (k in unique(FLOW$cohort)) {
  f <- FLOW[FLOW$cohort == k, ]
  log_msg(sprintf("07 FLOW %s: %d -> %d analysed", k, f$n[1], tail(f$n, 1)))
}

## Flow diagram, drawn as stacked boxes with the exclusion reason alongside
flow_plot <- function(k) {
  f <- FLOW[FLOW$cohort == k, ]
  f$y <- rev(seq_len(nrow(f)))
  f$lab <- sprintf("%s\nn = %d", f$step, f$n)
  ex <- f[!is.na(f$excluded) & f$excluded > 0, ]
  ggplot(f, aes(x = 0, y = y)) +
    geom_tile(width = 1.9, height = 0.66, fill = "grey95", colour = "grey40", linewidth = 0.3) +
    geom_text(aes(label = lab), size = 2.9, lineheight = 0.95) +
    { if (nrow(ex)) geom_text(data = ex, aes(x = 1.15, y = y + 0.5,
        label = sprintf("excluded n = %d", excluded)), size = 2.5, hjust = 0, colour = "grey30") } +
    geom_segment(data = f[-nrow(f), ], aes(x = 0, xend = 0, y = y - 0.33, yend = y - 0.67),
                 arrow = arrow(length = unit(0.06, "in")), linewidth = 0.3, colour = "grey40") +
    scale_x_continuous(limits = c(-1.1, 2.4)) +
    labs(title = k) + theme_void(base_size = 10) +
    theme(plot.title = element_text(size = 10, hjust = 0, face = "plain"))
}
for (k in unique(FLOW$cohort))
  save_fig(flow_plot(k), file.path(DIR_FIG, paste0("Figure_REMARK_flow_", gsub("-", "", k))),
           width = 4.6, height = 5.4)

## =====================================================================
## REMARK item 13: median follow-up among survivors (reverse Kaplan-Meier)
## =====================================================================
reverse_km <- function(time, event) {
  ## Censoring indicator reversed: the event of interest becomes "still
  ## under observation", so the estimator describes observation time.
  fit <- survfit(Surv(time, 1 - event) ~ 1)
  med <- unname(summary(fit)$table["median"])
  q <- tryCatch(quantile(fit, probs = c(0.25, 0.75))$quantile, error = function(e) c(NA, NA))
  c(median = med, q25 = unname(q[1]), q75 = unname(q[2]))
}
FU <- do.call(rbind, lapply(levels(G4$cohort), function(k) {
  d <- G4[G4$cohort == k, ]
  rk <- reverse_km(d$OS.time, d$OS.event)
  naive <- median(d$OS.time, na.rm = TRUE)
  alive <- d$OS.time[d$OS.event == 0]
  data.frame(cohort = k, n = nrow(d), deaths = sum(d$OS.event),
             censored = sum(d$OS.event == 0),
             median_reverse_km_months = rk["median"] / 30.4375,
             q25_months = rk["q25"] / 30.4375, q75_months = rk["q75"] / 30.4375,
             naive_median_months = naive / 30.4375,
             median_among_censored_months = median(alive, na.rm = TRUE) / 30.4375,
             max_followup_months = max(d$OS.time, na.rm = TRUE) / 30.4375,
             row.names = NULL)
}))
print(FU)
for (i in seq_len(nrow(FU))) {
  r <- FU[i, ]
  log_msg(sprintf("07 FOLLOW-UP %s: reverse KM median %.1f months (IQR %.1f-%.1f); naive median of observed time %.1f; %d censored of %d",
                  r$cohort, r$median_reverse_km_months, r$q25_months, r$q75_months,
                  r$naive_median_months, r$censored, r$n))
}

## =====================================================================
## REMARK item 10: detectable effect size
## =====================================================================
## Schoenfeld's relation for a continuous covariate expressed per 1 SD.
## Reported in two directions: the hazard ratio detectable at 80% and 90%
## power, and the power the study had against the estimate reported by the
## cohort that did find a signal.
detectable_hr <- function(events, power, alpha = 0.05)
  exp((qnorm(1 - alpha / 2) + qnorm(power)) / sqrt(events))
power_for_hr <- function(events, hr, alpha = 0.05)
  pnorm(sqrt(events) * abs(log(hr)) - qnorm(1 - alpha / 2))

REF_HR <- 1.38   # Chen 2017 in CGGA-325, the only positive estimate in the stratum
POW <- do.call(rbind, lapply(levels(G4$cohort), function(k) {
  d <- G4[G4$cohort == k, ]; e <- sum(d$OS.event)
  data.frame(cohort = k, n = nrow(d), events = e,
             hr_80 = detectable_hr(e, 0.80), hr_90 = detectable_hr(e, 0.90),
             power_vs_ref = power_for_hr(e, REF_HR), row.names = NULL)
}))
## TCGA, from the primary analysis table when available
tcga <- tryCatch({
  t2 <- read.csv(file.path(DIR_RES, "table2_figure1_values.csv"))
  i <- grep("TCGA", t2$cohort)[1]
  data.frame(cohort = "TCGA-GBM", n = t2$n[i], events = t2$deaths[i],
             hr_80 = detectable_hr(t2$deaths[i], 0.80),
             hr_90 = detectable_hr(t2$deaths[i], 0.90),
             power_vs_ref = power_for_hr(t2$deaths[i], REF_HR))
}, error = function(e) NULL)
if (!is.null(tcga)) POW <- rbind(POW, tcga)
print(POW)
for (i in seq_len(nrow(POW))) {
  r <- POW[i, ]
  log_msg(sprintf("07 POWER %s: %d events; detectable HR >= %.2f at 80%% and >= %.2f at 90%% power; power against HR %.2f was %.0f%%",
                  r$cohort, r$events, r$hr_80, r$hr_90, REF_HR, 100 * r$power_vs_ref))
}

## ---- Tables -------------------------------------------------------------------
tf <- FLOW %>% select(Cohort = cohort, Step = step, `n remaining` = n,
                      `n excluded` = excluded)
tf$`n excluded` <- ifelse(is.na(tf$`n excluded`), "", as.character(tf$`n excluded`))
ftF <- flextable::flextable(tf) %>% ft_style() %>%
  flextable::merge_v(j = "Cohort") %>%
  flextable::add_footer_lines(values = c(
    "Sequence of exclusions applied to each cohort, reconstructed from the source clinical and expression files with the filters used in the analysis, so that the counts reported here and those in the Results derive from the same operations.",
    "Recurrent and secondary tumours were excluded because the WHO 2021 definition of glioblastoma applies to primary IDH-wildtype disease and because their inclusion in an earlier version of this work was an error.",
    "WHO, World Health Organization."))
flextable::save_as_docx(ftF, path = file.path(DIR_TAB, "Table_REMARK11_flow.docx"))

tfu <- data.frame(
  Cohort = FU$cohort, `Patients (deaths)` = sprintf("%d (%d)", FU$n, FU$deaths),
  `Median follow-up, reverse Kaplan-Meier, months (IQR)` =
    sprintf("%.1f (%.1f, %.1f)", FU$median_reverse_km_months, FU$q25_months, FU$q75_months),
  `Median observed time, months` = sprintf("%.1f", FU$naive_median_months),
  `Maximum follow-up, months` = sprintf("%.1f", FU$max_followup_months),
  check.names = FALSE)
ftFU <- flextable::flextable(tfu) %>% ft_style() %>%
  flextable::add_footer_lines(values = c(
    "Median follow-up was estimated by the reverse Kaplan-Meier method, in which the censoring indicator is reversed so that the estimator describes the duration of observation rather than survival. The median of the observed times is given for comparison and is systematically shorter, because patients who die cease to contribute follow-up.",
    "IQR, interquartile range."))
flextable::save_as_docx(ftFU, path = file.path(DIR_TAB, "Table_REMARK13_followup.docx"))

tp <- data.frame(
  Cohort = POW$cohort, `Patients (events)` = sprintf("%d (%d)", POW$n, POW$events),
  `Detectable HR, 80% power` = sprintf("%.2f", POW$hr_80),
  `Detectable HR, 90% power` = sprintf("%.2f", POW$hr_90),
  `Power against HR 1.38` = sprintf("%.0f%%", 100 * POW$power_vs_ref),
  check.names = FALSE)
ftP <- flextable::flextable(tp) %>% ft_style() %>%
  flextable::add_footer_lines(values = c(
    "Detectable effect sizes were obtained from Schoenfeld's relation between the number of events and the log hazard ratio, for a continuous covariate expressed per 1 standard deviation, at a two-sided alpha of 0.05. They describe what each cohort was able to exclude rather than what it observed.",
    "The final column gives the power each cohort had against a hazard ratio of 1.38, the estimate obtained for the Chen 2017 signature in CGGA-325 and the only positive result in this stratum.",
    "HR, hazard ratio."))
flextable::save_as_docx(ftP, path = file.path(DIR_TAB, "Table_REMARK10_power.docx"))

## ---- Persist -------------------------------------------------------------------
saveRDS(list(flow = FLOW, followup = FU, power = POW, ref_hr = REF_HR),
        file.path(DIR_RES, "remark_items.rds"))
write.csv(FLOW, file.path(DIR_RES, "remark11_flow.csv"), row.names = FALSE)
write.csv(FU,   file.path(DIR_RES, "remark13_followup.csv"), row.names = FALSE)
write.csv(POW,  file.path(DIR_RES, "remark10_power.csv"), row.names = FALSE)

nums <- c("", "--- Layer 07: REMARK items 10, 11 and 13 ---",
  sprintf("FLOW %s | %s: n=%d%s", FLOW$cohort, FLOW$step, FLOW$n,
          ifelse(is.na(FLOW$excluded), "", sprintf(" (excluded %d)", FLOW$excluded))),
  "",
  sprintf("FOLLOW-UP %s: reverse Kaplan-Meier median %.1f months (IQR %.1f-%.1f); median observed time %.1f months; maximum %.1f months; %d of %d censored",
          FU$cohort, FU$median_reverse_km_months, FU$q25_months, FU$q75_months,
          FU$naive_median_months, FU$max_followup_months, FU$censored, FU$n),
  "",
  sprintf("POWER %s: %d events; detectable HR >= %.2f (80%%) and >= %.2f (90%%); power against HR %.2f = %.0f%%",
          POW$cohort, POW$events, POW$hr_80, POW$hr_90, REF_HR, 100 * POW$power_vs_ref))
cat(nums, sep = "\n")
write(nums, file = file.path(DIR_OUT, "manuscript_numbers.txt"), append = TRUE)
log_msg("Layer 07 complete: REMARK tables and flow diagrams exported to Sep_06/.")
