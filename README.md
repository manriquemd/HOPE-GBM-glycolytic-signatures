# Glycolytic gene signatures, IDH status and survival in glioma

Reproducible analysis pipeline for the manuscript *Glycolytic gene signatures, IDH status and survival in glioma: a three-cohort prognostic study* (Manrique-Guzmán S, Gutiérrez-Salmeán G, Albores-Méndez EM, Castro-Marín M, Uribe-Figueroa LI).

## What the pipeline does

The analysis (pipeline version 3.2, 2026-09-20) runs in ordered layers with a fixed random seed (1024):

| Step | File | Content |
|---|---|---|
| 00 to 06, 07 | `HOPE-GBM_pipeline_v3.2.Rmd` | Environment and seed; TCGA-GBM acquisition (GDC, STAR-Counts) and CGGA mRNAseq_693 / mRNAseq_325 import; selection of primary IDH-wildtype glioblastoma; DESeq2 VST and log2(RSEM + 1); HALLMARK_GLYCOLYSIS GSVA score and five-gene panel; Cox models adjusted for age and MGMT, proportional hazards diagnostics and sensitivity analyses; Verhaak subtypes; Chen 2017 signature with and without subtype adjustment; export of tables and figures |
| 06b | `scripts/layer06b_idh.R` | Chen signature and IDH status in TCGA |
| 06c | `scripts/layer06c_cgga_idh.R` | Nested adjustment for IDH status and WHO grade in the full CGGA glioma cohorts |
| 06d | `scripts/layer06d_generality.R` | Zhang 2021 signature; IDH-mutant stratum and grade |
| 06e | `scripts/layer06e_discordance.R` | Composition of the two CGGA IDH-wildtype grade IV strata |
| 06f | `scripts/layer06f_mgmt.R` | MGMT adjustment and score by methylation status |
| 06g | `scripts/layer06g_random_signatures.R` | Expression-matched random-signature null and sign-flip null |
| 06h | `scripts/layer06h_null_calibration.R` | Calibration of the null (dispersion, independent draws) |
| 06i | `scripts/layer06i_structure.R` | First principal component of expression and its removal |
| 07 | `scripts/layer07_remark.R` | REMARK items: patient flow, follow-up, power |

Run the R Markdown first; each script in `scripts/` is then sourced from the project folder in the order shown (`setwd(<project>); source("scripts/<file>")`), reading the objects saved by the previous steps. Every script fixes the seed at 1024.

The random-signature null draws 1,000 ten-gene sets per cohort and signature, matching each signature gene by mean and standard deviation of expression, excluding the signature's own genes, and fits the same Cox model to each set. The sign-flip null keeps the published genes and permutes the signs of their coefficients, separating the contribution of gene identity from that of weight direction. Both can be applied to any published transcriptomic prognostic signature.

## Data

No patient-level data are distributed here. All inputs are public:

* TCGA-GBM: NCI Genomic Data Commons (https://portal.gdc.cancer.gov), retrieved with TCGAbiolinks.
* CGGA mRNAseq_693 and mRNAseq_325 expression and clinical files: http://www.cgga.org.cn (download manually).
* MSigDB HALLMARK gene sets: msigdbr.

## Running

```r
Sys.setenv(HOPE_GBM_DIR = "/path/to/project",   # will hold data_raw/, data_processed/, results/, logs/
           CGGA_DIR     = "/path/to/CGGA")      # CGGA mRNAseq_693 / mRNAseq_325 files
rmarkdown::render("HOPE-GBM_pipeline_v3.2.Rmd")
```

Outputs are written to `Sep_06/Tables` and `Sep_06/Figures` under the project folder, with every number cited in the manuscript in `Sep_06/manuscript_numbers.txt`. Session information is printed at the end of the rendered report.

Software used for the published results: R 4.5.2, TCGAbiolinks 2.38.0, DESeq2 1.50.2, GSVA 2.4.8, msigdbr 26.1.0 (MSigDB 2026.1.Hs), survival 3.8, survminer 0.5.2, gtsummary 2.5.0.

## Citation

Please cite the article once published. A versioned release of this repository will be archived on Zenodo with a DOI.

## Contact

Salvador Manrique-Guzmán, MD, PhD (ORCID 0000-0001-5792-6363).

## License

MIT (see `LICENSE`).
