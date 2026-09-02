# CXCL13-associated Tpex-like immune niche in HCC

This repository contains the analysis scripts and processed outputs supporting:

> Spatial and clinical features of a CXCL13-associated immune niche containing progenitor-exhausted CD8 T cells in hepatocellular carcinoma

## Study components

- GSE149614: single-cell RNA-seq definition of T/NK states and CXCL13-expressing populations.
- GSE238264: Visium scoring and patient-level spatial summaries.
- TCGA-LIHC: overall-survival analyses.
- GSE140901, GSE279750, and GSE215011: exploratory clinical response-associated analyses.
- GSE202069: FAD-inferred subgroup analysis only; it is not a clinical response cohort.

## Reproduction

1. Clone this repository and create the directory structure described in `data/README.md`.
2. Download the public source files from GEO or GDC.
3. Set `HCC_PROJECT_ROOT` to the repository root.
4. Run scripts in numerical order. Scripts 07-09 generate the single-cell figure panels; scripts 10-16 reproduce the spatial, survival, normalization, FAD-inferred, and external-cohort analyses.

```bash
export HCC_PROJECT_ROOT=/path/to/CXCL13-Tpex-HCC-niche
Rscript scripts/01_import_scrnaseq.R
```

The `processed_data` directory contains numerical values underlying the final analyses. Raw public data are not redistributed.

## Interpretation boundary

The CXCL13-associated exhaustion score is a co-expression score. CXCR3 is included as a T-cell trafficking and effector-state marker; the score does not establish a receptor-mediated mechanism. Spatial spot-level analyses are descriptive. GSE202069 labels are FAD-inferred rather than clinical response labels.

## Citation

See `CITATION.cff`. Repository URL: https://github.com/liukairui621/CXCL13-Tpex-HCC-niche

## License

Code is released under the MIT License. Public source datasets remain subject to their original repository terms.
