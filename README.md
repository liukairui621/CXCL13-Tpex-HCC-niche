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
4. Run scripts in numerical order. Scripts 07-09 generate the single-cell analytical panels; scripts 10-16 reproduce the spatial, survival, normalization, FAD-inferred, and external-cohort analyses.
5. Run script 17 after the analytical outputs are available. It reproduces the exact submission geometry and 600-dpi TIFF export for Fig 1, Fig 6, and S6 Fig.

```bash
export HCC_PROJECT_ROOT=/path/to/CXCL13-Tpex-HCC-niche
Rscript scripts/01_import_scrnaseq.R
python scripts/17_build_submission_layouts.py
```

The `processed_data` directory contains numerical values underlying the final analyses. Raw public data are not redistributed.

Final submission layouts are versioned in `submission_figures`. Fig 6 and S6 Fig are rebuilt directly from the audited numerical tables. Fig 1 is assembled from lossless, version-controlled panel exports in `figure_assets`; their scientific source panels remain reproducible through script 07. This separation prevents plotting-library or manual-composition differences from being mistaken for analytical differences.

For the same typography as the archived submission files, install Arial or set `HCC_FIGURE_FONT` to an available metrically compatible font before running script 17. The build manifest records the actual software versions, font, statistics, TIFF properties, and SHA-256 checksums.

For a byte-stable TIFF rebuild, use Python 3.12.10 and install the pinned figure-rendering environment with `pip install -r requirements-figures-lock.txt`. PDF metadata include creation timestamps and are therefore not expected to have stable file hashes.

Patient-level spatial and small clinical-cohort comparisons use two-sided exact Wilcoxon rank-sum tests. Spot-level spatial tests are descriptive and use the asymptotic approximation because spots are not independent patient-level replicates.

## Interpretation boundary

The CXCL13-associated exhaustion score is a co-expression score. CXCR3 is included as a T-cell trafficking and effector-state marker; the score does not establish a receptor-mediated mechanism. Spatial spot-level analyses are descriptive. GSE202069 labels are FAD-inferred rather than clinical response labels.

## Citation

See `CITATION.cff`. Repository URL: https://github.com/liukairui621/CXCL13-Tpex-HCC-niche

## License

Code is released under the MIT License. Public source datasets remain subject to their original repository terms.
