# Public data layout

Download source files from the indicated public repositories and organize them as follows:

```
data/
  GSE149614/
    GSE149614_HCC.scRNAseq.S71915.count.txt.gz
    GSE149614_HCC.metadata.updated.txt.gz
  GSE238264/
    HCC1R/ ... HCC7NR/       # Space Ranger directories
  TCGA_LIHC/
    TCGA_LIHC_expression.gz
    TCGA_LIHC_survival.txt
    LIHC_clinicalMatrix.txt
  GSE140901/
    GSE140901_series_matrix.txt.gz
    rcc_files/
  GSE202069/
    GSE202069_tpm.txt.gz
  external_validation/
    raw/
      GSE215011_gene_description_human_samples.txt.gz
      GSE279750/
```

GSE279750 and GSE215011 response labels are read from deposited metadata. GSE202069 is used only for FAD-inferred subgrouping.
