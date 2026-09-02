#!/usr/bin/env python3
from pathlib import Path
import os
PROJECT_ROOT = Path(os.getenv("HCC_PROJECT_ROOT", ".")).resolve()
for subdir in ("raw", "processed", "tables", "stats", "panels", "logs"):
    (PROJECT_ROOT / "results" / "external_validation" / subdir).mkdir(parents=True, exist_ok=True)
"""Prepare expression matrices for GSE215011 and GSE279750."""

import gzip
import os
import sys

# --- GSE215011 ---
print("=== Processing GSE215011 ===")
infile = PROJECT_ROOT / "results/external_validation/raw/GSE215011_gene_description_human_samples.txt.gz"
outfile = PROJECT_ROOT / "results/external_validation/processed/GSE215011_expression_matrix.tsv"

# Sample metadata from SOFT
sample_meta = {
    "A18_440T": {"gsm": "GSM6619959", "response": "Responder"},
    "A19_174T": {"gsm": "GSM6619960", "response": "Responder"},
    "A19_171T": {"gsm": "GSM6619961", "response": "Responder"},
    "A19_121T": {"gsm": "GSM6619962", "response": "Responder"},
    "A19_1T":   {"gsm": "GSM6619963", "response": "Responder"},
    "A19_16T":  {"gsm": "GSM6619964", "response": "Nonresponder"},
    "A11_385T": {"gsm": "GSM6619965", "response": "Nonresponder"},
    "A17_152T": {"gsm": "GSM6619966", "response": "Nonresponder"},
    "A12_949T": {"gsm": "GSM6619967", "response": "Nonresponder"},
    "A16_557T": {"gsm": "GSM6619968", "response": "Nonresponder"},
}

# Read expression file
genes = {}  # gene_name -> {sample: fpkm}
sample_cols = []
with gzip.open(infile, "rt") as f:
    header = f.readline().strip().split("\t")
    # Columns: gene_id, then sample columns, then gene_name, gene_chr, etc.
    # Find sample columns (between gene_id and gene_name)
    gene_name_idx = header.index("gene_name")
    sample_cols = header[1:gene_name_idx]  # sample names
    print(f"Sample columns: {sample_cols}")

    for line in f:
        fields = line.strip().split("\t")
        gene_name = fields[gene_name_idx]
        if not gene_name or gene_name == "-":
            continue
        vals = {}
        for i, sc in enumerate(sample_cols):
            try:
                vals[sc] = float(fields[1 + i])
            except (ValueError, IndexError):
                vals[sc] = 0.0

        if gene_name not in genes:
            genes[gene_name] = vals
        else:
            # Duplicate gene name: keep max expression
            for sc in sample_cols:
                genes[gene_name][sc] = max(genes[gene_name].get(sc, 0), vals.get(sc, 0))

print(f"Unique genes: {len(genes)}, Samples: {len(sample_cols)}")

# Write expression matrix
with open(outfile, "w") as f:
    f.write("gene\t" + "\t".join(sample_cols) + "\n")
    for gene in sorted(genes.keys()):
        vals = [str(genes[gene].get(sc, 0)) for sc in sample_cols]
        f.write(gene + "\t" + "\t".join(vals) + "\n")
print(f"Written: {outfile}")

# Write metadata
meta_out = PROJECT_ROOT / "results/external_validation/tables/GSE215011_sample_metadata.csv"
with open(meta_out, "w") as f:
    f.write("sample_id,gsm_accession,response,tissue,disease,treatment,label_source\n")
    for sid in sample_cols:
        m = sample_meta.get(sid, {})
        f.write(f"{sid},{m.get('gsm','')},{m.get('response','')},tumor,HCC,Nivolumab,GEO_metadata_group_field\n")
print(f"Written: {meta_out}")


# --- GSE279750 ---
print("\n=== Processing GSE279750 ===")
import openpyxl

raw_dir = PROJECT_ROOT / "results/external_validation/raw/GSE279750"
xlsx_files = sorted([f for f in os.listdir(raw_dir) if f.endswith(".xlsx")])

# Sample metadata from SOFT
sample_meta_279 = {
    "GSM8579621": {"name": "NR1", "response": "Nonresponder"},
    "GSM8579622": {"name": "NR2", "response": "Nonresponder"},
    "GSM8579623": {"name": "NR3", "response": "Nonresponder"},
    "GSM8579624": {"name": "NR4", "response": "Nonresponder"},
    "GSM8579625": {"name": "R1", "response": "Responder"},
    "GSM8579626": {"name": "R2", "response": "Responder"},
    "GSM8579627": {"name": "R3", "response": "Responder"},
    "GSM8579628": {"name": "R4", "response": "Responder"},
    "GSM8579629": {"name": "R5", "response": "Responder"},
    "GSM8579630": {"name": "R6", "response": "Responder"},
}

all_data = {}  # sample_name -> {gene: fpkm}
sample_names_279 = []

for xf in xlsx_files:
    # Extract GSM and sample name from filename
    gsm = xf.split("_")[0]
    sname = sample_meta_279.get(gsm, {}).get("name", gsm)
    sample_names_279.append(sname)

    wb = openpyxl.load_workbook(os.path.join(raw_dir, xf), read_only=True)
    ws = wb.active
    gene_data = {}
    for row in ws.iter_rows(min_row=2, values_only=True):
        gene = row[0]
        fpkm = row[1] if row[1] is not None else 0
        if gene and isinstance(fpkm, (int, float)):
            if gene not in gene_data:
                gene_data[gene] = fpkm
            else:
                gene_data[gene] = max(gene_data[gene], fpkm)
    all_data[sname] = gene_data
    wb.close()
    print(f"  {xf}: {sname}, {len(gene_data)} genes")

# Union of all genes
all_genes_279 = sorted(set().union(*[set(d.keys()) for d in all_data.values()]))
print(f"Total unique genes: {len(all_genes_279)}, Samples: {len(sample_names_279)}")

# Write expression matrix
outfile_279 = PROJECT_ROOT / "results/external_validation/processed/GSE279750_expression_matrix.tsv"
with open(outfile_279, "w") as f:
    f.write("gene\t" + "\t".join(sample_names_279) + "\n")
    for gene in all_genes_279:
        vals = [str(all_data[s].get(gene, 0)) for s in sample_names_279]
        f.write(gene + "\t" + "\t".join(vals) + "\n")
print(f"Written: {outfile_279}")

# Write metadata
meta_out_279 = PROJECT_ROOT / "results/external_validation/tables/GSE279750_sample_metadata.csv"
with open(meta_out_279, "w") as f:
    f.write("sample_id,gsm_accession,response,tissue,disease,treatment,label_source,timepoint_status\n")
    for xf in xlsx_files:
        gsm = xf.split("_")[0]
        m = sample_meta_279.get(gsm, {})
        sname = m.get("name", gsm)
        resp = m.get("response", "")
        # Timepoint: metadata says "HCC samples from 10 patients receiving anti-PD-L1 therapy"
        # No explicit "pre-treatment" or "post-treatment" stated
        f.write(f"{sname},{gsm},{resp},HCC_tumor,HCC,anti-PD-L1,GEO_metadata_patient_field,not_explicitly_stated\n")
print(f"Written: {meta_out_279}")

print("\nAll expression matrices prepared.")
