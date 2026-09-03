# Submission figures

This directory contains the versioned PLOS ONE submission layouts for Fig 1, Fig 5, Fig 6, S5 Fig, and S6 Fig. Rebuild them with:

```bash
python scripts/17_build_submission_layouts.py
python scripts/18_build_fig5_s5_layouts.py
```

The scripts also write local editable PDF exports. PDFs are ignored by Git because their embedded creation timestamps are not byte-deterministic. The versioned TIFF files are deterministic. `layout_build_manifest.json` and `fig5_s5_build_manifest.json` record dimensions, compression, checksums, software versions, and displayed statistics.
