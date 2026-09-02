# Submission figures

This directory contains the versioned PLOS ONE submission layouts for Fig 1, Fig 6, and S6 Fig. Rebuild them with:

```bash
python scripts/17_build_submission_layouts.py
```

The script also writes local editable PDF exports for Fig 6 and S6 Fig. PDFs are ignored by Git because their embedded creation timestamps are not byte-deterministic. The versioned TIFF files are deterministic, and `layout_build_manifest.json` records their dimensions, compression, checksums, software versions, and displayed statistics.
