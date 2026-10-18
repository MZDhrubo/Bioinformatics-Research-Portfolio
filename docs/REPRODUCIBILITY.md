# Reproducibility notes

The command-line workflow records a SHA-256 checksum for the input, alignment files, and tree files, together with analysis settings and basic alignment QC metrics. Logs are written next to the run outputs.

For each research run, keep a private copy of:

- The exact input file and its source/provenance record
- The script version or Git commit used
- Python and external-tool versions
- Full and trimmed alignments
- IQ-TREE reports, trees and console logs
- The generated run manifest and QC summary

The public demo contains synthetic data only. Do not replace it with unpublished data in a public branch.
