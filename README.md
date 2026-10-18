# Protein Sequence Phylogeny Workflow Portfolio

A public-safe coding portfolio demonstrating a reproducible protein-sequence analysis workflow. It includes reusable Python code, input validation, alignment quality checks, command-line integration with MAFFT, ClipKIT and IQ-TREE, automated tests, and a synthetic demonstration dataset.

> **Research-data boundary:** This public portfolio intentionally excludes original research FASTA files, unpublished alignments and trees, project-specific audit spreadsheets, result tables, private logs, and manuscript interpretations. The example FASTA contains synthetic sequences only. Do not add unpublished or third-party data without explicit permission.

## What this demonstrates

- Python scripting and command-line interfaces
- FASTA input validation and duplicate-ID checks
- Reproducible external-tool orchestration with logged commands
- Alignment quality-control summaries
- SHA-256 checksums and run manifests
- Automated unit tests and GitHub Actions
- Research integrity: separation of public code from private data and results

## Workflow overview

1. Validate a protein FASTA file.
2. Align sequences with MAFFT.
3. Create a smart-gap-trimmed alignment with ClipKIT.
4. Infer maximum-likelihood trees for the full and trimmed alignments using IQ-TREE.
5. Save logs, alignment QC summaries, tree files, reports, and a manifest to a local results directory.

The pipeline is a reusable framework. Model choice, trimming policy, bootstrap settings, and biological interpretation must be justified for the actual study; a pipeline run does not by itself establish orthology or any biological conclusion.

## Requirements

- Python 3.10 or later
- MAFFT
- ClipKIT
- IQ-TREE (the executable may be named `iqtree3`, `iqtree2`, or `iqtree`)

Install Python packages:

```bash
python -m venv .venv
# Windows PowerShell:
.venv\Scripts\Activate.ps1
# macOS/Linux:
# source .venv/bin/activate
python -m pip install -r requirements.txt
```

Install MAFFT, ClipKIT and IQ-TREE from their official project sources and make their executables available on `PATH`. This portfolio does not automatically download or execute third-party installers.

## Run the synthetic example

```bash
python scripts/run_pipeline.py --input examples/synthetic_proteins.fasta --outdir results/demo --model JTT --bootstrap 100 --threads 2
```

The synthetic example is for software demonstration only, not biological inference. A successful run requires the external tools above. If the synthetic dataset is too small or unsuitable for a particular IQ-TREE release, use it only for input/QC testing and follow the tool's diagnostic output.

## Use with private research data

Keep original project inputs in `data/private/` on your own computer. That directory is ignored by Git. For example:

```bash
python scripts/run_pipeline.py --input data/private/your_protein_sequences.fasta --outdir results/private_run --model JTT --bootstrap 1000 --threads 4
```

**Before pushing to GitHub:** inspect `git status`, run the privacy checklist in `docs/RESEARCH_DATA_POLICY.md`, and confirm that no raw data, unpublished results, accessions, private logs, or credentials are staged. The pipeline's `results/` directory is also ignored by Git.

## Run tests

```bash
python -m pip install -r requirements-dev.txt
python -m pytest -q
```

Tests validate core parsing and QC functions. They do not claim that MAFFT, ClipKIT, or IQ-TREE have been run successfully in every environment; those tools must be installed and their run logs checked separately.

## Reproducibility and limitations

- Input and output files are hashed with SHA-256 in the local run manifest.
- Tool versions should be recorded with each real analysis and reported in the manuscript.
- The fixed substitution model and trimming strategy are configurable choices, not universally optimal defaults.
- Bootstrap support is not the probability that a biological hypothesis is true.
- Tree proximity alone does not prove orthology; sequence completeness, annotation, genomic context and synteny may be needed.
- The public repository is not the authoritative archive for the unpublished study. Keep the original analysis and its data in a controlled private location.

## Responsible portfolio use

Review `CONTRIBUTIONS.md` before presenting this repository to an admissions committee. Describe your own work, code review, testing, and any AI-assisted development accurately. Be prepared to explain each function and each analytical choice yourself.
