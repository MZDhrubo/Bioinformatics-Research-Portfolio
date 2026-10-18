# Research data publication checklist

This public repository is deliberately separated from unpublished research inputs and results.

## Do not commit unless cleared for public release

- Original FASTA/FASTQ or other sequence datasets supplied for the unpublished study
- Alignments, trees, candidate placement tables, audit workbooks, final results, and manuscript-specific interpretation
- Private reference lists or accession mappings where release is not permitted
- Logs that expose private file paths, usernames, tokens, server locations, or unpublished results
- Credentials, API keys, passwords, private keys, or `.env` files
- Data downloaded under a license or agreement that does not permit redistribution

## Before every push

1. Run `git status --short` and inspect every staged filename.
2. Review staged text files for sequence strings, accession identifiers, result values, local paths, and credentials.
3. Inspect notebook cells and outputs, not just notebook filenames.
4. Confirm permission to redistribute every non-synthetic dataset and third-party file.
5. Keep private inputs under `data/private/` and results under `results/`; both are ignored by Git.
6. If something sensitive was committed, deleting it in a later commit does not remove it from history. Revoke exposed credentials and follow the repository's history-cleaning procedure.

This checklist reduces accidental disclosure but is not a substitute for institutional, funder, journal, or data-provider requirements.
