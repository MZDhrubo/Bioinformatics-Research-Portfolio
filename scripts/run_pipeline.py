#!/usr/bin/env python3
"""Run a reproducible protein alignment and phylogeny workflow.

The script expects MAFFT, ClipKIT, and IQ-TREE to be installed separately.
It does not download tools or upload data. Keep unpublished inputs/results local.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

from Bio import SeqIO

VALID_AA = set("ACDEFGHIKLMNPQRSTVWYBXZJUO*-?.")
MISSING_AA = set("-?.X")


def sha256_file(path: Path) -> str:
    """Return the SHA-256 hex digest of a file without loading it all in memory."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_and_validate_fasta(path: Path, min_taxa: int = 3):
    """Read a protein FASTA and reject empty input, duplicate IDs, or invalid symbols."""
    if not path.is_file():
        raise ValueError(f"Input FASTA does not exist: {path}")
    records = list(SeqIO.parse(str(path), "fasta"))
    if len(records) < min_taxa:
        raise ValueError(f"Expected at least {min_taxa} sequences; found {len(records)}")
    ids = [record.id for record in records]
    if len(ids) != len(set(ids)):
        raise ValueError("Duplicate FASTA sequence IDs found")
    for record in records:
        sequence = str(record.seq).upper()
        if not sequence:
            raise ValueError(f"Empty sequence: {record.id}")
        invalid = sorted(set(sequence) - VALID_AA)
        if invalid:
            raise ValueError(f"Invalid protein symbols in {record.id}: {''.join(invalid)}")
    return records


def alignment_summary(path: Path) -> dict:
    """Summarize a protein alignment; all sequences must have equal length."""
    records = list(SeqIO.parse(str(path), "fasta"))
    if not records:
        raise ValueError(f"No sequences found in alignment: {path}")
    lengths = {len(record.seq) for record in records}
    if len(lengths) != 1:
        raise ValueError(f"Sequences have unequal lengths in alignment: {sorted(lengths)}")
    site_count = next(iter(lengths))
    if site_count == 0:
        raise ValueError("Alignment has zero sites")
    sequences = [str(record.seq).upper() for record in records]
    occupancy = []
    variable_sites = 0
    informative_sites = 0
    for column_index in range(site_count):
        observed = [seq[column_index] for seq in sequences if seq[column_index] not in MISSING_AA]
        occupancy.append(len(observed) / len(sequences))
        counts: dict[str, int] = {}
        for residue in observed:
            counts[residue] = counts.get(residue, 0) + 1
        if len(counts) > 1:
            variable_sites += 1
            if sum(count >= 2 for count in counts.values()) >= 2:
                informative_sites += 1
    return {
        "taxa": len(records),
        "sites": site_count,
        "gap_or_missing_fraction": round(1 - (sum(occupancy) / site_count), 6),
        "columns_occupancy_lt_25pct": sum(value < 0.25 for value in occupancy),
        "columns_occupancy_lt_50pct": sum(value < 0.50 for value in occupancy),
        "variable_sites": variable_sites,
        "parsimony_informative_sites": informative_sites,
        "sequence_ids_unique": len({record.id for record in records}) == len(records),
    }


def resolve_executable(candidates: Iterable[str]) -> str:
    """Find an installed executable from a list of common command names."""
    for candidate in candidates:
        located = shutil.which(candidate)
        if located:
            return located
    raise FileNotFoundError(
        "Could not find any of these executables on PATH: " + ", ".join(candidates)
    )


def run_logged(command: list[str], log_path: Path, stdout_path: Path | None = None) -> None:
    """Run a command and preserve its console output in a log file."""
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w", encoding="utf-8", errors="replace") as log_handle:
        if stdout_path is None:
            result = subprocess.run(
                command, stdout=log_handle, stderr=subprocess.STDOUT,
                text=True, errors="replace", check=False,
            )
        else:
            stdout_path.parent.mkdir(parents=True, exist_ok=True)
            with stdout_path.open("w", encoding="utf-8", errors="replace") as output_handle:
                result = subprocess.run(
                    command, stdout=output_handle, stderr=log_handle,
                    text=True, errors="replace", check=False,
                )
    if result.returncode != 0:
        raise RuntimeError(
            f"Command failed with exit code {result.returncode}. "
            f"See log: {log_path}\nCommand: {' '.join(command)}"
        )


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def run_pipeline(input_path: Path, outdir: Path, model: str, bootstrap: int, threads: int, seed: int) -> dict:
    """Run MAFFT, ClipKIT smart-gap trimming, and two IQ-TREE ML analyses."""
    records = read_and_validate_fasta(input_path, min_taxa=3)
    mafft = resolve_executable(["mafft", "mafft.bat"])
    clipkit = resolve_executable(["clipkit"])
    iqtree = resolve_executable(["iqtree3", "iqtree2", "iqtree"])

    outdir = outdir.resolve()
    full_dir = outdir / "alignments_full"
    trimmed_dir = outdir / "alignments_trimmed"
    tree_dir = outdir / "trees"
    log_dir = outdir / "logs"
    for directory in (full_dir, trimmed_dir, tree_dir, log_dir):
        directory.mkdir(parents=True, exist_ok=True)

    full_alignment = full_dir / "mafft_full.fasta"
    trimmed_alignment = trimmed_dir / "clipkit_smart_gap.fasta"
    run_logged(
        [mafft, "--auto", "--thread", str(threads), str(input_path.resolve())],
        log_dir / "mafft.log", stdout_path=full_alignment,
    )
    full_qc = alignment_summary(full_alignment)
    if full_qc["taxa"] != len(records):
        raise RuntimeError("MAFFT output taxon count differs from the input")

    run_logged(
        [clipkit, str(full_alignment), "-m", "smart-gap", "-o", str(trimmed_alignment)],
        log_dir / "clipkit.log",
    )
    trimmed_qc = alignment_summary(trimmed_alignment)
    if trimmed_qc["taxa"] != len(records):
        raise RuntimeError("ClipKIT output taxon count differs from the input")

    tree_records = []
    for label, alignment, qc in (
        ("full", full_alignment, full_qc),
        ("clipkit_smart_gap", trimmed_alignment, trimmed_qc),
    ):
        prefix = tree_dir / f"{label}.{model}.bootstrap{bootstrap}"
        run_logged(
            [iqtree, "-s", str(alignment), "-m", model, "-b", str(bootstrap),
             "-nt", str(threads), "-seed", str(seed), "-pre", str(prefix)],
            log_dir / f"iqtree_{label}.log",
        )
        treefile = Path(str(prefix) + ".treefile")
        report = Path(str(prefix) + ".iqtree")
        if not treefile.is_file() or not report.is_file():
            raise RuntimeError(
                f"Expected IQ-TREE output missing for {label}. Inspect {log_dir / f'iqtree_{label}.log'}"
            )
        tree_records.append({
            "alignment_variant": label,
            "taxa": qc["taxa"],
            "alignment_sites": qc["sites"],
            "alignment_file": str(alignment.relative_to(outdir)),
            "alignment_sha256": sha256_file(alignment),
            "tree_file": str(treefile.relative_to(outdir)),
            "tree_sha256": sha256_file(treefile),
            "iqtree_report": str(report.relative_to(outdir)),
        })

    qc_rows = [
        {"alignment_variant": "full_mafft", "alignment_file": str(full_alignment.relative_to(outdir)), **full_qc},
        {"alignment_variant": "clipkit_smart_gap", "alignment_file": str(trimmed_alignment.relative_to(outdir)), **trimmed_qc},
    ]
    write_csv(outdir / "alignment_qc_summary.csv", qc_rows)
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "input_filename": input_path.name,
        "input_sha256": sha256_file(input_path.resolve()),
        "input_taxa": len(records),
        "settings": {
            "aligner": "MAFFT --auto",
            "trimming": "ClipKIT smart-gap; full alignment also analyzed",
            "tree_method": "IQ-TREE maximum likelihood",
            "model": model,
            "bootstrap_replicates": bootstrap,
            "threads": threads,
            "random_seed": seed,
        },
        "tool_paths": {"mafft": mafft, "clipkit": clipkit, "iqtree": iqtree},
        "analyses": tree_records,
        "note": "Tree proximity alone does not establish orthology or biological causation.",
    }
    (outdir / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    return manifest


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path, help="Protein FASTA input file")
    parser.add_argument("--outdir", type=Path, default=Path("results/run"), help="Local output directory")
    parser.add_argument("--model", default="JTT", help="IQ-TREE amino-acid model (default: JTT)")
    parser.add_argument("--bootstrap", type=int, default=1000, help="Standard bootstrap replicates (default: 1000)")
    parser.add_argument("--threads", type=int, default=4, help="Threads per tool (default: 4)")
    parser.add_argument("--seed", type=int, default=20261018, help="Random seed for reproducibility")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if args.bootstrap < 1:
        print("ERROR: --bootstrap must be at least 1", file=sys.stderr)
        return 2
    if args.threads < 1:
        print("ERROR: --threads must be at least 1", file=sys.stderr)
        return 2
    try:
        manifest = run_pipeline(args.input, args.outdir, args.model, args.bootstrap, args.threads, args.seed)
    except Exception as exc:  # CLI boundary: report a useful failure instead of a traceback by default.
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    print(f"Finished {len(manifest['analyses'])} analyses. Local outputs: {args.outdir.resolve()}")
    print("Review the QC summary, manifest, IQ-TREE reports, and logs before interpreting results.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
