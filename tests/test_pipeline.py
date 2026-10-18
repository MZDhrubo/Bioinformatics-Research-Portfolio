from pathlib import Path

import pytest

from scripts.run_pipeline import alignment_summary, read_and_validate_fasta, sha256_file


def write_fasta(path: Path, content: str) -> Path:
    path.write_text(content, encoding="utf-8")
    return path


def test_valid_fasta_reads_unique_protein_records(tmp_path):
    path = write_fasta(
        tmp_path / "valid.fasta",
        ">A\nACDEFG\n>B\nACDEFA\n>C\nACDEYA\n",
    )
    records = read_and_validate_fasta(path)
    assert len(records) == 3
    assert [record.id for record in records] == ["A", "B", "C"]


def test_duplicate_ids_are_rejected(tmp_path):
    path = write_fasta(tmp_path / "duplicate.fasta", ">A\nACDE\n>A\nACDF\n>C\nACDG\n")
    with pytest.raises(ValueError, match="Duplicate"):
        read_and_validate_fasta(path)


def test_invalid_residue_is_rejected(tmp_path):
    path = write_fasta(tmp_path / "invalid.fasta", ">A\nACD!\n>B\nACDE\n>C\nACDF\n")
    with pytest.raises(ValueError, match="Invalid protein symbols"):
        read_and_validate_fasta(path)


def test_too_few_sequences_are_rejected(tmp_path):
    path = write_fasta(tmp_path / "small.fasta", ">A\nACDE\n>B\nACDF\n")
    with pytest.raises(ValueError, match="at least 3"):
        read_and_validate_fasta(path)


def test_alignment_summary_counts_sites_and_taxa(tmp_path):
    path = write_fasta(
        tmp_path / "aligned.fasta",
        ">A\nAC-D\n>B\nAT-D\n>C\nATGD\n>D\nACGD\n",
    )
    summary = alignment_summary(path)
    assert summary["taxa"] == 4
    assert summary["sites"] == 4
    assert summary["sequence_ids_unique"] is True
    assert summary["variable_sites"] >= 1
    assert 0 <= summary["gap_or_missing_fraction"] <= 1


def test_alignment_summary_rejects_unequal_lengths(tmp_path):
    path = write_fasta(tmp_path / "bad_alignment.fasta", ">A\nACD\n>B\nACDE\n>C\nACDF\n")
    with pytest.raises(ValueError, match="unequal lengths"):
        alignment_summary(path)


def test_sha256_is_stable(tmp_path):
    path = write_fasta(tmp_path / "small.txt", "portfolio test\n")
    assert sha256_file(path) == sha256_file(path)
    assert len(sha256_file(path)) == 64
