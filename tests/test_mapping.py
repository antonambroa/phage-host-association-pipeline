from pathlib import Path

import pytest

from phage_host_pipeline.mapping import build_sequence_index


def test_build_sequence_index_keeps_full_accession(tmp_path: Path):
    fasta = tmp_path / "host.fna"
    fasta.write_text(">NZ_CP059690.1 example\nACGT\n>NC_000001.2 example\nGGCC\n")
    rows = build_sequence_index([{"genome_id": "GCF_demo", "fasta_path": str(fasta)}])
    assert [row["sequence_id"] for row in rows] == ["NZ_CP059690.1", "NC_000001.2"]
    assert all(row["genome_id"] == "GCF_demo" for row in rows)


def test_duplicate_sequence_ids_across_genomes_fail(tmp_path: Path):
    a = tmp_path / "a.fna"
    b = tmp_path / "b.fna"
    a.write_text(">same\nAAAA\n")
    b.write_text(">same\nCCCC\n")
    with pytest.raises(ValueError, match="duplicate sequence id"):
        build_sequence_index([
            {"genome_id": "a", "fasta_path": str(a)},
            {"genome_id": "b", "fasta_path": str(b)},
        ])
