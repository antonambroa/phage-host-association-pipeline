from pathlib import Path

import pytest

from phage_host_pipeline.features import codon_usage_from_cds_fasta, gc_fraction


def test_gc_fraction_ignores_ambiguous_bases():
    assert gc_fraction("GCGCATATNN") == pytest.approx(0.5)
    assert gc_fraction("NNNN") == 0.0


def test_codon_usage_from_cds(tmp_path: Path):
    cds = tmp_path / "cds.fna"
    cds.write_text(">cds1\nAAAGGGAAA\n>cds2\nAAATTT\n")
    usage = codon_usage_from_cds_fasta(cds)
    assert sum(usage.values()) == pytest.approx(1.0)
    assert usage["AAA"] == pytest.approx(3 / 5)
    assert usage["GGG"] == pytest.approx(1 / 5)
    assert usage["TTT"] == pytest.approx(1 / 5)
