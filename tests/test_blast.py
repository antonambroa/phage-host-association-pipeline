from pathlib import Path

import pytest

from phage_host_pipeline.blast import aggregate_hits, hit_passes, read_blast_tsv


ROW1 = "q\th\t97.5\t600\t15\t0\t1\t600\t2\t601\t1e-80\t900\t1000\t5000000\n"
ROW2 = "q\th\t92.0\t520\t40\t1\t20\t539\t900\t1419\t1e-40\t500\t1000\t5000000\n"


def test_parse_filter_and_aggregate(tmp_path: Path):
    path = tmp_path / "hits.tsv"
    path.write_text(ROW1 + ROW2)
    hits = list(read_blast_tsv(path))
    assert len(hits) == 2
    assert hit_passes(hits[0], min_alignment_length=500, min_identity=95)
    assert not hit_passes(hits[1], min_alignment_length=500, min_identity=95)
    summary = aggregate_hits(hits)
    assert summary["hit_count"] == 2
    assert summary["max_pident"] == pytest.approx(97.5)
    assert summary["max_alignment_length"] == 600
    assert summary["min_evalue"] == pytest.approx(1e-80)
    assert summary["best_query_coverage"] == pytest.approx(0.6)


def test_bad_column_count_is_reported(tmp_path: Path):
    path = tmp_path / "bad.tsv"
    path.write_text("a\tb\tc\n")
    with pytest.raises(ValueError, match="expected 14 BLAST columns"):
        list(read_blast_tsv(path))
