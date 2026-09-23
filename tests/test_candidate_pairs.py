from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path


def test_hits_on_multiple_host_contigs_become_one_genome_pair(tmp_path: Path):
    blast_dir = tmp_path / "blast"
    blast_dir.mkdir()
    (blast_dir / "phage1.blast.tsv").write_text(
        "q1\tcontigA\t99\t600\t0\t0\t1\t600\t1\t600\t1e-50\t800\t1000\t1000\n"
        "q1\tcontigB\t98\t550\t0\t0\t10\t559\t5\t554\t1e-40\t700\t1000\t900\n"
    )
    index = tmp_path / "host_index.tsv"
    index.write_text(
        "sequence_id\tgenome_id\tfasta_path\n"
        "contigA\thost1\t/path/host1.fna\n"
        "contigB\thost1\t/path/host1.fna\n"
    )
    output = tmp_path / "pairs.tsv"

    subprocess.run(
        [
            sys.executable,
            "scripts/build_candidate_pairs.py",
            "--blast-dir", str(blast_dir),
            "--host-sequence-index", str(index),
            "--output", str(output),
        ],
        check=True,
        cwd=Path(__file__).resolve().parents[1],
    )

    with output.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 1
    assert rows[0]["host_id"] == "host1"
    assert rows[0]["host_sequence_count"] == "2"
    assert rows[0]["host_sequence_ids"] == "contigA;contigB"
    assert rows[0]["hit_count"] == "2"
