#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from phage_host_pipeline.blast import BLAST_OUTFMT
from phage_host_pipeline.mapping import read_manifest


def select_chunk(rows: list[dict[str, str]], chunk_index: int, chunk_count: int) -> list[dict[str, str]]:
    if chunk_count < 1:
        raise ValueError("chunk_count must be >= 1")
    if not 0 <= chunk_index < chunk_count:
        raise ValueError("chunk_index must satisfy 0 <= chunk_index < chunk_count")
    return [row for i, row in enumerate(rows) if i % chunk_count == chunk_index]


def main() -> None:
    p = argparse.ArgumentParser(description="Run blastn for a chunk of phage genomes against a host database.")
    p.add_argument("--phage-manifest", required=True)
    p.add_argument("--db-prefix", required=True)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--evalue", type=float, default=1e-5)
    p.add_argument("--chunk-count", type=int, default=1)
    p.add_argument("--chunk-index", type=int, default=None)
    p.add_argument("--max-target-seqs", type=int, default=None)
    args = p.parse_args()

    chunk_index = args.chunk_index
    if chunk_index is None:
        chunk_index = int(os.environ.get("SLURM_ARRAY_TASK_ID", "0"))

    rows = select_chunk(read_manifest(args.phage_manifest), chunk_index, args.chunk_count)
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    print(f"chunk={chunk_index}/{args.chunk_count} phage_genomes={len(rows)}")

    for row in rows:
        phage_id = row["genome_id"]
        output = outdir / f"{phage_id}.blast.tsv"
        if output.exists():
            print(f"skip existing: {output}")
            continue
        temporary_output = output.with_name(f"{output.name}.{os.getpid()}.tmp")
        temporary_output.unlink(missing_ok=True)
        cmd = [
            "blastn",
            "-task", "blastn",
            "-query", row["fasta_path"],
            "-db", args.db_prefix,
            "-evalue", str(args.evalue),
            "-outfmt", BLAST_OUTFMT,
            "-num_threads", str(args.threads),
            "-out", str(temporary_output),
        ]
        if args.max_target_seqs is not None:
            cmd.extend(["-max_target_seqs", str(args.max_target_seqs)])
        try:
            subprocess.run(cmd, check=True)
            temporary_output.replace(output)
        finally:
            temporary_output.unlink(missing_ok=True)
        print(f"phage={phage_id} output={output}")


if __name__ == "__main__":
    main()
