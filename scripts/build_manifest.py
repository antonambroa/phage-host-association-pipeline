#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from phage_host_pipeline.fasta import discover_fastas


def genome_id_from_path(path: Path) -> str:
    name = path.name
    for suffix in (".fna.gz", ".fasta.gz", ".fa.gz", ".fna", ".fasta", ".fa"):
        if name.lower().endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


def main() -> None:
    p = argparse.ArgumentParser(description="Build a two-column genome FASTA manifest.")
    p.add_argument("--root", required=True, help="Directory containing FASTA files")
    p.add_argument("--output", required=True)
    args = p.parse_args()

    fastas = discover_fastas(args.root)
    if not fastas:
        raise SystemExit(f"no FASTA files found under {args.root}")

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=["genome_id", "fasta_path"])
        writer.writeheader()
        seen_ids = set()
        for path in fastas:
            genome_id = genome_id_from_path(path)
            if genome_id in seen_ids:
                raise ValueError(f"duplicate genome_id {genome_id!r}; rename files or provide a curated manifest")
            seen_ids.add(genome_id)
            writer.writerow({"genome_id": genome_id, "fasta_path": str(path.resolve())})
    print(f"genomes={len(fastas)} output={out}")


if __name__ == "__main__":
    main()
