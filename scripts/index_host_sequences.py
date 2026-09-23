#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from phage_host_pipeline.mapping import build_sequence_index, read_manifest


def main() -> None:
    p = argparse.ArgumentParser(description="Index host FASTA record IDs to genome IDs and source files.")
    p.add_argument("--manifest", required=True)
    p.add_argument("--output", required=True)
    args = p.parse_args()

    rows = build_sequence_index(read_manifest(args.manifest))
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=["sequence_id", "genome_id", "fasta_path"])
        writer.writeheader()
        writer.writerows(rows)
    print(f"sequence_ids={len(rows)} output={out}")


if __name__ == "__main__":
    main()
