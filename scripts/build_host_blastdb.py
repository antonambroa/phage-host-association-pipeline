#!/usr/bin/env python3
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from phage_host_pipeline.fasta import open_text
from phage_host_pipeline.mapping import read_manifest


def concatenate(manifest: str, output: Path) -> int:
    rows = read_manifest(manifest)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w") as out:
        for row in rows:
            with open_text(row["fasta_path"]) as src:
                shutil.copyfileobj(src, out)
    return len(rows)


def main() -> None:
    p = argparse.ArgumentParser(description="Concatenate host FASTAs and build a nucleotide BLAST database.")
    p.add_argument("--manifest", required=True)
    p.add_argument("--combined-fasta", required=True)
    p.add_argument("--db-prefix", required=True)
    p.add_argument("--skip-makeblastdb", action="store_true", help="Only create the concatenated FASTA")
    args = p.parse_args()

    combined = Path(args.combined_fasta)
    n = concatenate(args.manifest, combined)
    print(f"host_genomes={n} combined_fasta={combined}")

    if args.skip_makeblastdb:
        return
    cmd = ["makeblastdb", "-in", str(combined), "-dbtype", "nucl", "-out", args.db_prefix, "-parse_seqids"]
    subprocess.run(cmd, check=True)
    print(f"blast_db={args.db_prefix}")


if __name__ == "__main__":
    main()
