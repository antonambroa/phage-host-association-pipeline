#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import shutil
import subprocess
import tempfile
from pathlib import Path

import sys
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from phage_host_pipeline.mapping import read_manifest


def materialize_if_gz(path: str, directory: Path) -> Path:
    source = Path(path)
    if source.suffix != ".gz":
        return source
    target = directory / source.stem
    with gzip.open(source, "rb") as src, target.open("wb") as dst:
        shutil.copyfileobj(src, dst)
    return target


def main() -> None:
    p = argparse.ArgumentParser(description="Predict host CDS once per genome with Prodigal and cache the FASTA output.")
    p.add_argument("--host-manifest", required=True)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--prodigal-mode", choices=["single", "meta"], default="single")
    args = p.parse_args()

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    for row in read_manifest(args.host_manifest):
        output = outdir / f"{row['genome_id']}.cds.fna"
        if output.exists():
            print(f"skip existing: {output}")
            continue
        with tempfile.TemporaryDirectory(prefix="phage-host-prodigal-") as tmp:
            tmpdir = Path(tmp)
            input_path = materialize_if_gz(row["fasta_path"], tmpdir)
            temporary_output = tmpdir / "predicted.cds.fna"
            subprocess.run(
                ["prodigal", "-i", str(input_path), "-d", str(temporary_output), "-p", args.prodigal_mode, "-q"],
                check=True,
            )
            shutil.move(str(temporary_output), output)
        print(f"host={row['genome_id']} cds={output}")


if __name__ == "__main__":
    main()
