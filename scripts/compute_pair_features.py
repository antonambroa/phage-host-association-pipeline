#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from phage_host_pipeline.features import codon_usage_from_cds_fasta, genome_basic_features
from phage_host_pipeline.mapping import read_manifest


def main() -> None:
    p = argparse.ArgumentParser(description="Add simple genome-level features to candidate pair tables.")
    p.add_argument("--pairs", required=True)
    p.add_argument("--phage-manifest", required=True)
    p.add_argument("--output", required=True)
    p.add_argument(
        "--host-cds-dir",
        default=None,
        help="Optional directory containing cached <host_id>.cds.fna files. Host codon usage is serialized as JSON.",
    )
    args = p.parse_args()

    phage_paths = {row["genome_id"]: row["fasta_path"] for row in read_manifest(args.phage_manifest)}
    genome_cache: dict[str, dict[str, float | int]] = {}
    codon_cache: dict[str, str] = {}

    with Path(args.pairs).open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise SystemExit("candidate-pair table contains no rows")

    fields = list(rows[0].keys()) + ["host_length", "phage_length", "host_gc", "phage_gc", "gc_abs_diff"]
    if args.host_cds_dir:
        fields.append("host_codon_usage_json")

    for row in rows:
        host_path = row["host_fasta"]
        phage_path = phage_paths.get(row["phage_id"])
        if not phage_path:
            raise ValueError(f"phage {row['phage_id']!r} is absent from the phage manifest")

        if host_path not in genome_cache:
            genome_cache[host_path] = genome_basic_features(host_path)
        if phage_path not in genome_cache:
            genome_cache[phage_path] = genome_basic_features(phage_path)
        host_features = genome_cache[host_path]
        phage_features = genome_cache[phage_path]
        row["host_length"] = str(host_features["genome_length"])
        row["phage_length"] = str(phage_features["genome_length"])
        row["host_gc"] = f"{host_features['gc_fraction']:.6f}"
        row["phage_gc"] = f"{phage_features['gc_fraction']:.6f}"
        row["gc_abs_diff"] = f"{abs(float(host_features['gc_fraction']) - float(phage_features['gc_fraction'])):.6f}"

        if args.host_cds_dir:
            host_id = row["host_id"]
            cds_path = Path(args.host_cds_dir) / f"{host_id}.cds.fna"
            if cds_path.exists():
                codon_cache.setdefault(host_id, json.dumps(codon_usage_from_cds_fasta(cds_path), sort_keys=True))
                row["host_codon_usage_json"] = codon_cache[host_id]
            else:
                row["host_codon_usage_json"] = ""

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    print(f"rows={len(rows)} output={out}")


if __name__ == "__main__":
    main()
