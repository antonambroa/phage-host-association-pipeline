#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from phage_host_pipeline.blast import aggregate_hits, format_float, hit_passes, read_blast_tsv
from phage_host_pipeline.mapping import read_sequence_index


def main() -> None:
    p = argparse.ArgumentParser(
        description="Aggregate BLAST evidence into candidate phage-host associations. "
        "Candidates are sequence-search evidence, not validated host assignments."
    )
    p.add_argument("--blast-dir", required=True)
    p.add_argument("--host-sequence-index", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--min-alignment-length", type=int, default=500)
    p.add_argument("--min-identity", type=float, default=0.0)
    p.add_argument("--max-evalue", type=float, default=1e-5)
    p.add_argument("--min-query-coverage", type=float, default=0.0)
    args = p.parse_args()

    host_index = read_sequence_index(args.host_sequence_index)
    grouped: dict[tuple[str, str, str], list] = defaultdict(list)
    subject_ids: dict[tuple[str, str, str], set[str]] = defaultdict(set)
    unmatched = set()

    for blast_path in sorted(Path(args.blast_dir).glob("*.blast.tsv")):
        phage_id = blast_path.name.removesuffix(".blast.tsv")
        for hit in read_blast_tsv(blast_path):
            if not hit_passes(
                hit,
                min_alignment_length=args.min_alignment_length,
                min_identity=args.min_identity,
                max_evalue=args.max_evalue,
                min_query_coverage=args.min_query_coverage,
            ):
                continue
            host = host_index.get(hit.sseqid)
            if host is None:
                unmatched.add(hit.sseqid)
                continue
            key = (phage_id, host["genome_id"], host["fasta_path"])
            grouped[key].append(hit)
            subject_ids[key].add(hit.sseqid)

    fields = [
        "phage_id", "host_id", "host_fasta", "host_sequence_count", "host_sequence_ids",
        "hit_count", "max_pident", "max_alignment_length", "min_evalue",
        "max_bitscore", "best_query_coverage", "best_subject_coverage",
    ]
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        for key in sorted(grouped):
            phage_id, host_id, host_fasta = key
            summary = aggregate_hits(grouped[key])
            ids = sorted(subject_ids[key])
            writer.writerow(
                {
                    "phage_id": phage_id,
                    "host_id": host_id,
                    "host_fasta": host_fasta,
                    "host_sequence_count": len(ids),
                    "host_sequence_ids": ";".join(ids),
                    **{k: format_float(v) if isinstance(v, float) else v for k, v in summary.items()},
                }
            )

    print(f"candidate_pairs={len(grouped)} output={out}")
    if unmatched:
        print(f"unmapped_subject_ids={len(unmatched)}", file=sys.stderr)


if __name__ == "__main__":
    main()
