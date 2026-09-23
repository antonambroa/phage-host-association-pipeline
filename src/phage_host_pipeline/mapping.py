from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterable

from .fasta import iter_fasta


def build_sequence_index(manifest_rows: Iterable[dict[str, str]]) -> list[dict[str, str]]:
    """Map every FASTA record id to its genome id and source FASTA.

    BLAST subject identifiers refer to sequence records, while downstream analyses
    usually work at assembly/genome level.  This index keeps that translation
    explicit and avoids regex-based accession guessing.
    """
    seen: dict[str, tuple[str, str]] = {}
    output: list[dict[str, str]] = []

    for row in manifest_rows:
        genome_id = row["genome_id"]
        fasta_path = row["fasta_path"]
        for record in iter_fasta(fasta_path):
            seqid = record.id
            previous = seen.get(seqid)
            current = (genome_id, fasta_path)
            if previous and previous != current:
                raise ValueError(
                    f"duplicate sequence id {seqid!r} occurs in multiple genomes: "
                    f"{previous[0]!r} and {genome_id!r}"
                )
            if previous:
                continue
            seen[seqid] = current
            output.append(
                {
                    "sequence_id": seqid,
                    "genome_id": genome_id,
                    "fasta_path": fasta_path,
                }
            )
    return output


def read_manifest(path: str | Path) -> list[dict[str, str]]:
    with Path(path).open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {"genome_id", "fasta_path"}
    if not rows:
        return []
    missing = required.difference(rows[0])
    if missing:
        raise ValueError(f"manifest missing required columns: {sorted(missing)}")
    return rows


def read_sequence_index(path: str | Path) -> dict[str, dict[str, str]]:
    with Path(path).open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        return {}
    required = {"sequence_id", "genome_id", "fasta_path"}
    missing = required.difference(rows[0])
    if missing:
        raise ValueError(f"sequence index missing required columns: {sorted(missing)}")
    return {row["sequence_id"]: row for row in rows}
