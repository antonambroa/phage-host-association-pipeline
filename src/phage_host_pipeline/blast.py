from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator

BLAST_FIELDS = (
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "qlen",
    "slen",
)

BLAST_OUTFMT = "6 " + " ".join(BLAST_FIELDS)


@dataclass(frozen=True)
class BlastHit:
    qseqid: str
    sseqid: str
    pident: float
    length: int
    mismatch: int
    gapopen: int
    qstart: int
    qend: int
    sstart: int
    send: int
    evalue: float
    bitscore: float
    qlen: int
    slen: int

    @property
    def query_coverage(self) -> float:
        return self.length / self.qlen if self.qlen else 0.0

    @property
    def subject_coverage(self) -> float:
        return self.length / self.slen if self.slen else 0.0


def parse_blast_row(row: list[str]) -> BlastHit:
    if len(row) != len(BLAST_FIELDS):
        raise ValueError(f"expected {len(BLAST_FIELDS)} BLAST columns, got {len(row)}")
    return BlastHit(
        qseqid=row[0],
        sseqid=row[1],
        pident=float(row[2]),
        length=int(row[3]),
        mismatch=int(row[4]),
        gapopen=int(row[5]),
        qstart=int(row[6]),
        qend=int(row[7]),
        sstart=int(row[8]),
        send=int(row[9]),
        evalue=float(row[10]),
        bitscore=float(row[11]),
        qlen=int(row[12]),
        slen=int(row[13]),
    )


def read_blast_tsv(path: str | Path) -> Iterator[BlastHit]:
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        return
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for line_number, row in enumerate(reader, start=1):
            if not row:
                continue
            try:
                yield parse_blast_row(row)
            except (TypeError, ValueError) as exc:
                raise ValueError(f"{path}:{line_number}: {exc}") from exc


def hit_passes(
    hit: BlastHit,
    *,
    min_alignment_length: int = 500,
    min_identity: float = 0.0,
    max_evalue: float = 1e-5,
    min_query_coverage: float = 0.0,
) -> bool:
    return (
        hit.length >= min_alignment_length
        and hit.pident >= min_identity
        and hit.evalue <= max_evalue
        and hit.query_coverage >= min_query_coverage
    )


def aggregate_hits(hits: Iterable[BlastHit]) -> dict[str, float | int]:
    hits = list(hits)
    if not hits:
        raise ValueError("cannot aggregate an empty hit set")
    return {
        "hit_count": len(hits),
        "max_pident": max(hit.pident for hit in hits),
        "max_alignment_length": max(hit.length for hit in hits),
        "min_evalue": min(hit.evalue for hit in hits),
        "max_bitscore": max(hit.bitscore for hit in hits),
        "best_query_coverage": max(hit.query_coverage for hit in hits),
        "best_subject_coverage": max(hit.subject_coverage for hit in hits),
    }


def format_float(value: float) -> str:
    if value == 0:
        return "0"
    if abs(value) < 1e-3 or abs(value) >= 1e4:
        return f"{value:.6e}"
    if math.isclose(value, round(value)):
        return str(int(round(value)))
    return f"{value:.6f}".rstrip("0").rstrip(".")
