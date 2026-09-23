from __future__ import annotations

import gzip
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator, TextIO

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

FASTA_SUFFIXES = (".fa", ".fasta", ".fna", ".fa.gz", ".fasta.gz", ".fna.gz")


def is_fasta_path(path: Path) -> bool:
    name = path.name.lower()
    return any(name.endswith(suffix) for suffix in FASTA_SUFFIXES)


@contextmanager
def open_text(path: str | Path) -> Iterator[TextIO]:
    path = Path(path)
    if path.suffix == ".gz":
        with gzip.open(path, "rt") as handle:
            yield handle
    else:
        with path.open("rt") as handle:
            yield handle


def iter_fasta(path: str | Path) -> Iterator[SeqRecord]:
    with open_text(path) as handle:
        yield from SeqIO.parse(handle, "fasta")


def sequence_ids(path: str | Path) -> list[str]:
    return [record.id for record in iter_fasta(path)]


def genome_sequence(path: str | Path) -> str:
    """Concatenate all records in one assembly FASTA into a single uppercase string."""
    return "".join(str(record.seq).upper() for record in iter_fasta(path))


def discover_fastas(root: str | Path) -> list[Path]:
    root = Path(root)
    return sorted(path for path in root.rglob("*") if path.is_file() and is_fasta_path(path))
