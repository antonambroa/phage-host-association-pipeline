from __future__ import annotations

import itertools
import math
from collections import Counter
from pathlib import Path

from .fasta import genome_sequence, iter_fasta

DNA_BASES = "ACGT"
CODONS = tuple("".join(parts) for parts in itertools.product(DNA_BASES, repeat=3))


def gc_fraction(sequence: str) -> float:
    sequence = sequence.upper()
    canonical = [base for base in sequence if base in DNA_BASES]
    if not canonical:
        return 0.0
    gc = sum(base in {"G", "C"} for base in canonical)
    return gc / len(canonical)


def genome_basic_features(fasta_path: str | Path) -> dict[str, float | int]:
    seq = genome_sequence(fasta_path)
    canonical_length = sum(base in DNA_BASES for base in seq)
    return {
        "genome_length": len(seq),
        "canonical_bases": canonical_length,
        "gc_fraction": gc_fraction(seq),
    }


def codon_usage_from_cds_fasta(cds_fasta: str | Path) -> dict[str, float]:
    counts: Counter[str] = Counter()
    total = 0
    for record in iter_fasta(cds_fasta):
        seq = str(record.seq).upper()
        for offset in range(0, len(seq) - 2, 3):
            codon = seq[offset : offset + 3]
            if len(codon) == 3 and all(base in DNA_BASES for base in codon):
                counts[codon] += 1
                total += 1
    if total == 0:
        return {codon: 0.0 for codon in CODONS}
    return {codon: counts[codon] / total for codon in CODONS}


def cosine_similarity(a: dict[str, float], b: dict[str, float]) -> float:
    keys = sorted(set(a) | set(b))
    dot = sum(a.get(key, 0.0) * b.get(key, 0.0) for key in keys)
    norm_a = math.sqrt(sum(a.get(key, 0.0) ** 2 for key in keys))
    norm_b = math.sqrt(sum(b.get(key, 0.0) ** 2 for key in keys))
    if norm_a == 0.0 or norm_b == 0.0:
        return 0.0
    return dot / (norm_a * norm_b)
