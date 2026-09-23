import gzip
from pathlib import Path

from phage_host_pipeline.fasta import genome_sequence, sequence_ids


def test_multirecord_fasta_and_gzip(tmp_path: Path):
    text = ">a\nACGT\n>b\nGGCC\n"
    plain = tmp_path / "x.fna"
    plain.write_text(text)
    gz = tmp_path / "x.fna.gz"
    with gzip.open(gz, "wt") as handle:
        handle.write(text)

    assert sequence_ids(plain) == ["a", "b"]
    assert sequence_ids(gz) == ["a", "b"]
    assert genome_sequence(plain) == "ACGTGGCC"
    assert genome_sequence(gz) == "ACGTGGCC"
