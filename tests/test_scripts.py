"""End-to-end checks for the command-line scripts."""

from __future__ import annotations

import csv
import os
import runpy
import subprocess
import sys
from pathlib import Path

import pytest

from phage_host_pipeline.fasta import sequence_ids

ROOT = Path(__file__).resolve().parents[1]


def run_script(name: str, *args: str, check: bool = True) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, f"scripts/{name}", *args],
        cwd=ROOT,
        check=check,
        capture_output=True,
        text=True,
    )


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_readme_example_end_to_end(tmp_path: Path):
    """The tiny synthetic example behaves as the README says."""
    run_script("build_manifest.py", "--root", "examples/hosts", "--output", str(tmp_path / "hosts.tsv"))
    run_script("build_manifest.py", "--root", "examples/phages", "--output", str(tmp_path / "phages.tsv"))
    run_script(
        "index_host_sequences.py",
        "--manifest", str(tmp_path / "hosts.tsv"),
        "--output", str(tmp_path / "host_index.tsv"),
    )
    run_script(
        "build_candidate_pairs.py",
        "--blast-dir", "examples/blast",
        "--host-sequence-index", str(tmp_path / "host_index.tsv"),
        "--output", str(tmp_path / "pairs.tsv"),
        "--min-alignment-length", "500",
    )
    run_script(
        "compute_pair_features.py",
        "--pairs", str(tmp_path / "pairs.tsv"),
        "--phage-manifest", str(tmp_path / "phages.tsv"),
        "--output", str(tmp_path / "features.tsv"),
    )

    rows = read_tsv(tmp_path / "features.tsv")
    assert [(row["phage_id"], row["host_id"]) for row in rows] == [("phage_demo", "host_alpha")]
    assert rows[0]["hit_count"] == "2"
    assert rows[0]["host_length"] == "1200"
    assert rows[0]["phage_length"] == "1000"


def test_concatenation_adds_missing_newline(tmp_path: Path):
    hosts = tmp_path / "hosts"
    hosts.mkdir()
    (hosts / "h1.fna").write_text(">h1\nACGT")  # no trailing newline
    (hosts / "h2.fna").write_text(">h2\nGGCC\n")
    run_script("build_manifest.py", "--root", str(hosts), "--output", str(tmp_path / "hosts.tsv"))

    combined = tmp_path / "combined.fna"
    run_script(
        "build_host_blastdb.py",
        "--manifest", str(tmp_path / "hosts.tsv"),
        "--combined-fasta", str(combined),
        "--db-prefix", str(tmp_path / "db"),
        "--skip-makeblastdb",
    )

    assert sequence_ids(combined) == ["h1", "h2"]
    assert not list(tmp_path.glob("*.tmp"))


def write_unmapped_case(tmp_path: Path) -> tuple[Path, Path]:
    blast_dir = tmp_path / "blast"
    blast_dir.mkdir()
    (blast_dir / "phage1.blast.tsv").write_text(
        "q1\tcontigA\t99\t600\t0\t0\t1\t600\t1\t600\t1e-50\t800\t1000\t1000\n"
        "q1\tref|contigB|\t98\t550\t0\t0\t10\t559\t5\t554\t1e-40\t700\t1000\t900\n"
    )
    index = tmp_path / "host_index.tsv"
    index.write_text(
        "sequence_id\tgenome_id\tfasta_path\n"
        "contigA\thost1\t/path/host1.fna\n"
        "contigB\thost1\t/path/host1.fna\n"
    )
    return blast_dir, index


def test_unmapped_subject_ids_stop_the_run(tmp_path: Path):
    blast_dir, index = write_unmapped_case(tmp_path)
    output = tmp_path / "pairs.tsv"
    result = run_script(
        "build_candidate_pairs.py",
        "--blast-dir", str(blast_dir),
        "--host-sequence-index", str(index),
        "--output", str(output),
        check=False,
    )
    assert result.returncode != 0
    assert "ref|contigB|" in result.stderr
    assert not output.exists()


def test_allow_unmapped_skips_them(tmp_path: Path):
    blast_dir, index = write_unmapped_case(tmp_path)
    output = tmp_path / "pairs.tsv"
    result = run_script(
        "build_candidate_pairs.py",
        "--blast-dir", str(blast_dir),
        "--host-sequence-index", str(index),
        "--output", str(output),
        "--allow-unmapped",
    )
    assert "unmapped_subject_ids=1" in result.stderr
    rows = read_tsv(output)
    assert len(rows) == 1
    assert rows[0]["hit_count"] == "1"


def fake_prodigal(bin_dir: Path, *, succeed: bool) -> None:
    """Stand-in for prodigal: writes some CDS output to the -d path, then exits."""
    bin_dir.mkdir(exist_ok=True)
    script = bin_dir / "prodigal"
    script.write_text(
        "#!/usr/bin/env python3\n"
        "import sys\n"
        "args = sys.argv[1:]\n"
        "out = args[args.index('-d') + 1]\n"
        "open(out, 'w').write('>cds1\\nATGAAATAA\\n')\n"
        f"sys.exit({0 if succeed else 1})\n"
    )
    script.chmod(0o755)


def run_predict_host_cds(tmp_path: Path, *, succeed: bool) -> tuple[subprocess.CompletedProcess, Path]:
    fake_prodigal(tmp_path / "bin", succeed=succeed)
    host = tmp_path / "host1.fna"
    host.write_text(">h\nACGTACGT\n")
    manifest = tmp_path / "hosts.tsv"
    manifest.write_text(f"genome_id\tfasta_path\nhost1\t{host}\n")
    outdir = tmp_path / "cds"
    env = {**os.environ, "PATH": f"{tmp_path / 'bin'}{os.pathsep}{os.environ.get('PATH', '')}"}
    result = subprocess.run(
        [sys.executable, "scripts/predict_host_cds.py", "--host-manifest", str(manifest), "--output-dir", str(outdir)],
        cwd=ROOT,
        env=env,
        capture_output=True,
        text=True,
    )
    return result, outdir


def test_predict_host_cds_promotes_output_only_on_success(tmp_path: Path):
    result, outdir = run_predict_host_cds(tmp_path, succeed=True)
    assert result.returncode == 0, result.stderr
    assert (outdir / "host1.cds.fna").read_text().startswith(">cds1")
    assert not list(outdir.glob("*.tmp"))


def test_failed_prodigal_leaves_no_output(tmp_path: Path):
    result, outdir = run_predict_host_cds(tmp_path, succeed=False)
    assert result.returncode != 0
    assert not (outdir / "host1.cds.fna").exists()
    assert not list(outdir.glob("*.tmp"))


def test_host_codon_usage_is_computed_once_per_host(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    import phage_host_pipeline.features as features

    calls: list[Path] = []
    original = features.codon_usage_from_cds_fasta

    def counting(path):
        calls.append(Path(path))
        return original(path)

    monkeypatch.setattr(features, "codon_usage_from_cds_fasta", counting)

    host = tmp_path / "host.fna"
    host.write_text(">h\nACGTACGT\n")
    phage = tmp_path / "phage.fna"
    phage.write_text(">p\nGGCC\n")
    cds_dir = tmp_path / "cds"
    cds_dir.mkdir()
    (cds_dir / "host1.cds.fna").write_text(">c1\nAAAGGG\n")
    phage_manifest = tmp_path / "phages.tsv"
    phage_manifest.write_text(f"genome_id\tfasta_path\np1\t{phage}\np2\t{phage}\np3\t{phage}\n")
    pairs = tmp_path / "pairs.tsv"
    pairs.write_text(
        "phage_id\thost_id\thost_fasta\n"
        + "".join(f"p{i}\thost1\t{host}\n" for i in (1, 2, 3))
    )
    output = tmp_path / "features.tsv"

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "compute_pair_features.py",
            "--pairs", str(pairs),
            "--phage-manifest", str(phage_manifest),
            "--host-cds-dir", str(cds_dir),
            "--output", str(output),
        ],
    )
    runpy.run_path(str(ROOT / "scripts" / "compute_pair_features.py"), run_name="__main__")

    assert len(calls) == 1
    rows = read_tsv(output)
    assert len(rows) == 3
    assert len({row["host_codon_usage_json"] for row in rows}) == 1
