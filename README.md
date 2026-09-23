# phage-host-association-pipeline

[![tests](https://github.com/antonambroa/phage-host-association-pipeline/actions/workflows/tests.yml/badge.svg)](https://github.com/antonambroa/phage-host-association-pipeline/actions/workflows/tests.yml)

A small, inspectable workflow for turning phage-vs-bacterial sequence-search results into **candidate phage-host association tables** and simple genome-level features.

The code grew out of an exploratory HPC workflow built around RefSeq genomes, BLAST and Slurm. The public version focuses on the parts that are reusable: explicit FASTA manifests, stable sequence-to-genome mapping, chunked BLAST execution, preservation of search evidence, configurable candidate filters and cached feature extraction.

> A BLAST hit is sequence-similarity evidence, not proof that a bacteriophage infects a particular host. Shared mobile elements, prophage sequence and other homologous regions can all produce strong matches. The output of this repository is therefore described as **candidate association evidence** and should be combined with independent biological or computational evidence before host assignment.

## Workflow

```text
host FASTAs ──> manifest ──> sequence-id index ──> BLAST database
                                                   ▲
                                                   │
phage FASTAs ─> manifest ──> chunked blastn ───────┘
                              │
                              v
                         BLAST evidence
                              │
                    configurable filtering
                              │
                              v
                   candidate association table
                              │
                              v
                    simple genome features
```

## Why this version differs from the exploratory scripts

The original working directory was useful for testing the idea, but several choices did not scale well or were difficult to reproduce. This version deliberately changes them:

- no cluster-specific absolute paths or personal directory names;
- BLAST filtering is separated from BLAST execution, so raw evidence is retained;
- thresholds are explicit command-line parameters rather than hidden constants;
- BLAST subject IDs are mapped from actual FASTA record IDs rather than guessed with regular expressions;
- multi-record host assemblies are indexed completely, not only from the first FASTA header;
- temporary files use isolated temporary directories instead of shared names such as `temp_host.fna`;
- host CDS prediction is cached once per host instead of being rerun for every phage-host pair;
- pair tables contain paths and derived features rather than duplicating entire genomes inside JSON documents;
- subprocess failures are propagated instead of silently continuing with incomplete output;
- BLAST and Prodigal outputs are promoted to their final paths only after successful completion, so interrupted jobs are not mistaken for finished work.

## Requirements

Python 3.10+ and Biopython are required. BLAST+ is needed for database construction/search. Prodigal is optional and is used only if host codon-usage features are requested.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[dev]'
pytest -q
```

## 1. Build FASTA manifests

A manifest gives each genome a stable ID and records its FASTA path.

```bash
python scripts/build_manifest.py \
  --root data/host_genomes \
  --output work/hosts.tsv

python scripts/build_manifest.py \
  --root data/phage_genomes \
  --output work/phages.tsv
```

For downloaded assemblies, rename or construct the manifest IDs according to the assembly accession you want to track. The workflow does not assume a particular filesystem layout.

## 2. Index host sequence IDs

BLAST returns subject sequence identifiers, which are not necessarily the same as assembly identifiers. The index is built directly from every FASTA record header.

Taxonomy is deliberately not guessed from free-text FASTA descriptions. If species or taxonomic identifiers are needed, join them from a trusted assembly/metadata table using the genome ID.

```bash
python scripts/index_host_sequences.py \
  --manifest work/hosts.tsv \
  --output work/host_sequence_index.tsv
```

## 3. Build the host BLAST database

```bash
python scripts/build_host_blastdb.py \
  --manifest work/hosts.tsv \
  --combined-fasta work/hosts.fna \
  --db-prefix work/blast/hosts
```

## 4. Run phage-vs-host BLAST

For a local/single-process run:

```bash
python scripts/run_blast.py \
  --phage-manifest work/phages.tsv \
  --db-prefix work/blast/hosts \
  --output-dir results/blast \
  --threads 4
```

For HPC use, `slurm/blast_array.sbatch` shows a generic Slurm-array pattern. Cluster-specific accounts, partitions and module names are intentionally omitted.

The BLAST output retains `qlen` and `slen` in addition to the standard alignment fields so candidate filtering can inspect coverage rather than relying only on raw aligned length.

## 5. Build candidate associations

The following reproduces the exploratory minimum-alignment-length rule while making it explicit:

```bash
python scripts/build_candidate_pairs.py \
  --blast-dir results/blast \
  --host-sequence-index work/host_sequence_index.tsv \
  --output results/candidate_pairs.tsv \
  --min-alignment-length 500 \
  --max-evalue 1e-5
```

Additional constraints can be added without rerunning BLAST, for example:

```bash
  --min-identity 95 \
  --min-query-coverage 0.1
```

Each output row represents one phage–host genome pair. Hits across multiple contigs from the same host assembly are aggregated into that pair while the contributing subject sequence IDs are retained. The row also preserves summary evidence such as hit count, maximum identity/alignment length, minimum E-value, best bit score and best query/subject coverage.

## 6. Add simple pair-level genome features

```bash
python scripts/compute_pair_features.py \
  --pairs results/candidate_pairs.tsv \
  --phage-manifest work/phages.tsv \
  --output results/candidate_pairs.features.tsv
```

This adds host/phage genome lengths, GC fractions and absolute GC difference. Full genome sequences are not copied into the result table.

The exploratory workflow also calculated host codon usage from Prodigal-predicted CDS. To keep that computation reproducible and avoid predicting the same host repeatedly, CDS files can first be cached once per host:

```bash
python scripts/predict_host_cds.py \
  --host-manifest work/hosts.tsv \
  --output-dir work/host_cds

python scripts/compute_pair_features.py \
  --pairs results/candidate_pairs.tsv \
  --phage-manifest work/phages.tsv \
  --host-cds-dir work/host_cds \
  --output results/candidate_pairs.features.tsv
```

Host codon frequencies are then included as a JSON column. They are retained as descriptive features rather than presented as evidence of host assignment by themselves.

## Tiny synthetic example

The repository contains tiny synthetic FASTAs and a synthetic BLAST table for testing the data-processing steps without downloading biological datasets.

```bash
python scripts/build_manifest.py --root examples/hosts --output /tmp/hosts.tsv
python scripts/build_manifest.py --root examples/phages --output /tmp/phages.tsv
python scripts/index_host_sequences.py --manifest /tmp/hosts.tsv --output /tmp/host_index.tsv
python scripts/build_candidate_pairs.py \
  --blast-dir examples/blast \
  --host-sequence-index /tmp/host_index.tsv \
  --output /tmp/candidate_pairs.tsv \
  --min-alignment-length 500
python scripts/compute_pair_features.py \
  --pairs /tmp/candidate_pairs.tsv \
  --phage-manifest /tmp/phages.tsv \
  --output /tmp/candidate_pairs.features.tsv
```

With the example data, the 300 bp beta-host hit is filtered out and the two alpha-host HSPs are retained as one candidate pair with `hit_count=2`.

## Data and reproducibility

No RefSeq genome collection, BLAST database, production result table or private infrastructure configuration is distributed here. The repository contains only tiny synthetic records used for tests and examples.

The default `500 bp` threshold is preserved from the original exploratory workflow for reproducibility; it is **not** presented as a validated biological decision rule. Appropriate thresholds depend on the biological question, reference database and validation strategy.

## License

No license has been assigned yet.
