# CheckU

CheckU evaluates bacterial and archaeal genomes with the UNI56 universal single-copy marker set. The program reads amino acid FASTA files or nucleotide assemblies, calls genes with Pyrodigal when needed, and scores markers with PyHMMER. Results include completeness, contamination, and per-marker hit tables.

## Requirements

- Linux x86_64 with enough CPU and RAM for HMMER searches
- FASTA inputs in plain or gzip form (`.faa`, `.fa`, `.fna`, and friends)

## Installation

### Option 1: pip (PyPI)

```bash
pip install checku
```

### Option 2: Pixi (development)

```bash
pixi install
```

## Quick Check

```bash
checku --help
```

If you are running from the repository with Pixi:

```bash
pixi run python -m checku --help
```

You should see the command line help without errors.

## Input Rules

- Provide either a single FASTA file or a directory of FASTA files.
- Protein files are used as-is. Nucleotide files trigger Pyrodigal gene calls.
- Compressed files (`.gz`) are supported; they are unpacked into the run workspace.

## Running The Pipeline

If you are running from the repository with Pixi, replace `checku` below with `pixi run python -m checku`.

### Pipeline Overview

The diagram below shows the main stages executed by CheckU.

```mermaid
graph TD
    A([Start run]) --> B[Collect FASTA inputs from file or directory]
    B --> C[Materialize gzipped files under `work/` when needed]
    C --> D{Detect sequence type}
    D -->|Protein| E[Use supplied protein FASTA]
    D -->|Nucleotide| F[Predict proteins with Pyrodigal]
    F --> E
    E --> G[Search UNI56 HMMs with pyhmmer]
    G --> H[Aggregate marker hits and completeness statistics]
    H --> I[Write `checku_summary.tsv`]
    H --> J[Write `details/checku_presence.tsv`]
    H --> K[Write raw hit tables in `details/hits/`]
    H --> L[Update checkpoint data and logs]
    H -.-> M[Optional: delete predicted proteins when `--clean-intermediate`]
    I --> N([Pipeline complete])
    J --> N
    K --> N
    L --> N
    M --> N
```

### Single Proteome

```bash
checku run \
  data/test_genomes/faa/IMGI2140918011.faa \
  --output-dir tmp/proteome_example \
  --cpus 4
```

### Directory Of Proteomes

```bash
checku run \
  data/test_genomes/faa \
  --output-dir tmp/proteome_batch \
  --cpus 8
```

### Single Assembly

```bash
checku run \
  data/test_genomes/fna/IMG2140918011.fna \
  --output-dir tmp/assembly_example \
  --cpus 4 \
  --clean-intermediate
```

Use `--clean-intermediate` if you do not need the predicted protein FASTA after the run.

## Custom Marker Sets

- The default marker file ships with CheckU (UNI56).
- Point `--hmm` to a different GA-calibrated `.hmm` file or to a directory that holds `.hmm` or `.hmm.gz` profiles.
- Every profile must define GA cutoffs. The run stops early if a profile is missing them or if names are duplicated.

Example:

```bash
checku run \
  /path/to/genomes \
  --hmm /path/to/custom_markers.hmm \
  --output-dir tmp/custom_markers \
  --cpus 8
```

## Outputs

All outputs live in the chosen `--output-dir`.

- `checku_summary.tsv` — per-genome summary with completeness, contamination, duplicate counts, and Pyrodigal gene statistics.
- `details/checku_presence.tsv` — marker presence/absence matrix.
- `details/hits/*.tsv` — raw pyhmmer hits with domain scores.
- `checkpoint/checku_checkpoint.json` — resume data for interrupted runs.
- `logs/checku.log` — timestamps, command line, and status messages.

## Resume And Logging

- Runs resume automatically when `--resume` is left on (default).
- Use `--no-resume` to start fresh; the older checkpoint is copied aside.
- Increase `--log-level` to `DEBUG` when you need extra detail.

## Verification Step

Small test data sets are stored under `data/test_genomes/`. After installation you can confirm the pipeline by running:

```bash
checku run data/test_genomes/faa --output-dir tmp/test_run --cpus 2
```

The command should finish without errors and produce the summary and presence tables described above.
