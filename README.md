# CheckU

CheckU evaluates bacterial and archaeal genomes with the UNI56 universal single-copy marker set. The program reads amino acid FASTA files or nucleotide assemblies, calls genes with Pyrodigal when needed, and scores markers with PyHMMER. Results include raw completeness, calibrated completeness, contamination, and per-marker hit tables.

## Requirements

- FASTA inputs in plain or gzip form (`.faa`, `.fa`, `.fna`, and friends)

## Installation (Recommended)

Make sure you have [`Pixi`](https://pixi.prefix.dev/latest/) installed:

```bash
curl -fsSL https://pixi.sh/install.sh | sh
```

Install `CheckU` with `Pixi`:

```bash
pixi global install \
  -c conda-forge \
  -c bioconda \
  -c https://repo.prefix.dev/astrogenomics \
  checku
```

### Quick test

Small test data sets ship with `CheckU`. After installation you can confirm the pipeline by running:

```bash
checku test
```

See the **Expected Results** section below for the expected output tables.

### Alternative: pip (PyPI)

```bash
pip install checku
```

### Developer install (Pixi)

If you want to download the code and develop locally:

```bash
git clone https://github.com/juanvillada/checku
cd checku
pixi install
```

## Quick check

```bash
checku --help
```

If you are running from the repository with `Pixi`:

```bash
pixi run python -m checku --help
```

You should see the command line help without errors.

## Input rules

- Provide either a single FASTA file or a directory of FASTA files.
- Protein files are used as-is. Nucleotide files trigger Pyrodigal gene calls.
- Compressed files (`.gz`) are supported; they are unpacked into the run workspace.

## Running the pipeline

If you are running from the repository with `Pixi`, replace `checku` below with `pixi run python -m checku`.

The examples below use the bundled test data from a source checkout. Replace the
paths with your own FASTA inputs, or run `checku test` after installation.

### Pipeline overview

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

### Single proteome

```bash
checku run \
  checku/data/test_genomes/faa/IMGI2140918011.faa \
  --output-dir tmp/proteome_example \
  --cpus 4
```

### Directory of proteomes

```bash
checku run \
  checku/data/test_genomes/faa \
  --output-dir tmp/proteome_batch \
  --cpus 8
```

### Single assembly

```bash
checku run \
  checku/data/test_genomes/fna/IMG2140918011.fna \
  --output-dir tmp/assembly_example \
  --cpus 4 \
  --clean-intermediate
```

Use `--clean-intermediate` if you do not need the predicted protein FASTA after the run.

## CheckU-Cal output

Every `checku_summary.tsv` now contains both:

- `completeness` — the raw UNI56 marker fraction (`markers_detected / 56 * 100`)
- `completeness_calibrated` — a bundled residual-corrected estimate (`CheckU-Cal`)

The bundled calibration table is used automatically. If you have genome-level metadata available, you can refine the lookup with:

- `--calibration-metadata <metadata.tsv>` — optional TSV/CSV with `genome_id` plus either `taxonomy_group` or GTDB fields such as `classification_gtdbtk`, `gtdbtk_domain`, and `gtdbtk_phylum`
- `--calibration-table <table.tsv>` — override the bundled manuscript-derived calibration table

When GTDB-style phylum metadata are provided, CheckU-Cal now uses exact phylum-specific residual tables where the shredded benchmark supports them. If an exact phylum/bin combination is too sparse, CheckU falls back to a coarse phylum grouping and then to domain-wide/global residual tables. This keeps the software-facing calibration more specific than the manuscript display groups while retaining safe fallbacks for poorly represented taxa.

### Calibration metadata template

The `genome_id` column must match the CheckU genome identifier, which is derived from the input filename stem. For example, `bin_001.fna.gz` becomes `bin_001`.

Minimal example:

```tsv
genome_id	gtdbtk_domain	gtdbtk_phylum	classification_gtdbtk
bin_001	d__Bacteria	p__Planctomycetota	d__Bacteria;p__Planctomycetota;c__Planctomycetes;o__Planctomycetales;f__Planctomycetaceae;g__Planctomyces;s__
bin_002	d__Archaea	p__Thermoproteota	d__Archaea;p__Thermoproteota;c__Nitrososphaeria;o__Nitrososphaerales;f__Nitrososphaeraceae;g__Nitrososphaera;s__
bin_003	d__Bacteria	p__Bacillota_A	d__Bacteria;p__Bacillota_A;c__Clostridia;o__Eubacteriales;f__Lachnospiraceae;g__Roseburia;s__
```

If you already maintain a curated calibration stratum, you can also provide `taxonomy_group` directly:

```tsv
genome_id	domain	taxonomy_group
bin_001	Bacteria	Planctomycetota
bin_002	Archaea	Thermoproteota
bin_003	Bacteria	Bacillota_A
```

Accepted metadata columns:

- genome id: `genome_id`, `record_id`, `taxon_oid`, `genome`, `name`
- domain: `domain`, `gtdbtk_domain`, `gtdb_domain`
- taxonomy group: `taxonomy_group`
- phylum: `gtdbtk_phylum`, `gtdbtk__phylum`, `gtdb_phylum`, `phylum_gtdb`, `phylum`
- full GTDB classification: `classification_gtdbtk`, `classification`, `gtdb_classification`, `gtdb_taxonomy`

Resolution order:

1. Use `taxonomy_group` directly if provided.
2. Otherwise derive the exact phylum from the phylum or GTDB classification fields.
3. If that exact phylum is not represented in the calibration table for the relevant completeness bin, fall back to the coarser group.
4. If taxonomy is absent or unsupported, fall back to domain-wide and then global calibration.

Example run with metadata:

```bash
checku run \
  /path/to/mags \
  --calibration-metadata metadata.tsv \
  --output-dir tmp/checku_with_calibration \
  --cpus 8
```

The summary and presence tables also report the calibration provenance:

- `calibration_domain`
- `calibration_taxonomy_group`
- `calibration_checku_bin`
- `calibration_n_train`

## Custom marker sets

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

- `checku_summary.tsv` — per-genome summary with raw completeness, calibrated completeness, contamination, duplicate counts, calibration provenance, and Pyrodigal gene statistics.
- `details/checku_presence.tsv` — marker presence/absence matrix.
- `details/hits/*.tsv` — raw pyhmmer hits with domain scores.
- `checkpoint/checku_checkpoint.json` — resume data for interrupted runs.
- `logs/checku.log` — timestamps, command line, and status messages.
- Output tables and logs record input/output locations using absolute paths for reproducibility.

## Resume and logging

- Runs resume automatically when `--resume` is left on (default).
- Use `--no-resume` to start fresh; the older checkpoint is copied aside.
- Increase `--log-level` to `DEBUG` when you need extra detail.

## Verification step

Small test data sets ship with `CheckU`. After installation you can confirm the pipeline by running:

```bash
checku test
```

The command should finish without errors and produce the summary and presence tables described above.

If you are running from the repository with `Pixi`:
```bash
pixi run python -m checku test
```

### Expected results (Bundled test data)

The tables below summarize the expected `checku_summary.tsv` values for the bundled FAA and FNA test sets.
Absolute paths (input/protein columns in the real table) are omitted for privacy.

The real output table also includes `calibration_domain`, `calibration_taxonomy_group`,
`calibration_checku_bin`, and `calibration_n_train`.

FAA (protein inputs):

| genome_id | markers_detected | completeness | completeness_calibrated | duplicated_markers | contamination |
| --- | --- | --- | --- | --- | --- |
| IMGI2140918011 | 55 | 98.21 | 89.20 | 0 | 0.00 |
| IMGI2645727657 | 56 | 100.00 | 89.09 | 0 | 0.00 |
| IMGI651324087 | 56 | 100.00 | 89.09 | 0 | 0.00 |
| IMGM3300027739_BIN74 | 36 | 64.29 | 57.07 | 0 | 0.00 |
| SCISO2808607008 | 55 | 98.21 | 89.20 | 1 | 1.79 |
| SDISOGCA_003484685.1 | 47 | 83.93 | 70.34 | 1 | 1.79 |
| SHISO2654587767 | 55 | 98.21 | 89.20 | 1 | 1.79 |
| SLISOGCF_900639865.1 | 56 | 100.00 | 89.09 | 1 | 1.79 |
| SRISO640427127 | 52 | 92.86 | 81.67 | 0 | 0.00 |
| SXGCA_000019745.1 | 55 | 98.21 | 89.20 | 0 | 0.00 |
| SXGCA_902860225.1_Azoamicus_ciliaticola | 51 | 91.07 | 79.88 | 0 | 0.00 |
| SXISO642555114 | 54 | 96.43 | 87.42 | 1 | 1.79 |

FNA (nucleotide inputs with Pyrodigal):

| genome_id | markers_detected | completeness | completeness_calibrated | duplicated_markers | contamination | pyrodigal_genes | pyrodigal_contigs |
| --- | --- | --- | --- | --- | --- | --- | --- |
| IMG2140918011 | 56 | 100.00 | 89.09 | 0 | 0.00 | 2974 | 78 |
| IMG2645727657 | 56 | 100.00 | 89.09 | 0 | 0.00 | 1516 | 1 |
| IMG2645727657_HALF | 46 | 82.14 | 68.55 | 0 | 0.00 | 821 | 1 |
| IMG651324087 | 56 | 100.00 | 89.09 | 0 | 0.00 | 2572 | 73 |
