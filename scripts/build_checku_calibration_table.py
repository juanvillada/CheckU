#!/usr/bin/env python3
"""Build the packaged CheckU-Cal residual table from the shredded benchmark."""

from __future__ import annotations

import argparse
import csv
import importlib.util
import sys
from collections import defaultdict
from pathlib import Path
from statistics import mean, median

CALIBRATION_MODULE_PATH = (
    Path(__file__).resolve().parents[1] / "checku" / "calibration.py"
)
CALIBRATION_SPEC = importlib.util.spec_from_file_location(
    "checku_calibration_module", CALIBRATION_MODULE_PATH
)
if CALIBRATION_SPEC is None or CALIBRATION_SPEC.loader is None:
    raise RuntimeError(f"Unable to load calibration helpers from {CALIBRATION_MODULE_PATH}")
calibration = importlib.util.module_from_spec(CALIBRATION_SPEC)
sys.modules[CALIBRATION_SPEC.name] = calibration
CALIBRATION_SPEC.loader.exec_module(calibration)

SUPPORTED_FINE_PHYLA = calibration.SUPPORTED_FINE_PHYLA
canonical_domain = calibration.canonical_domain
canonical_phylum = calibration.canonical_phylum
checku_bin = calibration.checku_bin
coarse_phylum_group = calibration.coarse_phylum_group

MIN_STRATUM_SIZE = 20


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build a software-facing CheckU calibration table with exact-phylum "
            "rows where benchmark support exists, plus coarse and global backoffs."
        )
    )
    parser.add_argument(
        "--benchmark-tsv",
        required=True,
        help="Path to labelled_benchmark_master.tsv.",
    )
    parser.add_argument(
        "--output",
        default=str(Path(__file__).resolve().parents[1] / "checku" / "data" / "checku_calibration.tsv"),
        help="Output TSV path.",
    )
    return parser.parse_args()


def as_float(raw: str | None) -> float | None:
    if raw is None:
        return None
    text = str(raw).strip()
    if not text or text.lower() == "nan":
        return None
    try:
        return float(text)
    except ValueError:
        return None


def row_taxonomy_labels(row: dict[str, str]) -> list[str]:
    labels: list[str] = []
    phylum = canonical_phylum(row.get("gtdbtk_phylum") or row.get("phylum"))
    if phylum and phylum in SUPPORTED_FINE_PHYLA:
        labels.append(phylum)

    coarse = coarse_phylum_group(phylum)
    if coarse and coarse not in labels:
        labels.append(coarse)

    if not labels:
        labels.append("Other")
    return labels


def fit_table(rows: list[dict[str, str]]) -> dict[tuple[str, str, str], list[float]]:
    residuals_by_key: dict[tuple[str, str, str], list[float]] = defaultdict(list)
    for row in rows:
        if row.get("record_type") != "shredded":
            continue
        raw_checku = as_float(row.get("checku_completeness_pct"))
        realized = as_float(row.get("realized_completeness_pct"))
        if raw_checku is None or realized is None:
            continue

        domain = canonical_domain(row.get("domain"))
        if domain not in {"Bacteria", "Archaea"}:
            continue

        residual = realized - raw_checku
        bin_label = checku_bin(raw_checku)

        for taxonomy_group in row_taxonomy_labels(row):
            residuals_by_key[(domain, taxonomy_group, bin_label)].append(residual)
        residuals_by_key[(domain, "all_taxa", bin_label)].append(residual)
        residuals_by_key[("all_domains", "all_taxa", bin_label)].append(residual)
        residuals_by_key[("all_domains", "all_taxa", "all_checku")].append(residual)
    return residuals_by_key


def write_table(output_path: Path, residuals_by_key: dict[tuple[str, str, str], list[float]]) -> tuple[int, int]:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    rows_written = 0
    exact_taxa: set[str] = set()
    with open(output_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "domain",
                "taxonomy_group",
                "checku_bin",
                "n_samples",
                "median_residual_pct",
                "mean_residual_pct",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for key in sorted(residuals_by_key):
            residuals = residuals_by_key[key]
            if len(residuals) < MIN_STRATUM_SIZE and key != ("all_domains", "all_taxa", "all_checku"):
                continue
            writer.writerow(
                {
                    "domain": key[0],
                    "taxonomy_group": key[1],
                    "checku_bin": key[2],
                    "n_samples": len(residuals),
                    "median_residual_pct": f"{median(residuals):.2f}",
                    "mean_residual_pct": f"{mean(residuals):.2f}",
                }
            )
            rows_written += 1
            if key[1] in SUPPORTED_FINE_PHYLA:
                exact_taxa.add(key[1])
    return rows_written, len(exact_taxa)


def main() -> None:
    args = parse_args()
    benchmark_path = Path(args.benchmark_tsv)
    output_path = Path(args.output)

    with open(benchmark_path, newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    residuals = fit_table(rows)
    rows_written, exact_taxa = write_table(output_path, residuals)
    print(
        f"Wrote {rows_written} calibration rows to {output_path} "
        f"covering {exact_taxa} exact phyla plus coarse/global backoffs."
    )


if __name__ == "__main__":
    main()
