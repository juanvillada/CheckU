#!/usr/bin/env python3
"""
CheckU: UNI56 marker detection and completeness reporting.

This utility scans input proteome FASTA files (or nucleotide assemblies converted
via Pyrodigal) against the UNI56 marker set using pyhmmer. It produces summary
tables capturing marker presence, completeness, and duplication counts per
genome, alongside detailed hit reports and JSON checkpoints for reproducibility.
"""

from __future__ import annotations

import csv
import gzip
import json
import logging
import os
import shutil
import sys
import shlex
from collections import Counter
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import typer

try:
    import pyhmmer
except ImportError as exc:  # pragma: no cover - import guard
    typer.secho(
        "Error: pyhmmer is required to run CheckU. "
        "Install it with 'pip install pyhmmer'.",
        fg=typer.colors.RED,
        err=True,
    )
    raise typer.Exit(1) from exc


app = typer.Typer(
    help="CheckU marker completeness profiling tool.",
    add_completion=False,
    context_settings={"allow_extra_args": True},
)

VERSION = "0.1.7"


def _data_roots() -> List[Path]:
    roots: List[Path] = []
    env_root = os.environ.get("CHECKU_DATA_DIR")
    if env_root:
        roots.append(Path(env_root))
    module_dir = Path(__file__).resolve().parent
    roots.append(module_dir / "data")
    roots.append(module_dir.parent / "data")

    unique: List[Path] = []
    seen: set[str] = set()
    for root in roots:
        key = str(root)
        if key not in seen:
            seen.add(key)
            unique.append(root)
    return unique


def _data_candidates(*parts: str) -> List[Path]:
    return [root.joinpath(*parts) for root in _data_roots()]


def resolve_data_path(*parts: str) -> Path:
    candidates = _data_candidates(*parts)
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[0]


def require_data_path(*parts: str) -> Path:
    candidates = _data_candidates(*parts)
    for candidate in candidates:
        if candidate.exists():
            return candidate
    locations = ", ".join(str(path) for path in candidates)
    raise FileNotFoundError(
        f"Data resource not found: {Path(*parts)}. Looked in: {locations}"
    )


DEFAULT_HMM_PATH = resolve_data_path("uni56.hmm")
CHECKPOINT_VERSION = 1

SUPPORTED_NUCLEOTIDE_SUFFIXES = {
    ".fa",
    ".fas",
    ".fasta",
    ".fna",
    ".fnn",
    ".ffn",
    ".fa.gz",
    ".fasta.gz",
    ".fna.gz",
}
SUPPORTED_PROTEIN_SUFFIXES = {
    ".faa",
    ".faa.gz",
    ".aa",
    ".aa.fasta",
    ".pep",
    ".pep.fasta",
    ".protein.faa",
}


def configure_logging(log_level: str, log_dir: Path) -> logging.Logger:
    """Initialise console and file logging."""
    log_dir = log_dir.expanduser().absolute()
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / "checku.log"

    numeric_level = getattr(logging, log_level.upper(), logging.INFO)
    formatter = logging.Formatter(
        fmt="%(asctime)s %(levelname)s %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    logger = logging.getLogger("CheckU")
    logger.setLevel(numeric_level)
    logger.handlers.clear()

    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)
    stream_handler.setLevel(numeric_level)

    file_handler = logging.FileHandler(log_file, encoding="utf-8")
    file_handler.setFormatter(formatter)
    file_handler.setLevel(logging.DEBUG)

    logger.addHandler(stream_handler)
    logger.addHandler(file_handler)
    logger.debug("Logging initialised at level %s", logging.getLevelName(numeric_level))
    return logger


def timestamp() -> str:
    """Return current UTC timestamp."""
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def derive_genome_id(path: Path) -> str:
    """Create a stable genome identifier from the input filename."""
    name = path.name
    if name.endswith(".gz"):
        name = name[:-3]
    stem = Path(name).stem
    return stem.replace(" ", "_")


def detect_sequence_type(fasta_path: Path) -> str:
    """Guess whether the FASTA file contains nucleotide or amino-acid sequences."""
    protein_letters = set("ABCDEFGHIKLMNPQRSTVWXYZ*")  # exclude J/O/U (rare)
    nucleotide_letters = set("ACGTRYKMSWBDHVN-.U")

    observed_protein = 0
    observed_nucleotide = 0

    with open(fasta_path, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if line.startswith(">"):
                continue
            seq = line.strip().upper()
            if not seq:
                continue
            observed_protein += sum(ch in protein_letters for ch in seq)
            observed_nucleotide += sum(ch in nucleotide_letters for ch in seq)
            if observed_protein + observed_nucleotide >= 2000:
                break

    if observed_protein == 0 and observed_nucleotide == 0:
        return "unknown"

    # If any residue outside nucleotide alphabet is seen, treat as protein.
    if observed_protein > observed_nucleotide:
        return "protein"

    return "nucleotide"


def parse_fasta(path: Path) -> Iterable[Tuple[str, str]]:
    """Yield (header, sequence) tuples from a FASTA file."""
    header: Optional[str] = None
    chunks: List[str] = []

    with open(path, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:].strip() or "sequence"
                chunks = []
            else:
                chunks.append(line.strip())
        if header is not None:
            yield header, "".join(chunks)


@dataclass
class MarkerSummary:
    """Summary of hits for a single marker."""

    hit_count: int
    unique_sequences: List[str]
    best_score: float
    best_evalue: float

    def to_dict(self) -> Dict[str, object]:
        return {
            "hit_count": self.hit_count,
            "unique_sequences": self.unique_sequences,
            "best_score": self.best_score,
            "best_evalue": self.best_evalue,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, object]) -> "MarkerSummary":
        return cls(
            hit_count=int(data.get("hit_count", 0)),
            unique_sequences=list(data.get("unique_sequences", [])),
            best_score=float(data.get("best_score", 0.0)),
            best_evalue=float(data.get("best_evalue", 1.0)),
        )


@dataclass
class GenomeResult:
    """Per-genome summary of marker detection."""

    genome_id: str
    input_path: str
    input_type: str
    proteins_path: str
    markers_total: int
    markers_detected: int
    completeness: float
    duplicated_markers: int
    contamination: float
    marker_summaries: Dict[str, MarkerSummary] = field(default_factory=dict)
    raw_hits_path: Optional[str] = None
    pyrodigal: Optional[Dict[str, object]] = None

    def to_dict(self) -> Dict[str, object]:
        return {
            "genome_id": self.genome_id,
            "input_path": self.input_path,
            "input_type": self.input_type,
            "proteins_path": self.proteins_path,
            "markers_total": self.markers_total,
            "markers_detected": self.markers_detected,
            "completeness": self.completeness,
            "duplicated_markers": self.duplicated_markers,
            "contamination": self.contamination,
            "marker_summaries": {
                key: summary.to_dict() for key, summary in self.marker_summaries.items()
            },
            "raw_hits_path": self.raw_hits_path,
            "pyrodigal": self.pyrodigal,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, object]) -> "GenomeResult":
        marker_data = data.get("marker_summaries", {}) or {}
        marker_summaries = {
            key: MarkerSummary.from_dict(value) for key, value in marker_data.items()
        }
        return cls(
            genome_id=str(data.get("genome_id")),
            input_path=str(data.get("input_path")),
            input_type=str(data.get("input_type", "unknown")),
            proteins_path=str(data.get("proteins_path")),
            markers_total=int(data.get("markers_total", 0)),
            markers_detected=int(data.get("markers_detected", 0)),
            completeness=float(data.get("completeness", 0.0)),
            duplicated_markers=int(data.get("duplicated_markers", 0)),
            contamination=float(data.get("contamination", 0.0)),
            marker_summaries=marker_summaries,
            raw_hits_path=data.get("raw_hits_path"),
            pyrodigal=data.get("pyrodigal"),
        )


class Checkpoint:
    """JSON checkpoint to support resumable runs."""

    def __init__(self, path: Path, resume: bool, logger: logging.Logger):
        self.path = path
        self.resume = resume
        self.logger = logger
        self.records: Dict[str, GenomeResult] = {}
        self.failures: List[Dict[str, object]] = []
        self.created_at = timestamp()
        self.updated_at = self.created_at

        if self.path.exists():
            if resume:
                self._load()
            else:
                backup = self.path.with_suffix(
                    f"{self.path.suffix}.bak-{datetime.utcnow().strftime('%Y%m%d%H%M%S')}"
                )
                shutil.copy2(self.path, backup)
                self.logger.info(
                    "Existing checkpoint moved to %s (resume disabled).", backup
                )
                self._write()  # reset file
        else:
            self._write()

    def _load(self) -> None:
        with open(self.path, "r", encoding="utf-8") as handle:
            data = json.load(handle)
        if data.get("version") != CHECKPOINT_VERSION:
            self.logger.warning(
                "Checkpoint version mismatch (found %s). Starting fresh.",
                data.get("version"),
            )
            self._write()
            return
        self.created_at = data.get("created_at", self.created_at)
        self.updated_at = data.get("updated_at", self.updated_at)
        self.records = {
            genome_id: GenomeResult.from_dict(record)
            for genome_id, record in (data.get("results", {}) or {}).items()
        }
        self.failures = list(data.get("failures", []) or [])
        self.logger.info(
            "Loaded checkpoint with %d completed genomes.", len(self.records)
        )

    def _write(self) -> None:
        payload = {
            "version": CHECKPOINT_VERSION,
            "created_at": self.created_at,
            "updated_at": self.updated_at,
            "results": {key: record.to_dict() for key, record in self.records.items()},
            "failures": self.failures,
        }
        self.path.parent.mkdir(parents=True, exist_ok=True)
        with open(self.path, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)

    def record_success(self, result: GenomeResult) -> None:
        self.records[result.genome_id] = result
        self.updated_at = timestamp()
        self._write()

    def record_failure(self, genome_id: str, input_path: Path, error: str) -> None:
        failure = {
            "genome_id": genome_id,
            "input_path": str(input_path),
            "error": error,
            "timestamp": timestamp(),
        }
        self.failures.append(failure)
        self.updated_at = timestamp()
        self._write()


class PyrodigalRunner:
    """Facade around Pyrodigal gene prediction."""

    def __init__(
        self,
        meta_mode: bool,
        translation_table: Optional[int],
        logger: logging.Logger,
    ):
        self.meta_mode = meta_mode
        self.translation_table = translation_table
        self.logger = logger

        try:
            import pyrodigal  # type: ignore
        except ImportError as exc:  # pragma: no cover - import guard
            raise RuntimeError(
                "pyrodigal is required for nucleotide FASTA inputs. "
                "Install it with 'pip install pyrodigal'."
            ) from exc

        self.pyrodigal = pyrodigal
        self.finder = pyrodigal.GeneFinder(meta=meta_mode, mask=True)

    def translate(self, fasta_path: Path, output_path: Path) -> Dict[str, object]:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        genes_total = 0
        contigs = 0

        with open(output_path, "w", encoding="utf-8") as protein_handle:
            for header, sequence in parse_fasta(fasta_path):
                contigs += 1
                if not sequence:
                    continue
                genes = self.finder.find_genes(sequence)
                if not genes:
                    continue
                genes_total += len(genes)
                genes.write_translations(
                    protein_handle,
                    sequence_id=header,
                    translation_table=self.translation_table,
                    include_stop=False,
                    full_id=True,
                )

        self.logger.info(
            "Pyrodigal predicted %d proteins across %d contigs (%s).",
            genes_total,
            contigs,
            fasta_path,
        )

        return {
            "meta_mode": self.meta_mode,
            "translation_table": self.translation_table,
            "contigs": contigs,
            "genes": genes_total,
            "proteins_path": str(output_path),
        }


class PyhmmerRunner:
    """Perform pyhmmer searches against UNI56 HMMs."""

    def __init__(self, hmm_path: Path, cpus: int, logger: logging.Logger):
        self.hmm_path = hmm_path
        self.cpus = max(1, cpus)
        self.logger = logger
        self.hmms, self.hmm_source_kind = self._load_hmms()
        self.marker_label = self._derive_marker_label()
        if not self.hmms:
            origin = (
                "directory" if self.hmm_source_kind == "directory" else "file"
            )
            raise ValueError(
                f"No HMM profiles detected in {origin} {self.hmm_path}. "
                "Provide GA-calibrated markers to continue."
            )
        self._validate_gathering_thresholds()
        self._validate_unique_names()
        self.marker_order = [hmm.name for hmm in self.hmms]
        self.marker_count = len(self.marker_order)
        self.logger.debug(
            "Loaded %d HMM profiles from %s %s.",
            len(self.hmms),
            "directory" if self.hmm_source_kind == "directory" else "file",
            self.hmm_path,
        )

    def _load_hmms(self) -> Tuple[List[pyhmmer.plan7.HMM], str]:
        if not self.hmm_path.exists():
            raise FileNotFoundError(f"HMM source not found: {self.hmm_path}")

        if self.hmm_path.is_dir():
            hmm_files = sorted(
                p
                for p in self.hmm_path.rglob("*")
                if p.is_file()
                and p.stat().st_size > 0
                and (
                    p.name.lower().endswith(".hmm")
                    or p.name.lower().endswith(".hmm.gz")
                )
            )
            if not hmm_files:
                raise ValueError(
                    "No HMM files found in directory "
                    f"{self.hmm_path}. Expected *.hmm or *.hmm.gz files."
                )
            hmms: List[pyhmmer.plan7.HMM] = []
            for hmm_file in hmm_files:
                with pyhmmer.plan7.HMMFile(str(hmm_file)) as handle:
                    hmms.extend(list(handle))
            return hmms, "directory"

        if self.hmm_path.is_file():
            with pyhmmer.plan7.HMMFile(str(self.hmm_path)) as hmm_file:
                hmms = list(hmm_file)
            return hmms, "file"

        raise ValueError(
            f"HMM source must be a file or directory, not {self.hmm_path}."
        )

    def _derive_marker_label(self) -> str:
        try:
            resolved = self.hmm_path.resolve()
            label = resolved.name
        except OSError:
            resolved = None
            label = ""
        if not label or label in {".", ""}:
            label = self.hmm_path.name
        if (not label or label in {".", ""}) and resolved is not None:
            label = resolved.parent.name if resolved.parent.name else str(resolved)
        if not label:
            label = str(self.hmm_path)
        return label

    def _validate_gathering_thresholds(self) -> None:
        missing = [
            hmm.name
            for hmm in self.hmms
            if not hmm.cutoffs.gathering_available()
        ]
        if missing:
            sample = ", ".join(missing[:5])
            if len(missing) > 5:
                sample += ", ..."
            raise ValueError(
                "All marker HMMs must define GA thresholds. "
                f"The following models are missing GA values: {sample}"
            )

    def _validate_unique_names(self) -> None:
        names = [hmm.name for hmm in self.hmms]
        duplicates = [name for name, count in Counter(names).items() if count > 1]
        if duplicates:
            sample = ", ".join(sorted(duplicates)[:5])
            if len(duplicates) > 5:
                sample += ", ..."
            raise ValueError(
                "Marker HMM names must be unique across the marker set. "
                f"Duplicate model names detected: {sample}"
            )

    def search(
        self,
        proteins_fasta: Path,
    ) -> Tuple[Dict[str, MarkerSummary], List[Dict[str, object]]]:
        with pyhmmer.easel.SequenceFile(str(proteins_fasta), digital=True) as seq_file:
            sequences = list(seq_file)

        if not sequences:
            self.logger.warning(
                "No sequences found in %s; skipping pyhmmer search.", proteins_fasta
            )
            return {}, []

        raw_hits: List[Dict[str, object]] = []
        marker_hits: Dict[str, MarkerSummary] = {}

        for hits in pyhmmer.hmmsearch(
            self.hmms,
            sequences,
            cpus=self.cpus,
            bit_cutoffs="gathering",
        ):
            marker = hits.query.name
            best_per_sequence: Dict[str, Dict[str, object]] = {}

            for hit in hits:
                seq_name = hit.name
                if not hit.domains:
                    continue

                best_domain = max(hit.domains, key=lambda d: d.score)
                alignment = best_domain.alignment
                record = {
                    "sequence": seq_name,
                    "full_score": float(hit.score),
                    "full_evalue": float(hit.evalue),
                    "full_bias": float(hit.bias),
                    "domain_score": float(best_domain.score),
                    "domain_evalue": float(best_domain.i_evalue),
                    "domain_bias": float(best_domain.bias),
                    "hmm_from": int(alignment.hmm_from) if alignment else None,
                    "hmm_to": int(alignment.hmm_to) if alignment else None,
                    "ali_from": int(alignment.target_from) if alignment else None,
                    "ali_to": int(alignment.target_to) if alignment else None,
                    "env_from": int(best_domain.env_from),
                    "env_to": int(best_domain.env_to),
                }

                previous = best_per_sequence.get(seq_name)
                if previous is None or record["domain_score"] > previous["domain_score"]:
                    best_per_sequence[seq_name] = record

            if best_per_sequence:
                scores = [entry["domain_score"] for entry in best_per_sequence.values()]
                evalues = [entry["domain_evalue"] for entry in best_per_sequence.values()]
                marker_hits[marker] = MarkerSummary(
                    hit_count=len(best_per_sequence),
                    unique_sequences=sorted(best_per_sequence.keys()),
                    best_score=float(max(scores)),
                    best_evalue=float(min(evalues)),
                )
                for record in best_per_sequence.values():
                    enriched = dict(record)
                    enriched["marker"] = marker
                    raw_hits.append(enriched)

        return marker_hits, raw_hits


@dataclass
class RunConfig:
    """Configuration for the CheckU pipeline."""

    input_path: Path
    output_dir: Path
    command_line: str
    hmm_path: Path
    cpus: int
    resume: bool
    fail_fast: bool
    meta_mode: bool
    translation_table: Optional[int]
    keep_intermediate: bool
    log_level: str


class CheckUPipeline:
    """End-to-end controller for CheckU marker detection."""

    def __init__(self, config: RunConfig):
        self.config = config
        self.logger = configure_logging(config.log_level, config.output_dir / "logs")
        self.logger.info("Starting CheckU pipeline (version %s).", VERSION)
        self.logger.info("Input path: %s", config.input_path)
        self.logger.info("Output directory: %s", config.output_dir)
        config.output_dir.mkdir(parents=True, exist_ok=True)
        if config.command_line:
            self.logger.info("CLI command: %s", config.command_line)

        checkpoint_path = config.output_dir / "checkpoint" / "checku_checkpoint.json"
        self.checkpoint = Checkpoint(checkpoint_path, config.resume, self.logger)
        self.work_dir = config.output_dir / "work"
        self.details_dir = config.output_dir / "details"
        self.details_dir.mkdir(parents=True, exist_ok=True)
        self.summary_path = config.output_dir / "checku_summary.tsv"

        self.pyhmmer_runner = PyhmmerRunner(config.hmm_path, config.cpus, self.logger)
        self.marker_order = self.pyhmmer_runner.marker_order
        self.logger.info(
            "Marker source (%s): %s (%s) (%d profiles with GA thresholds).",
            self.pyhmmer_runner.hmm_source_kind,
            self.pyhmmer_runner.marker_label,
            config.hmm_path,
            self.pyhmmer_runner.marker_count,
        )

        self.pyrodigal_runner: Optional[PyrodigalRunner] = None

    def _ensure_pyrodigal_runner_if_needed(self) -> None:
        if self.pyrodigal_runner is None:
            self.pyrodigal_runner = PyrodigalRunner(
                meta_mode=self.config.meta_mode,
                translation_table=self.config.translation_table,
                logger=self.logger,
            )

    def _materialize_input(self, path: Path, genome_id: str) -> Path:
        if path.suffix == ".gz":
            target = self.work_dir / genome_id / path.with_suffix("").name
            target.parent.mkdir(parents=True, exist_ok=True)
            if not target.exists():
                self.logger.debug("Decompressing %s to %s.", path, target)
                with gzip.open(path, "rt", encoding="utf-8", errors="ignore") as src:
                    with open(target, "w", encoding="utf-8") as dst:
                        shutil.copyfileobj(src, dst)
            return target

        return path

    def _collect_inputs(self) -> List[Path]:
        path = self.config.input_path
        if not path.exists():
            raise FileNotFoundError(f"Input path not found: {path}")

        if path.is_file():
            return [path]

        collected: List[Path] = []
        for suffixes in (SUPPORTED_PROTEIN_SUFFIXES, SUPPORTED_NUCLEOTIDE_SUFFIXES):
            for suffix in suffixes:
                collected.extend(path.rglob(f"*{suffix}"))

        unique = sorted(set(collected))
        if not unique:
            raise FileNotFoundError(
                f"No FASTA files detected under {path}. "
                "Supported suffixes include .faa, .fa, .fna, and gzipped equivalents."
            )
        return unique

    def _write_raw_hits(self, genome_id: str, raw_hits: List[Dict[str, object]]) -> Path:
        hits_dir = self.details_dir / "hits"
        hits_dir.mkdir(parents=True, exist_ok=True)
        output_path = hits_dir / f"{genome_id}_checku_hits.tsv"
        fieldnames = [
            "marker",
            "sequence",
            "full_score",
            "full_evalue",
            "full_bias",
            "domain_score",
            "domain_evalue",
            "domain_bias",
            "hmm_from",
            "hmm_to",
            "ali_from",
            "ali_to",
            "env_from",
            "env_to",
        ]
        with open(output_path, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            for record in raw_hits:
                writer.writerow(record)
        return output_path

    def _build_presence_table(self) -> List[Dict[str, object]]:
        rows: List[Dict[str, object]] = []
        marker_label = self.pyhmmer_runner.marker_label
        for genome_id, result in sorted(self.checkpoint.records.items()):
            row: Dict[str, object] = {
                "genome_id": genome_id,
                "markers_name": marker_label,
                "markers_total": result.markers_total,
                "markers_detected": result.markers_detected,
                "completeness": round(result.completeness, 2),
                "duplicated_markers": result.duplicated_markers,
                "contamination": round(result.contamination, 2),
            }
            for marker in self.marker_order:
                summary = result.marker_summaries.get(marker)
                row[marker] = 1 if summary and summary.hit_count > 0 else 0
            rows.append(row)
        return rows

    def _export_presence_table(self) -> None:
        presence_rows = self._build_presence_table()
        if not presence_rows:
            self.logger.warning(
                "No successful genomes to report; skipping presence table export."
            )
            return

        presence_path = self.details_dir / "checku_presence.tsv"
        fieldnames = ["genome_id", "markers_name"] + self.marker_order + [
            "markers_total",
            "markers_detected",
            "completeness",
            "duplicated_markers",
            "contamination",
        ]

        with open(presence_path, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            for row in presence_rows:
                writer.writerow(row)

        with open(self.summary_path, "w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(
                [
                    "genome_id",
                    "input_path",
                    "input_type",
                    "proteins_path",
                    "markers_name",
                    "markers_total",
                    "markers_detected",
                    "completeness",
                    "duplicated_markers",
                    "contamination",
                    "pyrodigal_genes",
                    "pyrodigal_contigs",
                ]
            )
            marker_label = self.pyhmmer_runner.marker_label
            for genome_id, result in sorted(self.checkpoint.records.items()):
                pyrodigal_info = result.pyrodigal or {}
                writer.writerow(
                    [
                        genome_id,
                        result.input_path,
                        result.input_type,
                        result.proteins_path,
                        marker_label,
                        result.markers_total,
                        result.markers_detected,
                        round(result.completeness, 2),
                        result.duplicated_markers,
                        round(result.contamination, 2),
                        pyrodigal_info.get("genes", ""),
                        pyrodigal_info.get("contigs", ""),
                    ]
                )

        self.logger.info(
            "Presence matrix written to %s and summary to %s.",
            presence_path,
            self.summary_path,
        )

    def _process_file(self, fasta_path: Path) -> None:
        genome_id = derive_genome_id(fasta_path)
        if genome_id in self.checkpoint.records:
            self.logger.info(
                "Skipping %s (already processed; use --no-resume to overwrite).",
                fasta_path,
            )
            return

        self.logger.info("Processing %s (genome id: %s).", fasta_path, genome_id)

        try:
            materialized = self._materialize_input(fasta_path, genome_id)
            sequence_type = detect_sequence_type(materialized)
            self.logger.debug(
                "Detected %s sequences for %s.", sequence_type, materialized
            )

            if sequence_type == "unknown":
                size = materialized.stat().st_size if materialized.exists() else None
                if size == 0:
                    raise ValueError(
                        f"{fasta_path} appears empty; expected FASTA sequences."
                    )
                raise ValueError(
                    f"Unable to determine sequence type for {fasta_path}. "
                    "Ensure the file contains nucleotide or amino-acid FASTA records."
                )

            proteins_path = materialized
            pyrodigal_info: Optional[Dict[str, object]] = None

            if sequence_type != "protein":
                if self.pyrodigal_runner is None:
                    try:
                        self._ensure_pyrodigal_runner_if_needed()
                    except RuntimeError as exc:
                        raise RuntimeError(
                            "Pyrodigal is required to process nucleotide FASTA inputs."
                        ) from exc
                proteins_path = (
                    self.work_dir / genome_id / f"{genome_id}_predicted_proteins.faa"
                )
                pyrodigal_info = self.pyrodigal_runner.translate(
                    materialized, proteins_path
                )
                if pyrodigal_info["genes"] == 0:
                    self.logger.warning(
                        "Pyrodigal predicted zero proteins for %s.", fasta_path
                    )

            marker_hits, raw_hits = self.pyhmmer_runner.search(proteins_path)
            marker_total = self.pyhmmer_runner.marker_count
            markers_detected = len(
                [summary for summary in marker_hits.values() if summary.hit_count > 0]
            )
            duplicated_markers = sum(
                1 for summary in marker_hits.values() if summary.hit_count > 1
            )
            completeness = (
                (markers_detected / marker_total) * 100
                if self.marker_order
                else 0.0
            )
            contamination = (
                (duplicated_markers / marker_total) * 100
                if self.marker_order
                else 0.0
            )

            raw_hits_path: Optional[str] = None
            if raw_hits:
                hits_path = self._write_raw_hits(genome_id, raw_hits)
                raw_hits_path = str(hits_path)

            result = GenomeResult(
                genome_id=genome_id,
                input_path=str(fasta_path),
                input_type=sequence_type,
                proteins_path=str(proteins_path),
                markers_total=marker_total,
                markers_detected=markers_detected,
                completeness=round(completeness, 2),
                duplicated_markers=duplicated_markers,
                contamination=round(contamination, 2),
                marker_summaries=marker_hits,
                raw_hits_path=raw_hits_path,
                pyrodigal=pyrodigal_info,
            )

            self.checkpoint.record_success(result)
            self.logger.info(
                "%s: %d/%d markers detected (%.2f%% completeness, %d duplicates, %.2f%% contamination).",
                genome_id,
                markers_detected,
                marker_total,
                result.completeness,
                duplicated_markers,
                result.contamination,
            )

            if not self.config.keep_intermediate and proteins_path != materialized:
                try:
                    proteins_path.unlink(missing_ok=True)
                except OSError as exc:
                    self.logger.warning(
                        "Failed to remove temporary protein file %s: %s",
                        proteins_path,
                        exc,
                    )

        except Exception as exc:
            self.logger.exception("Processing failed for %s: %s", fasta_path, exc)
            self.checkpoint.record_failure(genome_id, fasta_path, str(exc))
            if self.config.fail_fast:
                raise

    def run(self) -> None:
        self.logger.debug("Collecting input FASTA files.")
        inputs = self._collect_inputs()
        self.logger.info("Found %d FASTA files to evaluate.", len(inputs))

        for fasta_path in inputs:
            self._process_file(fasta_path)

        self._export_presence_table()

        if self.checkpoint.failures:
            self.logger.warning(
                "Completed with %d failures. Checkpoint: %s",
                len(self.checkpoint.failures),
                self.checkpoint.path,
            )
        else:
            self.logger.info("Pipeline finished without errors.")


def version_callback(value: bool) -> None:
    if value:
        typer.echo(f"checku version {VERSION}")
        raise typer.Exit()


def _run_pipeline(
    *,
    input_path: Path,
    output_dir: Path,
    hmm_path: Path,
    cpus: int,
    resume: bool,
    fail_fast: bool,
    meta_mode: bool,
    translation_table: Optional[int],
    keep_intermediate: bool,
    log_level: str,
    command_line: Optional[str] = None,
) -> None:
    if command_line is None:
        command_line = " ".join(shlex.quote(arg) for arg in sys.argv)
    input_path = input_path.expanduser().absolute()
    output_dir = output_dir.expanduser().absolute()
    hmm_path = hmm_path.expanduser().absolute()
    config = RunConfig(
        input_path=input_path,
        output_dir=output_dir,
        command_line=command_line,
        hmm_path=hmm_path,
        cpus=cpus,
        resume=resume,
        fail_fast=fail_fast,
        meta_mode=meta_mode,
        translation_table=translation_table,
        keep_intermediate=keep_intermediate,
        log_level=log_level,
    )

    pipeline = CheckUPipeline(config)
    pipeline.run()


@app.callback(invoke_without_command=True)
def main(
    ctx: typer.Context,
    output_dir: Optional[Path] = typer.Option(
        None,
        "--output-dir",
        "-o",
        help="Directory for all outputs. Defaults to output_CheckU_<DATETIME>.",
    ),
    hmm_path: Path = typer.Option(
        DEFAULT_HMM_PATH,
        "--hmm",
        help=(
            "Path to the UNI56 marker file or a directory of GA-calibrated HMMs "
            "(defaults to bundled copy)."
        ),
    ),
    cpus: int = typer.Option(
        1,
        "--cpus",
        "-c",
        min=1,
        help="Number of CPU cores to hand over to pyhmmer.",
    ),
    resume: bool = typer.Option(
        True,
        "--resume/--no-resume",
        help="Resume from existing checkpoint when available.",
    ),
    fail_fast: bool = typer.Option(
        False,
        "--fail-fast/--no-fail-fast",
        help="Stop immediately when a genome fails instead of continuing.",
    ),
    meta_mode: bool = typer.Option(
        True,
        "--meta/--single",
        help="Configure Pyrodigal to run in metagenomic mode.",
    ),
    translation_table: Optional[int] = typer.Option(
        None,
        "--translation-table",
        "-t",
        help="Override translation table for Pyrodigal predictions.",
    ),
    keep_intermediate: bool = typer.Option(
        True,
        "--keep-intermediate/--clean-intermediate",
        help="Retain intermediate protein FASTA files generated by Pyrodigal.",
    ),
    log_level: str = typer.Option(
        "INFO",
        "--log-level",
        "-l",
        help="Logging level (DEBUG, INFO, WARNING, ERROR).",
    ),
    version: bool = typer.Option(
        None,
        "--version",
        callback=version_callback,
        is_eager=True,
        help="Print version information and exit.",
    ),
) -> None:
    """Execute the CheckU pipeline."""
    if ctx.invoked_subcommand is not None:
        return
    if not ctx.args:
        raise typer.BadParameter("Missing INPUT_PATH.")
    if len(ctx.args) > 1:
        raise typer.BadParameter("Only one INPUT_PATH is supported.")
    input_path = Path(ctx.args[0])
    if output_dir is None:
        output_dir = Path(
            f"output_CheckU_{datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S')}"
        )
    _run_pipeline(
        input_path=input_path,
        output_dir=output_dir,
        hmm_path=hmm_path,
        cpus=cpus,
        resume=resume,
        fail_fast=fail_fast,
        meta_mode=meta_mode,
        translation_table=translation_table,
        keep_intermediate=keep_intermediate,
        log_level=log_level,
    )


@app.command()
def run(
    input_path: Path = typer.Argument(..., help="FASTA file or directory to process."),
    output_dir: Optional[Path] = typer.Option(
        None,
        "--output-dir",
        "-o",
        help="Directory for all outputs. Defaults to output_CheckU_<DATETIME>.",
    ),
    hmm_path: Path = typer.Option(
        DEFAULT_HMM_PATH,
        "--hmm",
        help=(
            "Path to the UNI56 marker file or a directory of GA-calibrated HMMs "
            "(defaults to bundled copy)."
        ),
    ),
    cpus: int = typer.Option(
        1,
        "--cpus",
        "-c",
        min=1,
        help="Number of CPU cores to hand over to pyhmmer.",
    ),
    resume: bool = typer.Option(
        True,
        "--resume/--no-resume",
        help="Resume from existing checkpoint when available.",
    ),
    fail_fast: bool = typer.Option(
        False,
        "--fail-fast/--no-fail-fast",
        help="Stop immediately when a genome fails instead of continuing.",
    ),
    meta_mode: bool = typer.Option(
        True,
        "--meta/--single",
        help="Configure Pyrodigal to run in metagenomic mode.",
    ),
    translation_table: Optional[int] = typer.Option(
        None,
        "--translation-table",
        "-t",
        help="Override translation table for Pyrodigal predictions.",
    ),
    keep_intermediate: bool = typer.Option(
        True,
        "--keep-intermediate/--clean-intermediate",
        help="Retain intermediate protein FASTA files generated by Pyrodigal.",
    ),
    log_level: str = typer.Option(
        "INFO",
        "--log-level",
        "-l",
        help="Logging level (DEBUG, INFO, WARNING, ERROR).",
    ),
) -> None:
    """Execute the CheckU pipeline."""
    if output_dir is None:
        output_dir = Path(
            f"output_CheckU_{datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S')}"
        )
    _run_pipeline(
        input_path=input_path,
        output_dir=output_dir,
        hmm_path=hmm_path,
        cpus=cpus,
        resume=resume,
        fail_fast=fail_fast,
        meta_mode=meta_mode,
        translation_table=translation_table,
        keep_intermediate=keep_intermediate,
        log_level=log_level,
    )


@app.command()
def test(
    output_dir: Optional[Path] = typer.Option(
        None,
        "--output-dir",
        "-o",
        help=(
            "Directory for test outputs (FAA/FNA subfolders will be created). "
            "Defaults to output_test_CheckU_<DATETIME>."
        ),
    ),
    hmm_path: Path = typer.Option(
        DEFAULT_HMM_PATH,
        "--hmm",
        help=(
            "Path to the UNI56 marker file or a directory of GA-calibrated HMMs "
            "(defaults to bundled copy)."
        ),
    ),
    cpus: int = typer.Option(
        1,
        "--cpus",
        "-c",
        min=1,
        help="Number of CPU cores to hand over to pyhmmer.",
    ),
    resume: bool = typer.Option(
        True,
        "--resume/--no-resume",
        help="Resume from existing checkpoint when available.",
    ),
    fail_fast: bool = typer.Option(
        False,
        "--fail-fast/--no-fail-fast",
        help="Stop immediately when a genome fails instead of continuing.",
    ),
    meta_mode: bool = typer.Option(
        True,
        "--meta/--single",
        help="Configure Pyrodigal to run in metagenomic mode.",
    ),
    translation_table: Optional[int] = typer.Option(
        None,
        "--translation-table",
        "-t",
        help="Override translation table for Pyrodigal predictions.",
    ),
    keep_intermediate: bool = typer.Option(
        True,
        "--keep-intermediate/--clean-intermediate",
        help="Retain intermediate protein FASTA files generated by Pyrodigal.",
    ),
    log_level: str = typer.Option(
        "INFO",
        "--log-level",
        "-l",
        help="Logging level (DEBUG, INFO, WARNING, ERROR).",
    ),
) -> None:
    """Run CheckU against bundled FAA and FNA test genomes."""
    test_root = require_data_path("test_genomes")
    if not test_root.is_dir():
        raise NotADirectoryError(f"Test dataset path is not a directory: {test_root}")
    if output_dir is None:
        output_dir = Path(
            f"output_test_CheckU_{datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S')}"
        )
    cases = [
        ("faa", test_root / "faa"),
        ("fna", test_root / "fna"),
    ]
    for label, case_path in cases:
        if not case_path.exists():
            raise FileNotFoundError(f"Test dataset not found: {case_path}")
        _run_pipeline(
            input_path=case_path,
            output_dir=output_dir / label,
            hmm_path=hmm_path,
            cpus=cpus,
            resume=resume,
            fail_fast=fail_fast,
            meta_mode=meta_mode,
            translation_table=translation_table,
            keep_intermediate=keep_intermediate,
            log_level=log_level,
            command_line=" ".join(shlex.quote(arg) for arg in sys.argv),
        )


if __name__ == "__main__":
    app()
