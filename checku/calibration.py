"""CheckU calibration helpers."""

from __future__ import annotations

import csv
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional, Tuple


SUPPORTED_FINE_PHYLA = frozenset(
    {
        "Acidobacteriota",
        "Actinomycetota",
        "Aquificota",
        "Armatimonadota",
        "Bacillota",
        "Bacillota_A",
        "Bacillota_B",
        "Bacillota_C",
        "Bacillota_D",
        "Bacillota_E",
        "Bacillota_F",
        "Bacillota_G",
        "Bacillota_I",
        "Bacteroidota",
        "Bdellovibrionota",
        "Caldisericota",
        "Calditrichota",
        "Campylobacterota",
        "Campylobacterota_A",
        "Chlamydiota",
        "Chloroflexota",
        "Chrysiogenota",
        "Coprothermobacterota",
        "Cyanobacteriota",
        "Deferribacterota",
        "Deinococcota",
        "Dependentiae",
        "Desulfobacterota",
        "Dictyoglomota",
        "Elusimicrobiota",
        "Fibrobacterota",
        "Fusobacteriota",
        "Gemmatimonadota",
        "Halobacteriota",
        "Methanobacteriota",
        "Methanobacteriota_A",
        "Methanobacteriota_B",
        "Methylomirabilota",
        "Micrarchaeota",
        "Myxococcota",
        "Nanoarchaeota",
        "Nitrospinota",
        "Nitrospirota",
        "Nitrospirota_A",
        "Patescibacteria",
        "Planctomycetota",
        "Pseudomonadota",
        "Spirochaetota",
        "Synergistota",
        "Thermoplasmatota",
        "Thermoproteota",
        "Thermotogota",
        "Verrucomicrobiota",
    }
)

COARSE_PHYLA = frozenset(
    {
        "Pseudomonadota",
        "Actinomycetota",
        "Bacteroidota",
        "Bacillota",
        "Cyanobacteriota",
        "Desulfobacterota",
        "Halobacteriota",
        "Thermoproteota",
    }
)

ARCHAEAL_PHYLA = frozenset(
    {
        "Halobacteriota",
        "Methanobacteriota",
        "Methanobacteriota_A",
        "Methanobacteriota_B",
        "Micrarchaeota",
        "Nanoarchaeota",
        "Thermoplasmatota",
        "Thermoproteota",
    }
)

PHYLUM_ALIASES = {
    "Proteobacteria": "Pseudomonadota",
    "Actinobacteria": "Actinomycetota",
    "Bacteroidetes": "Bacteroidota",
    "Firmicutes": "Bacillota",
    "Firmicutes_A": "Bacillota_A",
    "Firmicutes_B": "Bacillota_B",
    "Firmicutes_C": "Bacillota_C",
    "Firmicutes_D": "Bacillota_D",
    "Firmicutes_E": "Bacillota_E",
    "Firmicutes_F": "Bacillota_F",
    "Firmicutes_G": "Bacillota_G",
    "Cyanobacteria": "Cyanobacteriota",
    "Verrucomicrobia": "Verrucomicrobiota",
    "Planctomycetes": "Planctomycetota",
    "Spirochaetes": "Spirochaetota",
    "Chloroflexi": "Chloroflexota",
    "Chlamydiae": "Chlamydiota",
    "Deferribacteres": "Deferribacterota",
    "Deinococcus-Thermus": "Deinococcota",
    "Fusobacteria": "Fusobacteriota",
    "Aquificae": "Aquificota",
}


def clamp_pct(value: float | None) -> float | None:
    if value is None:
        return None
    return max(0.0, min(100.0, value))


def checku_bin(value: float | None, step: int = 5) -> str:
    if value is None:
        return ""
    if value >= 100.0:
        return "[95,100]"
    start = int((value // step) * step)
    if start < 0:
        start = 0
    end = start + step
    return f"[{start},{end})"


def strip_rank_prefix(raw: str | None, prefix: str) -> str:
    if raw is None:
        return ""
    value = raw.strip()
    if value.startswith(prefix):
        return value[len(prefix) :]
    return value


def canonical_phylum(raw_phylum: str | None) -> str:
    clean = strip_rank_prefix(raw_phylum, "p__")
    if not clean:
        return ""
    return PHYLUM_ALIASES.get(clean, clean)


def normalize_taxonomy_label(raw_label: str | None) -> str:
    clean = strip_rank_prefix(raw_label, "p__")
    if not clean:
        return ""
    if clean in {"Other", "all_taxa", "all_domains"}:
        return clean
    return PHYLUM_ALIASES.get(clean, clean)


def infer_domain_from_phylum(raw_phylum: str | None) -> str:
    clean = canonical_phylum(raw_phylum)
    if not clean:
        return ""
    if clean in ARCHAEAL_PHYLA:
        return "Archaea"
    if clean in SUPPORTED_FINE_PHYLA or clean in COARSE_PHYLA:
        return "Bacteria"
    return ""


def canonical_domain(raw: str | None) -> str:
    if raw is None:
        return ""
    value = raw.strip()
    if value in {"d__Bacteria", "Bacteria"}:
        return "Bacteria"
    if value in {"d__Archaea", "Archaea"}:
        return "Archaea"
    return value


def gtdb_parts(classification: str | None) -> Dict[str, str]:
    parts = {
        "domain": "",
        "phylum": "",
        "class": "",
        "order": "",
        "family": "",
        "genus": "",
        "species": "",
    }
    if not classification:
        return parts
    tokens = classification.split(";")
    keys = ["domain", "phylum", "class", "order", "family", "genus", "species"]
    for key, token in zip(keys, tokens):
        parts[key] = token
    return parts


def phylum_group(raw_phylum: str | None) -> str:
    clean = canonical_phylum(raw_phylum)
    if not clean:
        return "Other"
    if clean in SUPPORTED_FINE_PHYLA:
        return clean
    if clean.startswith("Bacillota"):
        return "Bacillota"
    if clean in COARSE_PHYLA:
        return clean
    return "Other"


def coarse_phylum_group(raw_phylum: str | None) -> str:
    clean = canonical_phylum(raw_phylum)
    if not clean:
        return "Other"
    if clean.startswith("Bacillota"):
        return "Bacillota"
    if clean in COARSE_PHYLA:
        return clean
    return "Other"


def detect_delimiter(path: Path) -> str:
    with open(path, "r", encoding="utf-8") as handle:
        first_line = handle.readline()
    return "\t" if "\t" in first_line else ","


def first_nonempty(row: Dict[str, str], candidates: Tuple[str, ...]) -> str:
    for key in candidates:
        value = row.get(key, "")
        if value and value.strip():
            return value.strip()
    return ""


@dataclass(frozen=True)
class CalibrationChoice:
    domain: str
    taxonomy_group: str
    checku_bin: str
    n_samples: int
    median_residual_pct: float
    mean_residual_pct: float


class CalibrationTable:
    """Load and apply packaged CheckU residual calibration tables."""

    def __init__(self, path: Path, logger: logging.Logger):
        self.path = path
        self.logger = logger
        self.records: Dict[Tuple[str, str, str], CalibrationChoice] = {}
        self._load()

    def _load(self) -> None:
        if not self.path.exists():
            raise FileNotFoundError(f"Calibration table not found: {self.path}")

        delimiter = detect_delimiter(self.path)
        with open(self.path, "r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter=delimiter)
            for row in reader:
                if row.get("fold") and row["fold"] != "all":
                    continue
                domain = row.get("domain", "").strip()
                taxonomy_group = row.get("taxonomy_group", "").strip() or "Other"
                checku_group = row.get("checku_bin", "").strip()
                if not checku_group:
                    continue
                key = (domain, taxonomy_group, checku_group)
                self.records[key] = CalibrationChoice(
                    domain=domain,
                    taxonomy_group=taxonomy_group,
                    checku_bin=checku_group,
                    n_samples=int(float(row.get("n_samples", "0") or 0)),
                    median_residual_pct=float(row.get("median_residual_pct", "0") or 0.0),
                    mean_residual_pct=float(row.get("mean_residual_pct", "0") or 0.0),
                )

        self.logger.info(
            "Loaded %d calibration strata from %s.", len(self.records), self.path
        )

    @staticmethod
    def _taxonomy_candidates(metadata: Dict[str, str]) -> list[str]:
        candidates: list[str] = []
        seen: set[str] = set()

        raw_candidates = [
            metadata.get("taxonomy_group", ""),
            metadata.get("phylum", ""),
            gtdb_parts(metadata.get("classification", "")).get("phylum", ""),
        ]
        for raw in raw_candidates:
            label = normalize_taxonomy_label(raw)
            if not label:
                continue
            for candidate in (label, coarse_phylum_group(label)):
                if candidate and candidate not in seen:
                    seen.add(candidate)
                    candidates.append(candidate)
        return candidates

    @classmethod
    def _backoffs(
        cls, raw_checku: float, metadata: Dict[str, str]
    ) -> list[Tuple[str, str, str]]:
        domain = canonical_domain(metadata.get("domain"))
        if not domain:
            domain = infer_domain_from_phylum(
                metadata.get("phylum") or metadata.get("taxonomy_group")
            )
        checku_group = checku_bin(raw_checku)

        candidates: list[Tuple[str, str, str]] = []
        seen: set[Tuple[str, str, str]] = set()

        for taxonomy_group in cls._taxonomy_candidates(metadata):
            key = (domain, taxonomy_group, checku_group)
            if key not in seen:
                seen.add(key)
                candidates.append(key)

        for key in (
            (domain, "all_taxa", checku_group),
            ("all_domains", "all_taxa", checku_group),
            ("all_domains", "all_taxa", "all_checku"),
        ):
            if key not in seen:
                seen.add(key)
                candidates.append(key)

        return candidates

    def apply(
        self, raw_checku: float | None, metadata: Optional[Dict[str, str]] = None
    ) -> Tuple[float | None, Optional[Dict[str, object]]]:
        if raw_checku is None:
            return None, None

        metadata = metadata or {}

        chosen: Optional[CalibrationChoice] = None
        chosen_key: Optional[Tuple[str, str, str]] = None
        for candidate in self._backoffs(raw_checku, metadata):
            chosen = self.records.get(candidate)
            if chosen is not None:
                chosen_key = candidate
                break

        if chosen is None or chosen_key is None:
            return raw_checku, None

        corrected = clamp_pct(raw_checku + chosen.median_residual_pct)
        return corrected, {
            "domain": metadata.get("domain", ""),
            "taxonomy_group": metadata.get("taxonomy_group", ""),
            "key": chosen_key,
            "n_samples": chosen.n_samples,
            "median_residual_pct": chosen.median_residual_pct,
            "mean_residual_pct": chosen.mean_residual_pct,
        }


class CalibrationMetadata:
    """Optional genome-level metadata used to refine CheckU-Cal lookups."""

    GENOME_ID_COLUMNS = ("genome_id", "record_id", "taxon_oid", "genome", "name")
    DOMAIN_COLUMNS = ("domain", "gtdbtk_domain", "gtdb_domain")
    TAXONOMY_GROUP_COLUMNS = ("taxonomy_group",)
    PHYLUM_COLUMNS = ("gtdbtk_phylum", "gtdbtk__phylum", "gtdb_phylum", "phylum_gtdb", "phylum")
    CLASSIFICATION_COLUMNS = (
        "classification_gtdbtk",
        "classification",
        "gtdb_classification",
        "gtdb_taxonomy",
    )

    def __init__(self, path: Optional[Path], logger: logging.Logger):
        self.path = path
        self.logger = logger
        self.records: Dict[str, Dict[str, str]] = {}
        if path is not None:
            self._load(path)

    def _load(self, path: Path) -> None:
        if not path.exists():
            raise FileNotFoundError(f"Calibration metadata file not found: {path}")

        delimiter = detect_delimiter(path)
        with open(path, "r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter=delimiter)
            for row in reader:
                genome_id = first_nonempty(row, self.GENOME_ID_COLUMNS)
                if not genome_id:
                    continue

                classification = first_nonempty(row, self.CLASSIFICATION_COLUMNS)
                parts = gtdb_parts(classification)
                phylum = canonical_phylum(
                    first_nonempty(row, self.PHYLUM_COLUMNS) or parts.get("phylum", "")
                )
                domain = canonical_domain(
                    first_nonempty(row, self.DOMAIN_COLUMNS) or parts.get("domain", "")
                )
                if not domain:
                    domain = infer_domain_from_phylum(phylum)

                taxonomy_group = normalize_taxonomy_label(
                    first_nonempty(row, self.TAXONOMY_GROUP_COLUMNS)
                )
                if not taxonomy_group:
                    taxonomy_group = phylum_group(phylum)

                self.records[genome_id] = {
                    "domain": domain,
                    "taxonomy_group": taxonomy_group or "Other",
                    "phylum": phylum,
                    "classification": classification,
                }

        self.logger.info(
            "Loaded calibration metadata for %d genomes from %s.",
            len(self.records),
            path,
        )

    def lookup(self, genome_id: str) -> Dict[str, str]:
        return self.records.get(genome_id, {})
