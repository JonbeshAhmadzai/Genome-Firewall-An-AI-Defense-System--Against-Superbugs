"""Resolve a configured pathogen/antibiotic pair to its pipeline backend."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re

from targets_config import pathogen_spec, target_spec


def slugify(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", value.lower()).strip("_")


@dataclass(frozen=True)
class PipelineRoute:
    species: str
    antibiotic: str
    kind: str
    reader: str
    amrfinder_organism: str | None
    model_path: Path
    molecular_targets: tuple[str, ...]
    annotation_markers: dict[str, tuple[str, ...]]
    enabled: bool

    @property
    def supported(self) -> bool:
        """Whether a local reader backend exists for this route."""

        return self.reader == "amrfinderplus" and self.kind == "bacterium" and bool(self.molecular_targets)

    def unsupported_reason(self) -> str:
        if self.supported:
            return ""
        if self.kind == "bacterium" and not self.molecular_targets:
            return f"No molecular target genes are configured for {self.species} / {self.antibiotic}."
        return f"No implemented {self.reader!r} reader for {self.kind} pathogen {self.species}."


def resolve_route(
    species: str,
    antibiotic: str,
    *,
    model_root: Path = Path("models/cohort30"),
) -> PipelineRoute:
    """Resolve all runtime choices from ``targets_config.py``."""

    pathogen = pathogen_spec(species)
    target = target_spec(species, antibiotic)
    return PipelineRoute(
        species=species,
        antibiotic=antibiotic,
        kind=str(pathogen.get("kind", "unknown")),
        reader=str(pathogen.get("reader", "unconfigured")),
        amrfinder_organism=pathogen.get("amrfinder_organism"),
        model_path=model_root / f"{slugify(species)}__{slugify(antibiotic)}.joblib",
        molecular_targets=tuple(target["molecular_targets"]),
        annotation_markers=dict(target["annotation_markers"]),
        enabled=bool(target.get("enabled")) and bool(pathogen.get("enabled")),
    )
