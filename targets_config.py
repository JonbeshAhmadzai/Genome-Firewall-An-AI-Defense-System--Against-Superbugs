"""Single configuration registry for Genome Firewall coverage.

Adding a bacterial species or antibiotic starts here. Configuration makes a
target discoverable; the pipeline still requires sufficient labelled data,
grouped evaluation, calibration, and human review before it is trusted.
"""

from __future__ import annotations

from typing import Any


# ``reader`` is the feature-extraction backend. AMRFinderPlus is implemented
# today; ``virus_placeholder`` reserves a clear extension point without
# pretending that viral support already exists.
PATHOGEN_CONFIG: dict[str, dict[str, Any]] = {
    "Escherichia coli": {
        "kind": "bacterium",
        "taxon_id": 562,
        "reader": "amrfinderplus",
        "amrfinder_organism": "Escherichia",
        "enabled": True,
    },
    "Klebsiella pneumoniae": {
        "kind": "bacterium",
        "taxon_id": 573,
        "reader": "amrfinderplus",
        "amrfinder_organism": "Klebsiella",
        "enabled": False,
    },
    "Staphylococcus aureus": {
        "kind": "bacterium",
        "taxon_id": 1280,
        "reader": "amrfinderplus",
        "amrfinder_organism": "Staphylococcus",
        "enabled": False,
    },
    "Pseudomonas aeruginosa": {
        "kind": "bacterium",
        "taxon_id": 287,
        "reader": "amrfinderplus",
        "amrfinder_organism": "Pseudomonas_aeruginosa",
        "enabled": False,
    },
    "Acinetobacter baumannii": {
        "kind": "bacterium",
        "taxon_id": 470,
        "reader": "amrfinderplus",
        "amrfinder_organism": "Acinetobacter_baumannii",
        "enabled": False,
    },
    # Register a virus here once a validated virus-specific reader/database
    # exists. It is intentionally disabled and cannot enter the AMRFinder path.
    "Influenza A virus": {
        "kind": "virus",
        "taxon_id": 11320,
        "reader": "virus_placeholder",
        "enabled": False,
    },
}

# BV-BRC IDs remain available to existing fetchers and are derived from the
# registry so a new organism only needs one configuration entry.
TAXON_IDS = {species: int(spec["taxon_id"]) for species, spec in PATHOGEN_CONFIG.items()}

# This is only the preliminary phenotype-data gate. Genetic-group counts and
# held-out generalization are additional gates later in the pipeline.
GATE = {
    "min_total": 200,
    "min_per_class": 50,
    "min_minority_frac": 0.15,
}

ANTIBIOTIC_ALIASES = {
    "ciprofloxacin": "ciprofloxacin",
    "cipro": "ciprofloxacin",
    "ampicillin": "ampicillin",
    "amp": "ampicillin",
    "gentamicin": "gentamicin",
    "gentamycin": "gentamicin",
    "ceftriaxone": "ceftriaxone",
    "cro": "ceftriaxone",
    "meropenem": "meropenem",
    "mero": "meropenem",
    "oxacillin": "oxacillin",
    "oxa": "oxacillin",
    "oseltamivir": "oseltamivir",
    "tamiflu": "oseltamivir",
}

# Curated target-gene hints used by the deterministic compatibility gate. The
# gate is conservative: a missing or uncertain target never permits a
# "likely to work" result.
MOLECULAR_TARGETS: dict[str, tuple[str, ...]] = {
    "ciprofloxacin": ("gyrA", "parC"),
    "ampicillin": ("ftsI", "mrdA", "pbpA", "dacB"),
    "gentamicin": ("rpsL", "rpsJ", "rrs"),
    "ceftriaxone": ("ftsI", "mrdA", "pbpA", "dacB"),
    "meropenem": ("ftsI", "mrdA", "pbpA"),
    "oxacillin": ("mecA", "pbp2a"),
}

# General genome-annotation markers used to verify that target machinery is
# represented. These are compatibility markers, not resistance determinants.
TARGET_ANNOTATION_MARKERS: dict[str, dict[str, tuple[str, ...]]] = {
    "ciprofloxacin": {
        "genes": ("gyrA", "gyrB", "parC", "parE"),
        "product_keywords": ("DNA gyrase subunit", "DNA topoisomerase IV subunit"),
    },
    "ampicillin": {
        "genes": ("ftsI", "mrdA", "mrcA", "mrcB", "dacA", "dacB", "pbpA"),
        "product_keywords": ("penicillin-binding protein", "peptidoglycan transpeptidase"),
    },
    "ceftriaxone": {
        "genes": ("ftsI", "mrdA", "mrcA", "mrcB", "dacA", "dacB", "pbpA"),
        "product_keywords": ("penicillin-binding protein", "peptidoglycan transpeptidase"),
    },
    "gentamicin": {
        "genes": ("rrs", "rrsA", "rrsB", "rrsC", "rpsL"),
        "product_keywords": ("16S ribosomal RNA", "30S ribosomal protein S12"),
    },
    "meropenem": {
        "genes": ("ftsI", "mrdA", "pbpA"),
        "product_keywords": ("penicillin-binding protein",),
    },
    "oxacillin": {
        "genes": ("mecA", "mecC", "pbp2a"),
        "product_keywords": ("penicillin-binding protein 2a",),
    },
}

# One registry entry controls whether the pair is collected/selected now.
# Gentamicin and future species are defined but disabled until their data pass
# the gates. To add coverage, add one entry here and the relevant target hints.
TARGET_DEFINITIONS: dict[tuple[str, str], dict[str, Any]] = {
    ("Escherichia coli", "ciprofloxacin"): {"enabled": True},
    ("Escherichia coli", "ampicillin"): {"enabled": True},
    ("Escherichia coli", "ceftriaxone"): {"enabled": True},
    ("Escherichia coli", "gentamicin"): {"enabled": False},
    ("Klebsiella pneumoniae", "meropenem"): {"enabled": False},
    ("Staphylococcus aureus", "oxacillin"): {"enabled": False},
    ("Pseudomonas aeruginosa", "meropenem"): {"enabled": False},
    ("Acinetobacter baumannii", "meropenem"): {"enabled": False},
    ("Influenza A virus", "oseltamivir"): {"enabled": False},
}

CANDIDATE_TARGETS: tuple[tuple[str, str], ...] = tuple(TARGET_DEFINITIONS)
ACTIVE_TARGETS: tuple[tuple[str, str], ...] = tuple(
    pair
    for pair, spec in TARGET_DEFINITIONS.items()
    if bool(spec.get("enabled"))
    and bool(PATHOGEN_CONFIG.get(pair[0], {}).get("enabled"))
    and PATHOGEN_CONFIG.get(pair[0], {}).get("kind") == "bacterium"
)

# Backwards-compatible names used by existing fetchers and Blueprint 1.
TARGETS = CANDIDATE_TARGETS
FUTURE_TARGETS = tuple(pair for pair in CANDIDATE_TARGETS if pair not in ACTIVE_TARGETS)


def pathogen_spec(species: str) -> dict[str, Any]:
    """Return a configured pathogen or raise an actionable error."""

    try:
        return PATHOGEN_CONFIG[species]
    except KeyError as error:
        raise ValueError(f"Species is not configured: {species}. Add it to PATHOGEN_CONFIG.") from error


def target_spec(species: str, antibiotic: str) -> dict[str, Any]:
    """Return pair configuration, including species-independent target hints."""

    if (species, antibiotic) not in TARGET_DEFINITIONS:
        raise ValueError(f"Target pair is not configured: {species} / {antibiotic}.")
    return {
        **TARGET_DEFINITIONS[(species, antibiotic)],
        "species": species,
        "antibiotic": antibiotic,
        "molecular_targets": MOLECULAR_TARGETS.get(antibiotic, ()),
        "annotation_markers": TARGET_ANNOTATION_MARKERS.get(antibiotic, {"genes": (), "product_keywords": ()}),
    }


def enabled_species(*, kind: str | None = None) -> tuple[str, ...]:
    return tuple(
        species
        for species, spec in PATHOGEN_CONFIG.items()
        if spec.get("enabled") and (kind is None or spec.get("kind") == kind)
    )


def target_pairs_for_species(species: str) -> tuple[tuple[str, str], ...]:
    """Return all configured antibiotic pairs for one species."""

    return tuple(pair for pair in TARGET_DEFINITIONS if pair[0] == species)


def validate_config() -> list[str]:
    """Return actionable configuration errors before an expensive pipeline run."""

    errors: list[str] = []
    for species, spec in PATHOGEN_CONFIG.items():
        for field in ("kind", "reader", "taxon_id"):
            if field not in spec:
                errors.append(f"PATHOGEN_CONFIG[{species!r}] is missing {field!r}")
    for species, antibiotic in TARGET_DEFINITIONS:
        if species not in PATHOGEN_CONFIG:
            errors.append(f"Target pair {species}/{antibiotic} references an unknown species")
        if antibiotic not in ANTIBIOTIC_ALIASES:
            errors.append(f"Target pair {species}/{antibiotic} has no antibiotic alias")
        if TARGET_DEFINITIONS[(species, antibiotic)].get("enabled") and not MOLECULAR_TARGETS.get(antibiotic):
            errors.append(f"Enabled target {species}/{antibiotic} has no molecular target genes")
    return errors


# Confidence below this level is abstained as no-call.
TRUST = {
    "min_confidence": 0.70,
    "reliability_bins": 10,
}
