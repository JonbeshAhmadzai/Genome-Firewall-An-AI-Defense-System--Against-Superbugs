"""Configuration for candidate Genome Firewall coverage.

Candidate pairs are audited before they can become validated product coverage.
Nothing in this file implies that a model is clinically or scientifically approved.
"""

from __future__ import annotations


CANDIDATE_TARGETS: tuple[tuple[str, str], ...] = (
    ("Escherichia coli", "ciprofloxacin"),
    ("Escherichia coli", "ampicillin"),
    ("Escherichia coli", "gentamicin"),
    ("Escherichia coli", "ceftriaxone"),
)

# Backwards-compatible name used by Blueprint 1. These are candidates until the
# data, grouped evaluation, calibration, and human-review gates all pass.
TARGETS = CANDIDATE_TARGETS

FUTURE_TARGETS: tuple[tuple[str, str], ...] = (
    ("Klebsiella pneumoniae", "meropenem"),
    ("Staphylococcus aureus", "oxacillin"),
    ("Pseudomonas aeruginosa", "meropenem"),
    ("Acinetobacter baumannii", "meropenem"),
)

TAXON_IDS = {
    "Escherichia coli": 562,
    "Klebsiella pneumoniae": 573,
    "Staphylococcus aureus": 1280,
    "Pseudomonas aeruginosa": 287,
    "Acinetobacter baumannii": 470,
}

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
}

# Curated target-gene hints used by the deterministic compatibility gate. The
# gate is conservative: a missing or uncertain target never permits a
# "likely to work" result. These names are matching hints for annotations, not
# phenotype labels and not a claim that any one gene fully explains activity.
MOLECULAR_TARGETS: dict[str, tuple[str, ...]] = {
    "ciprofloxacin": ("gyrA", "parC"),
    "ampicillin": ("ftsI", "mrdA", "pbpA", "dacB"),
    "gentamicin": ("rpsL", "rpsJ", "rrs"),
    "ceftriaxone": ("ftsI", "mrdA", "pbpA", "dacB"),
}

# General genome-annotation markers used to verify that the molecular target
# machinery is represented. These are compatibility markers, not resistance
# determinants and never override the phenotype model.
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
}

# Confidence below this level is abstained as no-call. This is a starting
# policy; final thresholds are selected on the calibration split only.
TRUST = {
    "min_confidence": 0.70,
    "reliability_bins": 10,
}
