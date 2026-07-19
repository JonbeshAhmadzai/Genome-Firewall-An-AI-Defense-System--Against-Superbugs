"""Build an auditable, lower-dimensional AMRFinderPlus training matrix.

The raw AMRFinderPlus matrix has one binary column for every detected gene or
mutation allele.  That is useful evidence for a report, but it is a poor
training representation when there are many rare or co-occurring markers and
few labelled genomes.  This script creates a *new* feature matrix for model
training while deliberately leaving the raw matrix and evidence tables intact.

The script is intentionally implemented with the Python standard library.  It
can therefore run in a minimal environment and makes its transformations easy
to audit.  It reads the repository's existing CSV format; no AMRFinderPlus
database, network access, or model fitting is required.

Default behaviour
-----------------
* Use genomes shared by the feature and label tables to decide whether a
  feature varies in the training cohort.
* Ignore source markers that are constant on those training genomes.
* Do not emit one-off source-marker columns.  A rare marker may still
  contribute to a broader mechanism/family feature if several rare variants
  together occur often enough.
* Keep only AMR classes relevant to at least one configured project drug:
  beta-lactam, aminoglycoside, phenicol, efflux, and multidrug.
* Replace individual marker columns with binary, interpretable aggregates:
  AMR class, the AMRFinder class/subclass label, and conservative
  gene/mutation family.
* Add compact numeric burden counts.  Two are restricted to the configured
  drug classes; two more span every AMR class (including off-target ones such
  as quinolone or sulfonamide) to capture the multidrug-resistant-lineage
  co-resistance signal that the drug-scope filter otherwise discards.  These
  cross-class counts are the challenge's "statistical association only" tier and
  are kept separate from the interpretable mechanism features.
* Drop engineered features that are constant, below the prevalence threshold,
  or exact duplicates of a more specific retained feature.

Why this script does NOT use outcome labels to remove "unimportant" features
-------------------------------------------------------------------------
Selecting columns by their association with resistant/susceptible labels on
the complete dataset would leak held-out-label information into evaluation.
The script instead writes a label-association *audit* CSV for exploration.  If
we later add supervised selection (for example elastic-net logistic
regression), it must be fitted inside every grouped training fold in the model
pipeline, never once before the train/test split.

The output feature names begin with ``refined_`` so that the existing
``scripts/train_models.py`` command can consume them with
``--feature-prefix refined_``.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]

# These are the union of ``relevant_amr_classes`` in configs/antibiotics.yaml
# at the time this script was written.  They are command-line configurable so
# adding a drug does not require changing the source code.
DEFAULT_RELEVANT_CLASSES = "BETA-LACTAM,AMINOGLYCOSIDE,PHENICOL,EFFLUX,MULTIDRUG"


@dataclass
class EvidenceMetadata:
    """Curated metadata attached to one raw ``amr_`` source column."""

    classes: set[str] = field(default_factory=set)
    subclasses: set[tuple[str, str]] = field(default_factory=set)
    element_symbols: set[str] = field(default_factory=set)
    subtypes: set[str] = field(default_factory=set)


@dataclass
class CandidateFeature:
    """An engineered feature before final prevalence/deduplication filtering."""

    name: str
    kind: str
    source_features: tuple[str, ...]
    description: str
    # A lower value wins when two candidates have identical training profiles.
    # We prefer the most specific biological representation over a broad total.
    priority: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create an auditable reduced AMRFinderPlus feature matrix for training."
    )
    parser.add_argument(
        "--features",
        type=Path,
        default=ROOT / "data" / "interim" / "ecoli_amrfinder_features.csv",
        help="Raw AMRFinderPlus binary feature CSV (one row per genome).",
    )
    parser.add_argument(
        "--labels",
        type=Path,
        default=ROOT / "data" / "processed" / "ecoli_cohort_labels.csv",
        help="Lab-label CSV. Its genome IDs define the training cohort for prevalence filtering.",
    )
    parser.add_argument(
        "--evidence",
        type=Path,
        default=ROOT / "data" / "interim" / "ecoli_amrfinder_evidence.csv",
        help="AMRFinderPlus evidence CSV containing feature/class/subclass metadata.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "data" / "interim" / "ecoli_amrfinder_refined_features.csv",
        help="Output wide matrix containing genome_id and refined_* feature columns.",
    )
    parser.add_argument(
        "--catalog-out",
        type=Path,
        default=ROOT / "data" / "interim" / "ecoli_amrfinder_refinement_catalog.csv",
        help="Audit CSV describing each raw marker and its contribution or removal reason.",
    )
    parser.add_argument(
        "--summary-out",
        type=Path,
        default=ROOT / "data" / "interim" / "ecoli_amrfinder_refinement_summary.csv",
        help="Audit CSV describing every engineered candidate and final selection decision.",
    )
    parser.add_argument(
        "--label-audit-out",
        type=Path,
        default=ROOT / "data" / "interim" / "ecoli_amrfinder_refined_label_audit.csv",
        help="Exploratory per-drug label association audit; it is never used for filtering.",
    )
    parser.add_argument(
        "--manifest-out",
        type=Path,
        default=ROOT / "data" / "interim" / "ecoli_amrfinder_refinement_manifest.json",
        help="JSON record of the retained feature definitions and run parameters.",
    )
    parser.add_argument(
        "--min-prevalence",
        type=int,
        default=2,
        help=(
            "Minimum number of labelled training genomes in which an engineered "
            "feature must be non-zero. Default 2 removes singleton features."
        ),
    )
    parser.add_argument(
        "--relevant-classes",
        default=DEFAULT_RELEVANT_CLASSES,
        help=(
            "Comma-separated AMR classes allowed in the predictive matrix. "
            f"Default: {DEFAULT_RELEVANT_CLASSES}"
        ),
    )
    parser.add_argument(
        "--include-stable-markers",
        action="store_true",
        help=(
            "Also include individual raw marker columns that pass all structural "
            "filters. Disabled by default to favour lower-dimensional mechanism features."
        ),
    )
    return parser.parse_args()


def read_csv(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    """Read a CSV with clear errors instead of silently accepting a bad input."""

    if not path.exists():
        raise FileNotFoundError(f"Input file does not exist: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ValueError(f"CSV has no header row: {path}")
        rows = list(reader)
        return rows, list(reader.fieldnames)


def require_columns(path: Path, columns: Iterable[str], available: Iterable[str]) -> None:
    """Fail early when a supposedly compatible input table has changed shape."""

    available_set = set(available)
    missing = [column for column in columns if column not in available_set]
    if missing:
        raise ValueError(f"{path} is missing required columns: {', '.join(missing)}")


def parse_binary(value: str, genome_id: str, feature: str) -> int:
    """Accept normal CSV representations of binary values and reject everything else."""

    text = (value or "").strip()
    if text in {"0", "0.0"}:
        return 0
    if text in {"1", "1.0"}:
        return 1
    raise ValueError(
        f"Expected binary 0/1 value for genome {genome_id!r}, feature {feature!r}; got {value!r}"
    )


def slug(value: str) -> str:
    """Create deterministic, CSV/model-safe lower-case identifiers."""

    value = str(value or "").strip().lower()
    value = re.sub(r"[^a-z0-9]+", "_", value)
    return value.strip("_") or "unknown"


def class_tokens(value: str) -> set[str]:
    """Split composite AMR classes such as ``AMINOGLYCOSIDE/QUINOLONE``."""

    return {part.strip().upper() for part in str(value or "").split("/") if part.strip()}


def source_base(feature: str) -> str:
    """Remove the repository's generated ``amr_`` prefix and type suffix."""

    base = feature.removeprefix("amr_")
    for suffix in ("_point_disrupt", "_point", "_amr"):
        if base.endswith(suffix):
            return base[: -len(suffix)]
    return base


def derive_family(feature: str, metadata: EvidenceMetadata) -> str:
    """Return a conservative, human-readable gene or mutation family key.

    This deliberately avoids an overly clever biological parser.  Where an
    allele family is obvious from common AMRFinder naming (for example
    blaCTX-M-15 -> blaCTX-M), it is grouped.  Otherwise the complete normalized
    marker name is retained as its own family rather than guessing incorrectly.
    A future project-specific curated family map can replace this heuristic.
    """

    base = source_base(feature)
    subtype_is_point = any(value.upper() == "POINT" for value in metadata.subtypes)
    if subtype_is_point or feature.endswith("_point") or feature.endswith("_point_disrupt"):
        # AMRFinder point-marker names begin with the affected gene, such as
        # gyrA_S83L.  Aggregating by target gene avoids one column per codon/allele.
        return f"{base.split('_', 1)[0]}_point_mutation"

    # These rules cover common allele-rich resistance families in the current
    # matrix.  The expression is intentionally small and transparent.
    family_rules = (
        (r"^(blactx_m)(?:_.*)?$", r"\1"),
        (r"^(blatem)(?:_.*)?$", r"\1"),
        (r"^(blashv)(?:_.*)?$", r"\1"),
        (r"^(blacmy)(?:_.*)?$", r"\1"),
        (r"^(blaoxa)(?:_.*)?$", r"\1"),
        (r"^(aac_\d+)(?:_.*)?$", r"\1"),
        (r"^(aph_\d+)(?:_.*)?$", r"\1"),
        (r"^(ant_\d+)(?:_.*)?$", r"\1"),
        (r"^(dfra)(?:_?\d+.*)?$", r"\1"),
        (r"^(cata|catb)(?:_?\d+.*)?$", r"\1"),
        (r"^(qnrb|qnrs)(?:_?\d+.*)?$", r"\1"),
        (r"^(mcr)(?:_?\d+.*)?$", r"\1"),
    )
    for pattern, replacement in family_rules:
        if re.match(pattern, base):
            return re.sub(pattern, replacement, base)
    return base


def load_evidence_metadata(path: Path) -> dict[str, EvidenceMetadata]:
    """Combine all AMRFinder hit rows belonging to the same feature column."""

    evidence_rows, columns = read_csv(path)
    require_columns(path, ["feature", "class", "subclass", "element_symbol", "subtype"], columns)
    metadata: dict[str, EvidenceMetadata] = defaultdict(EvidenceMetadata)

    for row in evidence_rows:
        feature = (row["feature"] or "").strip()
        if not feature:
            continue
        entry = metadata[feature]
        entry.classes.update(class_tokens(row["class"]))
        # Retain class/subclass pairs because a subclass name can be ambiguous
        # when used without its parent AMR class.
        for amr_class in class_tokens(row["class"]):
            subclass = str(row["subclass"] or "").strip().upper()
            if subclass and subclass != amr_class:
                entry.subclasses.add((amr_class, subclass))
        if row["element_symbol"]:
            entry.element_symbols.add(str(row["element_symbol"]).strip())
        if row["subtype"]:
            entry.subtypes.add(str(row["subtype"]).strip())
    return dict(metadata)


def write_csv(path: Path, columns: list[str], rows: Iterable[dict[str, object]]) -> None:
    """Write deterministic CSVs, creating only the requested output directories."""

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def profile_key(values: list[int]) -> tuple[int, ...]:
    """Use the training-genome feature profile to identify exact duplicate columns."""

    return tuple(values)


def label_association_rows(
    labels: list[dict[str, str]],
    retained: list[CandidateFeature],
    feature_values: dict[str, dict[str, int]],
) -> list[dict[str, object]]:
    """Make a transparent association audit without influencing feature selection.

    The audit can help a researcher inspect whether a mechanism appears more
    often among resistant or susceptible isolates for a given drug.  It is not
    a causal claim and must not be used to pre-filter a held-out evaluation set.
    """

    result: list[dict[str, object]] = []
    for antibiotic in sorted({row["antibiotic"] for row in labels if row.get("antibiotic")}):
        drug_rows = [
            row
            for row in labels
            if row.get("antibiotic") == antibiotic and row.get("label") in {"resistant", "susceptible"}
        ]
        for candidate in retained:
            resistant_present = susceptible_present = 0
            resistant_absent = susceptible_absent = 0
            for row in drug_rows:
                value = feature_values[candidate.name].get(row["genome_id"])
                # Ignore a label whose genome is not present in the feature table.
                if value is None:
                    continue
                if value:
                    if row["label"] == "resistant":
                        resistant_present += 1
                    else:
                        susceptible_present += 1
                elif row["label"] == "resistant":
                    resistant_absent += 1
                else:
                    susceptible_absent += 1

            present_total = resistant_present + susceptible_present
            absent_total = resistant_absent + susceptible_absent
            result.append(
                {
                    "antibiotic": antibiotic,
                    "feature": candidate.name,
                    "feature_kind": candidate.kind,
                    "resistant_present": resistant_present,
                    "susceptible_present": susceptible_present,
                    "resistant_absent": resistant_absent,
                    "susceptible_absent": susceptible_absent,
                    "resistant_rate_when_present": (
                        round(resistant_present / present_total, 6) if present_total else ""
                    ),
                    "resistant_rate_when_absent": (
                        round(resistant_absent / absent_total, 6) if absent_total else ""
                    ),
                    "warning": "Exploratory only; do not use this full-cohort audit for pre-split feature selection.",
                }
            )
    return result


def main() -> None:
    args = parse_args()
    if args.min_prevalence < 2:
        raise SystemExit("--min-prevalence must be at least 2; singleton features must be removed.")

    relevant_classes = {
        value.strip().upper() for value in args.relevant_classes.split(",") if value.strip()
    }
    if not relevant_classes:
        raise SystemExit("--relevant-classes must contain at least one AMR class.")

    raw_rows, all_columns = read_csv(args.features)
    label_rows, label_columns = read_csv(args.labels)
    require_columns(args.features, ["genome_id"], all_columns)
    require_columns(args.labels, ["genome_id", "antibiotic", "label"], label_columns)
    if not raw_rows:
        raise SystemExit(f"No feature rows found in {args.features}")

    # Only AMRFinderPlus feature columns are candidates.  A future version can
    # explicitly accept a different prefix, but silently mixing metadata into
    # the model input would be unsafe.
    source_features = [column for column in all_columns if column.startswith("amr_")]
    if not source_features:
        raise SystemExit(f"No columns beginning with 'amr_' found in {args.features}")

    raw_by_id: dict[str, dict[str, int]] = {}
    for row in raw_rows:
        genome_id = str(row["genome_id"] or "").strip()
        if not genome_id:
            raise ValueError(f"Blank genome_id in {args.features}")
        if genome_id in raw_by_id:
            raise ValueError(f"Duplicate genome_id {genome_id!r} in {args.features}")
        raw_by_id[genome_id] = {
            feature: parse_binary(row.get(feature, ""), genome_id, feature) for feature in source_features
        }

    labelled_genome_ids = {
        str(row["genome_id"] or "").strip() for row in label_rows if str(row["genome_id"] or "").strip()
    }
    training_ids = sorted(set(raw_by_id).intersection(labelled_genome_ids))
    if not training_ids:
        raise SystemExit("No genome IDs overlap between the feature and label tables.")

    metadata_by_feature = load_evidence_metadata(args.evidence)
    source_prevalence = {
        feature: sum(raw_by_id[genome_id][feature] for genome_id in training_ids)
        for feature in source_features
    }
    training_n = len(training_ids)

    # Classify raw marker columns before deriving aggregates.  The catalogue
    # records this decision so a researcher can inspect every exclusion.
    source_status: dict[str, dict[str, object]] = {}
    eligible_for_aggregation: list[str] = []
    for feature in source_features:
        metadata = metadata_by_feature.get(feature)
        prevalence = source_prevalence[feature]
        classes = metadata.classes if metadata else set()
        matching_classes = sorted(classes.intersection(relevant_classes))
        if not metadata:
            status, reason = "excluded", "missing_evidence_metadata"
        elif not matching_classes:
            status, reason = "excluded", "class_not_relevant_to_configured_drugs"
        elif prevalence == 0:
            status, reason = "excluded", "absent_from_labelled_training_genomes"
        elif prevalence == training_n:
            status, reason = "excluded", "constant_on_labelled_training_genomes"
        else:
            # A singleton source marker can still contribute to an aggregate.
            # It will never be emitted by itself unless --include-stable-markers
            # is set and it passes the final prevalence check (which it cannot at
            # the default threshold of two).
            status = "aggregate_candidate"
            reason = "eligible_for_biological_aggregation"
            eligible_for_aggregation.append(feature)
        source_status[feature] = {
            "status": status,
            "reason": reason,
            "matching_classes": matching_classes,
        }

    # Co-resistance burden set: every non-constant marker that carries any AMR
    # class annotation, NOT only the configured drug classes.  Individually these
    # off-target markers (quinolone, sulfonamide, tetracycline, ...) are dropped
    # from the interpretable matrix, but their *total count* is a legitimate
    # statistical signal: multidrug-resistant lineages accumulate determinants
    # across many classes at once, so overall burden correlates with resistance
    # to the target drugs. This is the challenge's "statistical association only"
    # evidence tier, kept separate from the mechanism features above so the two
    # are never confused.  Constant-present intrinsic genes are excluded (they
    # would only add a fixed offset); markers absent on the training cohort are
    # kept in the definition because they can still fire for a future genome.
    coresistance_features = sorted(
        feature
        for feature in source_features
        if source_prevalence[feature] < training_n
        and metadata_by_feature.get(feature)
        and metadata_by_feature[feature].classes
    )

    # Build groups from all non-constant, drug-relevant markers.  Boolean OR is
    # intentional: multiple alleles in one family represent presence of that
    # mechanism, rather than pretending that a genome with two copies has twice
    # the biological effect.
    groups: dict[tuple[str, str], set[str]] = defaultdict(set)
    for feature in eligible_for_aggregation:
        metadata = metadata_by_feature[feature]
        for amr_class in sorted(metadata.classes.intersection(relevant_classes)):
            groups[("class", amr_class)].add(feature)
            for parent_class, subclass in sorted(metadata.subclasses):
                if parent_class == amr_class:
                    # ``Subclass`` is AMRFinder's own drug-spectrum label (for
                    # example BETA-LACTAM/CEPHALOSPORIN).  Preserve that wording
                    # instead of overstating it as a mechanistic annotation.
                    groups[("subclass", f"{amr_class}__{subclass}")].add(feature)
            family = derive_family(feature, metadata)
            groups[("family", f"{amr_class}__{family}")].add(feature)

    candidates: list[CandidateFeature] = []
    kind_priority = {"family": 0, "subclass": 1, "class": 2, "summary": 3, "marker": 4}
    for (kind, group_name), members in sorted(groups.items()):
        safe_group_name = slug(group_name)
        candidates.append(
            CandidateFeature(
                name=f"refined_{kind}_{safe_group_name}",
                kind=kind,
                source_features=tuple(sorted(members)),
                description=f"Binary presence of one or more {kind} markers in {group_name}.",
                priority=kind_priority[kind],
            )
        )

    # Two compact burden features capture information shared across mechanisms
    # without adding one raw column per gene.  They are numeric counts, not
    # binary indicators; logistic regression can consume them directly.
    if eligible_for_aggregation:
        candidates.extend(
            [
                CandidateFeature(
                    name="refined_summary_relevant_marker_count",
                    kind="summary",
                    source_features=tuple(sorted(eligible_for_aggregation)),
                    description="Count of non-constant, drug-relevant raw AMR markers present.",
                    priority=kind_priority["summary"],
                ),
                CandidateFeature(
                    name="refined_summary_relevant_class_count",
                    kind="summary",
                    source_features=tuple(sorted(eligible_for_aggregation)),
                    description="Count of configured AMR classes with at least one marker present.",
                    priority=kind_priority["summary"],
                ),
            ]
        )

    # Cross-class co-resistance burden.  These span every AMR class, capturing
    # the multidrug-resistant-lineage signal that the drug-scope filter above
    # deliberately excludes from the mechanism features.  They are the only
    # features permitted to draw on off-target markers, and only as aggregate
    # counts, never as identifiable per-drug evidence.
    if coresistance_features:
        candidates.extend(
            [
                CandidateFeature(
                    name="refined_summary_total_marker_count",
                    kind="summary",
                    source_features=tuple(coresistance_features),
                    description="Count of present AMR markers across all classes (co-resistance burden).",
                    priority=kind_priority["summary"],
                ),
                CandidateFeature(
                    name="refined_summary_total_class_count",
                    kind="summary",
                    source_features=tuple(coresistance_features),
                    description="Count of distinct AMR classes (any class) with at least one marker present.",
                    priority=kind_priority["summary"],
                ),
            ]
        )

    if args.include_stable_markers:
        # This flag is useful once a larger cohort exists.  It enables an
        # explicit comparison between a compact mechanism model and one that
        # retains sufficiently common individual markers.
        for feature in eligible_for_aggregation:
            candidates.append(
                CandidateFeature(
                    name=f"refined_marker_{feature.removeprefix('amr_')}",
                    kind="marker",
                    source_features=(feature,),
                    description=f"Individual stable AMRFinderPlus marker copied from {feature}.",
                    priority=kind_priority["marker"],
                )
            )

    # Calculate candidate values for every supplied genome.  Prevalence and
    # removal decisions below use only ``training_ids`` so the three unlabelled
    # current rows do not influence the training feature schema.
    candidate_values: dict[str, dict[str, int]] = {}
    candidate_summary: list[dict[str, object]] = []
    for candidate in candidates:
        values: dict[str, int] = {}
        for genome_id, feature_row in raw_by_id.items():
            member_values = [feature_row[feature] for feature in candidate.source_features]
            if candidate.name in {
                "refined_summary_relevant_marker_count",
                "refined_summary_total_marker_count",
            }:
                values[genome_id] = sum(member_values)
            elif candidate.name in {
                "refined_summary_relevant_class_count",
                "refined_summary_total_class_count",
            }:
                # The relevant variant restricts to configured drug classes; the
                # total variant counts every AMR class present in the genome.
                restrict = (
                    relevant_classes
                    if candidate.name == "refined_summary_relevant_class_count"
                    else None
                )
                present_classes = set()
                for feature in candidate.source_features:
                    if not feature_row[feature]:
                        continue
                    feature_classes = metadata_by_feature[feature].classes
                    present_classes.update(
                        feature_classes.intersection(restrict) if restrict else feature_classes
                    )
                values[genome_id] = len(present_classes)
            else:
                values[genome_id] = int(any(member_values))
        candidate_values[candidate.name] = values

    # Keep only features with variation and at least the requested number of
    # positive training genomes.  For count features, prevalence means non-zero
    # prevalence; a count that is always zero or always non-zero is still not a
    # useful standalone predictor in this cohort.
    prelim_retained: list[CandidateFeature] = []
    for candidate in candidates:
        training_values = [candidate_values[candidate.name][genome_id] for genome_id in training_ids]
        prevalence = sum(value > 0 for value in training_values)
        unique_values = len(set(training_values))
        if unique_values < 2:
            status, reason = "removed", "constant_on_labelled_training_genomes"
        elif prevalence < args.min_prevalence:
            status, reason = "removed", f"prevalence_below_{args.min_prevalence}"
        else:
            status, reason = "candidate_retained_before_duplicate_check", "varies_and_meets_prevalence"
            prelim_retained.append(candidate)
        candidate_summary.append(
            {
                "feature": candidate.name,
                "feature_kind": candidate.kind,
                "source_feature_count": len(candidate.source_features),
                "source_features": ";".join(candidate.source_features),
                "training_prevalence": prevalence,
                "training_prevalence_rate": round(prevalence / training_n, 6),
                "training_unique_values": unique_values,
                "status": status,
                "reason": reason,
                "description": candidate.description,
            }
        )

    # A family feature can occasionally have the same training profile as its
    # parent subclass or class feature.  Retaining both would double-count the
    # same evidence.  Keep the more specific/lower-priority-number candidate,
    # with deterministic name ordering as a final tie-breaker.
    retained: list[CandidateFeature] = []
    seen_profiles: dict[tuple[int, ...], CandidateFeature] = {}
    summary_by_name = {row["feature"]: row for row in candidate_summary}
    for candidate in sorted(prelim_retained, key=lambda item: (item.priority, item.name)):
        training_profile = profile_key(
            [candidate_values[candidate.name][genome_id] for genome_id in training_ids]
        )
        existing = seen_profiles.get(training_profile)
        if existing is None:
            seen_profiles[training_profile] = candidate
            retained.append(candidate)
        else:
            entry = summary_by_name[candidate.name]
            entry["status"] = "removed"
            entry["reason"] = f"exact_duplicate_of_{existing.name}"

    retained_names = {candidate.name for candidate in retained}
    for row in candidate_summary:
        if row["feature"] in retained_names:
            row["status"] = "retained"
            row["reason"] = "varies_and_meets_prevalence_after_duplicate_check"

    if not retained:
        raise SystemExit(
            "No refined features survived. Lower --min-prevalence or check the evidence/label cohort overlap."
        )

    # Record exactly how every original marker was treated.  This makes a
    # reduced matrix reviewable rather than an opaque black-box transformation.
    source_to_retained = defaultdict(list)
    for candidate in retained:
        for feature in candidate.source_features:
            source_to_retained[feature].append(candidate.name)
    catalog_rows: list[dict[str, object]] = []
    for feature in sorted(source_features):
        metadata = metadata_by_feature.get(feature)
        initial = source_status[feature]
        contributions = sorted(source_to_retained.get(feature, []))
        if contributions:
            final_status = "contributes_to_retained_aggregate"
            final_reason = ";".join(contributions)
        else:
            final_status = str(initial["status"])
            final_reason = str(initial["reason"])
        catalog_rows.append(
            {
                "source_feature": feature,
                "training_prevalence": source_prevalence[feature],
                "training_prevalence_rate": round(source_prevalence[feature] / training_n, 6),
                "amr_classes": ";".join(sorted(metadata.classes)) if metadata else "",
                "amr_subclasses": (
                    ";".join(f"{amr_class}/{subclass}" for amr_class, subclass in sorted(metadata.subclasses))
                    if metadata
                    else ""
                ),
                "element_symbols": ";".join(sorted(metadata.element_symbols)) if metadata else "",
                "derived_family": derive_family(feature, metadata) if metadata else "",
                "matching_configured_classes": ";".join(initial["matching_classes"]),
                "final_status": final_status,
                "final_reason_or_retained_features": final_reason,
            }
        )

    # Keep all input genome rows in the refined matrix.  Existing training code
    # already inner-joins this table with labels, while retaining extra rows
    # helps inspect the transformation and supports future cohort expansion.
    refined_rows: list[dict[str, object]] = []
    for genome_id in [str(row["genome_id"]).strip() for row in raw_rows]:
        output_row: dict[str, object] = {"genome_id": genome_id}
        for candidate in retained:
            output_row[candidate.name] = candidate_values[candidate.name][genome_id]
        refined_rows.append(output_row)

    output_columns = ["genome_id", *[candidate.name for candidate in retained]]
    write_csv(args.output, output_columns, refined_rows)
    write_csv(
        args.catalog_out,
        [
            "source_feature",
            "training_prevalence",
            "training_prevalence_rate",
            "amr_classes",
            "amr_subclasses",
            "element_symbols",
            "derived_family",
            "matching_configured_classes",
            "final_status",
            "final_reason_or_retained_features",
        ],
        catalog_rows,
    )
    write_csv(
        args.summary_out,
        [
            "feature",
            "feature_kind",
            "source_feature_count",
            "source_features",
            "training_prevalence",
            "training_prevalence_rate",
            "training_unique_values",
            "status",
            "reason",
            "description",
        ],
        sorted(candidate_summary, key=lambda row: str(row["feature"])),
    )
    write_csv(
        args.label_audit_out,
        [
            "antibiotic",
            "feature",
            "feature_kind",
            "resistant_present",
            "susceptible_present",
            "resistant_absent",
            "susceptible_absent",
            "resistant_rate_when_present",
            "resistant_rate_when_absent",
            "warning",
        ],
        label_association_rows(label_rows, retained, candidate_values),
    )

    # The manifest is deliberately sufficient to inspect exactly what went into
    # a model artifact.  A future prediction-time transformer can consume this
    # file to apply the same retained definitions to a new AMRFinder matrix.
    manifest = {
        "feature_type": "amrfinderplus_refined_biological_aggregates",
        "input_features": str(args.features),
        "input_labels": str(args.labels),
        "input_evidence": str(args.evidence),
        "labelled_training_genomes": training_ids,
        "labelled_training_genome_count": training_n,
        "input_genome_count": len(raw_by_id),
        "input_source_feature_count": len(source_features),
        "retained_feature_count": len(retained),
        "min_prevalence": args.min_prevalence,
        "relevant_amr_classes": sorted(relevant_classes),
        "include_stable_markers": args.include_stable_markers,
        "label_leakage_note": (
            "Outcome labels define the training cohort and an exploratory audit only; "
            "they do not determine retained features. Any supervised selection must run inside each training fold."
        ),
        "retained_features": [
            {
                "name": candidate.name,
                "kind": candidate.kind,
                "source_features": list(candidate.source_features),
                "description": candidate.description,
            }
            for candidate in retained
        ],
    }
    args.manifest_out.parent.mkdir(parents=True, exist_ok=True)
    args.manifest_out.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    removed_source_count = sum(
        1 for feature in source_features if not source_to_retained.get(feature)
    )
    print(f"Input genomes: {len(raw_by_id)}")
    print(f"Labelled training genomes used for filtering: {training_n}")
    print(f"Input AMRFinder source features: {len(source_features)}")
    print(f"Retained refined features: {len(retained)}")
    print(f"Raw markers with no retained contribution: {removed_source_count}")
    print(f"Wrote refined matrix: {args.output}")
    print(f"Wrote source-marker catalogue: {args.catalog_out}")
    print(f"Wrote engineered-feature summary: {args.summary_out}")
    print(f"Wrote exploratory label audit: {args.label_audit_out}")
    print(f"Wrote refinement manifest: {args.manifest_out}")


if __name__ == "__main__":
    main()
