from __future__ import annotations

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd

from data_fetchers.fetch_bvbrc import (
    collapse_labels,
    evaluate_gate,
    extract_observations,
    resolve_columns,
    separate_conflicts,
)
from data_fetchers.select_genomes import select
from src.predictor.grouped_split import grouped_three_way_split
from src.predictor.predict import predict_target
from src.predictor.train import train_target
from src.predictor.target_gate import apply_target_gate, target_presence_from_features
from src.report.clinician_report import LAB_DISCLAIMER, build_report
from src.genome_reader.fasta_qc import fasta_stats
from src.genome_reader.grouping import assign_homology_groups
from src.genome_reader.build_features import add_evidence_features, parse_amrfinder
from src.genome_reader.target_presence import build_target_presence, find_target_markers
from src.genome_reader.run_amrfinder import run_amrfinder
from src.trust.calibration import apply_no_call, reliability_table
from src.trust.evidence import classify_evidence
from src.validation.scorecard import scorecard
from src.pipeline.router import resolve_route
from targets_config import ACTIVE_TARGETS, validate_config


class PipelineContractTests(unittest.TestCase):
    def test_config_registry_routes_bacteria_and_rejects_unimplemented_viruses(self) -> None:
        self.assertEqual(validate_config(), [])
        self.assertIn(("Escherichia coli", "ciprofloxacin"), ACTIVE_TARGETS)
        bacterial = resolve_route("Escherichia coli", "ciprofloxacin")
        self.assertTrue(bacterial.supported)
        virus = resolve_route("Influenza A virus", "oseltamivir")
        self.assertFalse(virus.supported)
        self.assertIn("reader", virus.unsupported_reason())

    def test_small_cohort_selector_balances_every_passing_target(self) -> None:
        rows = []
        for antibiotic in ("ciprofloxacin", "ampicillin", "ceftriaxone"):
            for phenotype in ("resistant", "susceptible"):
                rows.extend(
                    {
                        "genome_id": f"{antibiotic}-{phenotype}-{index}",
                        "species": "Escherichia coli",
                        "antibiotic": antibiotic,
                        "phenotype": phenotype,
                    }
                    for index in range(20)
                )
        selected, manifest = select(pd.DataFrame(rows), max_genomes=60, min_per_class=10, seed=42)
        counts = selected.groupby(["antibiotic", "phenotype"]).size()
        self.assertTrue(all(counts.loc[pair] >= 10 for pair in counts.index))
        self.assertEqual(manifest["selected_genomes"], 60)

    def test_calibrated_training_writes_model_receipt(self) -> None:
        rows = [
            {
                "genome_id": f"g{index}",
                "species": "Escherichia coli",
                "antibiotic": "ciprofloxacin",
                "phenotype": "resistant" if index % 2 else "susceptible",
                "y": index % 2,
                "homology_group_id": f"group-{index}",
                "gene:gyrA": index % 2,
                "gene:other": (index // 2) % 2,
            }
            for index in range(90)
        ]
        with TemporaryDirectory() as directory:
            receipt = train_target(
                pd.DataFrame(rows),
                output_dir=Path(directory),
                species="Escherichia coli",
                antibiotic="ciprofloxacin",
            )
            self.assertTrue(Path(receipt["model_path"]).exists())
            self.assertTrue(Path(receipt["split_path"]).exists())
            self.assertEqual(receipt["split"]["rows"]["train"] + receipt["split"]["rows"]["calibration"] + receipt["split"]["rows"]["test"], 90)

    def test_amrfinder_features_exclude_non_amr_plus_hits(self) -> None:
        with TemporaryDirectory() as directory:
            path = Path(directory) / "g1.tsv"
            pd.DataFrame(
                [
                    {
                        "Element symbol": "blaTEM-1",
                        "Element name": "beta-lactamase",
                        "Type": "AMR",
                        "Subtype": "AMR",
                        "Class": "BETA-LACTAM",
                        "% Identity to reference": "100.0",
                    },
                    {
                        "Element symbol": "stressGene",
                        "Element name": "stress response",
                        "Type": "STRESS",
                        "Subtype": "STRESS",
                        "Class": "STRESS",
                        "% Identity to reference": "99.0",
                    },
                ]
            ).to_csv(path, sep="\t", index=False)
            _, features, evidence = parse_amrfinder(path)
            self.assertEqual(features, {"gene:blaTEM-1"})
            self.assertEqual(evidence["amr_class"].tolist(), ["BETA-LACTAM"])

    def test_evidence_feature_engineering_adds_counts(self) -> None:
        matrix = pd.DataFrame({"genome_id": ["g1"], "gene:bla": [1]})
        evidence = pd.DataFrame(
            {
                "genome_id": ["g1", "g1"],
                "feature": ["gene:bla", "gene:gyrA_S83L"],
                "gene_symbol": ["bla", "gyrA_S83L"],
                "amr_class": ["BETA-LACTAM", "QUINOLONE"],
                "amr_subclass": ["BETA-LACTAM", "QUINOLONE"],
                "subtype": ["AMR", "POINT"],
            }
        )
        engineered = add_evidence_features(matrix, evidence)
        self.assertEqual(int(engineered.loc[0, "amr_hit_count"]), 2)
        self.assertEqual(int(engineered.loc[0, "amr_class_count"]), 2)
        self.assertEqual(int(engineered.loc[0, "amr_point_mutation_count"]), 1)

    def test_amrfinder_command_options_are_validated(self) -> None:
        with self.assertRaises(ValueError):
            run_amrfinder(Path("missing.fna"), Path("out.tsv"), threads=0)

    def test_lab_audit_preserves_conflicts_and_drops_them_from_labels(self) -> None:
        raw = pd.DataFrame(
            {
                "genome_id": ["g1", "g1", "g2", "g3"],
                "genome_name": ["Escherichia coli x"] * 4,
                "species": ["Escherichia coli"] * 4,
                "antibiotic": ["CIPRO", "ciprofloxacin", "amp", "ciprofloxacin"],
                "phenotype": ["Resistant", "Susceptible", "Susceptible", "Resistant"],
                "evidence": ["Phenotype", "Phenotype", "Predicted", "Phenotype"],
            }
        )
        columns = resolve_columns(raw.columns)
        observations = extract_observations(
            raw,
            columns,
            (("Escherichia coli", "ciprofloxacin"), ("Escherichia coli", "ampicillin")),
        )
        clean, conflicts = separate_conflicts(observations)
        labels = collapse_labels(clean)
        self.assertEqual(len(conflicts), 2)
        self.assertEqual(labels["genome_id"].tolist(), ["g3"])
        gate = evaluate_gate(labels, (("Escherichia coli", "ciprofloxacin"),))
        self.assertFalse(bool(gate.loc[0, "preliminary_gate_pass"]))

    def test_prediction_and_report_keep_target_gate_and_disclaimer(self) -> None:
        class DummyModel:
            def predict_proba(self, frame: pd.DataFrame):
                return [[0.1, 0.9] for _ in range(len(frame))]

        predictions = predict_target(
            {"model": DummyModel(), "feature_columns": ["gene:gyrA"]},
            pd.DataFrame({"genome_id": ["g1"], "gene:gyrA": [1]}),
            target_present=None,
            target_names=("gyrA",),
        )
        self.assertEqual(predictions.loc[0, "prediction"], "no-call")
        self.assertIn(LAB_DISCLAIMER, build_report(predictions, species="Escherichia coli", antibiotic="ciprofloxacin"))

    def test_fasta_qc_and_duplicate_grouping_are_reproducible(self) -> None:
        with TemporaryDirectory() as directory:
            root = Path(directory)
            sequence = "ACGT" * 80
            first = root / "g1.fna"
            second = root / "g2.fna"
            third = root / "g3.fna"
            first.write_text(f">g1\n{sequence}\n")
            second.write_text(f">g2\n{sequence}\n")
            third.write_text(f">g3\n{'ACGT' * 40 + 'G' * 160}\n")
            stats = fasta_stats(first)
            self.assertEqual(stats["total_bases"], 320)
            groups = assign_homology_groups([first, second, third])
            self.assertEqual(groups["g1"], groups["g2"])
            self.assertNotEqual(groups["g1"], groups["g3"])

    def test_target_gate_never_allows_work_without_target(self) -> None:
        decision = apply_target_gate("likely to work", None, target_names=("gyrA", "parC"))
        self.assertEqual(decision.prediction, "no-call")
        self.assertEqual(decision.target_status, "uncertain")

        decision = apply_target_gate("likely to work", False, target_names=("gyrA",))
        self.assertEqual(decision.prediction, "no-call")
        self.assertEqual(decision.target_status, "absent")

    def test_target_presence_empty_features_is_unknown(self) -> None:
        self.assertIsNone(target_presence_from_features({}, ("gyrA",)))
        self.assertIs(target_presence_from_features({"gyrA": True}, ("gyrA",)), True)
        self.assertIs(target_presence_from_features({"parC": True}, ("gyrA",)), False)

    def test_annotation_target_presence_is_tristate_and_marker_based(self) -> None:
        annotations = pd.DataFrame(
            {
                "gene": ["gyrA"],
                "product": ["DNA gyrase subunit A"],
            }
        )
        markers = {"genes": ("gyrA",), "product_keywords": ("DNA gyrase subunit",)}
        self.assertEqual(find_target_markers(annotations, markers), ["gene:gyra", "product:dna gyrase subunit"])
        frame = build_target_presence(
            ["g1", "g2"],
            {"g1": annotations},
            {"ciprofloxacin": markers},
            {"ciprofloxacin": ("gyrA",)},
        )
        self.assertEqual(frame.loc[frame.genome_id.eq("g1"), "target_status"].item(), "present")
        self.assertIsNone(frame.loc[frame.genome_id.eq("g2"), "target_present"].item())

    def test_grouped_split_has_no_group_overlap(self) -> None:
        frame = pd.DataFrame(
            {
                "genome_id": [f"g{i}" for i in range(18)],
                "homology_group_id": [f"group-{i // 2}" for i in range(18)],
                "y": [i % 2 for i in range(18)],
            }
        )
        assignment = grouped_three_way_split(frame, random_state=7)
        groups = {
            split: set(assignment.frame.loc[assignment.frame["split"].eq(split), "homology_group_id"])
            for split in ("train", "calibration", "test")
        }
        self.assertFalse(groups["train"] & groups["calibration"])
        self.assertFalse(groups["train"] & groups["test"])
        self.assertFalse(groups["calibration"] & groups["test"])

    def test_no_call_and_reliability_outputs(self) -> None:
        decisions = apply_no_call([0.95, 0.52, 0.10], min_confidence=0.70)
        self.assertEqual(decisions["prediction"].tolist(), ["likely to fail", "no-call", "likely to work"])
        reliability, brier = reliability_table([1, 0, 0], [0.95, 0.52, 0.10])
        self.assertEqual(len(reliability), 10)
        self.assertGreaterEqual(brier, 0)
        self.assertLessEqual(brier, 1)

    def test_evidence_types_do_not_call_association_causal(self) -> None:
        self.assertEqual(classify_evidence(["gene:gyrA"], ["gene:gyrA"]), "known_resistance_determinant")
        self.assertEqual(classify_evidence(["gene:unknown"], ["gene:gyrA"]), "statistical_association")
        self.assertEqual(classify_evidence([], ["gene:gyrA"]), "no_known_signal")

    def test_scorecard_contains_required_metrics(self) -> None:
        result = scorecard([0, 1, 0, 1], [0.10, 0.90, 0.40, 0.60])
        for key in (
            "balanced_accuracy",
            "resistant_recall",
            "susceptible_recall",
            "f1",
            "auroc",
            "pr_auc",
            "brier",
            "no_call_rate",
            "called_accuracy",
        ):
            self.assertIn(key, result)


if __name__ == "__main__":
    unittest.main()
