import json
import unittest

from variant_pipeline.alignments import parse_sam
from variant_pipeline.caller import call_candidates, fisher_exact_two_sided
from variant_pipeline.config import load_manifest
from variant_pipeline.io import read_fasta
from variant_pipeline.models import Parameters


class PipelineEvidenceTest(unittest.TestCase):
    def setUp(self):
        self.config = load_manifest("examples/manifest.json")
        self.reference = read_fasta(self.config.reference_fasta)["synthetic-region"]
        self.alignments = parse_sam(self.config.reads_sam, self.config.parameters)
        self.candidates = {
            candidate.variant_id: candidate
            for candidate in call_candidates(
                self.alignments, self.reference, self.config.parameters
            )
        }

    def test_known_substitution_and_indels_have_stable_identity(self):
        self.assertIn("synthetic-region:50:C>G", self.candidates)
        self.assertIn("synthetic-region:8:TACG>T", self.candidates)
        self.assertIn("synthetic-region:58:A>ACGT", self.candidates)

    def test_overlapping_pair_counts_once(self):
        self.assertEqual(
            self.candidates["synthetic-region:58:A>ACGT"].evidence["supports"]["templates"],
            3,
        )

    def test_strand_bias_keeps_substitution_in_review(self):
        candidate = self.candidates["synthetic-region:50:C>G"]
        self.assertEqual(candidate.status, "review")
        self.assertLess(candidate.evidence["coverage"]["strand_bias_fisher_p"], 0.05)

    def test_low_quality_evidence_is_separated(self):
        candidate = self.candidates["synthetic-region:50:C>G"]
        self.assertEqual(candidate.evidence["low_quality_templates"], 2)
        self.assertEqual(candidate.evidence["opposes"]["reference_templates"], 6)

    def test_insufficient_coverage_is_reported(self):
        candidate = self.candidates["synthetic-region:60:G>T"]
        self.assertEqual(candidate.status, "review")
        self.assertIsNotNone(candidate.evidence["coverage_note"])
        self.assertEqual(candidate.evidence["coverage"]["template_depth"], 4)

    def test_balanced_high_coverage_variant_passes(self):
        candidate = self.candidates["synthetic-region:65:A>T"]
        self.assertEqual(candidate.status, "candidate")
        self.assertEqual(candidate.evidence["supports"]["templates"], 6)
        self.assertEqual(candidate.evidence["opposes"]["reference_templates"], 6)

    def test_manifest_carries_sample_and_reference_context(self):
        self.assertEqual(self.config.sample_id, "synthetic-demo-001")
        self.assertEqual(self.config.reference_version, "demo-v1")

    def test_fisher_exact(self):
        self.assertLess(fisher_exact_two_sided(1, 3, 6, 0), 0.05)
        self.assertEqual(fisher_exact_two_sided(0, 0, 0, 0), 1.0)


if __name__ == "__main__":
    unittest.main()
