import json

from variant_pipeline.association import (
    AssociationRelation,
    MATCH_INSUFFICIENT,
    PENDING,
    RELATION_INCONSISTENT,
    RESOLVED,
    reconcile_associations,
)
from variant_pipeline.caller import call_variants
from variant_pipeline.config import Parameters
from variant_pipeline.fasta import Reference
from variant_pipeline.normalize import normalize
from variant_pipeline.samio import Fragment, Observation


PARAMS = Parameters(min_coverage_for_calling=1, min_alt_count=1)
REF = Reference("c", "ACGTACGT")


def _fragment(qname, source, allele_specs, strand="F"):
    fragment = Fragment(qname=qname, strand=strand, source_id=source,
                        source_name=source, batch_id=source)
    for pos, alt in allele_specs:
        variant = normalize("c", pos, REF.sequence[pos - 1], alt, REF)
        fragment.alt[variant.canonical_key] = Observation(
            variant=variant, base_quals=(40,), strand=strand, low_quality=False)
        fragment.covered.append((pos, pos))
    return fragment


def _relation(a, b, source_a="s1", source_b="s1", **kwargs):
    return AssociationRelation(
        association_id="unit-1", molecule_id="mol-1", a=a, b=b,
        source_id_a=source_a, source_id_b=source_b, **kwargs)


def test_late_endpoint_is_retained_as_pending_without_dropping_arrived_end():
    first = _fragment("end-a", "s1", [(4, "A")])
    relation = _relation("s1::end-a", "s1::end-b")
    result = reconcile_associations([first], [relation])[0]
    assert result.status == PENDING
    assert result.present_fragment_ids == ["s1::end-a"]
    assert result.missing_fragment_ids == ["s1::end-b"]
    assert result.single_fragment_candidates == {
        "s1::end-a": ["c:g.4T>A"]
    }
    assert result.linked_candidates == ["c:g.4T>A"]


def test_same_qname_different_source_is_relation_inconsistent_not_auto_paired():
    fragments = [
        _fragment("same", "s1", [(4, "A")]),
        _fragment("same", "s2", [(4, "A")]),
    ]
    relation = _relation("s1::same", "s1::same")
    # The explicit endpoint IDs are different records; source mismatch is
    # reported instead of silently joining the imported same-name fragment.
    relation = AssociationRelation(
        "unit-1", "mol-1", "s1::same", "s2::same", "s1", "s1")
    result = reconcile_associations(fragments, [relation])[0]
    assert result.status == RELATION_INCONSISTENT
    assert result.issues[0]["code"] == RELATION_INCONSISTENT


def test_orientation_distance_contradiction_distinguishes_basis_from_no_match():
    fragments = [
        _fragment("a", "s1", [(4, "A")], strand="F"),
        _fragment("b", "s1", [(4, "G")], strand="F"),
    ]
    relation = _relation("s1::a", "s1::b", orientation="opposite")
    result = reconcile_associations(fragments, [relation])[0]
    assert result.status == RELATION_INCONSISTENT
    assert any(i["code"] == RELATION_INCONSISTENT for i in result.issues)

    fragments[1].alt.clear()
    fragments[1].covered = []
    result = reconcile_associations(fragments, [relation])[0]
    assert result.status == MATCH_INSUFFICIENT


def test_association_excludes_single_fragment_candidate_and_counts_one_unit():
    fragments = [
        _fragment("a", "s1", [(4, "A"), (8, "A")], strand="F"),
        _fragment("b", "s1", [(4, "G"), (8, "A")], strand="R"),
    ]
    relation = _relation(
        "s1::a", "s1::b", orientation="opposite", min_distance=4, max_distance=4)
    result = call_variants(fragments, {"c": REF}, PARAMS, [relation])
    by_id = {c.variant_id: c for c in result.candidates}
    assert by_id["c:g.4T>A"].alt_fragments == 1
    assert by_id["c:g.8T>A"].alt_fragments == 1
    assert by_id["c:g.8T>A"].independent_alt_fragments == 1
    assert by_id["c:g.8T>A"].alt_fragment_ids == ["unit-1"]
    assert by_id["c:g.4T>G"].alt_fragments == 0
    assert by_id["c:g.4T>G"].alt_fragments == 0
    unit = result.evidence_units[0]
    assert unit["status"] == RESOLVED
    assert unit["single_fragment_candidates"]["s1::a"] == [
        "c:g.4T>A", "c:g.8T>A"]
    assert unit["linked_candidates"] == ["c:g.4T>A", "c:g.8T>A"]
    assert unit["excluded_candidates"] == ["c:g.4T>G"]


def test_withdrawn_association_restores_original_single_fragment_ambiguity():
    fragments = [
        _fragment("a", "s1", [(4, "G")], strand="F"),
        _fragment("b", "s1", [(4, "G")], strand="R"),
    ]
    relation = _relation("s1::a", "s1::b", orientation="same",
                         status="withdrawn")
    result = call_variants(fragments, {"c": REF}, PARAMS, [relation])
    candidate = result.candidates[0]
    assert candidate.alt_fragments == 2
    assert result.evidence_units[0]["status"] == "withdrawn"


def test_revised_association_restores_original_single_fragment_ambiguity():
    fragments = [
        _fragment("a", "s1", [(4, "A"), (8, "A")], "F"),
        _fragment("b", "s1", [(4, "G"), (8, "A")], "R"),
    ]
    relation = _relation(
        "s1::a", "s1::b", orientation="opposite", min_distance=4, max_distance=4,
        status="revised")
    result = call_variants(fragments, {"c": REF}, PARAMS, [relation])
    by_id = {c.variant_id: c for c in result.candidates}
    assert by_id["c:g.4T>A"].alt_fragments == 1
    assert by_id["c:g.4T>G"].alt_fragments == 1
    assert result.evidence_units[0]["status"] == "revised"
