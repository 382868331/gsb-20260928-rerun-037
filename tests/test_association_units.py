import json

from variant_pipeline.association import (
    AssociationRelation,
    INTRA_MOLECULE_CONFLICT,
    MATCH_INSUFFICIENT,
    MOLECULE_CONFLICT,
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
    # Different alleles at the same site are an intra-molecule conflict:
    # orientation cannot pick a winner and neither end is dropped.
    assert result.status == MOLECULE_CONFLICT
    assert result.conflicting_candidates == ["c:g.4T>A", "c:g.4T>G"]
    assert result.linked_candidates == []
    assert result.excluded_candidates == []
    conflict = result.conflicting_combinations[0]
    assert conflict["code"] == INTRA_MOLECULE_CONFLICT
    assert conflict["a_variant_id"] == "c:g.4T>A"
    assert conflict["b_variant_id"] == "c:g.4T>G"
    assert conflict["pos"] == 4 and conflict["end"] == 4
    assert any(i["code"] == INTRA_MOLECULE_CONFLICT for i in result.issues)

    fragments[1].alt.clear()
    fragments[1].covered = []
    result = reconcile_associations(fragments, [relation])[0]
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
    # The same-site T>A vs T>G disagreement is an intra-molecule conflict:
    # both candidates remain, neither receives support, the site abstains
    # once for the molecule (no import-order winner).
    assert by_id["c:g.4T>A"].alt_fragments == 0
    assert by_id["c:g.4T>A"].conflicting_overlap_fragments == 1
    assert by_id["c:g.4T>G"].alt_fragments == 0
    assert by_id["c:g.4T>G"].conflicting_overlap_fragments == 1
    assert by_id["c:g.4T>A"].other_alt_fragments == 0
    # The shared 8T>A x 8T>A pair has distance 0, violating min_distance=4,
    # so the relation excludes it; it stays recorded (not silently merged).
    assert by_id["c:g.8T>A"].alt_fragments == 0
    assert by_id["c:g.8T>A"].other_alt_fragments == 2
    unit = result.evidence_units[0]
    assert unit["status"] == MOLECULE_CONFLICT
    assert unit["single_fragment_candidates"]["s1::a"] == [
        "c:g.4T>A", "c:g.8T>A"]
    assert unit["linked_candidates"] == []
    assert unit["conflicting_candidates"] == ["c:g.4T>A", "c:g.4T>G"]
    assert unit["excluded_candidates"] == ["c:g.8T>A"]
    assert unit["valid_combinations"] == []
    assert unit["excluded_combinations"] == [{
        "s1::a": "c:g.8T>A", "s1::b": "c:g.8T>A",
        "chrom": "c", "distance": "0"}]


def test_same_site_conflict_counts_as_one_molecule_and_retains_both_sides():
    fragments = [
        _fragment("a", "s1", [(4, "A")], strand="F"),
        _fragment("b", "s1", [(4, "G")], strand="R"),
    ]
    relation = _relation("s1::a", "s1::b", orientation="opposite")
    result = call_variants(fragments, {"c": REF}, PARAMS, [relation])
    by_id = {c.variant_id: c for c in result.candidates}
    unit = result.evidence_units[0]
    assert unit["independent_evidence_units"] == 1
    assert unit["status"] == MOLECULE_CONFLICT
    for vid in ("c:g.4T>A", "c:g.4T>G"):
        candidate = by_id[vid]
        assert candidate.alt_fragments == 0
        assert candidate.independent_alt_fragments == 0
        assert candidate.conflicting_overlap_fragments == 1
        assert candidate.depth == 1
    # Both single-fragment candidates are still visible.
    assert set(unit["single_fragment_candidates"]["s1::a"]) == {"c:g.4T>A"}
    assert set(unit["single_fragment_candidates"]["s1::b"]) == {"c:g.4T>G"}


def test_compatible_non_overlapping_candidates_still_resolve():
    fragments = [
        _fragment("a", "s1", [(4, "A")], strand="F"),
        _fragment("b", "s1", [(8, "A")], strand="R"),
    ]
    relation = _relation(
        "s1::a", "s1::b", orientation="opposite", min_distance=4, max_distance=4)
    result = call_variants(fragments, {"c": REF}, PARAMS, [relation])
    by_id = {c.variant_id: c for c in result.candidates}
    assert by_id["c:g.4T>A"].alt_fragments == 1
    assert by_id["c:g.8T>A"].alt_fragments == 1
    unit = result.evidence_units[0]
    assert unit["status"] == RESOLVED
    assert unit["linked_candidates"] == ["c:g.4T>A", "c:g.8T>A"]
    assert unit["conflicting_candidates"] == []


def test_same_allele_at_overlapping_site_is_not_a_conflict():
    fragments = [
        _fragment("a", "s1", [(4, "A")], strand="F"),
        _fragment("b", "s1", [(4, "A")], strand="R"),
    ]
    relation = _relation("s1::a", "s1::b", orientation="opposite")
    result = call_variants(fragments, {"c": REF}, PARAMS, [relation])
    by_id = {c.variant_id: c for c in result.candidates}
    candidate = by_id["c:g.4T>A"]
    assert candidate.alt_fragments == 1
    assert candidate.independent_alt_fragments == 1
    assert candidate.conflicting_overlap_fragments == 0
    unit = result.evidence_units[0]
    assert unit["status"] == RESOLVED


def _indel_fragment(qname, source, pos, ref, alt, strand="F"):
    fragment = Fragment(qname=qname, strand=strand, source_id=source,
                        source_name=source, batch_id=source)
    variant = normalize("c", pos, ref, alt, REF)
    fragment.alt[variant.canonical_key] = Observation(
        variant=variant, base_quals=(40,), strand=strand, low_quality=False)
    fragment.covered.append((variant.pos, variant.pos + len(variant.ref) - 1))
    return fragment


def test_overlapping_insertion_and_deletion_conflict_retains_both():
    # REF = ACGTACGT (c:g.4 = T, c:g.5 = A).
    # Insertion between pos 4 and 5: anchor ref "T", alt "TA".
    end_a = _indel_fragment("a", "s1", 4, "T", "TA", strand="F")
    ins_id = next(iter(end_a.alt.values())).variant.vid
    # Deletion of ref bases 4-5: anchor "TA", alt "T".
    end_b = _indel_fragment("b", "s1", 4, "TA", "T", strand="R")
    del_id = next(iter(end_b.alt.values())).variant.vid
    relation = _relation("s1::a", "s1::b", orientation="opposite")
    result = call_variants([end_a, end_b], {"c": REF}, PARAMS, [relation])
    by_id = {c.variant_id: c for c in result.candidates}
    assert by_id[ins_id].alt_fragments == 0
    assert by_id[del_id].alt_fragments == 0
    assert by_id[ins_id].conflicting_overlap_fragments == 1
    assert by_id[del_id].conflicting_overlap_fragments == 1
    unit = result.evidence_units[0]
    assert unit["status"] == MOLECULE_CONFLICT
    assert unit["conflicting_candidates"] == sorted([ins_id, del_id])
    assert unit["excluded_candidates"] == []
    conflict = unit["conflicting_combinations"][0]
    assert conflict["a_variant_id"] == ins_id
    assert conflict["b_variant_id"] == del_id


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


def test_withdrawn_conflicting_association_restores_single_fragment_support():
    fragments = [
        _fragment("a", "s1", [(4, "A")], strand="F"),
        _fragment("b", "s1", [(4, "G")], strand="R"),
    ]
    relation = _relation("s1::a", "s1::b", orientation="opposite",
                         status="withdrawn")
    result = call_variants(fragments, {"c": REF}, PARAMS, [relation])
    by_id = {c.variant_id: c for c in result.candidates}
    # Without the (withdrawn) relation, each end is its own fragment again:
    # each allele gets one independent support, no conflict abstention.
    assert by_id["c:g.4T>A"].alt_fragments == 1
    assert by_id["c:g.4T>G"].alt_fragments == 1
    assert by_id["c:g.4T>A"].conflicting_overlap_fragments == 0
    unit = result.evidence_units[0]
    assert unit["status"] == "withdrawn"
    assert unit["linked_candidates"] == ["c:g.4T>A", "c:g.4T>G"]
    assert any(i["code"] == "association_withdrawn" for i in unit["issues"])


def test_late_endpoint_conflict_is_resolved_when_it_arrives_as_one_molecule():
    first = _fragment("end-a", "s1", [(4, "A")])
    relation = _relation("s1::end-a", "s1::end-b", orientation="opposite")
    pending = reconcile_associations([first], [relation])[0]
    assert pending.status == PENDING
    assert pending.independent_evidence_units == 1
    second = _fragment("end-b", "s1", [(4, "G")], strand="R")
    result = call_variants([first, second], {"c": REF}, PARAMS, [relation])
    by_id = {c.variant_id: c for c in result.candidates}
    unit = result.evidence_units[0]
    # Once the late end arrives, both fragments collapse into one molecule
    # with an intra-molecule conflict: depth 1, no doubled voting.
    assert unit["status"] == MOLECULE_CONFLICT
    assert unit["independent_evidence_units"] == 1
    assert by_id["c:g.4T>A"].depth == 1
    assert by_id["c:g.4T>G"].depth == 1
    assert by_id["c:g.4T>A"].alt_fragments == 0
    assert by_id["c:g.4T>G"].alt_fragments == 0


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


def test_same_display_name_different_sources_resolve_by_source_and_id():
    # Endpoints are located by (source, fragment id) jointly: the same
    # display name from two sources is kept as two distinct fragments and
    # the relation resolves identically regardless of import order.
    def make():
        return [
            _fragment("dup", "s1", [(4, "A")], "F"),
            _fragment("dup", "s2", [(8, "A")], "R"),
        ]
    relation = AssociationRelation(
        "unit-1", "mol-1", "dup", "dup", "s1", "s2",
        orientation="opposite", min_distance=4, max_distance=4)
    outcomes = []
    for fragments in (make(), make()[::-1]):
        result = call_variants(fragments, {"c": REF}, PARAMS, [relation])
        unit = result.evidence_units[0]
        by_id = {c.variant_id: c for c in result.candidates}
        outcomes.append((
            unit["status"], unit["valid_combinations"],
            unit["independent_evidence_units"],
            by_id["c:g.4T>A"].alt_fragments, by_id["c:g.8T>A"].alt_fragments,
            by_id["c:g.4T>A"].depth, by_id["c:g.8T>A"].depth))
    assert outcomes[0] == outcomes[1]
    status, combinations, units, alt_a, alt_b, depth_a, depth_b = outcomes[0]
    assert status == RESOLVED
    assert combinations == [{
        "s1::dup": "c:g.4T>A", "s2::dup": "c:g.8T>A",
        "chrom": "c", "distance": "4"}]
    assert (units, alt_a, alt_b, depth_a, depth_b) == (1, 1, 1, 1, 1)
    # Neither source is collapsed as a duplicate import; diagnostics carry
    # the source-qualified identity of each endpoint.
    assert unit["duplicate_fragment_ids"] == []
    assert not any(i["code"] == "duplicate_import" for i in unit["issues"])
    assert unit["fragment_ids"] == ["s1::dup", "s2::dup"]
    assert unit["present_fragment_ids"] == ["s1::dup", "s2::dup"]
    assert set(unit["single_fragment_candidates"]) == {"s1::dup", "s2::dup"}


def test_same_display_name_same_source_reimport_is_one_duplicate_unit():
    # A genuine repeated import (same source, same id) still collapses into
    # one evidence unit and is reported as a duplicate import.
    fragments = [
        _fragment("dup", "s1", [(4, "A")], "F"),
        _fragment("dup", "s1", [(4, "A")], "F"),
        _fragment("other", "s1", [(8, "A")], "R"),
    ]
    relation = _relation(
        "dup", "other", orientation="opposite", min_distance=4, max_distance=4)
    outcomes = []
    for order in (fragments, fragments[::-1]):
        result = reconcile_associations(order, [relation])[0]
        outcomes.append((result.status, result.valid_combinations,
                         result.duplicate_fragment_ids))
    assert outcomes[0] == outcomes[1]
    result = reconcile_associations(fragments, [relation])[0]
    assert result.status == RESOLVED
    assert result.duplicate_fragment_ids == ["s1::dup"]
    assert any(i["code"] == "duplicate_import" for i in result.issues)
    assert result.independent_evidence_units == 1


def test_missing_endpoint_reports_source_qualified_identity():
    present = _fragment("dup", "s1", [(4, "A")], "F")
    relation = AssociationRelation(
        "unit-1", "mol-1", "dup", "dup", "s1", "s2")
    result = reconcile_associations([present], [relation])[0]
    assert result.status == PENDING
    assert result.present_fragment_ids == ["s1::dup"]
    # The missing end is reported with its own source, not merged with the
    # same-named fragment that did arrive from the other source.
    assert result.missing_fragment_ids == ["s2::dup"]
    issue = next(i for i in result.issues if i["code"] == PENDING)
    assert issue["fragment_id"] == "s2::dup"


def test_cross_source_same_name_conflicting_alleles_abstain_consistently():
    def make():
        return [
            _fragment("dup", "s1", [(4, "A")], "F"),
            _fragment("dup", "s2", [(4, "G")], "R"),
        ]
    relation = AssociationRelation(
        "unit-1", "mol-1", "dup", "dup", "s1", "s2", orientation="opposite")
    outcomes = []
    for fragments in (make(), make()[::-1]):
        result = call_variants(fragments, {"c": REF}, PARAMS, [relation])
        unit = result.evidence_units[0]
        by_id = {c.variant_id: c for c in result.candidates}
        outcomes.append((
            unit["status"], unit["conflicting_candidates"],
            unit["conflicting_combinations"],
            unit["independent_evidence_units"],
            by_id["c:g.4T>A"].alt_fragments, by_id["c:g.4T>G"].alt_fragments,
            by_id["c:g.4T>A"].conflicting_overlap_fragments,
            by_id["c:g.4T>G"].conflicting_overlap_fragments,
            by_id["c:g.4T>A"].depth, by_id["c:g.4T>G"].depth))
    assert outcomes[0] == outcomes[1]
    (status, conflicting, combinations, units, alt_a, alt_g,
     conf_a, conf_g, depth_a, depth_g) = outcomes[0]
    assert status == MOLECULE_CONFLICT
    assert conflicting == ["c:g.4T>A", "c:g.4T>G"]
    assert combinations[0]["a_variant_id"] == "c:g.4T>A"
    assert combinations[0]["b_variant_id"] == "c:g.4T>G"
    assert (units, alt_a, alt_g, conf_a, conf_g, depth_a, depth_g) == (
        1, 0, 0, 1, 1, 1, 1)
