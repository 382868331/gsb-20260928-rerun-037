"""Tests for multi-source compatible combination search.

These cases are written independently from the association test helpers'
relation layout so that the search contract is checked on its own terms:

* required co-occurrence and mutual exclusion across relations;
* source + fragment identity matching (no same-name substitution);
* a fragment shared by several relations is one entity / one source;
* repeated imports never add entities or independent sources;
* missing endpoints stay unknown rather than failing closed or being filled;
* ranking covers the most relations, then uses the fewest independent sources;
* infeasible searches report the conflicting constraints and endpoints;
* search never confirms a relation; adoption is explicit and keeps sources;
* results are invariant to relation / fragment / constraint ordering.
"""

from variant_pipeline.association import AssociationRelation
from variant_pipeline.combination import (
    CombinationConstraint,
    FEASIBLE,
    FORBIDDEN,
    INFEASIBLE,
    REQUIRED,
    UNKNOWN,
    adopt_solution,
    load_combination_constraints,
    search_combinations,
)
from variant_pipeline.caller import call_variants

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import test_association_units as base


PARAMS = base.PARAMS
REF = base.REF
_fragment = base._fragment


def _relation(aid, a, b, source_a="s1", source_b="s1", **kwargs):
    return AssociationRelation(
        aid, f"mol-{aid}", a, b, source_a, source_b, **kwargs)


def _constraint(cid, kind, relation_ids):
    return CombinationConstraint(cid, kind, tuple(relation_ids))


def _three_fragments():
    return [
        _fragment("x", "s1", [(4, "A")], strand="F"),
        _fragment("y", "s1", [(8, "A")], strand="R"),
        _fragment("z", "s2", [(6, "C")], strand="R"),
    ]


def _two_relations(fragments=None):
    fragments = fragments or _three_fragments()
    r1 = _relation("r1", "s1::x", "s1::y", orientation="opposite")
    r2 = _relation("r2", "s1::x", "s2::z", "s1", "s2",
                   orientation="opposite")
    return fragments, r1, r2


def test_required_relations_must_cooccur_and_share_fragment_counts_once():
    fragments, r1, r2 = _two_relations()
    result = search_combinations(
        fragments, [r1, r2],
        [_constraint("must", REQUIRED, ["r1", "r2"])])
    assert result.status == FEASIBLE
    best = result.solutions[0]
    assert best.covered_relations == ["r1", "r2"]
    # s1::x is shared by both relations but is a single entity/source.
    assert best.fragment_entities == ["s1::x", "s1::y", "s2::z"]
    assert best.independent_sources == ["s1", "s2"]
    report = result.constraint_reports[0]
    assert report["state"] == "satisfied"
    assert report["endpoints"] == ["s1::x", "s1::y", "s2::z"]


def test_solutions_rank_most_relations_then_fewest_sources():
    fragments, r1, r2 = _two_relations()
    result = search_combinations(fragments, [r1, r2], [])
    # Every solution covers as many relations as possible; single-relation
    # solutions exist but rank below the two-relation one.
    cover_counts = [len(s.covered_relations) for s in result.solutions]
    assert cover_counts[0] == max(cover_counts) == 2
    assert 1 in cover_counts
    # Tie-break prefers the plan using fewer independent sources.
    top_tier = [s for s in result.solutions
                if len(s.covered_relations) == cover_counts[0]]
    assert len(top_tier[0].independent_sources) <= len(
        top_tier[-1].independent_sources)


def test_shared_fragment_must_keep_one_candidate_across_relations():
    # x carries two non-overlapping alleles; pairing it with different
    # relations forces one consistent x-candidate per joint assignment.
    fragments = [
        _fragment("x", "s1", [(4, "A"), (8, "A")], strand="F"),
        _fragment("y", "s1", [(4, "A")], strand="R"),
        _fragment("w", "s1", [(8, "A")], strand="R"),
    ]
    r1 = _relation("r1", "s1::x", "s1::y", orientation="opposite")
    r2 = _relation("r2", "s1::x", "s1::w", orientation="opposite")
    result = search_combinations(
        fragments, [r1, r2],
        [_constraint("must", REQUIRED, ["r1", "r2"])])
    full = [s for s in result.solutions if s.covered_relations == ["r1", "r2"]]
    assert full, "the required joint assignment should exist"
    for solution in full:
        assert (solution.assignments["r1"]["s1::x"]
                == solution.assignments["r2"]["s1::x"])


def test_forbidden_relations_never_cooccur_in_one_solution():
    fragments, r1, r2 = _two_relations()
    result = search_combinations(
        fragments, [r1, r2],
        [_constraint("fx", FORBIDDEN, ["r1", "r2"])])
    assert result.status == FEASIBLE
    for solution in result.solutions:
        assert not ({"r1", "r2"} <= set(solution.covered_relations))
    assert result.solutions[0].covered_relations in (["r1"], ["r2"])
    assert result.constraint_reports[0]["state"] == "satisfied"


def test_forbidden_required_when_all_must_hold_is_infeasible_with_endpoints():
    fragments, r1, r2 = _two_relations()
    constraints = [
        _constraint("must", REQUIRED, ["r1", "r2"]),
        _constraint("fx", FORBIDDEN, ["r1", "r2"]),
    ]
    result = search_combinations(fragments, [r1, r2], constraints)
    assert result.status == INFEASIBLE
    assert result.solutions == []
    codes = {c["code"] for c in result.conflicts}
    assert "forbidden_relations_cooccur" in codes
    conflict = next(c for c in result.conflicts
                    if c["code"] == "forbidden_relations_cooccur")
    assert set(conflict["relation_id"].split(",")) == {"r1", "r2"}
    assert set(conflict["endpoints"]) == {"s1::x", "s1::y", "s2::z"}


def test_required_missing_endpoint_is_unknown_not_infeasible():
    fragments = [_fragment("x", "s1", [(4, "A")], strand="F")]
    r1 = _relation("r1", "s1::x", "s1::y", orientation="opposite")
    result = search_combinations(
        fragments, [r1], [_constraint("must", REQUIRED, ["r1"])])
    assert result.status == UNKNOWN
    assert result.unknown_relations == ["r1"]
    # A searchable plan for the arrived relation is still retained.
    assert result.solutions[0].covered_relations == []


def test_same_named_fragment_from_other_source_cannot_substitute():
    # The relation expects s1::q but only s2::q exists; the endpoint stays
    # missing (unknown) and the s2 same-name fragment is never adopted.
    fragments = [
        _fragment("x", "s1", [(4, "A")], strand="F"),
        _fragment("q", "s2", [(4, "A")], strand="R"),
    ]
    r1 = _relation("r1", "s1::x", "s1::q", "s1", "s1",
                   orientation="opposite")
    result = search_combinations(
        fragments, [r1], [_constraint("must", REQUIRED, ["r1"])])
    assert result.status == UNKNOWN
    assert result.unknown_relations == ["r1"]
    entities = {ep for s in result.solutions for ep in s.fragment_entities}
    sources = {src for s in result.solutions for src in s.independent_sources}
    assert "s2::q" not in entities
    assert "s2" not in sources


def test_unknown_relation_id_is_infeasible_without_name_substitution():
    fragments, r1, _r2 = _two_relations()
    result = search_combinations(
        fragments, [r1],
        [_constraint("must", REQUIRED, ["r1", "ghost"])])
    assert result.status == INFEASIBLE
    assert result.unsatisfiable_relations == ["ghost"]
    assert result.conflicts[0]["relation_id"] == "ghost"


def test_duplicate_import_does_not_add_entity_or_source():
    fragments, r1, r2 = _two_relations()
    fragments.append(_fragment("x", "s1", [(4, "A")], strand="F"))
    result = search_combinations(
        fragments, [r1, r2],
        [_constraint("must", REQUIRED, ["r1", "r2"])])
    best = result.solutions[0]
    assert best.fragment_entities == ["s1::x", "s1::y", "s2::z"]
    assert best.independent_sources == ["s1", "s2"]


def test_search_does_not_confirm_relation_and_adoption_is_explicit():
    fragments, r1, r2 = _two_relations()
    result = search_combinations(
        fragments, [r1, r2],
        [_constraint("must", REQUIRED, ["r1", "r2"])])
    assert result.adopted is None
    adopted = adopt_solution(result, 0, rationale="researcher confirms")
    assert adopted["covered_relations"] == ["r1", "r2"]
    assert adopted["independent_sources"] == ["s1", "s2"]
    assert adopted["rationale"] == "researcher confirms"
    # Adopting does not mutate the search result itself.
    assert result.adopted is None


def test_cannot_adopt_an_infeasible_search():
    fragments, r1, _r2 = _two_relations()
    result = search_combinations(
        fragments, [r1],
        [_constraint("must", REQUIRED, ["ghost"])])
    assert result.status == INFEASIBLE
    try:
        adopt_solution(result)
    except ValueError:
        pass
    else:
        raise AssertionError("adoption of infeasible search must fail")


def test_result_invariant_to_input_permutation():
    fragments, r1, r2 = _two_relations()
    constraints = [
        _constraint("must", REQUIRED, ["r1", "r2"]),
        _constraint("fx", FORBIDDEN, ["r1"]),  # single-member exclusion is a no-op
    ]
    baseline = search_combinations(fragments, [r1, r2], constraints).to_dict()
    variants = [
        search_combinations(fragments[::-1], [r2, r1], constraints[::-1]).to_dict(),
        search_combinations(fragments[::-1], [r1, r2], constraints).to_dict(),
        search_combinations(fragments, [r2, r1], constraints[::-1]).to_dict(),
    ]
    for other in variants:
        assert other == baseline


def test_constraint_loader_roundtrip(tmp_path):
    path = tmp_path / "constraints.json"
    path.write_text(
        '{"combination_constraints": ['
        '{"constraint_id": "c1", "kind": "required", "relations": ["r1"]},'
        '{"constraint_id": "c2", "kind": "forbidden", "relations": ["r1", "r2"]}'
        "]}", encoding="utf-8")
    loaded = load_combination_constraints(path)
    assert [c.constraint_id for c in loaded] == ["c1", "c2"]
    assert loaded[1].relations == ("r1", "r2")


def test_pipeline_exposes_combination_search_in_call_result():
    fragments, r1, r2 = _two_relations()
    result = call_variants(
        fragments, {"c": REF}, PARAMS, [r1, r2],
        [_constraint("must", REQUIRED, ["r1", "r2"])])
    assert result.combination_search is not None
    assert result.combination_search["status"] == FEASIBLE
    best = result.combination_search["solutions"][0]
    assert best["covered_relations"] == ["r1", "r2"]
    assert best["independent_sources"] == ["s1", "s2"]
