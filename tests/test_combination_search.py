import json

import pytest

from variant_pipeline.association import (
    AssociationRelation,
    reconcile_associations,
)
from variant_pipeline.combination import (
    MAX_COVERAGE,
    MIN_SOURCES,
    CombinationConstraint,
    adopt_combination,
    constraints_fingerprint,
    load_combination_adoption,
    load_combination_constraints,
    search_combinations,
)
from variant_pipeline.caller import call_variants
from variant_pipeline.config import (
    Parameters,
    ReferenceConfig,
    RunConfig,
    SampleConfig,
)
from variant_pipeline.fasta import Reference
from variant_pipeline.normalize import normalize
from variant_pipeline.report import build_payload, write_json
from variant_pipeline.samio import Fragment, Observation
from variant_pipeline.store import RunStore


PARAMS = Parameters(min_coverage_for_calling=1, min_alt_count=1)
REF = Reference("c", "ACGTACGT")


def _unit(aid, status="resolved", present=(), missing=()):
    return {
        "association_id": aid,
        "status": status,
        "fragment_ids": list(present) + list(missing),
        "present_fragment_ids": list(present),
        "missing_fragment_ids": list(missing),
        "relation_id": "rule-1",
        "relation_version": "1",
    }


def _constraint(cid, kind, ids):
    return CombinationConstraint(cid, kind, tuple(ids))


def test_co_required_group_is_all_or_none_and_exclusive_is_respected():
    units = [
        _unit("u1", present=["s1::a", "s1::b"]),
        _unit("u2", present=["s2::c", "s2::d"]),
        _unit("u3", present=["s3::e", "s3::f"]),
    ]
    constraints = [
        _constraint("co", "co_required", ["u1", "u2"]),
        _constraint("me", "mutually_exclusive", ["u2", "u3"]),
    ]
    result = search_combinations(units, constraints)
    assert result.status == "ok"
    plan = result.plans[MAX_COVERAGE]
    # u3 cannot join because it is mutually exclusive with u2; the
    # co-required pair stays together.
    assert plan["association_ids"] == ["u1", "u2"]
    assert plan["covered_relation_count"] == 2
    assert plan["status"] == "candidate"
    # Covering the co-required pair needs two sources; the mutually
    # exclusive u3 alone covers one relation with a single source.
    lean = result.plans[MIN_SOURCES]
    assert lean["association_ids"] == ["u3"]
    assert lean["independent_source_count"] == 1
    assert result.comparison["same_plan"] is False


def test_objectives_can_select_different_plans_for_comparison():
    units = [
        _unit("A", present=["s1::a", "s2::b"]),
        _unit("B", present=["s2::c", "s3::d"]),
        _unit("C", present=["s1::e", "s1::f"]),
    ]
    constraints = [
        _constraint("me1", "mutually_exclusive", ["A", "C"]),
        _constraint("me2", "mutually_exclusive", ["B", "C"]),
    ]
    result = search_combinations(units, constraints)
    assert result.comparison["same_plan"] is False
    wide = result.plans[MAX_COVERAGE]
    lean = result.plans[MIN_SOURCES]
    assert wide["association_ids"] == ["A", "B"]
    assert wide["independent_source_count"] == 3
    assert lean["association_ids"] == ["C"]
    assert lean["independent_source_count"] == 1


def test_missing_endpoint_stays_unknown_and_blocks_co_required_group():
    units = [
        _unit("u1", present=["s1::a", "s1::b"]),
        _unit("u2", status="pending_endpoint",
              present=["s1::c"], missing=["s1::late"]),
    ]
    constraints = [_constraint("co", "co_required", ["u1", "u2"])]
    result = search_combinations(units, constraints)
    assert result.status == "ok"
    assert [u["association_id"] for u in result.unknown_relations] == ["u2"]
    assert result.unknown_relations[0]["missing_fragment_ids"] == ["s1::late"]
    assert result.blocked_constraints[0]["blocked_by_association_ids"] == ["u2"]
    # The pending relation is never confirmable, so the co-required group
    # cannot be covered.
    assert result.plans[MAX_COVERAGE]["association_ids"] == []


def test_shared_fragment_counts_once_and_sources_are_distinct():
    units = [
        _unit("u1", present=["s1::x", "s1::y"]),
        _unit("u2", present=["s1::x", "s2::z"]),
    ]
    result = search_combinations(units, [])
    plan = result.plans[MAX_COVERAGE]
    assert plan["association_ids"] == ["u1", "u2"]
    # s1::x is shared by both relations but is one entity.
    assert plan["entity_fragment_ids"] == ["s1::x", "s1::y", "s2::z"]
    assert plan["entity_count"] == 3
    assert plan["independent_sources"] == ["s1", "s2"]
    assert plan["independent_source_count"] == 2


def test_conflicting_constraints_report_endpoints_and_no_solution():
    units = [
        _unit("u1", present=["s1::a", "s1::b"]),
        _unit("u2", present=["s2::c"], missing=["s2::d"]),
    ]
    constraints = [
        _constraint("co", "co_required", ["u1", "u2"]),
        _constraint("me", "mutually_exclusive", ["u1", "u2"]),
    ]
    result = search_combinations(units, constraints)
    assert result.status == "no_solution"
    assert result.plans == {}
    conflict = result.conflicts[0]
    assert conflict["co_required_constraint"] == "co"
    assert conflict["mutually_exclusive_constraint"] == "me"
    assert conflict["shared_association_ids"] == ["u1", "u2"]
    assert conflict["endpoints"]["u1"]["present_fragment_ids"] == [
        "s1::a", "s1::b"]
    assert conflict["endpoints"]["u2"]["missing_fragment_ids"] == ["s2::d"]


def test_input_permutation_yields_equivalent_results():
    units = [
        _unit("u1", present=["s1::a"]),
        _unit("u2", present=["s2::b"]),
        _unit("u3", present=["s3::c"]),
        _unit("u4", status="pending_endpoint", present=["s1::d"],
              missing=["s1::e"]),
    ]
    constraints = [
        _constraint("co", "co_required", ["u1", "u2"]),
        _constraint("me", "mutually_exclusive", ["u2", "u3"]),
    ]
    baseline = search_combinations(units, constraints).to_dict()
    permuted = search_combinations(
        list(reversed(units)), list(reversed(constraints))).to_dict()
    assert baseline == permuted
    assert (constraints_fingerprint(constraints)
            == constraints_fingerprint(list(reversed(constraints))))


def test_unknown_constraint_association_id_is_rejected():
    with pytest.raises(ValueError, match="unknown association ids"):
        search_combinations(
            [_unit("u1", present=["s1::a"])],
            [_constraint("co", "co_required", ["u1", "nope"])])
    with pytest.raises(ValueError, match="unknown constraint kind"):
        CombinationConstraint("bad", "sometimes", ("u1", "u2"))


def _fragment(qname, source, specs, strand="F"):
    frag = Fragment(qname=qname, strand=strand, source_id=source,
                    source_name=source, batch_id=source)
    for pos, alt in specs:
        variant = normalize("c", pos, REF.sequence[pos - 1], alt, REF)
        frag.alt[variant.canonical_key] = Observation(
            variant=variant, base_quals=(40,), strand=strand, low_quality=False)
        frag.covered.append((pos, pos))
    return frag


def test_duplicate_import_adds_no_support_weight_in_search():
    fragments = [
        _fragment("a", "s1", [(4, "A")], "F"),
        _fragment("b", "s1", [(8, "A")], "R"),
        _fragment("a", "s1", [(4, "A")], "F"),  # repeated import, collapsed
    ]
    relations = [
        AssociationRelation("u1", "m1", "s1::a", "s1::b", "s1", "s1",
                            orientation="opposite",
                            min_distance=4, max_distance=4),
        AssociationRelation("u2", "m2", "s1::a", "s1::b", "s1", "s1",
                            orientation="opposite",
                            min_distance=4, max_distance=4),
    ]
    units = [u.to_dict() for u in reconcile_associations(fragments, relations)]
    result = search_combinations(units, [])
    plan = result.plans[MAX_COVERAGE]
    # Two relations share both fragments: entities counted once, and the
    # duplicate import contributes no extra source or entity.
    assert plan["entity_fragment_ids"] == ["s1::a", "s1::b"]
    assert plan["entity_count"] == 2
    assert plan["independent_sources"] == ["s1"]


def _publish_store(tmp_path, relations, fragments):
    store = RunStore(tmp_path / "store")
    result = call_variants(fragments, {"c": REF}, PARAMS, relations)
    cfg = RunConfig(
        sample=SampleConfig(sample_id="s"),
        reference=ReferenceConfig(
            fasta_path=str(tmp_path / "ref.fasta"), name="c", assembly="r1"),
        parameters=PARAMS,
    )
    (tmp_path / "ref.fasta").write_text(">c\nACGTACGT\n")
    reads = tmp_path / "reads.sam"
    reads.write_text("@HD\n")
    staging = tmp_path / "stage"
    staging.mkdir()
    write_json(build_payload(cfg, result), staging / "candidates.json")
    for name in ("candidates.tsv", "coverage.tsv", "report.txt"):
        (staging / name).write_text("")
    fp = store.fingerprint(config_dict=cfg.to_dict(), reads_path=reads,
                           fasta_path=tmp_path / "ref.fasta")
    info = store.publish(
        {n: staging / n for n in (
            "candidates.json", "candidates.tsv", "coverage.tsv", "report.txt")},
        config_dict=cfg.to_dict(), reads_path=reads,
        fasta_path=tmp_path / "ref.fasta", fingerprint=fp)
    return store, info


def test_adoption_records_selected_sources_and_keeps_base_version(tmp_path):
    fragments = [
        _fragment("a", "s1", [(4, "A")], "F"),
        _fragment("b", "s2", [(8, "A")], "R"),
    ]
    relations = [
        AssociationRelation("u1", "m1", "s1::a", "s2::b", "s1", "s2",
                            orientation="opposite",
                            min_distance=4, max_distance=4,
                            relation_id="rule-1"),
    ]
    store, info = _publish_store(tmp_path, relations, fragments)
    payload = json.loads(
        (info.path / "candidates.json").read_text(encoding="utf-8"))
    with pytest.raises(ValueError, match="repeats an association id"):
        _constraint("co", "co_required", ["u1", "u1"])

    # The search alone must not confirm anything.
    constraints = []
    result = search_combinations(payload["evidence_units"], constraints)
    plan = result.plans[MAX_COVERAGE]
    assert plan["status"] == "candidate"
    assert load_combination_adoption(store, info.version) is None

    record = adopt_combination(
        store, info.version, constraints, plan["plan_id"],
        rationale="reviewed")
    assert record["status"] == "adopted"
    assert record["selected_sources"] == ["s1", "s2"]
    assert record["independent_source_count"] == 2
    assert record["does_not_modify_base_version"] is True
    # Re-adopting archives the previous adoption instead of overwriting it.
    again = adopt_combination(store, info.version, constraints, plan["plan_id"])
    assert again["supersedes"] == record["adoption_id"]
    assert len(again["history"]) == 1
    assert load_combination_adoption(store, info.version)["adoption_id"] == \
        again["adoption_id"]
    # The adoption is part of the export audit trail.
    audit = store.audit_sources(info.version)
    assert audit["combination_adoption"]["plan_id"] == plan["plan_id"]
    with pytest.raises(ValueError, match="unknown plan"):
        adopt_combination(store, info.version, constraints, "plan-missing")


def test_load_constraints_from_json(tmp_path):
    path = tmp_path / "constraints.json"
    path.write_text(json.dumps({"constraints": [
        {"constraint_id": "co", "kind": "co_required",
         "association_ids": ["u1", "u2"]},
        {"constraint_id": "me", "kind": "mutually_exclusive",
         "association_ids": ["u2", "u3"]},
    ]}))
    constraints = load_combination_constraints(path)
    assert [c.constraint_id for c in constraints] == ["co", "me"]
    assert constraints[0].association_ids == ("u1", "u2")
