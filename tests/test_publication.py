from __future__ import annotations

from pathlib import Path

from tools.check_public_hygiene import tracked_paths
from validation.publication_audit import _is_irreducible_f3, audit_repository


ROOT = Path(__file__).resolve().parents[1]


def test_publication_is_complete_audited_and_frozen() -> None:
    report = audit_repository(ROOT)
    assert report["status"] == "passed"
    assert report["support_count"] == 420
    assert report["factor_orbit_branches"] == 2
    assert report["affine_correlation_shift_checks"] == 29_400
    assert report["four_shift_symbolic_core"] == "passed"
    assert report["contraction_rows_checked"] == 2
    assert report["mutation_boundaries_checked"] == [
        "factor_branch",
        "support_certificate",
        "affine_representative_and_stabilizer",
        "four_shift_terms",
        "fibre_to_phase_table",
        "contraction_conventions",
    ]


def test_publication_tree_is_small_and_canonical() -> None:
    paths = tracked_paths(ROOT)
    assert len(paths) <= 75
    assert "CITATION.cff" in paths
    assert "checksums.sha256" in paths
    assert [path for path in paths if path.endswith(".pdf")] == [
        "paper/cw105_36_nonexistence.pdf"
    ]
    assert [path for path in paths if path.endswith(".tex")] == [
        "paper/cw105_36_nonexistence.tex"
    ]
    assert [path for path in paths if path.endswith(".bib")] == [
        "paper/references.bib"
    ]
    assert [path for path in paths if path.endswith(".sha256")] == [
        "checksums.sha256"
    ]
    assert [path for path in paths if path.endswith(".md")] == ["README.md"]
    assert [path for path in paths if "mod3_weight12_support_masks" in path] == [
        "data/mod3_weight12_support_masks.txt"
    ]


def test_paper_has_exactly_the_focused_scope() -> None:
    source = (ROOT / "paper" / "cw105_36_nonexistence.tex").read_text(
        encoding="utf-8"
    )
    for required in (
        r"\CGW(35,12;6)",
        r"\CW(105,36)",
        "420 supports",
        "Factor-orbit lemma",
        "four-shift obstruction",
        "complete solution locus",
        r"C_{ah}(y)=C_h(x)",
        "Theorem~4.1",
        r"\gcd(35,36)=1",
        r"4\equiv 3^{10}\pmod {35}",
        r"B_1(X)",
        r"B_2(X)",
        r"R(t,s)",
        r"\sum_{s=0}^2\omega^{-s}R(t,s)",
    ):
        assert required in source


def test_publication_auditor_distinguishes_irreducible_f3_factors() -> None:
    assert _is_irreducible_f3((1, 1, 1, 1, 1))
    assert _is_irreducible_f3((1, 1, 1, 1, 1, 1, 1))
    assert not _is_irreducible_f3((2, 0, 1))
