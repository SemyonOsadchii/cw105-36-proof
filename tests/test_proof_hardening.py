from __future__ import annotations

from src.verify_fibre_character_bridge import (
    EXPECTED_CONTRACTION_ROWS,
    run_exact_checks,
)


def test_symbolic_fibre_character_bridge_is_exhaustive_and_exact() -> None:
    report = run_exact_checks()

    assert report["crt_points_checked"] == 105
    assert report["character_coefficients_checked"] == 9
    assert report["legal_columns_checked"] == 27
    assert report["active_character_states_checked"] == 6
    assert report["norm_one_minus_omega"] == 3


def test_both_multiplier_fixed_contraction_rows_are_handled_directly() -> None:
    report = run_exact_checks()

    assert len(EXPECTED_CONTRACTION_ROWS) == 2
    assert report["contraction_rows_checked"] == 2
    assert report["active_columns_by_row"] == [12, 12]
    assert report["forced_constant_columns_by_row"] == [4, 4]
    assert report["contraction_correlations_checked"] == 70
