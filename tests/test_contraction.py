from __future__ import annotations

from src.derive_icw import (
    affine_transform,
    canonical_affine_sign,
    enumerate_icw,
    enumerate_scalar_candidates_dfs,
    enumerate_scalar_candidates_mitm,
    integer_periodic_autocorrelation,
    multiplier_orbits,
    multiplier_order,
    normalized_b,
)

EXPECTED_ORBITS = (
    (0,),
    (1, 4, 16, 29, 11, 9),
    (2, 8, 32, 23, 22, 18),
    (3, 12, 13, 17, 33, 27),
    (5, 20, 10),
    (6, 24, 26, 34, 31, 19),
    (7, 28),
    (14, 21),
    (15, 25, 30),
)


def test_multiplier_four_orbits() -> None:
    assert multiplier_order() == 6
    assert multiplier_orbits() == EXPECTED_ORBITS


def test_two_scalar_enumerators_agree_exactly() -> None:
    sizes = tuple(len(orbit) for orbit in EXPECTED_ORBITS)
    dfs = enumerate_scalar_candidates_dfs(sizes)
    mitm = enumerate_scalar_candidates_mitm(sizes)
    assert dfs == mitm
    assert len(dfs) == 1434


def test_unique_icw_equivalence_class_is_reproduced() -> None:
    result = enumerate_icw()
    assert len(result.invariant_solutions) == 2
    assert len(result.equivalence_classes) == 1
    assert normalized_b() in result.invariant_solutions
    assert all(
        integer_periodic_autocorrelation(row) == (36,) + (0,) * 34
        for row in result.invariant_solutions
    )


def test_two_invariant_representatives_are_explicitly_equivalent() -> None:
    normalized = normalized_b()
    image = affine_transform(normalized, multiplier=3, shift=0, sign=1)
    assert {i for i, value in enumerate(image) if value == 3} == {15, 25, 30}
    assert canonical_affine_sign(normalized) == canonical_affine_sign(image)


def test_normalized_polynomial_directly() -> None:
    row = normalized_b()
    assert sum(row) == 6
    assert sum(value * value for value in row) == 36
    assert integer_periodic_autocorrelation(row) == (36,) + (0,) * 34
