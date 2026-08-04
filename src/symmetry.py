"""Exact surviving symmetries of the normalized CW(105,36) phase fibre.

The normalized contracted row in CRT q-coordinates is

    b'_0=-3,  b'_15=b'_25=b'_30=3.

Its affine stabilizer inside AGL(1,Z_35) has no translations and consists of
12 multipliers.  The phase equations additionally admit multiplication of all
nonzero phases by a sixth root and complex conjugation.  Everything here is
finite and exact; no floating-point roots of unity are used.
"""

from __future__ import annotations

import math
from collections import Counter
from dataclasses import dataclass
from typing import Iterable, Sequence

from .derive_icw import affine_transform
from .phase_reduction import (
    Q_MODULUS,
    SPECIAL_POSITIONS,
    normalized_b_grid,
)


@dataclass(frozen=True)
class SupportOrbitReport:
    allowed_positions: tuple[int, ...]
    stabilizer: tuple[int, ...]
    action_orbits: tuple[tuple[int, ...], ...]
    fixed_support_counts: tuple[tuple[int, int], ...]
    orbit_count: int


def units_modulo(n: int) -> tuple[int, ...]:
    return tuple(value for value in range(n) if math.gcd(value, n) == 1)


def affine_stabilizer_of_normalized_contraction() -> tuple[tuple[int, int], ...]:
    """Return all ``(multiplier, shift)`` fixing the normalized b' exactly."""
    row = normalized_b_grid()
    return tuple(
        (multiplier, shift)
        for multiplier in units_modulo(Q_MODULUS)
        for shift in range(Q_MODULUS)
        if affine_transform(row, multiplier, shift, 1) == row
    )


def contraction_multiplier_stabilizer() -> tuple[int, ...]:
    """Return the multiplier projection of the exact affine stabilizer."""
    stabilizer = affine_stabilizer_of_normalized_contraction()
    if any(shift != 0 for _, shift in stabilizer):
        raise AssertionError("unexpected nonzero translation in normalized stabilizer")
    return tuple(multiplier for multiplier, _ in stabilizer)


CONTRACTION_STABILIZER = contraction_multiplier_stabilizer()
ALLOWED_PHASE_POSITIONS = tuple(
    q for q in range(Q_MODULUS) if q not in SPECIAL_POSITIONS
)


def position_action_orbits(
    positions: Iterable[int] = ALLOWED_PHASE_POSITIONS,
    multipliers: Sequence[int] = CONTRACTION_STABILIZER,
) -> tuple[tuple[int, ...], ...]:
    """Return orbits of a position set under a multiplier group."""
    position_set = {int(q) % Q_MODULUS for q in positions}
    unseen = set(position_set)
    result: list[tuple[int, ...]] = []
    while unseen:
        start = min(unseen)
        orbit = tuple(
            sorted({(int(multiplier) * start) % Q_MODULUS for multiplier in multipliers})
        )
        if not set(orbit) <= position_set:
            raise ValueError("the supplied positions are not invariant under the multipliers")
        unseen.difference_update(orbit)
        result.append(orbit)
    return tuple(result)


def multiplier_cycles_on_allowed(multiplier: int) -> tuple[tuple[int, ...], ...]:
    """Cycle decomposition of one stabilizer multiplier on the 31 allowed q's."""
    h = int(multiplier) % Q_MODULUS
    if h not in CONTRACTION_STABILIZER:
        raise ValueError(f"{multiplier} is not in the normalized contraction stabilizer")
    unseen = set(ALLOWED_PHASE_POSITIONS)
    cycles: list[tuple[int, ...]] = []
    while unseen:
        start = min(unseen)
        cycle: list[int] = []
        current = start
        while current not in cycle:
            if current not in unseen:
                raise AssertionError("multiplier cycles overlap unexpectedly")
            cycle.append(current)
            unseen.remove(current)
            current = (h * current) % Q_MODULUS
        if current != start:
            raise AssertionError("cycle did not close at its start")
        cycles.append(tuple(cycle))
    return tuple(cycles)


def fixed_weight_support_count(multiplier: int, weight: int = 12) -> int:
    """Count invariant supports of a given size by a cycle-subset polynomial."""
    if not 0 <= weight <= len(ALLOWED_PHASE_POSITIONS):
        return 0
    coefficients = [0] * (weight + 1)
    coefficients[0] = 1
    for cycle in multiplier_cycles_on_allowed(multiplier):
        length = len(cycle)
        for degree in range(weight, length - 1, -1):
            coefficients[degree] += coefficients[degree - length]
    return coefficients[weight]


def support_orbit_report(weight: int = 12) -> SupportOrbitReport:
    """Apply Burnside's lemma to allowed supports of the specified weight."""
    fixed = tuple(
        (multiplier, fixed_weight_support_count(multiplier, weight))
        for multiplier in CONTRACTION_STABILIZER
    )
    numerator = sum(count for _, count in fixed)
    if numerator % len(CONTRACTION_STABILIZER):
        raise AssertionError("Burnside numerator is not divisible by group order")
    return SupportOrbitReport(
        allowed_positions=ALLOWED_PHASE_POSITIONS,
        stabilizer=CONTRACTION_STABILIZER,
        action_orbits=position_action_orbits(),
        fixed_support_counts=fixed,
        orbit_count=numerator // len(CONTRACTION_STABILIZER),
    )


def transform_support(support: Iterable[int], multiplier: int) -> tuple[int, ...]:
    h = int(multiplier) % Q_MODULUS
    if h not in CONTRACTION_STABILIZER:
        raise ValueError(f"{multiplier} is not in the normalized contraction stabilizer")
    values = {int(q) % Q_MODULUS for q in support}
    if values & SPECIAL_POSITIONS:
        raise ValueError("a normalized phase support cannot contain a special position")
    return tuple(sorted((h * q) % Q_MODULUS for q in values))


def canonical_support(support: Iterable[int]) -> tuple[int, ...]:
    raw = tuple(int(q) % Q_MODULUS for q in support)
    values = tuple(sorted(set(raw)))
    if len(values) != len(raw):
        raise ValueError("support contains duplicate positions")
    if set(values) & SPECIAL_POSITIONS:
        raise ValueError("a normalized phase support cannot contain a special position")
    return min(transform_support(values, h) for h in CONTRACTION_STABILIZER)


def transform_phases(
    phases: Sequence[int | None],
    *,
    multiplier: int = 1,
    global_phase: int = 0,
    conjugate: bool = False,
) -> tuple[int | None, ...]:
    """Apply a normalized-fibre phase symmetry.

    Old position ``q`` is sent to ``multiplier*q``.  A nonzero phase exponent
    ``e`` is sent to ``e+global_phase`` or ``-e+global_phase`` under
    conjugation, all modulo six.
    """
    if len(phases) != Q_MODULUS:
        raise ValueError(f"expected {Q_MODULUS} phase positions")
    h = int(multiplier) % Q_MODULUS
    if h not in CONTRACTION_STABILIZER:
        raise ValueError(f"{multiplier} is not in the normalized contraction stabilizer")
    eta = int(global_phase) % 6
    transformed: list[int | None] = [None] * Q_MODULUS
    for q, phase in enumerate(phases):
        target = (h * q) % Q_MODULUS
        if phase is None:
            transformed[target] = None
            continue
        exponent = int(phase)
        if not 0 <= exponent < 6:
            raise ValueError(f"invalid phase exponent at q={q}: {phase!r}")
        transformed[target] = ((-exponent if conjugate else exponent) + eta) % 6
    if any(transformed[q] is not None for q in SPECIAL_POSITIONS):
        raise AssertionError("a stabilizer action moved phase support into a special position")
    return tuple(transformed)


def phase_key(phases: Sequence[int | None]) -> tuple[int, ...]:
    """Comparable exact encoding: inactive=-1, active exponents=0,...,5."""
    return tuple(-1 if phase is None else int(phase) for phase in phases)


def canonical_phases(phases: Sequence[int | None]) -> tuple[int | None, ...]:
    """Canonicalize under H x (mu_6 semidirect conjugation), 144 actions."""
    candidates = (
        transform_phases(
            phases,
            multiplier=multiplier,
            global_phase=global_phase,
            conjugate=conjugate,
        )
        for multiplier in CONTRACTION_STABILIZER
        for global_phase in range(6)
        for conjugate in (False, True)
    )
    return min(candidates, key=phase_key)


def cycle_type(multiplier: int) -> tuple[tuple[int, int], ...]:
    """Return ``(cycle_length, multiplicity)`` pairs on allowed positions."""
    counts = Counter(len(cycle) for cycle in multiplier_cycles_on_allowed(multiplier))
    return tuple(sorted(counts.items()))
