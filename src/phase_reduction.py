"""Exact CRT and sixth-root phase reduction for lifting CW(105,36).

No floating-point arithmetic is used. Sixth roots are represented in the
Eisenstein integer basis (1, omega).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Sequence

from .eisenstein import (
    Eisenstein,
    ONE_MINUS_OMEGA,
    OMEGA,
    OMEGA2,
    SIXTH_ROOTS,
    ZERO,
)

N = 105
Q_MODULUS = 35
R_MODULUS = 3
SPECIAL_NEGATIVE = 0
SPECIAL_POSITIVE = (15, 25, 30)
SPECIAL_POSITIONS = frozenset((SPECIAL_NEGATIVE, *SPECIAL_POSITIVE))
ACTIVE_PHASE_COUNT = 12

# phase exponent k -> (a_{q,0},a_{q,1},a_{q,2})
PHASE_TO_COLUMN: dict[int, tuple[int, int, int]] = {
    0: (1, -1, 0),
    1: (1, 0, -1),
    2: (0, 1, -1),
    3: (-1, 1, 0),
    4: (-1, 0, 1),
    5: (0, -1, 1),
}
COLUMN_TO_PHASE = {column: exponent for exponent, column in PHASE_TO_COLUMN.items()}


@dataclass(frozen=True)
class PhaseVerificationResult:
    valid: bool
    active_count: int
    correlations: tuple[Eisenstein, ...]
    errors: tuple[str, ...]


def index_from_qr(q: int, r: int) -> int:
    return (3 * int(q) + 35 * int(r)) % N


def qr_from_index(index: int) -> tuple[int, int]:
    i = int(index) % N
    # 3^{-1}=12 mod 35; 35^{-1}=2 mod 3.
    return (12 * i) % Q_MODULUS, (2 * i) % R_MODULUS


def vector_to_grid(values: Sequence[int]) -> tuple[tuple[int, int, int], ...]:
    row = tuple(int(value) for value in values)
    if len(row) != N:
        raise ValueError(f"expected a length-{N} vector")
    grid = [[0] * R_MODULUS for _ in range(Q_MODULUS)]
    for index, value in enumerate(row):
        q, r = qr_from_index(index)
        grid[q][r] = value
    return tuple(tuple(column) for column in grid)  # type: ignore[return-value]


def grid_to_vector(grid: Sequence[Sequence[int]]) -> tuple[int, ...]:
    if len(grid) != Q_MODULUS or any(len(column) != R_MODULUS for column in grid):
        raise ValueError("expected a 35 by 3 grid")
    row = [0] * N
    for q, column in enumerate(grid):
        for r, value in enumerate(column):
            row[index_from_qr(q, r)] = int(value)
    return tuple(row)


def normalized_b_original() -> tuple[int, ...]:
    row = [0] * Q_MODULUS
    row[0] = -3
    for index in (5, 10, 20):
        row[index] = 3
    return tuple(row)


def normalized_b_grid() -> tuple[int, ...]:
    row = [0] * Q_MODULUS
    row[SPECIAL_NEGATIVE] = -3
    for q in SPECIAL_POSITIVE:
        row[q] = 3
    return tuple(row)


def contract_mod_35(values: Sequence[int]) -> tuple[int, ...]:
    row = tuple(int(value) for value in values)
    if len(row) != N:
        raise ValueError(f"expected a length-{N} vector")
    return tuple(row[j] + row[j + 35] + row[j + 70] for j in range(35))


def contract_grid(values: Sequence[int]) -> tuple[int, ...]:
    return tuple(sum(column) for column in vector_to_grid(values))


def reconstruct_from_phase(phases: Sequence[int | None]) -> tuple[int, ...]:
    """Reconstruct the ternary length-105 row from a length-35 phase sequence."""
    if len(phases) != Q_MODULUS:
        raise ValueError("expected 35 phase positions")
    grid: list[tuple[int, int, int]] = []
    for q, phase in enumerate(phases):
        if q == SPECIAL_NEGATIVE:
            if phase is not None:
                raise ValueError("the forced negative column must have zero phase")
            grid.append((-1, -1, -1))
        elif q in SPECIAL_POSITIVE:
            if phase is not None:
                raise ValueError("forced positive columns must have zero phase")
            grid.append((1, 1, 1))
        elif phase is None:
            grid.append((0, 0, 0))
        else:
            exponent = int(phase)
            if exponent not in PHASE_TO_COLUMN:
                raise ValueError(f"invalid sixth-root exponent at q={q}: {phase}")
            grid.append(PHASE_TO_COLUMN[exponent])
    return grid_to_vector(grid)


def extract_phase(values: Sequence[int]) -> tuple[int | None, ...]:
    """Invert ``reconstruct_from_phase`` for rows in the normalized lift fibre."""
    grid = vector_to_grid(values)
    phases: list[int | None] = []
    for q, column in enumerate(grid):
        if q == SPECIAL_NEGATIVE:
            if column != (-1, -1, -1):
                raise ValueError("q=0 is not the forced (-1,-1,-1) column")
            phases.append(None)
        elif q in SPECIAL_POSITIVE:
            if column != (1, 1, 1):
                raise ValueError(f"q={q} is not the forced (1,1,1) column")
            phases.append(None)
        elif column == (0, 0, 0):
            phases.append(None)
        elif column in COLUMN_TO_PHASE:
            phases.append(COLUMN_TO_PHASE[column])
        else:
            raise ValueError(f"column q={q} is not an allowed normalized lift column: {column}")
    return tuple(phases)


def phases_to_eisenstein(phases: Sequence[int | None]) -> tuple[Eisenstein, ...]:
    if len(phases) != Q_MODULUS:
        raise ValueError("expected 35 phase positions")
    result: list[Eisenstein] = []
    for phase in phases:
        if phase is None:
            result.append(ZERO)
        else:
            exponent = int(phase)
            if not 0 <= exponent < 6:
                raise ValueError(f"invalid sixth-root exponent: {phase}")
            result.append(SIXTH_ROOTS[exponent])
    return tuple(result)


def d_sequence_from_grid(grid: Sequence[Sequence[int]]) -> tuple[Eisenstein, ...]:
    if len(grid) != Q_MODULUS or any(len(column) != R_MODULUS for column in grid):
        raise ValueError("expected a 35 by 3 grid")
    result = []
    for a0, a1, a2 in grid:
        result.append(int(a0) + int(a1) * OMEGA + int(a2) * OMEGA2)
    return tuple(result)


def d_sequence_from_vector(values: Sequence[int]) -> tuple[Eisenstein, ...]:
    return d_sequence_from_grid(vector_to_grid(values))


def eisenstein_periodic_autocorrelation(
    values: Sequence[Eisenstein],
) -> tuple[Eisenstein, ...]:
    row = tuple(values)
    n = len(row)
    return tuple(
        sum(
            (row[q] * row[(q + shift) % n].conjugate() for q in range(n)),
            ZERO,
        )
        for shift in range(n)
    )


def phase_autocorrelation(phases: Sequence[int | None]) -> tuple[Eisenstein, ...]:
    return eisenstein_periodic_autocorrelation(phases_to_eisenstein(phases))


def verify_phase_sequence(phases: Sequence[int | None]) -> PhaseVerificationResult:
    errors: list[str] = []
    if len(phases) != Q_MODULUS:
        raise ValueError("expected 35 phase positions")

    for q in SPECIAL_POSITIONS:
        if phases[q] is not None:
            errors.append(f"special position q={q} must be zero")

    for q, phase in enumerate(phases):
        if phase is not None and (not isinstance(phase, int) or not 0 <= phase < 6):
            errors.append(f"position q={q} has invalid phase {phase!r}")

    active_count = sum(phase is not None for phase in phases)
    if active_count != ACTIVE_PHASE_COUNT:
        errors.append(f"active phase count is {active_count}, expected {ACTIVE_PHASE_COUNT}")

    correlations = phase_autocorrelation(phases)
    target = (Eisenstein(ACTIVE_PHASE_COUNT, 0),) + (ZERO,) * (Q_MODULUS - 1)
    if correlations != target:
        bad = [
            (shift, value.as_pair())
            for shift, value in enumerate(correlations)
            if value != target[shift]
        ]
        errors.append(f"phase autocorrelations fail at {bad}")

    return PhaseVerificationResult(
        valid=not errors,
        active_count=active_count,
        correlations=correlations,
        errors=tuple(errors),
    )


def grid_correlations(values: Sequence[int]) -> tuple[tuple[int, int, int], ...]:
    """Return C(t,s)=sum_{q,r} a[q,r] a[q+t,r+s]."""
    grid = vector_to_grid(values)
    result: list[tuple[int, int, int]] = []
    for t in range(Q_MODULUS):
        by_r_shift = []
        for s in range(R_MODULUS):
            by_r_shift.append(
                sum(
                    grid[q][r] * grid[(q + t) % Q_MODULUS][(r + s) % R_MODULUS]
                    for q in range(Q_MODULUS)
                    for r in range(R_MODULUS)
                )
            )
        result.append(tuple(by_r_shift))  # type: ignore[arg-type]
    return tuple(result)


def direct_correlations_via_grid(values: Sequence[int]) -> tuple[int, ...]:
    """Recover length-105 correlations through the CRT grid formulation."""
    correlations = grid_correlations(values)
    return tuple(correlations[(12 * shift) % 35][(2 * shift) % 3] for shift in range(105))


def fourier_pair_from_grid_correlation(correlation: Sequence[int]) -> Eisenstein:
    """C0+C1*omega^2+C2*omega in the (1,omega) basis."""
    if len(correlation) != 3:
        raise ValueError("expected three r-shift correlations")
    c0, c1, c2 = (int(value) for value in correlation)
    return Eisenstein(c0 - c1, c2 - c1)


def inverse_three_point_data(total: int, fourier: Eisenstein) -> tuple[int, int, int]:
    """Recover (C0,C1,C2) from S=C0+C1+C2 and T=C0+C1*w^2+C2*w."""
    p, q = fourier.as_pair()
    numerators = (total + 2 * p - q, total - p - q, total - p + 2 * q)
    if any(value % 3 for value in numerators):
        raise ValueError("the supplied data do not invert to integer correlations")
    return tuple(value // 3 for value in numerators)  # type: ignore[return-value]


def validate_phase_fourier_identity(phases: Sequence[int | None]) -> bool:
    """Check D=(1-w)u and the exact grid/Fourier correlation identities."""
    values = reconstruct_from_phase(phases)
    grid = vector_to_grid(values)
    d_values = d_sequence_from_grid(grid)
    u_values = phases_to_eisenstein(phases)
    if d_values != tuple(ONE_MINUS_OMEGA * value for value in u_values):
        return False

    grid_corr = grid_correlations(values)
    d_corr = eisenstein_periodic_autocorrelation(d_values)
    u_corr = eisenstein_periodic_autocorrelation(u_values)
    b_grid = contract_grid(values)
    b_corr = tuple(
        sum(b_grid[q] * b_grid[(q + t) % Q_MODULUS] for q in range(Q_MODULUS))
        for t in range(Q_MODULUS)
    )

    for t in range(Q_MODULUS):
        c = grid_corr[t]
        if sum(c) != b_corr[t]:
            return False
        if fourier_pair_from_grid_correlation(c) != d_corr[t]:
            return False
        if d_corr[t] != 3 * u_corr[t]:
            return False
        if inverse_three_point_data(b_corr[t], d_corr[t]) != c:
            return False
    return True


def admissible_phase_skeleton(
    active_positions: Iterable[int], exponents: Iterable[int]
) -> tuple[int | None, ...]:
    """Convenience constructor used by tests and small exact experiments."""
    positions = tuple(int(q) for q in active_positions)
    phases = tuple(int(k) for k in exponents)
    if len(positions) != len(phases):
        raise ValueError("positions and exponents must have equal length")
    result: list[int | None] = [None] * Q_MODULUS
    for q, exponent in zip(positions, phases, strict=True):
        if not 0 <= q < Q_MODULUS:
            raise ValueError(f"position outside Z_35: {q}")
        if q in SPECIAL_POSITIONS:
            raise ValueError(f"special position cannot carry a phase: {q}")
        if result[q] is not None:
            raise ValueError(f"duplicate position: {q}")
        if not 0 <= exponent < 6:
            raise ValueError(f"invalid phase exponent: {exponent}")
        result[q] = exponent
    return tuple(result)
