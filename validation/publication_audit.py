#!/usr/bin/env python3
"""Independent exact audit of the focused publication repository.

This module deliberately imports no project proof code.  It checks the small
machine-readable proof specification, the complete support certificate, the
coordinate conventions, and the publication manifest using only the Python
standard library.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path, PurePosixPath
from typing import Iterable, Sequence


N = 35
WEIGHT = 12
EXPECTED_ORIGINAL_CONTRACTION = {0: -3, 5: 3, 10: 3, 20: 3}
EXPECTED_GRID_CONTRACTION = {0: -3, 15: 3, 25: 3, 30: 3}
EXPECTED_ORIGINAL_CONTRACTIONS = (
    EXPECTED_ORIGINAL_CONTRACTION,
    EXPECTED_GRID_CONTRACTION,
)
EXPECTED_GRID_CONTRACTIONS = (
    EXPECTED_GRID_CONTRACTION,
    EXPECTED_ORIGINAL_CONTRACTION,
)
EXPECTED_MUTATION_BOUNDARIES = [
    "factor_branch",
    "support_certificate",
    "affine_representative_and_stabilizer",
    "four_shift_terms",
    "fibre_to_phase_table",
    "contraction_conventions",
]


class AuditError(ValueError):
    """Raised when a publication invariant fails."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise AuditError(message)


def _trim_f3(polynomial: Iterable[int]) -> tuple[int, ...]:
    values = [int(value) % 3 for value in polynomial]
    while len(values) > 1 and values[-1] == 0:
        values.pop()
    return tuple(values or [0])


def _multiply_f3(
    left: Sequence[int], right: Sequence[int]
) -> tuple[int, ...]:
    product = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            product[i + j] = (product[i + j] + int(a) * int(b)) % 3
    return _trim_f3(product)


def _monic_f3(polynomial: Sequence[int]) -> tuple[int, ...]:
    values = _trim_f3(polynomial)
    _require(values != (0,), "zero polynomial cannot be normalized")
    inverse = 1 if values[-1] == 1 else 2
    return tuple((inverse * value) % 3 for value in values)


def _divmod_f3(
    numerator: Sequence[int], denominator: Sequence[int]
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    remainder = list(_trim_f3(numerator))
    divisor = _trim_f3(denominator)
    _require(divisor != (0,), "polynomial division by zero")
    quotient = [0] * max(1, len(remainder) - len(divisor) + 1)
    inverse = 1 if divisor[-1] == 1 else 2
    while remainder != [0] and len(remainder) >= len(divisor):
        shift = len(remainder) - len(divisor)
        coefficient = remainder[-1] * inverse % 3
        quotient[shift] = coefficient
        for index, value in enumerate(divisor):
            remainder[shift + index] = (
                remainder[shift + index] - coefficient * value
            ) % 3
        while len(remainder) > 1 and remainder[-1] == 0:
            remainder.pop()
    return _trim_f3(quotient), _trim_f3(remainder)


def _gcd_f3(left: Sequence[int], right: Sequence[int]) -> tuple[int, ...]:
    a = _trim_f3(left)
    b = _trim_f3(right)
    while b != (0,):
        _, remainder = _divmod_f3(a, b)
        a, b = b, remainder
    return _monic_f3(a)


def _multiply_mod_f3(
    left: Sequence[int], right: Sequence[int], modulus: Sequence[int]
) -> tuple[int, ...]:
    return _divmod_f3(_multiply_f3(left, right), modulus)[1]


def _pow_mod_f3(
    base: Sequence[int], exponent: int, modulus: Sequence[int]
) -> tuple[int, ...]:
    _require(exponent >= 0, "negative polynomial exponent")
    result: tuple[int, ...] = (1,)
    power = _divmod_f3(base, modulus)[1]
    remaining = exponent
    while remaining:
        if remaining & 1:
            result = _multiply_mod_f3(result, power, modulus)
        power = _multiply_mod_f3(power, power, modulus)
        remaining >>= 1
    return result


def _subtract_f3(
    left: Sequence[int], right: Sequence[int]
) -> tuple[int, ...]:
    length = max(len(left), len(right))
    return _trim_f3(
        (int(left[index]) if index < len(left) else 0)
        - (int(right[index]) if index < len(right) else 0)
        for index in range(length)
    )


def _is_irreducible_f3(polynomial: Sequence[int]) -> bool:
    """Apply the exact finite-field Rabin irreducibility criterion."""
    factor = _monic_f3(polynomial)
    degree = len(factor) - 1
    if degree < 1:
        return False
    prime_divisors = {
        candidate
        for candidate in range(2, degree + 1)
        if degree % candidate == 0
        and all(candidate % divisor for divisor in range(2, int(math.sqrt(candidate)) + 1))
    }
    x_polynomial = (0, 1)
    for prime in prime_divisors:
        frobenius = _pow_mod_f3(
            x_polynomial, 3 ** (degree // prime), factor
        )
        if _gcd_f3(factor, _subtract_f3(frobenius, x_polynomial)) != (1,):
            return False
    final_frobenius = _pow_mod_f3(x_polynomial, 3**degree, factor)
    final_difference = _subtract_f3(final_frobenius, x_polynomial)
    return _divmod_f3(final_difference, factor)[1] == (0,)


def _periodic_f3_correlation(word: Sequence[int]) -> tuple[int, ...]:
    _require(len(word) == N, "explicit ternary word must have length 35")
    return tuple(
        sum(int(word[q]) * int(word[(q + shift) % N]) for q in range(N)) % 3
        for shift in range(N)
    )


def _support_from_mask(mask: int) -> tuple[int, ...]:
    _require(0 <= int(mask) < 1 << N, "support mask is outside 35 bits")
    return tuple(q for q in range(N) if (int(mask) >> q) & 1)


def _support_mask(points: Iterable[int]) -> int:
    mask = 0
    for q in points:
        _require(0 <= int(q) < N, "support coordinate is outside Z_35")
        mask |= 1 << int(q)
    return mask


def _affine_image(
    points: Sequence[int], multiplier: int, translation: int
) -> tuple[int, ...]:
    _require(math.gcd(int(multiplier), N) == 1, "affine multiplier is not a unit")
    return tuple(
        sorted({(int(multiplier) * int(q) + int(translation)) % N for q in points})
    )


def _eisenstein_multiply(
    left: tuple[int, int], right: tuple[int, int]
) -> tuple[int, int]:
    """Multiply a+b*omega using omega^2=-1-omega."""
    a, b = left
    c, d = right
    return a * c - b * d, a * d + b * c - b * d


def _ratio_monomial(pair: tuple[int, int]) -> tuple[int, ...]:
    numerator, denominator = pair
    result = [0] * N
    result[numerator] += 1
    result[denominator] -= 1
    return tuple(result)


def _add_monomials(
    left: Sequence[int], right: Sequence[int]
) -> tuple[int, ...]:
    return tuple(int(a) + int(b) for a, b in zip(left, right, strict=True))


def _subtract_monomials(
    left: Sequence[int], right: Sequence[int]
) -> tuple[int, ...]:
    return tuple(int(a) - int(b) for a, b in zip(left, right, strict=True))


def _verify_four_shift_symbolic_core(
    terms: dict[int, tuple[tuple[int, int], ...]],
) -> None:
    """Verify the short characteristic-two contradiction as formal ratios."""
    c2 = terms[2]
    c4 = terms[4]
    c8 = terms[8]
    c10 = terms[10]
    _require(
        (len(c2), len(c4), len(c8), len(c10)) == (4, 3, 2, 6),
        "four-shift term counts do not support the symbolic core",
    )

    p_representatives = tuple(_ratio_monomial(pair) for pair in c8)
    c10_monomials = tuple(_ratio_monomial(pair) for pair in c10)
    used_targets: set[int] = set()
    for term in c2:
        term_monomial = _ratio_monomial(term)
        matches = [
            index
            for p in p_representatives
            for index, target in enumerate(c10_monomials)
            if index not in used_targets
            and _add_monomials(p, term_monomial) == target
        ]
        _require(matches, "p*C2 does not match four distinct terms of C10")
        used_targets.add(matches[0])
    _require(len(used_targets) == 4, "p*C2 does not select four C10 terms")

    leftover = [
        c10_monomials[index]
        for index in range(len(c10_monomials))
        if index not in used_targets
    ]
    _require(len(leftover) == 2, "C10 does not leave exactly two terms")
    leftover_relation = _subtract_monomials(leftover[0], leftover[1])

    cancelling_pair: tuple[int, int] | None = None
    c4_monomials = tuple(_ratio_monomial(pair) for pair in c4)
    for first in range(3):
        for second in range(first + 1, 3):
            relation = _subtract_monomials(
                c4_monomials[first], c4_monomials[second]
            )
            if relation == leftover_relation or tuple(-x for x in relation) == leftover_relation:
                cancelling_pair = (first, second)
    _require(
        cancelling_pair is not None,
        "the two remaining C10 terms do not force a C4 cancellation",
    )
    remaining = set(range(3)) - set(cancelling_pair)
    _require(len(remaining) == 1, "C4 does not leave one nonzero ratio")


def _parse_certificate(path: Path) -> tuple[int, ...]:
    masks = tuple(
        int(line)
        for raw in path.read_text(encoding="ascii").splitlines()
        if (line := raw.strip()) and not line.startswith("#")
    )
    _require(masks == tuple(sorted(masks)), "support certificate is not sorted")
    _require(len(masks) == len(set(masks)) == 420, "support certificate must have 420 unique masks")
    for mask in masks:
        _require(len(_support_from_mask(mask)) == WEIGHT, "certificate mask does not have weight 12")
    return masks


def _parse_manifest(path: Path) -> dict[str, str]:
    entries: dict[str, str] = {}
    previous = ""
    for line in path.read_text(encoding="ascii").splitlines():
        parts = line.split("  ", 1)
        _require(len(parts) == 2, "malformed publication manifest line")
        digest, name = parts
        _require(len(digest) == 64 and all(c in "0123456789abcdef" for c in digest), "invalid SHA-256 digest")
        manifest_path = PurePosixPath(name)
        _require(
            name
            and "\\" not in name
            and not manifest_path.is_absolute()
            and ".." not in manifest_path.parts
            and str(manifest_path) == name,
            "manifest entry is not a safe normalized relative path",
        )
        _require(name > previous, "publication manifest entries are not strictly sorted")
        _require(name not in entries, "duplicate publication manifest entry")
        entries[name] = digest
        previous = name
    return entries


def audit_repository(root: Path) -> dict[str, object]:
    """Audit the canonical publication files and return a deterministic report."""
    root = Path(root)
    _require(root.is_dir(), f"publication repository does not exist: {root}")
    spec_path = root / "data" / "proof-spec.json"
    certificate_path = root / "data" / "mod3_weight12_support_masks.txt"
    manifest_path = root / "checksums.sha256"
    for required in (spec_path, certificate_path, manifest_path):
        _require(required.is_file(), f"missing publication file: {required.name}")

    spec = json.loads(spec_path.read_text(encoding="utf-8"))
    _require(spec.get("schema") == "cw105-proof-spec-v1", "unexpected proof-spec schema")

    support_spec = spec["support_classification"]
    canonical = tuple(int(q) for q in support_spec["canonical_support"])
    _require(len(canonical) == len(set(canonical)) == WEIGHT, "canonical support must contain 12 distinct coordinates")
    _require(canonical == tuple(sorted(canonical)), "canonical support must be sorted")
    _require(support_spec["orbit_size"] == 420, "declared affine orbit size is incorrect")
    _require(
        support_spec["coordinate_convention"] == "bit q is coefficient q in Z_35",
        "incorrect support serialization convention",
    )
    units = tuple(a for a in range(N) if math.gcd(a, N) == 1)
    orbit = {
        _affine_image(canonical, a, b)
        for a in units
        for b in range(N)
    }
    _require(len(orbit) == 420, "canonical affine orbit does not have size 420")
    _require(min(orbit) == canonical, "declared support is not the canonical orbit representative")
    stabilizer = tuple(
        (a, b)
        for a in units
        for b in range(N)
        if _affine_image(canonical, a, b) == canonical
    )
    declared_stabilizer = tuple(
        (int(pair[0]), int(pair[1])) for pair in support_spec["affine_stabilizer"]
    )
    _require(stabilizer == declared_stabilizer, "declared affine stabilizer is incorrect")

    masks = _parse_certificate(certificate_path)
    orbit_masks = tuple(sorted(_support_mask(points) for points in orbit))
    _require(masks == orbit_masks, "support certificate is not the complete canonical affine orbit")

    affine_correlation_shift_checks = 0
    canonical_set = set(canonical)
    for a in units:
        for b in range(N):
            image_set = set(_affine_image(canonical, a, b))
            for shift in range(N):
                mapped_pairs = {
                    ((a * q + b) % N, (a * (q + shift) + b) % N)
                    for q in canonical
                    if (q + shift) % N in canonical_set
                }
                image_shift = (a * shift) % N
                image_pairs = {
                    (q, (q + image_shift) % N)
                    for q in image_set
                    if (q + image_shift) % N in image_set
                }
                _require(
                    mapped_pairs == image_pairs,
                    "affine relabelling does not preserve a correlation term list",
                )
                affine_correlation_shift_checks += 1

    f3 = spec["characteristic_three"]
    self_factors = tuple(tuple(int(v) for v in factor) for factor in f3["self_reciprocal_factors"])
    reciprocal_pair = tuple(tuple(int(v) for v in factor) for factor in f3["reciprocal_pair"])
    _require(len(self_factors) == 3 and len(reciprocal_pair) == 2, "unexpected F3 factor-orbit structure")
    _require(
        all(
            _monic_f3(tuple(reversed(factor))) == _monic_f3(factor)
            for factor in self_factors
        ),
        "declared F3 fixed factor is not self-reciprocal",
    )
    _require(
        _monic_f3(tuple(reversed(reciprocal_pair[0])))
        == _monic_f3(reciprocal_pair[1]),
        "declared F3 pair is not reciprocal",
    )
    _require(
        all(_is_irreducible_f3(factor) for factor in (*self_factors, *reciprocal_pair)),
        "declared F3 factor is reducible",
    )
    _require(
        len(set((*self_factors, *reciprocal_pair))) == 5,
        "declared F3 factorization is not square-free",
    )
    product: tuple[int, ...] = (1,)
    for factor in (*self_factors, *reciprocal_pair):
        product = _multiply_f3(product, factor)
    _require(product == (2,) + (0,) * 34 + (1,), "declared factors do not multiply to X^35-1 over F3")

    branch_generators = []
    for reciprocal_factor in reciprocal_pair:
        branch_generator: tuple[int, ...] = (1,)
        for factor in (*self_factors, reciprocal_factor):
            branch_generator = _multiply_f3(branch_generator, factor)
        _require(
            len(branch_generator) - 1 == 23,
            "factor-orbit branch does not generate a dimension-12 cyclic code",
        )
        branch_generators.append(branch_generator)
    _require(
        len(set(branch_generators)) == 2,
        "factor-orbit partition does not give exactly two distinct branches",
    )

    branch_index = int(f3["explicit_branch_index"])
    _require(branch_index in (0, 1), "explicit factor branch index must be 0 or 1")
    generator: tuple[int, ...] = (1,)
    for factor in (*self_factors, reciprocal_pair[branch_index]):
        generator = _multiply_f3(generator, factor)
    declared_generator = tuple(int(v) for v in f3["explicit_generator"])
    _require(generator == declared_generator, "declared explicit factor branch generator is incorrect")

    message = tuple(int(v) for v in f3["explicit_message"])
    _require(len(message) == 12, "explicit ternary message must have length 12")
    word_polynomial = _multiply_f3(message, generator)
    word = word_polynomial + (0,) * (N - len(word_polynomial))
    declared_word = tuple(int(v) for v in f3["explicit_word"])
    _require(word == declared_word, "explicit ternary word is not message times branch generator")
    word_support = tuple(q for q, value in enumerate(word) if value % 3)
    _require(word_support == canonical, "explicit ternary word does not anchor the canonical support")
    _require(_periodic_f3_correlation(word) == (0,) * N, "explicit ternary word fails direct F3 autocorrelation")

    f4 = spec["characteristic_two"]
    core_shifts = tuple(int(value) for value in f4["core_shifts"])
    _require(core_shifts == (2, 4, 8, 10), "unexpected four-shift set")
    support_set = set(canonical)
    regenerated_terms: dict[int, tuple[tuple[int, int], ...]] = {}
    for shift in core_shifts:
        expected_terms = tuple(
            (q, (q + shift) % N)
            for q in canonical
            if (q + shift) % N in support_set
        )
        declared_terms = tuple(
            (int(pair[0]), int(pair[1]))
            for pair in f4["correlation_terms"][str(shift)]
        )
        _require(declared_terms == expected_terms, f"incorrect correlation terms for shift {shift}")
        regenerated_terms[shift] = expected_terms
    _verify_four_shift_symbolic_core(regenerated_terms)

    fibre = spec["fibre_phase_convention"]
    table = tuple(tuple(int(v) for v in row) for row in fibre["phase_to_column"])
    _require(len(table) == 6 and len(set(table)) == 6, "phase table must contain six distinct rows")
    sixth_roots = ((1, 0), (1, 1), (0, 1), (-1, 0), (-1, -1), (0, -1))
    for exponent, column in enumerate(table):
        _require(tuple(sorted(column)) == (-1, 0, 1), "active fibre column is not a permutation of (-1,0,1)")
        a0, a1, a2 = column
        fourier_value = (a0 - a2, a1 - a2)
        expected = _eisenstein_multiply((1, -1), sixth_roots[exponent])
        _require(fourier_value == expected, f"incorrect fibre-to-phase row at exponent {exponent}")

    crt = fibre["crt_index"]
    modulus = int(crt["modulus"])
    q_modulus = int(crt["q_modulus"])
    r_modulus = int(crt["r_modulus"])
    q_coefficient = int(crt["q_coefficient"])
    r_coefficient = int(crt["r_coefficient"])
    _require((modulus, q_modulus, r_modulus) == (105, 35, 3), "incorrect CRT moduli")
    _require(
        (q_coefficient, r_coefficient) == (3, 35),
        "incorrect CRT coordinate coefficients",
    )
    indices = {
        (q_coefficient * q + r_coefficient * r) % modulus
        for q in range(q_modulus)
        for r in range(r_modulus)
    }
    _require(indices == set(range(modulus)), "declared CRT index map is not a bijection")
    original = {int(key): int(value) for key, value in fibre["normalized_contraction_original"].items()}
    grid = {int(key): int(value) for key, value in fibre["normalized_contraction_grid"].items()}
    _require(original == EXPECTED_ORIGINAL_CONTRACTION, "incorrect original contraction convention")
    _require(grid == EXPECTED_GRID_CONTRACTION, "incorrect grid contraction convention")
    original_row = tuple(original.get(q, 0) for q in range(q_modulus))
    grid_row = tuple(grid.get(q, 0) for q in range(q_modulus))
    _require(
        grid_row == tuple(original_row[(q_coefficient * q) % q_modulus] for q in range(q_modulus)),
        "CRT coordinate map does not send the original contraction to the grid contraction",
    )
    original_rows = tuple(
        {int(key): int(value) for key, value in row.items()}
        for row in fibre["multiplier_fixed_contractions_original"]
    )
    grid_rows = tuple(
        {int(key): int(value) for key, value in row.items()}
        for row in fibre["multiplier_fixed_contractions_grid"]
    )
    _require(
        original_rows == EXPECTED_ORIGINAL_CONTRACTIONS,
        "incorrect list of multiplier-fixed contractions",
    )
    _require(
        grid_rows == EXPECTED_GRID_CONTRACTIONS,
        "incorrect CRT images of multiplier-fixed contractions",
    )
    for original_sparse, grid_sparse in zip(original_rows, grid_rows, strict=True):
        original_full = tuple(
            original_sparse.get(q, 0) for q in range(q_modulus)
        )
        grid_full = tuple(grid_sparse.get(q, 0) for q in range(q_modulus))
        _require(
            grid_full
            == tuple(
                original_full[(q_coefficient * q) % q_modulus]
                for q in range(q_modulus)
            ),
            "a multiplier-fixed contraction has the wrong CRT image",
        )
        _require(
            all(
                original_full[(4 * q) % q_modulus] == original_full[q]
                for q in range(q_modulus)
            ),
            "declared contraction is not fixed by multiplier 4",
        )
        correlation = tuple(
            sum(
                original_full[q] * original_full[(q + shift) % q_modulus]
                for q in range(q_modulus)
            )
            for shift in range(q_modulus)
        )
        _require(
            correlation == (36,) + (0,) * 34,
            "declared contraction fails direct periodic correlation",
        )
        forced_columns = sum(value != 0 for value in grid_full)
        _require(
            forced_columns == 4 and (36 - 3 * forced_columns) // 2 == 12,
            "declared contraction does not force four constant and twelve active columns",
        )

    manifest = _parse_manifest(manifest_path)
    for name, expected_digest in manifest.items():
        target = root / PurePosixPath(name)
        _require(target.is_file(), f"manifest target is missing: {name}")
        actual = hashlib.sha256(target.read_bytes()).hexdigest()
        _require(actual == expected_digest, f"publication checksum mismatch: {name}")

    return {
        "status": "passed",
        "support_count": len(masks),
        "affine_orbit_size": len(orbit),
        "affine_stabilizer": [list(pair) for pair in stabilizer],
        "affine_correlation_shift_checks": affine_correlation_shift_checks,
        "factor_orbit_branches": len(branch_generators),
        "explicit_f3_word_support": list(word_support),
        "core_shifts": list(core_shifts),
        "four_shift_symbolic_core": "passed",
        "contraction_rows_checked": len(original_rows),
        "manifest_entries": len(manifest),
        "mutation_boundaries_checked": EXPECTED_MUTATION_BOUNDARIES,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("repository", type=Path, nargs="?", default=Path("."))
    args = parser.parse_args()
    print(json.dumps(audit_repository(args.repository), sort_keys=True, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
