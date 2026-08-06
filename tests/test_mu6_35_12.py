from __future__ import annotations

import json
import shutil
import subprocess
from collections import Counter
from pathlib import Path

import pytest

from src.eisenstein_congruences import mod3_periodic_autocorrelation
from src.mu6_35_12 import (
    CANONICAL_SUPPORT,
    CORE_SHIFTS,
    affine_stabilizer,
    affine_support_orbit,
    enumerate_all_mod3_weight12_candidates_numpy,
    exhaustive_four_shift_core_check,
    mu6_35_12_nonexistence_proof,
    mu6_result_summary,
)


def read_masks(path: Path) -> tuple[int, ...]:
    return tuple(
        int(line)
        for line in path.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("#")
    )


def test_generalized_mod3_support_classification_is_complete() -> None:
    candidates_by_code = enumerate_all_mod3_weight12_candidates_numpy()
    assert tuple(len(code) for code in candidates_by_code) == (420, 420)

    supports_by_code = tuple(
        {candidate.support for candidate in code} for code in candidates_by_code
    )
    assert tuple(len(values) for values in supports_by_code) == (210, 210)
    assert supports_by_code[0].isdisjoint(supports_by_code[1])

    for candidates in candidates_by_code:
        multiplicities = Counter(candidate.support for candidate in candidates)
        assert set(multiplicities.values()) == {2}
        for candidate in candidates:
            assert len(candidate.support) == 12
            assert mod3_periodic_autocorrelation(candidate.sequence) == (0,) * 35

    supports = supports_by_code[0] | supports_by_code[1]
    orbit = set(affine_support_orbit(CANONICAL_SUPPORT))
    assert len(supports) == len(orbit) == 420
    assert supports == orbit
    assert min(orbit) == CANONICAL_SUPPORT
    assert affine_stabilizer(CANONICAL_SUPPORT) == ((1, 0), (29, 28))


def test_four_shift_hand_identity_is_exhaustively_regression_checked() -> None:
    result = exhaustive_four_shift_core_check()
    assert CORE_SHIFTS == (2, 4, 8, 10)
    assert result.normalized_assignments == 3**11 == 177_147
    assert result.c8_zero == 59_049
    assert result.c2_c8_zero == 15_309
    assert result.c2_c8_c10_zero == 5_103
    assert result.all_four_zero == 0
    assert result.algebraic_identity_failures == 0


def test_generalized_theorem_summary() -> None:
    result = mu6_35_12_nonexistence_proof()
    summary = mu6_result_summary(result)
    assert summary["uses_prescribed_zero_positions"] is False
    assert summary["matrix_form"] == "No circulant CGW(35,12;6) exists."
    assert summary["mod3_candidates_per_code"] == [420, 420]
    assert summary["mod3_supports_per_code"] == [210, 210]
    assert summary["mod3_support_count"] == 420
    assert summary["affine_orbit_count"] == 1
    assert summary["affine_orbit_size"] == 420
    assert summary["affine_stabilizer"] == [[1, 0], [29, 28]]
    assert summary["canonical_support"] == list(CANONICAL_SUPPORT)
    assert summary["core_counts"]["all_four_zero"] == 0
    assert summary["conclusion"] == "nonexistence"


@pytest.mark.skipif(shutil.which("g++") is None, reason="g++ is not installed")
def test_independent_cpp_generalized_checker_matches_python(tmp_path: Path) -> None:
    executable = tmp_path / "mu6_independent"
    subprocess.run(
        [
            "g++",
            "-std=c++20",
            "-O3",
            "-Wall",
            "-Wextra",
            "-Werror",
            "-pedantic",
            "src/mu6_35_12_independent.cpp",
            "-o",
            str(executable),
        ],
        check=True,
    )
    completed = subprocess.run(
        [str(executable)], check=True, capture_output=True, text=True
    )
    cpp = json.loads(completed.stdout)
    python = mu6_result_summary(mu6_35_12_nonexistence_proof())
    keys = (
        "uses_prescribed_zero_positions",
        "f3_factorization_valid",
        "mod3_candidates_per_code",
        "mod3_supports_per_code",
        "mod3_support_intersection_size",
        "mod3_support_count",
        "affine_orbit_size",
        "affine_stabilizer",
        "canonical_support",
        "core_shifts",
        "core_counts",
        "support_masks",
    )
    assert {key: cpp[key] for key in keys} == {key: python[key] for key in keys}


@pytest.mark.skipif(shutil.which("g++") is None, reason="g++ is not installed")
def test_direct_all_supports_mod2_enumerator_finds_no_assignment(tmp_path: Path) -> None:
    certificate = Path("data/mod3_weight12_support_masks.txt")
    masks = read_masks(certificate)
    python_masks = tuple(
        mu6_result_summary(mu6_35_12_nonexistence_proof())["support_masks"]
    )
    assert masks == python_masks
    assert len(masks) == len(set(masks)) == 420

    executable = tmp_path / "mu6_direct_all"
    subprocess.run(
        [
            "g++",
            "-std=c++20",
            "-O3",
            "-Wall",
            "-Wextra",
            "-Werror",
            "-pedantic",
            "src/mu6_35_12_direct_all.cpp",
            "-o",
            str(executable),
        ],
        check=True,
    )
    crlf_certificate = tmp_path / "support-masks-crlf.txt"
    crlf_certificate.write_bytes(
        b"\r\n".join(
            line.encode("ascii") for line in certificate.read_text().splitlines()
        )
        + b"\r\n"
    )
    completed = subprocess.run(
        [str(executable), str(crlf_certificate)],
        check=True,
        capture_output=True,
        text=True,
    )
    result = json.loads(completed.stdout)
    assert result["support_count"] == 420
    assert result["uses_prescribed_zero_positions"] is False
    assert result["global_phase_fixed"] is True
    assert result["assignments_per_support"] == 3**11
    assert result["assignments_checked"] == 420 * 3**11 == 74_401_740
    assert result["compatible_support_count"] == 0
    assert result["zero_correlation_assignments_fixed_phase"] == 0
    assert sum(result["first_failure_shift_counts"]) == result["assignments_checked"]
