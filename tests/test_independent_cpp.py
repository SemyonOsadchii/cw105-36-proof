from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path

import pytest

from src.derive_icw import enumerate_icw

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "src" / "derive_icw_independent.cpp"


def test_standalone_cpp_base7_sweep_matches_python_exactly(tmp_path: Path) -> None:
    compiler = shutil.which("g++")
    if compiler is None:
        pytest.skip("g++ is unavailable")

    binary = tmp_path / "derive_icw_independent"
    subprocess.run(
        [compiler, "-std=c++20", "-O3", "-DNDEBUG", str(SOURCE), "-o", str(binary)],
        check=True,
        timeout=60,
    )
    completed = subprocess.run(
        [str(binary)],
        check=True,
        text=True,
        capture_output=True,
        timeout=60,
    )
    cpp = json.loads(completed.stdout)
    python = enumerate_icw()

    assert cpp["total_assignments"] == 7**9
    assert cpp["scalar_candidate_count"] == 1434
    assert {tuple(row) for row in cpp["scalar_candidates"]} == set(
        python.scalar_candidates_dfs
    )
    assert cpp["invariant_solution_count"] == 2
    assert {tuple(row) for row in cpp["invariant_solutions"]} == set(
        python.invariant_solutions
    )
    assert cpp["equivalent_by_unit_3"] is True
