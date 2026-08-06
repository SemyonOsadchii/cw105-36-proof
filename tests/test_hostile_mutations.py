from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys

from validation.run_hostile_mutations import run_mutation_suite


ROOT = Path(__file__).resolve().parents[1]


def test_every_hostile_publication_mutation_is_rejected() -> None:
    report = run_mutation_suite(ROOT)
    assert report["status"] == "passed"
    assert report["case_count"] == 7
    assert report["cases"] == [
        {
            "name": "factor_branch",
            "caught": True,
            "error_contains": "factor branch generator",
        },
        {
            "name": "support_mask",
            "caught": True,
            "error_contains": "support certificate",
        },
        {
            "name": "affine_representative",
            "caught": True,
            "error_contains": "canonical",
        },
        {
            "name": "affine_stabilizer",
            "caught": True,
            "error_contains": "affine stabilizer",
        },
        {
            "name": "four_shift_coefficient",
            "caught": True,
            "error_contains": "correlation terms",
        },
        {
            "name": "fibre_to_phase_table",
            "caught": True,
            "error_contains": "phase",
        },
        {
            "name": "contraction_convention",
            "caught": True,
            "error_contains": "contraction",
        },
    ]


def test_hostile_mutation_suite_runs_as_a_standalone_cli() -> None:
    completed = subprocess.run(
        [sys.executable, "validation/run_hostile_mutations.py", "."],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    report = json.loads(completed.stdout)
    assert report["status"] == "passed"
    assert report["case_count"] == 7
