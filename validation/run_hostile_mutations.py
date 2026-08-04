#!/usr/bin/env python3
"""Run semantic hostile mutations against the publication audit."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import tempfile
from pathlib import Path
from typing import Callable

if __package__:
    from .publication_audit import AuditError, audit_repository
else:
    from publication_audit import AuditError, audit_repository


Mutation = Callable[[Path], None]


def _load_spec(root: Path) -> dict[str, object]:
    return json.loads((root / "data" / "proof-spec.json").read_text(encoding="utf-8"))


def _write_spec(root: Path, spec: dict[str, object]) -> None:
    (root / "data" / "proof-spec.json").write_text(
        json.dumps(spec, sort_keys=True, indent=2) + "\n",
        encoding="utf-8",
        newline="\n",
    )


def _manifest_names(root: Path) -> tuple[str, ...]:
    return tuple(
        line.split("  ", 1)[1]
        for line in (root / "checksums.sha256").read_text(encoding="ascii").splitlines()
    )


def _rewrite_manifest(root: Path) -> None:
    names = _manifest_names(root)
    (root / "checksums.sha256").write_text(
        "".join(
            f"{hashlib.sha256((root / name).read_bytes()).hexdigest()}  {name}\n"
            for name in names
        ),
        encoding="ascii",
        newline="\n",
    )


def _mutate_factor_branch(root: Path) -> None:
    spec = _load_spec(root)
    f3 = spec["characteristic_three"]
    assert isinstance(f3, dict)
    generator = f3["explicit_generator"]
    assert isinstance(generator, list)
    generator[0] = (int(generator[0]) + 1) % 3
    _write_spec(root, spec)


def _mutate_support_mask(root: Path) -> None:
    path = root / "data" / "mod3_weight12_support_masks.txt"
    lines = path.read_text(encoding="ascii").splitlines()
    comments = [line for line in lines if line.startswith("#")]
    masks = sorted(int(line) for line in lines if line and not line.startswith("#"))
    original = masks[0]
    present = [q for q in range(35) if (original >> q) & 1]
    absent = [q for q in range(35) if not ((original >> q) & 1)]
    replacement = None
    for removed in present:
        for added in absent:
            candidate = original ^ (1 << removed) ^ (1 << added)
            if candidate not in masks:
                replacement = candidate
                break
        if replacement is not None:
            break
    assert replacement is not None
    masks[0] = replacement
    masks.sort()
    path.write_text(
        "\n".join((*comments, *(str(mask) for mask in masks))) + "\n",
        encoding="ascii",
        newline="\n",
    )


def _mutate_affine_representative(root: Path) -> None:
    spec = _load_spec(root)
    classification = spec["support_classification"]
    assert isinstance(classification, dict)
    canonical = classification["canonical_support"]
    assert isinstance(canonical, list)
    canonical[-1] = 29
    _write_spec(root, spec)


def _mutate_affine_stabilizer(root: Path) -> None:
    spec = _load_spec(root)
    classification = spec["support_classification"]
    assert isinstance(classification, dict)
    stabilizer = classification["affine_stabilizer"]
    assert isinstance(stabilizer, list)
    assert isinstance(stabilizer[1], list)
    stabilizer[1][1] = 27
    _write_spec(root, spec)


def _mutate_four_shift_coefficient(root: Path) -> None:
    spec = _load_spec(root)
    characteristic_two = spec["characteristic_two"]
    assert isinstance(characteristic_two, dict)
    terms = characteristic_two["correlation_terms"]
    assert isinstance(terms, dict)
    shift_two = terms["2"]
    assert isinstance(shift_two, list)
    assert isinstance(shift_two[0], list)
    shift_two[0][1] = 3
    _write_spec(root, spec)


def _mutate_fibre_to_phase_table(root: Path) -> None:
    spec = _load_spec(root)
    fibre = spec["fibre_phase_convention"]
    assert isinstance(fibre, dict)
    table = fibre["phase_to_column"]
    assert isinstance(table, list)
    table[0] = [1, 0, -1]
    _write_spec(root, spec)


def _mutate_contraction_convention(root: Path) -> None:
    spec = _load_spec(root)
    fibre = spec["fibre_phase_convention"]
    assert isinstance(fibre, dict)
    contraction = fibre["normalized_contraction_original"]
    assert isinstance(contraction, dict)
    contraction["6"] = contraction.pop("5")
    _write_spec(root, spec)


MUTATIONS: tuple[tuple[str, str, Mutation], ...] = (
    ("factor_branch", "factor branch generator", _mutate_factor_branch),
    ("support_mask", "support certificate", _mutate_support_mask),
    ("affine_representative", "canonical", _mutate_affine_representative),
    ("affine_stabilizer", "affine stabilizer", _mutate_affine_stabilizer),
    ("four_shift_coefficient", "correlation terms", _mutate_four_shift_coefficient),
    ("fibre_to_phase_table", "phase", _mutate_fibre_to_phase_table),
    ("contraction_convention", "contraction", _mutate_contraction_convention),
)


def run_mutation_suite(repository: Path) -> dict[str, object]:
    repository = Path(repository).resolve()
    if not repository.is_dir():
        raise FileNotFoundError(repository)
    manifest_names = _manifest_names(repository)
    cases: list[dict[str, object]] = []
    with tempfile.TemporaryDirectory(prefix="cw105-hostile-mutations-") as raw_temp:
        temp_root = Path(raw_temp).resolve()
        for index, (name, expected_error, mutate) in enumerate(MUTATIONS):
            mutated = (temp_root / f"{index:02d}-{name}").resolve()
            if temp_root not in mutated.parents:
                raise RuntimeError("mutation path escaped its temporary root")
            mutated.mkdir()
            for relative in (*manifest_names, "checksums.sha256"):
                source = repository / relative
                target = mutated / relative
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(source, target)
            mutate(mutated)
            _rewrite_manifest(mutated)
            caught = False
            try:
                audit_repository(mutated)
            except AuditError as error:
                caught = expected_error in str(error)
            cases.append(
                {
                    "name": name,
                    "caught": caught,
                    "error_contains": expected_error,
                }
            )
    return {
        "status": "passed" if all(bool(case["caught"]) for case in cases) else "failed",
        "case_count": len(cases),
        "cases": cases,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("repository", type=Path, nargs="?", default=Path("."))
    args = parser.parse_args()
    report = run_mutation_suite(args.repository)
    print(json.dumps(report, sort_keys=True, indent=2))
    return 0 if report["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
