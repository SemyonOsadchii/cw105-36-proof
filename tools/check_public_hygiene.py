"""Reject private state, publication-irrelevant paths, and local machine data."""

from __future__ import annotations

import subprocess
from pathlib import Path


ROOT = Path(__file__).parents[1]
FORBIDDEN_TRACKED_PATHS = {
    "AGENTS.md",
    "CONTRIBUTING.md",
    "STATUS.md",
    "VALIDATION_REPORT.md",
    "WORKSPACE.md",
}
FORBIDDEN_TRACKED_PREFIXES = (
    ".codex/",
    "docs/superpowers/",
    "experiments/",
    "logs/",
    "proofs/",
    "release/",
    "results/",
    "validation/blind/",
)
FORBIDDEN_TEXT = (
    b"/mnt/c/" + b"Users/",
    b"C:\\Users\\",
    b"C:/" + b"Users/",
    b".config/" + b"superpowers",
    b"AppData/" + b"Local/Temp",
)
REQUIRED_IGNORE_RULES = {
    ".codex/",
    ".claude/",
    ".cursor/",
    ".continue/",
    ".aider*",
    "docs/superpowers/",
    "AGENTS.md",
    "CLAUDE.md",
    "GEMINI.md",
    "WORKSPACE.md",
}
TEXT_SUFFIXES = {
    ".bib",
    ".cpp",
    ".formula",
    ".g",
    ".json",
    ".log",
    ".md",
    ".opb",
    ".pbp",
    ".proof",
    ".py",
    ".sage",
    ".sh",
    ".sha256",
    ".tex",
    ".toml",
    ".txt",
    ".yaml",
    ".yml",
}


def tracked_paths(root: Path) -> list[str]:
    try:
        output = subprocess.check_output(
            ["git", "ls-files", "-z", "--cached", "--others", "--exclude-standard"],
            cwd=root,
            stderr=subprocess.DEVNULL,
        )
    except (FileNotFoundError, subprocess.CalledProcessError):
        manifest = root / "checksums.sha256"
        paths = [
            line.split("  ", 1)[1]
            for line in manifest.read_text(encoding="ascii").splitlines()
        ]
        return sorted((*paths, "checksums.sha256"))
    return sorted(
        {
            item.decode("utf-8").replace("\\", "/")
            for item in output.split(b"\0")
            if item and (root / item.decode("utf-8")).is_file()
        }
    )


def audit(root: Path = ROOT) -> list[str]:
    failures = []
    paths = tracked_paths(root)
    if len(paths) > 75:
        failures.append(f"tracked publication tree has {len(paths)} files; maximum is 75")
    for relative in paths:
        if relative in FORBIDDEN_TRACKED_PATHS or relative.startswith(
            FORBIDDEN_TRACKED_PREFIXES
        ):
            failures.append(f"internal path is tracked: {relative}")
            continue
        path = root / relative
        if path.suffix.lower() not in TEXT_SUFFIXES:
            continue
        content = path.read_bytes()
        for marker in FORBIDDEN_TEXT:
            if marker in content:
                failures.append(
                    f"private path marker {marker!r} in {relative}"
                )

    ignore_rules = {
        line.strip()
        for line in (root / ".gitignore").read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    }
    for rule in sorted(REQUIRED_IGNORE_RULES - ignore_rules):
        failures.append(f"missing .gitignore rule: {rule}")

    manifest_path = root / "checksums.sha256"
    if not manifest_path.is_file():
        failures.append("missing root checksums.sha256")
    else:
        manifest_names = {
            line.split("  ", 1)[1]
            for line in manifest_path.read_text(encoding="ascii").splitlines()
            if "  " in line
        }
        expected_names = set(paths) - {"checksums.sha256"}
        if manifest_names != expected_names:
            failures.append("checksums.sha256 does not cover exactly the publication tree")
    return failures


def main() -> int:
    failures = audit()
    if failures:
        for failure in failures:
            print(f"ERROR: {failure}")
        return 1
    print(f"public hygiene: PASS ({len(tracked_paths(ROOT))} tracked files)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
