#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
KEEP_ARTIFACTS=false
if [[ -n "${CW105_ARTIFACT_DIR:-}" ]]; then
  OUT="$CW105_ARTIFACT_DIR"
  KEEP_ARTIFACTS=true
  if [[ -e "$OUT" ]] && find "$OUT" -mindepth 1 -print -quit | grep -q .; then
    echo "Refusing to overwrite nonempty artifact directory: $OUT" >&2
    exit 2
  fi
  mkdir -p "$OUT"
else
  OUT="$(mktemp -d -t cw105-proof.XXXXXXXX)"
  trap 'rm -rf -- "$OUT"' EXIT
fi

for command in python g++ latexmk sha256sum cmp; do
  command -v "$command" >/dev/null || {
    echo "Missing required command: $command" >&2
    exit 2
  }
done

mkdir -p "$OUT/bin" "$OUT/logs" "$OUT/paper" "$OUT/results"
cd "$ROOT"

echo "[1/9] Primary Python support classification"
python -m src.mu6_35_12 \
  --output "$OUT/results/python.json" \
  --masks-output "$OUT/results/support-masks.txt" \
  > "$OUT/logs/python.stdout"
cmp "$OUT/results/support-masks.txt" data/mod3_weight12_support_masks.txt

echo "[2/9] Independent C++ support classification"
g++ -std=c++20 -O3 -Wall -Wextra -Werror -pedantic \
  src/mu6_35_12_independent.cpp -o "$OUT/bin/mu6-independent"
"$OUT/bin/mu6-independent" > "$OUT/results/cpp-independent.json"

echo "[3/9] Direct all-support phase check"
g++ -std=c++20 -O3 -Wall -Wextra -Werror -pedantic \
  src/mu6_35_12_direct_all.cpp -o "$OUT/bin/mu6-direct-all"
"$OUT/bin/mu6-direct-all" "$OUT/results/support-masks.txt" \
  > "$OUT/results/cpp-direct.json"

echo "[4/9] Exact reconciliation"
python validation/reconcile_mu6_35_12.py \
  --python "$OUT/results/python.json" \
  --cpp "$OUT/results/cpp-independent.json" \
  --direct "$OUT/results/cpp-direct.json" \
  --masks "$OUT/results/support-masks.txt" \
  --output "$OUT/results/reconciliation.json" \
  > "$OUT/logs/reconciliation.stdout"

echo "[5/9] Contraction and fibre-character bridge"
python -m src.derive_icw --output "$OUT/results/contraction-python.json" \
  > "$OUT/logs/contraction-python.stdout"
g++ -std=c++20 -O3 -Wall -Wextra -Werror -pedantic \
  src/derive_icw_independent.cpp -o "$OUT/bin/contraction-independent"
"$OUT/bin/contraction-independent" \
  > "$OUT/results/contraction-independent.json"
python src/verify_fibre_character_bridge.py \
  > "$OUT/results/fibre-character-bridge.txt"

echo "[6/9] Independent publication audit and hostile mutations"
python validation/publication_audit.py . \
  > "$OUT/results/publication-audit.json"
python validation/run_hostile_mutations.py . \
  > "$OUT/results/hostile-mutations.json"

echo "[7/9] Focused test suite"
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python -m pytest -q \
  | tee "$OUT/logs/pytest.txt"

echo "[8/9] Deterministic paper build"
cp paper/cw105_36_nonexistence.tex paper/references.bib "$OUT/paper/"
(
  cd "$OUT/paper"
  latexmk -pdf -interaction=nonstopmode -halt-on-error \
    cw105_36_nonexistence.tex > "$OUT/logs/latexmk.txt"
)
if grep -E -i \
  'LaTeX Warning|Overfull|Underfull|undefined|multiply defined' \
  "$OUT/paper/cw105_36_nonexistence.log"; then
  echo "Paper build contains a LaTeX warning." >&2
  exit 1
fi
cmp "$OUT/paper/cw105_36_nonexistence.pdf" \
  paper/cw105_36_nonexistence.pdf

echo "[9/9] Publication manifest"
sha256sum -c checksums.sha256 \
  | tee "$OUT/logs/checksums.txt"

echo "All exact publication checks passed."
if [[ "$KEEP_ARTIFACTS" == true ]]; then
  echo "Artifacts: $OUT"
fi
