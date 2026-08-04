# Nonexistence of $\mathrm{CGW}(35,12;6)$ and $\mathrm{CW}(105,36)$

This repository contains an exact computer-assisted proof of two statements:

1. no circulant generalized weighing matrix $\mathrm{CGW}(35,12;6)$ exists;
2. consequently, no circulant weighing matrix $\mathrm{CW}(105,36)$ exists.

The proof uses exact integer and finite-field arithmetic throughout.
The manuscript and computations are ready for specialist review, but the
result has not yet been peer reviewed.

## Paper

- [Focused manuscript (PDF)](paper/cw105_36_nonexistence.pdf)
- [LaTeX source](paper/cw105_36_nonexistence.tex)
- [Bibliography](paper/references.bib)

## Proof outline

1. Reducing a hypothetical sixth-root perfect sequence modulo
   $1-\omega$ places its support in one of two ternary cyclic $[35,12]$
   codes.
2. Complete enumeration of both codes produces 420 distinct weight-12
   supports.
3. Those supports form one affine orbit in $\mathbb Z/35\mathbb Z$, represented by
   $\{0,1,2,3,7,10,12,16,21,22,26,28\}$.
4. On that representative, four correlations over $\mathbb F_4$, at shifts
   $2,4,8,10$, give an explicit contradiction.
5. A separate C++ program redundantly checks all
   $420\cdot 3^{11}=74{,}401{,}740$ normalized phase assignments.
6. The multiplier-fixed contraction of a hypothetical $\mathrm{CW}(105,36)$ and its
   nontrivial three-fibre character would produce the excluded
   $\mathrm{CGW}(35,12;6)$.

## Reproduce

Requirements are Python 3.11 or newer, NumPy, pytest, a C++20 compiler,
`latexmk` with a standard LaTeX installation, and `sha256sum`.

From a fresh clone, install the Python package and run the complete local
reproduction:

```bash
python -m pip install -e '.[test]' && ./run_all.sh
```

The script uses a disposable output directory. It:

- regenerates the complete support family with Python;
- compares it byte-for-byte with the canonical 420-support certificate;
- runs the independent C++ support classifier;
- checks every support and normalized phase assignment with separate C++;
- reconciles the three exact outputs;
- reproduces both multiplier-fixed contraction rows in Python and C++;
- verifies the fibre-character bridge symbolically;
- runs the focused test and hostile-mutation suites;
- rebuilds the paper and compares it with the committed PDF;
- verifies the root SHA-256 manifest.

Set `CW105_ARTIFACT_DIR` to retain the generated JSON, logs, binaries, and
rebuilt paper:

```bash
CW105_ARTIFACT_DIR=/tmp/cw105-proof ./run_all.sh
```

## Independent validation

The primary implementation is
[`src/mu6_35_12.py`](src/mu6_35_12.py).
The standalone C++ classifier shares no project code with it, and the direct
C++ checker uses the generated support list but independently tests every
remaining phase assignment.

Two additional end-to-end implementations are retained as source:
[`validation/sage/validator.sage`](validation/sage/validator.sage) and
[`validation/gap/validator.g`](validation/gap/validator.g).
The single GitHub Actions workflow runs them in SageMath 10.9 and GAP 4.12.1,
then reconciles their independently generated support lists with the
canonical certificate.

## Trust boundary

- Every theorem-level decision uses exact arithmetic.
- The support certificate is an output of complete finite enumeration, not a
  search input to the primary classifier.
- The short $\mathbb F_4$ contradiction is checked symbolically and by exhaustive
  enumeration.
- The direct all-support checker is redundant validation, not a logical
  premise of the four-shift proof.
- Hostile tests mutate each factor branch, support mask, affine convention,
  correlation identity, fibre table, and contraction convention; every
  mutation must be rejected.
- Only complete exact computations enter the proof; floating-point decisions
  and incomplete searches do not.
- The imported multiplier and contraction results are identified and cited
  in the manuscript.
- Independent human review and peer review remain outside this repository's
  validation boundary.

## Repository map

| Path | Purpose |
|---|---|
| `paper/` | One focused manuscript, its source, and bibliography |
| `data/` | Canonical 420-support certificate and proof specification |
| `src/` | Primary and independent exact implementations |
| `tests/` | Theorem regressions, publication checks, hostile mutations |
| `validation/` | Standalone audit, reconciliation, SageMath, and GAP |
| `checksums.sha256` | One manifest covering the publication tree |
| `.github/workflows/ci.yml` | The complete Linux and CAS validation workflow |

## AI-assisted research disclosure

OpenAI GPT-5.6 Sol, accessed through ChatGPT and Codex, was used as a substantive interactive research and software-development assistant. It contributed candidate algebraic reductions, proof-step stress testing, exact-checking code, adversarial validation design, and manuscript organization and editing. All retained arguments, computations, code, and citations were reviewed and verified by the human author, who takes full responsibility for the work.

## Author and citation

Semyon Osadchii, independent researcher.

Contact: `Sposadchii@gmail.com` or `sp.osadchii@student.han.nl`.

Citation metadata is provided in [`CITATION.cff`](CITATION.cff).
