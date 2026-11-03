# Archive — retired material (kept, not deleted)

Nothing here is needed to run or understand the project; the notebooks replace all of it. It is kept as a record.

| Folder | What it is | Why it was retired |
|---|---|---|
| `sprint1_original/` | the first Sprint-1 script and figure (127-chip m-sequence, Gaussian filter, fixed delay) | superseded by notebook 01; the model was changed so the ZZB (which needs a random delay with a prior) applies — see notebook 00 §0.2 |
| `scripts_v1/` | `run_baseline.py`, `run_threshold_study.py`, `sanity_checks.py` and the first package `src/tde/` | consolidated: the run scripts became notebooks, the sanity checks moved into the notebook of the sprint they belong to, the package shrank to `sdet/` (4 modules) |
| `docs_v1/` | `theory.md`, `results_baseline.md`, `validation_report.md` | their content now lives in the notebooks. `results_baseline.md` describes an **older, pre-fix** version — do not quote its numbers. `validation_report.md` is the record of the independent validation (still accurate as a list of what was checked and fixed) |
| `validation_2026-10-07/` | scripts and outputs of the independent numerical check and code audit | evidence behind the validation report. The scripts' output paths point to a temporary folder that no longer exists; edit them before re-running |
| `results_v1_UNTRUSTED_smoke_tests/` | outputs of small (300–1000-trial) smoke tests; the full runs were stopped for low memory | **not valid results** — kept only so nothing is lost |

Changes since the archived code (all small): the pulse is now smoothly tapered over its last chip (removes a tiny
slope discontinuity that kept the ZZB from converging exactly to the CRLB); the DLL has an `s_curve` method (the
Sprint-2 deliverable); the threshold utilities live inside notebook 05.
