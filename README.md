# SDET Project 13 — Threshold effects in time-delay tracking: CRLB vs ZZB for DLL and batch MLE (BPSK)

We estimate the delay of a known BPSK signal in noise with a **batch maximum-likelihood estimator** and an
**early–late delay-locked loop (DLL)**, and compare them with the **Cramér–Rao lower bound (CRLB)** and the
**Ziv–Zakai bound (ZZB)**, focusing on the low-SNR *threshold* where errors suddenly become large.

**Start here:** open the notebooks in order. Each goes *math → implementation → experiment → result → meaning → why
it matters next*, and is saved **with its outputs**.

## Status at the 70% evaluation

| Sprint | Deliverable (`docs/sprint_plan.jpg`) | Where | Status |
|---|---|---|---|
| 0 | problem, model, notation, literature | `notebooks/00_problem_and_model.ipynb` | ✅ |
| 1 | BPSK signal, batch MLE, error-vs-SNR + CRLB | `notebooks/01_sprint1_batch_mle.ipynb` | ✅ |
| 2 | early–late DLL: S-curve + tracking run | `notebooks/02_sprint2_dll.ipynb` | ✅ |
| 3 | CRLB derived, matched to data | `notebooks/03_sprint3_crlb.ipynb` | ✅ |
| 4 | ZZB derived + intuition (+ numerical check) | `notebooks/04_sprint4_zzb.ipynb` | ✅ |
| 5 | four-curve comparison, DLL variance theory, DLL failure modes | `notebooks/05_beyond_sprint4.ipynb` | ✅ |
| 6 | threshold located precisely (CIs, 5 pilots), DLL vs MLE | `notebooks/06_sprint6_threshold.ipynb` | ✅ |
| 7–8 | research extension, report, slides | `docs/research_proposal.md` | ⏳ **proposed, not started** |

## Results registry (every number traceable to an executed notebook)

| Result | Value | Evidence |
|---|---|---|
| MLE efficiency above threshold | MSE/CRLB = 0.997 (14–30 dB, 2000 trials/pt) | NB01, NB03 |
| CRLB bandwidth check | β² time vs frequency domain agree to 0.05 % | NB03 |
| CRLB pilot dependence | β² 1.77–4.27 over 400 random pilots | NB03 |
| ZZB limits | → Ta²/12 (1.332 vs 1.333) and → CRLB (0.997 at 40 dB) | NB04 |
| DLL pull-in range / locked jitter | ±1.80 chips / 1.04–1.09 × CRLB (by design of K) | NB02 |
| **Thresholds, pilot 2026** (MSE > 2·CRLB, 95 % CI) | ZZB 5.95 · MMSE 8.54 [8.35, 9.07] · MLE 9.52 [9.00, 9.66] · DLL 9.38 [8.72, 9.98] dB | NB06 |
| Gap decomposition | bound looseness 2.59 dB + MLE sub-optimality 0.85 dB [0.25, 1.24] | NB06 |
| Pilot spread of the MLE–ZZB gap | 2.7–5.4 dB over 5 pilots | NB06 |
| DLL vs MLE | MSE criterion: no difference; gross-error criterion: DLL 0.68 dB more robust [0.40, 0.99] | NB06 |
| DLL failures below threshold | mostly loss of lock (0.42 of 0.54 at 0 dB); MSE grows with time | NB05 |

**In short:** the MLE follows the CRLB down to ≈ 12 dB, then wrong-peak errors take over (threshold 9.5 dB). The ZZB
captures this qualitatively but is ≈ 2.6 dB optimistic even for the best possible estimator. The DLL matches the MLE's
MSE threshold, makes fewer gross errors near threshold (it sees 10× more data), and fails much worse below it
(loss of lock).

## Documents for the evaluation (`docs/`)
| File | What it is |
|---|---|
| `project_audit_70pct.md` | history, completed / partial / planned work, honest 70 % assessment, likely criticisms |
| `research_proposal.md` | literature position, 5 candidate directions compared, recommended proposal + roadmap |
| `presentation_guide.md` | 13 slides: headline, content, figure, what to say, likely questions; equation sheet |
| `question_bank.md` | ≈ 40 questions with answers and follow-ups |
| `literature_survey.md` | what is known vs open (Oct 13), with sources |
| `project_brief_notes.md`, `project_brief_screenshot.png`, `sprint_plan.jpg` | the original brief and plan |

## Repository layout
```
notebooks/   00–06: the whole project, in order (the experiments)
sdet/        the reusable code, imported by the notebooks
  model.py       Config (all parameters), pulse, BPSK waveform, noise, prior
  estimators.py  BatchMLE (+ MMSE), EarlyLateDLL (+ S-curve, pull-in, linear theory)
  bounds.py      crlb, zzb, autocorrelation
  style.py       plot colours
tests/       fast regression tests of sdet against analytic facts
figures/     all figures, named by sprint (s0_… to s6_…)
data/        s1_mle.npz — Sprint-1 errors reused by notebooks 03–04
docs/        evaluation documents (above)
archive/     earlier scripts, validation evidence, untrusted smoke tests (see archive/README.md)
```
Configuration: every parameter lives in `sdet.Config` (pilot length, samples/chip, roll-off, prior width, DLL
spacing/block length/gain, trials). Notebooks create `Config()` and change only what they state explicitly.

## Setup and running
```bash
pip install -r requirements.txt          # Python 3.13 used; numpy, scipy, matplotlib, jupyter, pytest
python -m pytest -q tests                # ~35 s, 7 tests
jupyter lab notebooks/                   # then Run All in each notebook
```
Run notebook 01 before 03 and 04 (they reuse `data/s1_mle.npz`). Runtimes: 00–04 ≈ 1–2 min each, 05 ≈ 10–15 min,
06 ≈ 20 min. All random seeds are fixed, so results reproduce exactly.

## Key references
Kay (1993) Vol. I §3.11, Ch. 7 · Ziv & Zakai (1969) · Chazan, Zakai & Ziv (1975) · Bell et al. (1997) · Mallat et al.
(2014) · Emmanuele et al. (ICASSP 2012) · Betz & Kolodziejski (2009) · Hasson & Bobrovsky (1995). Full lists:
notebook 00 §0.7 and `docs/literature_survey.md`.
