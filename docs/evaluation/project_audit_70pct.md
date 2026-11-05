# Project audit and 70% milestone assessment (2026-11-05)

**Sources used:** the repository (`github.com/nilavra17ghosh/TIME-DELAY-TRACKING`, single commit "Mid-analysis",
2026-11-04), the archived scripts and validation evidence, `docs/brief/project_brief_notes.md`, the official brief
screenshot, the sprint plan image, and our working sessions from 2026-10-06 to 2026-10-13.

**Missing (not available to this audit — please supply if they exist):** the course **evaluation rubric**, any
**previous report or presentation** submitted to the professors, and the team's own notes from Oct 13 – Nov 4.
**Audit finding:** the repository is identical to the Oct 13 state; **no technical work happened between Oct 13 and
Nov 4.** Notebook 06 (precise thresholds) and the documents in `docs/` dated Nov 5 were added for this evaluation.

---

## 1. Reconstructed history

**Problem (from the official brief, Project 13, mentor: Kiruba Selvakumar).** Estimate the delay of a BPSK signal with
a DLL (real time) and a batch MLE; obtain error variance vs SNR; compare with the CRLB; find the SNR where they diverge;
formulate, derive and simulate the ZZB; compare; explain the DLL–MLE difference.

**Research questions we formulated.**
1. Where does each estimator leave the CRLB, and why?
2. Does the ZZB predict that departure?
3. Why do the DLL and MLE differ?
4. (Added Oct 7–13) How much of the ZZB–MLE gap is bound looseness vs estimator sub-optimality? Is there a bound
   that describes the DLL's loss of lock?

**Key design decisions (and the alternatives rejected).**

| Decision | Rejected alternative | Reason |
|---|---|---|
| Raised-cosine pulse, analytic delay | rectangular chips; FFT fractional delay (the original Sprint-1 script) | CRLB = 0 for rectangular chips; analytic delay is exact |
| Random delay with uniform prior ($T_a=4$) | fixed delay (Sprint-1 script); $T_a=1$ | ZZB needs a prior; $T_a=1$ masks the threshold |
| SNR = batch $E/N_0$ | per-sample SNR | bounds depend on total energy |
| MSE metric | variance | estimators are biased below threshold; ZZB bounds MSE |
| DLL: acquisition + block replicas + offset calibration + windows that follow θ + same pilot repeated | full replica over truncated window; per-block normalisation; fixed windows; different pilot | each alternative produced a measured artefact (bias of 33 %/100 % of CRLB; noise amplification; DLL "beating" the CRLB) |
| Loop gain $K=2/17$ ($1/(2B_L)$ = 64 chips) | arbitrary $K$ | equal steady-state information with the batch |
| Two threshold criteria with bootstrap CIs | single MSE criterion, no CIs | MSE near threshold is outlier-dominated |
| Smoothly tapered pulse truncation | hard truncation | hard cut made the ZZB miss the CRLB limit by ~1 % |
| Notebooks + 4-module library | many scripts | understandability (your request, Oct 9) |

---

## 2. Status of every component

### 2.1 Completed (code + executed results + verified)
| Component | Evidence |
|---|---|
| Signal/noise/prior model and its checks | NB00 (Nyquist, energy, SNR calibration) |
| Batch MLE, error-vs-SNR curve, error distributions, the wrong-peak mechanism | NB01; MLE/CRLB = 0.997 above 14 dB |
| Early–late DLL: S-curve, pull-in (±1.80 chips), tracking runs, steady-state jitter (1.04–1.09 × CRLB) | NB02 |
| CRLB derivation, two independent checks (FFT 0.05 %), pilot dependence (β² 1.77–4.27), overlay, three failure reasons | NB03 |
| ZZB derivation (6 steps + both limits), intuition plots, numerical ZZB, ZZB below MLE within CI | NB04 |
| Four-curve comparison (+ MMSE, + DLL-valid bound), DLL theory (σ²‖G‖²), DLL failure decomposition, MSE-vs-time | NB05 |
| Precise thresholds with CIs, 5 pilots, DLL vs MLE significance | **NB06** (added for this evaluation) |
| Prior-width study ($T_a$ = 1, 4, 16) with CIs | **NB07** (added for this evaluation) |
| Independent validation (literature, re-implementation, audit) | docs/surveys/01_independent_validation_2026-10-07.md |
| Literature survey + novelty assessment | docs/surveys/02_literature_survey_2026-10-13.md, docs/evaluation/research_proposal.md §0 |
| Fast regression tests | tests/test_core.py (7 tests, ~35 s) |

### 2.2 Partially completed
| Component | What is missing |
|---|---|
| Explanation of the pilot-to-pilot spread of the ZZB–MLE gap | measured, not yet explained (per-code MIE prediction not done) |
| "Validation of the ZZB" (brief) | ZZB vs simulation done; a written discussion for the report not yet done |
| Report / slides (Sprint 7) | not started |

### 2.3 Planned but not implemented
* Method-of-interval-errors (MIE) prediction inside the project code (only in archived validation scripts).
* Any time-varying delay, Doppler, carrier phase, unknown data bits, multipath.
* The research extension (dynamic decomposition, docs/evaluation/research_proposal.md) — **entirely proposed**.

### 2.4 New, unexplored ideas
* Sequential/marginal large-error bound for code tracking (stretch goal of the proposal).
* Loss-of-lock risk score (D2). Code selection by threshold (D4). Robot UWB transfer (D5, judged forced).

---

## 3. Is this genuinely ~70%?

Against the **official brief**, every required element exists: DLL ✓, batch MLE ✓, variance(MSE)-vs-SNR ✓,
CRLB ✓, divergence SNR ✓ (with CIs after NB06), ZZB formulated/derived/simulated ✓, comparison ✓, DLL-vs-MLE
explanation ✓ (mechanisms + significance). Against the **sprint plan**, Sprints 1–6 of 7–8 are done; Sprint 7 (report,
slides) is not.

**Honest verdict: the *technical* work is ≈ 75–80 % of the brief; the *deliverable* work (report, polished
presentation, written discussion) is behind.** Claiming "70 %" is defensible and slightly conservative. What is *not*
defensible is presenting the research extension as anything but a proposal.

---

## 4. What a professor will likely criticise (and our answer)

| Likely criticism | Severity | Answer / mitigation |
|---|---|---|
| "This is all known theory — what is yours?" | high | Be upfront: course-level reproduction, carefully verified; own observations (MMSE vs ZZB, pilot spread, DLL failure split); proposal for a genuine extension (question bank F1–F4) |
| "Your DLL vs MLE comparison is unfair" | medium | Steady-state matched by design; threshold region not equal-information; DLL-valid 640-chip ZZB shown (D3) |
| "The ZZB is 3 dB off — is it wrong?" | medium | Limits verified; independent implementation; decomposition MMSE vs MLE; pilot dependence (E4) |
| "Threshold definition is arbitrary" | medium | Pre-registered, two criteria, CIs (E2) |
| "Only synthetic data" | low–medium | Required for ground truth and bound comparison (E6) |
| "No Doppler / second-order loop / real GNSS" | low | Scope of the brief; listed as extension (G4) |
| "Why did the model change after Sprint 1?" | low | ZZB needs a random delay; exact analytic delay; story unchanged (NB00 §0.2) |
| "Nothing changed in three weeks" | **real risk** | Show NB06 (new, precise, with CIs) and the proposal; do not hide the gap |

---

## 5. Minimum work before the evaluation (status)
1. ✅ Precise threshold study with confidence intervals — notebook 06.
2. ✅ Research proposal with literature positioning — docs/evaluation/research_proposal.md.
3. ✅ Presentation guide and question bank — docs/evaluation/presentation_guide.md, docs/evaluation/question_bank.md.
4. ⬜ **You:** read NB00 → NB01 → NB03 → NB04 → NB06 and the question bank sections A, C, F (≈ 3–4 hours).
5. ⬜ Optional, if time: build the slides from the presentation guide.
