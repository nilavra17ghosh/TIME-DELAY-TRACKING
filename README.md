# SDET Project 13 — Threshold effects in time-delay tracking: CRLB vs ZZB for DLL and batch MLE (BPSK)

**Start here.** Open the notebooks in order; each one goes *math → implementation → experiment → result →
interpretation → why it matters next*. All notebooks are saved **with their outputs**, so you can read them without
running anything.

## Status at the midterm checkpoint (Sprint 4)

| Sprint | Planned deliverable (`docs/sprint_plan.jpg`) | Where | Status |
|---|---|---|---|
| 0 | problem, model, notation | `notebooks/00_problem_and_model.ipynb` | ✅ |
| 1 | BPSK signal, batch MLE, first error-vs-SNR curve + CRLB | `notebooks/01_sprint1_batch_mle.ipynb` | ✅ |
| 2 | early–late DLL: S-curve + tracking run | `notebooks/02_sprint2_dll.ipynb` | ✅ |
| 3 | CRLB derived properly, matched to Sprint-1 data | `notebooks/03_sprint3_crlb.ipynb` | ✅ |
| 4 | ZZB derived for this case + intuition | `notebooks/04_sprint4_zzb.ipynb` | ✅ (+ numerical check, labelled) |
| 5–6 | four-curve plot, threshold, DLL vs MLE | `notebooks/05_beyond_sprint4.ipynb` | ⟶ **already completed / beyond current Sprint 4** |
| 7 | report + slides | — | upcoming |

### The results in four sentences
1. The batch MLE attains the CRLB above ≈ 12–14 dB (MSE/CRLB = 0.997 averaged over 14–30 dB).
2. Below that, rare **gross errors** (the estimator picks a wrong correlation peak) make the MSE break away from the
   CRLB; the CRLB cannot see this because it is a *local* bound that also ignores the prior.
3. The **ZZB**, derived here by turning estimation into binary detection between delays $h$ apart, reduces to the
   prior variance at low SNR and to the CRLB at high SNR, and bends away from the CRLB at ≈ 6 dB — about 3 dB before
   the MLE does (9.4 dB). Even the best possible (MMSE) estimator breaks away only ≈ 0.6 dB earlier than the MLE, so
   around the threshold the ZZB itself is optimistic by roughly 2–3 dB.
4. The DLL tracks with jitter ≈ 1.04–1.09 × the CRLB once locked (by design of its loop gain), pulls in from
   ±1.8 chips, and has an extra failure mode the batch MLE does not have: **loss of lock**.

### Honest assessment: is Sprint 4 enough for a midterm?
Sprints 1–4 alone give: two working estimators, a derived and verified CRLB, and a derived ZZB. That is a solid
*theory* checkpoint, but on its own Sprint 4 ends with a derivation and no comparison of the ZZB against data.
**The minimum addition that makes it substantial is the numerical ZZB overlaid on the Sprint-1 MLE results** — which
is already done and shown in notebook 04, §4.5 (labelled as Sprint-5 material). Nothing else needs to be added for
the midterm. Notebook 05 can be shown as "work already ahead of schedule", not as the checkpoint itself.

### What to tell your professor
> "Till Sprint 4 I have built the BPSK model, the batch MLE and the early–late DLL, derived the CRLB and the
> Ziv–Zakai bound for this exact model, and verified each piece against theory. The MLE follows the CRLB at high SNR
> and breaks away at about 9–10 dB because of wrong-peak errors; the ZZB, which I derived via binary detection between
> delays, captures that break-away and reduces to the CRLB and to the prior variance in the two limits. Next I will
> finish the four-curve comparison and locate the threshold precisely with confidence intervals — a first version of
> both already exists — and then explain the DLL vs MLE difference, where loss of lock looks like the open research
> question."

## Plan for the upcoming sprints
* **Sprint 5:** four-curve plot (CRLB / ZZB / MLE / DLL) — first version done (notebook 05 §C); finalise the
  presentation.
* **Sprint 6:** precise threshold with confidence intervals (notebook 05 §E has the method; run it with more trials),
  and the DLL-vs-MLE explanation (§D: acquisition vs loss of lock).
* **Sprint 7–8:** report and slides. Optional research angle: a bound for the DLL's loss-of-lock threshold.

## Layout
```
README.md            this file
notebooks/           00 … 05: the whole project, in order
sdet/                the only Python code (imported by the notebooks)
  model.py           signal, pulse, noise, prior            (Sprints 0–1)
  estimators.py      batch MLE (S1), DLL (S2)  [+ MMSE, DLL theory: beyond S4]
  bounds.py          CRLB (S3), ZZB (S4; numerics = S5)
  style.py           plot colours
figures/             every figure, saved by the notebooks (for slides)
data/                one small file of Sprint-1 results reused by Sprints 3–4
docs/                literature_survey.md (what is known vs open), project brief notes, sprint plan image
archive/             earlier scripts, docs, validation evidence (see archive/README.md)
steven-m-kay-vol1.pdf  reference textbook
```

## How to run
Requirements: Python 3 with numpy, scipy, matplotlib, jupyter. From the `notebooks/` folder open any notebook and
*Run All*. Run notebook 01 before 03 and 04 (they reuse its saved results). Runtimes: 00–04 about 1–2 min each,
05 about 10–15 min.

## Key references
Kay, *Fundamentals of Statistical Signal Processing Vol. I* (1993) §3.11, Ch. 7 · Ziv & Zakai, IEEE T-IT 1969 ·
Chazan, Zakai & Ziv, IEEE T-IT 1975 · Bell, Steinberg, Ephraim & Van Trees, IEEE T-IT 1997 · Kaplan & Hegarty,
*Understanding GPS/GNSS*, 3rd ed. (2017) · Betz & Kolodziejski, IEEE T-AES 2009. Full list in notebook 00, §0.7.
