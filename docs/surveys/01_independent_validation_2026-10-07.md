# Validation report: is the baseline valid? (2026-10-07)

Three independent checks were run on the baseline (`docs/theory.md`,
`docs/results_baseline.md`, `src/tde/`):

| Track | Method | Evidence |
|---|---|---|
| A. Literature | Each formula and claim checked against primary sources (IEEE papers, theses, lecture notes reproducing Kay) | sources cited below |
| B. Independent numerics | CRLB, ZZB and MLE re-implemented from scratch (no project code), 5,000–60,000-trial Monte Carlo, bootstrap CIs, 5 pilots, MIE theory | `validation/numeric/` |
| C. Code audit | Adversarial line-by-line review with targeted tests | `validation/audit/` |

(The scripts in `validation/` were copied from the session scratchpad; their
output paths still point there, so edit `OUT` paths before re-running.)

## Verdict

**The core is valid.** The model, the CRLB, the ZZB, the MLE and the DLL
implementation are correct, and every number reproduces. No estimator uses the
true delay. **Several interpretive claims are wrong or overstated** and must be
corrected before the report. None of them changes the main story (CRLB valid
above ≈ 12 dB; threshold ≈ 9 dB for the MLE; ZZB predicts its shape but not its
exact location).

---

## 1. Confirmed (independent evidence)

| Item | Project | Independent | Source |
|---|---|---|---|
| ZZB formula, ½ and 2(1−h/Ta) factors, P_min = Q(√(E/N0(1−ρ))) | — | correct | Bell et al. 1997 eq. form; Chazan–Zakai–Ziv 1975; Mishra et al. arXiv:1701.03898 |
| ZZB limits Ta²/12 and CRLB | — | correct analytically | Mishra et al.; own derivation |
| CRLB = 1/(2·E/N0·β²), E/N0 = Σs²/(2σ²) | — | correct | Kay Vol. 1 §3.11 (via Fowler EECE 522 notes) |
| β² of the pilot | 2.4484 | 2.4480 (continuous FFT) | B |
| ZZB, all 21 SNRs | summary.csv | max deviation 0.12 % | B |
| ZZB threshold | 5.95 dB | 5.95 dB (robust to n_h, sps, span, valley filling) | B |
| MLE MSE, 7 SNRs | summary.csv | all within MC error (\|z\| < 1.7) | B (independent MLE with golden-section search) |
| MLE curve vs theory | — | lies between two MIE approximations (union / Rice) for Ta = 4 and 16 | B |
| No ground truth used | — | DLL calibration at wrong sampling phases gives identical MSE (1.0785 vs 1.0785) | C |
| Calibration b_j needed | "≈100 % of CRLB bias" | without it MSE/CRLB = 2.07 at 30 dB | C |
| DLL loop formulas N_eff = (2−K)/K, B_L·T = K/(2(2−K)) | — | exact, and confirmed by simulation (16.08 vs 16) | A |
| Signs, windows, RNG independence, `threshold_snr`, E/N0 consistency MLE↔DLL | — | correct | C |

## 2. Problems found, and corrections

| # | Severity | Claim / issue | Finding | Correction |
|---|---|---|---|---|
| P1 | **Major** | "DLL and MLE use the same energy"; "DLL gap to ZZB 4.4 dB" | True only for steady-state jitter. The DLL processes 640 chips and recovers from wrong acquisitions during tracking (3.5 % of trials at 6 dB), so in the threshold region it uses ~10× the MLE's data. The 64-chip ZZB is **not a bound on the DLL**. | Remove the "gap to ZZB" for the DLL. Either compare against a ZZB for the DLL's data, or give both estimators equal total data. |
| P2 | **Major** | "DLL efficiency loss = ‖G‖²‖s′‖²/⟨s′,G⟩² by Cauchy–Schwarz" | Wrong reference waveform: over the 640-chip stream the ratio is 12.6, and a narrower loop (N_eff = 64) gives 0.26 × batch CRLB. DLL ≈ batch CRLB holds **by construction of B_L**, not as a bound. | Rewrite theory §5. |
| P3 | **Major** | Threshold values as single numbers; "DLL ≈ 1.2 dB worse" | 2,000 trials on a 2 dB grid are too coarse. The "highest 2×CRLB crossing" rule can be triggered by a single outlier. With 40k trials and a 0.5 dB grid: MLE **9.4 [9.3, 9.6] dB**; DLL **10.4 [9.6, 11.3] dB**. The DLL–MLE difference ≈ **1 ± 0.9 dB, not significant**. | Report CIs, use a finer grid and ≥2×10⁴ trials near threshold, or define the threshold by gross-error probability. |
| P4 | **Major** | "ZZB is ~3 dB optimistic" (general claim) | Pilot-dependent: 2.5–5.5 dB over 5 pilots. The ZZB bounds the **MMSE** estimator and the MLE is not MMSE, so the gap mixes bound looseness with MLE suboptimality. The literature reports the ZZB as quite tight for TDE (≈0.5 dB to MLE in Mallat et al. 2014; Xu–Baggeroer–Richmond 2004). | Qualify the claim; add a posterior-mean (MMSE) estimator to separate the two effects. |
| P5 | Major | "Gap shrinks with Ta: 3.0 → 1.8 dB" | Direction holds on average (mean over 5 pilots 3.4 → 2.0 dB) but not for every pilot. The 1,000-trial Ta = 16 run cannot resolve it (one outlier exceeds 2×CRLB). | Report as a trend with spread, using more trials. |
| P6 | Major | "MLE → Ta²/6 at low SNR (uniform guess)" | The match is coincidental. Under pure noise the argmax piles up at the window edges (~20 % of estimates at an endpoint), giving MSE → **3.17**. The −10 dB value (2.49) still carries signal information and is not a plateau. | Rewrite (results §3.4, theory §6). |
| P7 | Major | "Pull-in ≈ ±0.75 chip"; Fig 5 "separates the two DLL failure mechanisms" | The S-curve restores out to **~1.9 chips**. `acq_fail = |e|>0.5` counts recoverable acquisitions, so the categories overlap (acq + lock − total up to 0.18). | Define acquisition failure as \|e\| > pull-in range, and make the categories partition the total. |
| P8 | Major | "DLL MSE 70 chips² at −10 dB, unbounded random walk" | DLL MSE below threshold grows with stream length (40 / 69 / 114 chips² for 80 / 160 / 320 updates). The fixed block windows (not following θ) also shape it. | Report low-SNR DLL MSE versus observation time; let the windows follow θ. |
| P9 | Minor | "Per-block slopes ≈ 0 / negative", "naive theory 10× off" | True for an **earlier** design (random stream, per-block normalisation), not the final one (here k_j ∈ [4.7, 28.1]). | Reword as design history. |
| P10 | Minor | Wording | "you can't do worse than guessing the prior mean" should read "no estimator can do *better* than Ta²/12 at zero SNR". β is used for both the roll-off and the Gabor bandwidth. Theory §6 says failure MSE ≈ Ta²/12. | Fix the text; rename the roll-off to α. |
| P11 | Minor | "No threshold at Ta = 1" | It is masked by the prior (MSE ≤ ~0.17), not absent. Mallat et al. Part 2: the asymptotic threshold can be nearly independent of the prior for baseband signals. | Say "masked"; add the caveat. |
| P12 | Minor | Numerics and figures | Pulse truncation is t−k ∈ [−6, 7), not \|t\| ≤ 6 (negligible). The ZZB n_h = 2001 gives 0.996 instead of 1.005 at 40 dB. Valley filling is inactive here. The DLL linear-theory curve is drawn into the non-linear region. The constrained MLE beats the CRLB at prior edges (~1 % effect). Sanity check 3 is tautological. | Small code and doc fixes. |
| P13 | Minor | References | Weiss–Weinstein 1983 is *passive* TDE (random signals). Kaplan & Hegarty edition missing. | Add Betz & Kolodziejski 2009 (TAES 45(4)), Sadler & Kozick 2006, Ianniello 1982, Mallat et al. 2014, Holmes & Biederman 1978, Groh et al. 2012. |

## 3. Novelty of the proposed research directions (literature check)

| Direction | Status | Prior work |
|---|---|---|
| Explain the ZZB–MLE gap and its Ta dependence | **Largely done** | Ianniello–Weinstein–Weiss ICASSP 1983; Mallat et al. 2014 (Parts 1–2); Richmond IEEE TIT 2006; Xu et al. IEEE TSP 2004. Separating MMSE-vs-MLE from bound looseness for PN pilots is incremental. |
| Bound capturing DLL loss of lock | **Most open** | MTLL theory exists (Holmes & Biederman 1978; Groh et al. 2012); batch ZZB for code tracking (Côté, McGill thesis 2009); recursive Bayesian bounds in general (Tichavský et al. 1998; recursive WWB). No recursive ZZB/WWB for code-tracking loss of lock was found. |
| ZZB with unknown BPSK data bits | **Largely done** for OFDM | Graff & Humphreys arXiv:2407.19106 (2024); Decarli & Dardari ICASSP 2014. A single-carrier DS-BPSK version would be incremental. |

## 4. Recommended fix plan

1. **Docs:** apply P1, P2, P4, P6, P9–P11 and P13 (text only).
2. **Code:**
   - Threshold estimation with bootstrap CIs, a finer grid and more trials near threshold (P3, P5).
   - DLL failure classification from the measured pull-in range (P7).
   - DLL windows that follow θ, and low-SNR DLL MSE reported versus time (P8).
   - A posterior-mean (MMSE) estimator (P4).
   - Minor items (P12).
3. **Fair DLL comparison** (P1): add a "data-matched" ZZB (or an equal-total-data MLE) so the DLL is compared against a bound that actually applies to it.
4. Re-run and update `results_baseline.md`.
