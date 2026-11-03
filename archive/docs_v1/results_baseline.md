# Baseline results

Reproduce with:

    python sanity_checks.py     # 10 checks of the code against theory (~25 s)
    python run_baseline.py      # full experiment (~8 min); --quick for a smoke test

Setup: 64-chip known BPSK pilot, raised-cosine (β = 0.5), 4 samples/chip, AWGN,
prior τ ~ U(width T_a = 4 chips), 2000 Monte Carlo trials per SNR, E/N0 = −10…30 dB.
DLL: acquisition on a 0.5-chip grid, early–late d = 0.5, 4-chip updates,
first-order loop with 1/(2B_L) = 64 chips (same energy as the MLE batch).
Figures: `results/figures/`; numbers: `results/tables/summary.csv`.

## 1. Code verification (`sanity_checks.py`): all 10 pass

| Check | Result |
|---|---|
| SNR calibration | empirical E/N0 = 9.995 for target 10 |
| β² averaged over 60 pilots vs pulse | 2.832 vs 2.838 |
| MLE noise-free error | 0.05 % of CRLB at 30 dB |
| DLL noise-free error | ≈ 0 (1.6·10⁻¹⁰ chips²) |
| ZZB, low-SNR limit | 1.323 vs T_a²/12 = 1.333 |
| ZZB, high-SNR limit | ZZB/CRLB = 0.996 at 40 dB |
| MLE efficiency at 25 dB | MSE/CRLB = 0.986 |
| DLL linear theory at 25 dB | MSE/theory = 0.996 |

## 2. Main result (`fig4_mse_vs_snr.png`, `fig5_threshold.png`)

| Region (batch E/N0) | What happens |
|---|---|
| **> 14 dB: asymptotic** | MLE MSE/CRLB = 1.01 at 30 dB (efficient). DLL = 1.08, matching its exact linear theory (1.07). ZZB ≈ CRLB. |
| **≈ 6–12 dB: threshold** | Gross errors (\|e\| > 0.5 chip) appear: MLE 0.6 % at 10 dB, 8 % at 6 dB. MSE leaves the CRLB. |
| **< 2 dB: no information** | MLE MSE → 2.5 ≈ T_a²/6; DLL MSE grows to 70 chips² (random walk after loss of lock). |

**Threshold SNR** (MSE first exceeds 2 × CRLB, batch E/N0; per-chip = −18.1 dB):

| | Threshold | Gap to ZZB |
|---|---|---|
| ZZB | 5.95 dB | — |
| Batch MLE | 9.11 dB | 3.2 dB |
| DLL | 10.33 dB | 4.4 dB |

## 3. Interpretation

1. **The CRLB fails below ≈ 12 dB for both estimators.** It is a local bound
   (curvature of the likelihood at the true τ) and cannot see the wrong peaks
   that cause gross errors (fig 3: the error histogram goes from a narrow
   Gaussian at 18 dB to a broad distribution spanning the whole prior at 2 dB).

2. **The ZZB predicts the threshold *shape* but is ~3 dB optimistic about its
   location.** It correctly rises above the CRLB, peaks, and saturates at the
   prior variance. The 3 dB gap is expected: the ZZB is a bound built from
   *pairwise* (binary) tests between τ and τ + h. The MLE must beat *all*
   competing delays at once (an M-ary problem), so its error probability is
   larger than any single pairwise P_min. Closing this gap is a known open
   theme (see the research directions below).

3. **At very low SNR the CRLB is above both the ZZB and the prior variance**,
   i.e. it is not a valid description at all there (the estimators are biased,
   and a biased estimator is not bound by the CRLB).

4. **The MLE's low-SNR MSE is T_a²/6, not T_a²/12.** With no information the MLE
   picks a point essentially uniformly in the window, independent of τ:
   E[(U₁ − U₂)²] = 2·T_a²/12. The ZZB saturates at T_a²/12 because it bounds
   *every* estimator, including "always answer the prior mean". The MLE is not
   the MMSE estimator; a Bayesian (posterior-mean) estimator would approach the
   ZZB more closely at low SNR.

5. **DLL vs MLE.**
   * *High SNR: DLL is 0.3 dB worse (1.08 × CRLB).* By Cauchy–Schwarz,
     a linearised DLL equals the CRLB only if its equivalent weighting waveform
     G is proportional to s′. Early–late with d = 0.5 is a finite-difference
     approximation of s′, and the exponential loop memory is not a uniform
     window. The linear theory σ²‖G‖² predicts the loss exactly (1.07).
   * *Threshold: DLL is ≈ 1.2 dB worse.* Fig 5 (right) separates its two failure
     mechanisms: **wrong acquisition** (same wrong-peak mechanism as the MLE,
     on a coarser grid) and **loss of lock** (noise drives the error out of the
     S-curve's linear range during tracking, after which the loop drifts).
   * *Low SNR: DLL is far worse (70 vs 2.5 chips²).* The MLE's answer is
     always inside the prior window. The DLL has no such constraint; once
     unlocked it random-walks. Its MSE grows with observation time, so the DLL
     is **not** described by the batch ZZB below threshold. Loss of lock is a
     dynamics problem (mean time to lose lock), not an ambiguity problem.

6. **Prior width controls the threshold** (`fig6_prior_width.png`).
   ZZB thresholds: T_a = 1 → none, 4 → 5.95 dB, 16 → 11.56 dB. MLE thresholds:
   none / 8.97 / 13.36 dB. With T_a = 1 chip (≈ one correlation main lobe) there
   is no wrong peak to jump to, so no threshold exists. More delay uncertainty
   means more competing peaks and a higher threshold. Interestingly the
   ZZB-to-MLE gap *shrinks* (3.0 dB → 1.8 dB) as T_a grows.

## 4. Implementation pitfalls found (worth a paragraph in the report)

* Rectangular chips give an infinite-bandwidth CRLB of 0, so a band-limited pulse is required.
* The CRLB depends on the specific pilot pattern (β² ∈ [1.77, 4.23] over 64-chip
  pilots), so the DLL and MLE must use the same pilot. Otherwise the DLL appeared
  to "beat" the CRLB (MSE/CRLB = 0.90).
* DLL block-correlation bias: a hard-truncated correlation window gave a bias of
  33 % of the CRLB. Block-only replicas reduced it, and a periodic pilot then
  needed explicit code-pattern offset calibration (bias of 100 % of the CRLB
  without it).
* Normalising the discriminator per block is unstable (some block slopes ≈ 0).
* Simple "independent per-update noise" DLL theory was off by 10×; the exact
  equivalent-weighting analysis fixes it.

## 5. Research directions this baseline points to

1. **Close the 3 dB ZZB-to-MLE gap**: compare with the method of interval
   errors (MIE), which predicts MLE MSE through the threshold, and with tighter
   Bayesian bounds (Weiss–Weinstein). Why does the gap shrink with T_a?
2. **A bound for the DLL**: the batch ZZB does not capture loss of lock.
   Candidates are recursive Bayesian bounds (tracking WWB/BCRB) and mean time to lose lock.
3. **Design trade-offs**: the threshold and efficiency as functions of the
   early–late spacing d, the loop bandwidth B_L and the acquisition grid.
4. **Unknown data bits** (non-pilot BPSK): ZZB and MLE with the bits as
   nuisance parameters.
