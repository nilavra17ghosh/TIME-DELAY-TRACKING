# Theory behind the baseline (matched exactly to the code)

Every formula here is implemented in `src/tde/`. Units: time and delay in **chips**
(T_c = 1), MSE in **chips²**, SNR = **total batch E/N0** (64-chip pilot).
Revised 2026-10-07 after independent validation (`docs/validation_report.md`).

Notation: **α** = raised-cosine roll-off (0.5); **β²** = mean-square (Gabor)
bandwidth of the waveform. (Earlier drafts used β for both.)

---

## 1. Signal and observation model  (`signal_model.py`)

Known BPSK pilot `a_k ∈ {±1}`, k = 0..N_c−1 (N_c = 64), raised-cosine pulse p(t)
(roll-off α = 0.5, truncated symmetrically to |t| ≤ 6 chips):

    s(t) = Σ_k a_k p(t − k)

Sampled at f_s = 4 samples/chip, real baseband, coherent:

    r[n] = s(t_n − τ) + w[n],    w[n] ~ N(0, σ²) i.i.d.,    t_n = n / 4

**Why raised-cosine instead of rectangular chips?** The CRLB (§3) is
1/(2·E/N0·β²). A rectangular chip has infinite bandwidth, so β² = ∞ and the
CRLB = 0, which is meaningless. Real front-ends band-limit the signal, so a
band-limited pulse is both physical and makes the bound well-defined. p(t) is
evaluated analytically, so a fractional delay τ is applied exactly.

**SNR definition.** For a real signal sampled at or above its Nyquist rate,
with noise PSD N0/2 and σ² = N0·f_s/2,

    E/N0 = Σ_n s[n]² / (2σ²)     ⇒     σ² = Σ_n s[n]² / (2·E/N0)

Per-chip SNR: Ec/N0 = E/N0 − 10·log10(64) = E/N0 − 18.1 dB.

**Prior.** τ ~ U[τ_c − T_a/2, τ_c + T_a/2] (default T_a = 4 chips). The prior is
needed by the (Bayesian) ZZB and the MMSE estimator, and it defines the search
range of every estimator.

---

## 2. Batch estimators  (`estimators.BatchMLE`)

**MLE.** With white Gaussian noise,

    ln p(r; θ) = const − (1/2σ²) Σ_n (r[n] − s(t_n − θ))²

    τ̂_ML = argmax_θ  J(θ),     J(θ) = Σ_n r[n] s(t_n − θ) − ½ Σ_n s(t_n − θ)²

i.e. a matched-filter (correlator) search over the prior window. σ² does not
appear, so the MLE needs no SNR knowledge. Implementation: J on a 1/64-chip
grid, then a parabolic refinement (noise-free error ≤ 0.1 % of the CRLB at 30 dB).

**MMSE (posterior mean).** With the uniform prior, p(θ | r) ∝ exp(J(θ)/σ²) on
the window, and

    τ̂_MMSE = ∫ θ exp(J(θ)/σ²) dθ / ∫ exp(J(θ)/σ²) dθ

computed from the same J grid. It needs σ² (the SNR). This is the estimator
with the smallest possible MSE, so it is the one the ZZB (a bound on *every*
estimator) should be compared with. Its limits: CRLB at high SNR, and T_a²/12
(the prior mean) at zero SNR (both verified in `sanity_checks.py`).

---

## 3. Cramér–Rao lower bound  (`bounds.crlb`)

For r[n] = s(t_n − τ) + w[n] with Gaussian w (Kay Vol. 1, §3.11, range estimation):

    I(τ) = (1/σ²) Σ_n s'(t_n − τ)²

    var(τ̂) ≥ σ² / Σ_n s'(t_n−τ)²  =  1 / (2 · E/N0 · β²),
    β² = Σ_n s'[n]² / Σ_n s[n]²   ( = ∫(2πf)²|S(f)|²df / ∫|S(f)|²df )

Precision comes from the **slope** of the waveform. The discrete β² (2.4483)
matches the frequency-domain value (2.4471) computed independently.

**Pilot dependence.** Averaged over random pilots, β² = 2.832, which matches the
isolated pulse's 2.838. But across sampled 64-chip pilots (400 draws), β² ranged
from 1.77 to 4.27, because what matters is the bandwidth of the actual pattern:
more sign changes give a lower CRLB. Our pilot (seed 2026) has β² = 2.448.
Both estimators therefore use the **same** pilot.

**Limits of the CRLB.** (i) It bounds the variance of *unbiased* estimators.
Below threshold every estimator is biased, so the CRLB does not apply, and at
low SNR it even exceeds the prior variance. (ii) It is *local*: it uses only
the likelihood's curvature at the true τ, so it cannot see distant competing
peaks. (iii) Estimators constrained to the prior window can beat it near the
window edges (measured MSE/CRLB ≈ 0.6 for τ within 0.05 chip of an edge),
which is a small (~1 %) effect overall.

---

## 4. Ziv–Zakai bound  (`bounds.zzb`)

**Idea.** Turn estimation into binary detection. An estimator that is accurate
to within h/2 can decide between "τ = φ" and "τ = φ + h", so its error cannot be
smaller than the best detector allows.

**General form** (Ziv–Zakai 1969; Chazan–Zakai–Ziv 1975; Bell et al. 1997):

    MSE ≥ ½ ∫₀^∞ h · V{ ∫ (p(φ) + p(φ+h)) · P_min(φ, φ+h) dφ } dh,
    V{g}(h) = max_{η ≥ 0} g(h + η)   (valley filling)

**Our case.** With a uniform prior of width T_a, ∫(p(φ)+p(φ+h))dφ = 2(1 − h/T_a).
The binary test is between two equally likely known signals in white Gaussian
noise:

    P_min(h) = Q(d/2),  d² = ‖s_φ − s_{φ+h}‖²/σ² = 2·Σs²·(1 − ρ(h))/σ²
    ⇒  P_min(h) = Q( √( E/N0 · (1 − ρ(h)) ) )

    ZZB = ∫₀^{T_a} h · V{ (1 − h/T_a) · Q(√(E/N0·(1 − ρ(h)))) } dh

ρ(h) is the normalised autocorrelation of the pilot waveform. We treat it as
independent of the reference point φ, an approximation that holds to ~2·10⁻⁵
here. Valley filling is numerically inactive for this pilot (sidelobes
|ρ| ≤ 0.13; < 0.01 % effect).

**Limits (verified):**
* E/N0 → 0: Q(0) = ½ gives ZZB → T_a²/12. No estimator can do *better* than
  the prior variance without data, and the prior mean achieves it.
* E/N0 → ∞: 1 − ρ(h) ≈ β²h²/2 gives ZZB → 1/(2·E/N0·β²) = CRLB.

**Why it captures the threshold.** The CRLB only sees ρ(h) near h = 0. The
ZZB integrates P_min(h) over all separations up to T_a, so it rises when large
separations become confusable, which is the mechanism that produces gross errors.

**How tight is it?** The ZZB bounds the MMSE estimator, and is built from
*pairwise* tests while a real estimator must reject *all* competing delays at
once. For TDE the literature finds it fairly tight (Mallat et al. 2014 report
≈ 0.5 dB to the MLE threshold for their signals; Xu, Baggeroer & Richmond 2004
find it the closest of the Bayesian bounds to MLE performance). The gap is
problem- and pilot-dependent; see `run_threshold_study.py`.

---

## 5. Early–late DLL  (`estimators.EarlyLateDLL`)

**Acquisition.** The first code period is correlated with the pilot on a
0.5-chip grid over the prior window, and the maximum is taken.

**Tracking.** The pilot repeats periodically (like a GNSS PRN code). Every
M = 4 chips (update j):

    E_j = Σ_n r[n] c_j(t_n − θ_j + d/2),   L_j = Σ_n r[n] c_j(t_n − θ_j − d/2)
    D_j = (E_j − L_j − b_j) / k̄
    θ_{j+1} = θ_j − K · D_j

* c_j: the local replica containing only block j's chips. It is correlated
  over a window positioned by the DLL's **own** timing θ_j, as a receiver's
  integrate-and-dump follows its local code clock.
* b_j: the block's noise-free discriminator output at zero error. It is non-zero
  because c_j also correlates with neighbouring chips, and the pattern repeats
  every period (uncorrected: MSE/CRLB = 2.07 at 30 dB). It depends only on the
  known pilot: recomputing it at other sampling phases changes nothing measurable.
* k̄: mean discriminator slope (constant gain).
* d = 0.5 chip early–late spacing.
* **Pull-in range** (measured, noise-free): the loop converges from initial
  errors up to **1.80 chips**, so acquisition errors up to that are recovered.

**Loop bandwidth.** For a first-order loop, N_eff = (2 − K)/K and
B_L·T_I = K/(2(2 − K)) (Stephens & Thomas 1995; checked by simulation). We set
K = 2/(N_eff + 1) with N_eff = 16, so 1/(2 B_L) = 64 chips = the MLE batch.

**What this matching does and does not mean.** It equalises the *steady-state
tracking jitter*: in lock, the DLL's error is driven by roughly 64 chips' worth
of noise, like the MLE. But the DLL **observes 640 chips** (160 updates × 4), ten
times the batch energy. Two consequences:
1. In the threshold region the DLL can use the extra data. Tracking recovers
   from wrong acquisitions within the pull-in range, so it is *not* an
   equal-information comparison there.
2. The 64-chip CRLB and ZZB are **not lower bounds on the DLL**. The valid bounds
   are those for its 640-chip observation (`zzb_dll_data` in the results, about
   10 dB lower). The DLL sits well above them (the loop discards most of the
   information).

**Linearised DLL theory (exact for the linear loop).** Near lock,
E_j − L_j − b_j ≈ k_j ε_j + ⟨w, g_j⟩ with g_j = early − late replica, so

    ε_{j+1} = a_j ε_j − (K/k̄)⟨w, g_j⟩,   a_j = 1 − K k_j/k̄
    ε_final = (Π a_j) ε_0 − ⟨w, G⟩,   G = Σ_j (K/k̄) Π_{i>j} a_i · g_j
    var(ε) = σ² ‖G‖²

This matches simulation to within ~2 % in the locked region. (An earlier,
naive per-update formula failed under a previous design that normalised each
block by its own slope; with constant gain a naive formula is also close.)

**Why the DLL is about 1.07 × the batch CRLB.** This holds **by construction of
B_L**, not as a bound. A narrower loop (N_eff = 64) gives about 0.26 × the
batch CRLB, because it averages more data. The residual 7 % comes from the
early–late difference being only a finite-difference approximation of the
waveform derivative s′ (plus ~2 % because the periodic stream's per-period
information differs from the 64-chip batch).

---

## 6. The threshold: different mechanisms

| | Batch MLE / MMSE | DLL |
|---|---|---|
| Locked / local regime | MSE → CRLB | MSE → σ²‖G‖² (≈ 1.07 × batch CRLB by design) |
| Global failure | **wrong peak** within the prior | **acquisition outside the pull-in range** (> 1.8 chips) |
| Dynamic failure | none (one-shot) | **loss of lock**: noise drives ε beyond the S-curve, then the loop random-walks |
| MSE far below threshold | MLE: **3.1 chips²** on pure noise (argmax piles up at the window edges, ~20 % of estimates exactly on an edge), above the T_a²/6 = 2.67 of a uniform guess; MMSE: → T_a²/12 = 1.33 | **grows linearly with observation time** (random walk); no single value |
| Covered by which bound | 64-chip ZZB (the MMSE gets closest) | 640-chip ZZB is valid but loose; loss of lock needs dynamic analysis (mean time to lose lock; Holmes & Biederman 1978) |

**Threshold definitions** (fixed before running):
* *MSE criterion:* highest E/N0 where MSE/CRLB first exceeds 2 (3 dB).
* *Gross-error criterion:* highest E/N0 where P(|e| > 0.5 chip) first exceeds
  1 %. One outlier changes this by 1/n, whereas it can multiply the MSE, so the
  gross criterion is more robust.

Both are reported with 95 % bootstrap confidence intervals (`tde/analysis.py`).

**Role of the prior width.** With T_a = 1 chip (about one correlation main lobe)
the threshold is **masked** rather than absent: the prior clips all errors, so the
MSE stays below about 0.17 chips² and the CRLB at low SNR exceeds the prior
variance. A visible threshold needs T_a ≫ the correlation width. More delay
uncertainty means more competing peaks and a later threshold. (Mallat et al.
2014, Part 2, note that some threshold definitions are nearly independent of the
prior for baseband signals; which threshold is meant matters.)

---

## 7. Key references

* J. Ziv, M. Zakai, "Some lower bounds on signal parameter estimation," *IEEE Trans. Inf. Theory*, 15(3):386–391, 1969.
* D. Chazan, M. Zakai, J. Ziv, "Improved lower bounds on signal parameter estimation," *IEEE Trans. Inf. Theory*, 21(1):90–93, 1975.
* K. L. Bell, Y. Steinberg, Y. Ephraim, H. L. Van Trees, "Extended Ziv–Zakai lower bound for vector parameter estimation," *IEEE Trans. Inf. Theory*, 43(2):624–637, 1997.
* H. L. Van Trees, K. L. Bell (eds.), *Bayesian Bounds for Parameter Estimation and Nonlinear Filtering/Tracking*, Wiley-IEEE Press, 2007.
* B. M. Sadler, R. J. Kozick, "A survey of time delay estimation performance bounds," *IEEE SAM Workshop*, pp. 282–288, 2006.
* J. P. Ianniello, "Time delay estimation via cross-correlation in the presence of large estimation errors," *IEEE Trans. ASSP*, 30(6):998–1003, 1982.
* A. Mallat, S. Gezici, D. Dardari, C. Craeye, L. Vandendorpe, "Statistics of the MLE and approximate upper and lower bounds — Parts 1 & 2," 2014 (arXiv:1405.4645, arXiv:1405.4647).
* W. Xu, A. B. Baggeroer, C. D. Richmond, "Bayesian bounds for matched-field parameter estimation," *IEEE Trans. Signal Process.*, 52(12):3293–3305, 2004.
* C. D. Richmond, "Mean-squared error and threshold SNR prediction of maximum-likelihood signal parameter estimation with estimated colored noise covariances," *IEEE Trans. Inf. Theory*, 52(5):2146–2164, 2006.
* A. Weiss, E. Weinstein, "Fundamental limitations in passive time delay estimation — Part I: Narrow-band systems," *IEEE Trans. ASSP*, 31(2):472–486, 1983. *(Passive TDE with random signals; context only, not our known-pilot case.)*
* J. W. Betz, K. R. Kolodziejski, "Generalized theory of code tracking with an early-late discriminator, Parts I & II," *IEEE Trans. Aerosp. Electron. Syst.*, 45(4):1538–1564, 2009.
* S. A. Stephens, J. B. Thomas, "Controlled-root formulation for digital phase-locked loops," *IEEE Trans. Aerosp. Electron. Syst.*, 31(1):78–95, 1995.
* W. K. Holmes, L. Biederman, "Delay-lock-loop mean time to lose lock," *IEEE Trans. Commun.*, 26(11):1549–1557, 1978.
* E. D. Kaplan, C. J. Hegarty (eds.), *Understanding GPS/GNSS: Principles and Applications*, 3rd ed., Artech House, 2017.
* S. M. Kay, *Fundamentals of Statistical Signal Processing, Vol. I: Estimation Theory*, Prentice Hall, 1993 (§3.11 range estimation; Ch. 7 MLE).

Page numbers were verified from publisher or author listings where available
(see `docs/validation_report.md`); re-check the exact pages before submission.
