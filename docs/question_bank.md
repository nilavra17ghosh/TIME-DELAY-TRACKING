# Question bank for the 70% evaluation (with answers and follow-ups)

How to use: read the **short answer** aloud in your own words; know the **why** well enough to answer the
**follow-up**. Notebook references (NB0x §y) point to where the full derivation or evidence lives.
Threshold numbers come from the precise study in notebook 06 (pilot 2026 unless stated).

---

## A. Problem and concepts

**A1. In one sentence, what is your project?**
We estimate the delay of a known BPSK signal in noise with two estimators — a batch maximum-likelihood estimator and
an early–late delay-locked loop — and compare them with two lower bounds, the CRLB and the Ziv–Zakai bound, to
understand the *threshold*: the SNR below which errors suddenly become large.

**A2. What exactly is the threshold effect?**
Above threshold, errors are small and Gaussian — the estimate sits on the true correlation peak. Below it, noise
sometimes makes a *wrong* delay look best, and the estimator jumps chips away. The MSE then rises far faster than the
CRLB predicts. It is a change in the *type* of error (NB01 §1.4, §1.6).
* *Follow-up: "Is it a sharp transition?"* — No, it is a region a few dB wide; that is why we define it with a fixed
  criterion (MSE > 2·CRLB, or P(gross error) > 1 %) and give confidence intervals.

**A3. Why do you plot MSE and not variance?**
Below threshold the estimators are biased for a given delay (NB03 §3.5), and the ZZB is a bound on MSE. Variance would
hide the bias and could not be compared with the ZZB.

**A4. Why is the CRLB not valid at low SNR?** Three reasons (NB03 §3.6):
1. it bounds only *unbiased* estimators, and ours become biased;
2. it is *local*: it uses only the curvature of the likelihood at the true delay, so it cannot see wrong peaks;
3. it ignores the prior: below ≈ −8 dB it even exceeds the prior variance $T_a^2/12$, which no sensible bound should.
* *Follow-up: "So is the CRLB wrong?"* — No; it is a correct bound for a class of estimators our estimators leave at
  low SNR. It is not *informative* there.

**A5. Can any estimator beat the CRLB?** Yes, three ways, all seen in our data:
a biased estimator; an estimator constrained to the prior window near its edges (MSE/CRLB ≈ 0.63 there, NB03); and an
estimator using more data than the CRLB assumes (a narrow-band DLL uses 640 chips, the CRLB assumes 64).

**A6. Why does the ZZB need a prior, and is a uniform prior realistic?**
The ZZB is Bayesian: it bounds the MSE averaged over the prior. That is exactly why it stays finite at zero SNR
($T_a^2/12$). Uniform is the natural model after coarse acquisition ("the delay is somewhere in this window").
* *Follow-up: "How does the result depend on $T_a$?"* — With $T_a = 1$ chip the threshold is masked (no wrong peak to
  jump to); larger $T_a$ means more competing peaks and a later threshold (NB00 §0.2; literature: Mallat et al. 2014).

---

## B. Signal model and implementation choices

**B1. Why a raised-cosine pulse instead of rectangular chips?**
CRLB $=1/(2\,(E/N_0)\,\beta^2)$ with $\beta^2$ the mean-square bandwidth. Rectangular chips have infinite bandwidth
⇒ CRLB = 0, which is meaningless. Real receivers band-limit anyway.

**B2. Why 4 samples per chip?** The RC spectrum ends at $(1+\alpha)/2 = 0.75$ cycles/chip < Nyquist (2 cycles/chip),
so samples carry all information, and sums equal integrals (needed for $E/N_0=\sum s^2/(2\sigma^2)$).

**B3. Why did you taper the pulse?** A hard cut at ±6 chips left a tiny slope jump (≈0.005). With samples landing on
it, the sampled autocorrelation's curvature disagreed with β², so the ZZB did not converge to the CRLB. A smooth taper
over the last chip fixed it; the Nyquist property and β² (to 4 digits) are unchanged.

**B4. How is the delay applied? Interpolation?** No interpolation — the waveform is evaluated analytically at
$t_n-\tau$, so fractional delays are exact.

**B5. Why a known pilot, coherent, real baseband?** It is the cleanest setting in which the CRLB and ZZB have closed
forms; unknown data bits, carrier phase and Doppler are extensions (and ZZBs with unknown data exist for OFDM —
Graff & Humphreys 2024).

**B6. Why does the MLE need no SNR knowledge but the MMSE does?** In $\ln p = J/\sigma^2+$const, $\sigma^2$ scales the
whole function, so the arg-max is unchanged; the posterior mean weights by $\exp(J/\sigma^2)$, which depends on it.

**B7. Is your MLE exact?** It maximises $J$ on a 1/64-chip grid with parabolic refinement; the noise-free error is
0.02 % of the CRLB at 30 dB (NB01 §1.3), so it is effectively exact.

---

## C. Mathematics (expect "derive this on the board")

**C1. Derive the MLE.** Gaussian likelihood → expand the square → $\theta$-independent $\sum r^2$ drops →
maximise $J(\theta)=\sum r[n]s(t_n-\theta)-\frac12\sum s(t_n-\theta)^2$ (NB01 §1.1).

**C2. Derive the CRLB.** Score $=\frac1{\sigma^2}\sum w[n]\mu_n'$ (zero mean = regularity); Fisher information
$=\frac1{\sigma^2}\sum\mu_n'^2$ with $\mu_n'=-s'(t_n-\tau)$; substitute $\sigma^2=\sum s^2/(2E/N_0)$
⇒ $1/(2(E/N_0)\beta^2)$ (NB03 §3.1).
* *Follow-up: "Why is $\beta^2$ the Gabor bandwidth?"* — Parseval: differentiation multiplies the spectrum by $j2\pi f$.

**C3. Prove $\mathbb E[X^2]=\int_0^\infty 2x\,P(X\ge x)\,dx$ for $X\ge0$.**
$X^2=\int_0^X 2x\,dx=\int_0^\infty 2x\,\mathbf 1\{X\ge x\}dx$; take expectations and swap (Tonelli, non-negative
integrand).

**C4. Why is any estimator's tail probability lower-bounded by $P_{\min}$?** Use the estimator as a detector
between $\tau=\varphi$ and $\tau=\varphi+h$ (decide $H_1$ if $\hat\tau\ge\varphi+h/2$). Its errors are exactly
$\{e\ge h/2\}$ under $H_0$ and $\{e<-h/2\}$ under $H_1$; no detector beats the optimal one (NB04 §4.1, step 3).

**C5. Why $P_{\min}=Q(d/2)$?** Equal priors, two known signals in white Gaussian noise: the LRT statistic
$T=\langle\mathbf r,s_1-s_0\rangle$ is Gaussian with means $\|s_1-s_0\|^2$ apart and std $\sigma\|s_1-s_0\|$;
threshold half-way ⇒ $Q(\|s_1-s_0\|/2\sigma)$; and $\|s_1-s_0\|^2=2\sum s^2(1-\rho(h))$ ⇒ $Q(\sqrt{(E/N_0)(1-\rho)})$.

**C6. Show the ZZB → CRLB at high SNR.** Small $h$ only: $1-\rho(h)\approx\frac12\beta^2h^2$ (since
$\rho''(0)=-\beta^2$); $\int_0^\infty hQ(ah)dh=1/(4a^2)$ (parts); $a^2=\frac12(E/N_0)\beta^2$ ⇒ $1/(2(E/N_0)\beta^2)$.
* *Follow-up: "Why is $\rho''(0)=-\beta^2$?"* — $\frac{d^2}{dh^2}\sum s(t)s(t-h)\big|_0=\sum s\,s''=-\sum s'^2$
  (summation by parts for a band-limited, Nyquist-sampled signal).

**C7. Show the ZZB → $T_a^2/12$ at zero SNR.** $Q(0)=\frac12$, $\int_0^{T_a}\frac h2(1-\frac h{T_a})dh=\frac{T_a^2}{4}-\frac{T_a^2}6=\frac{T_a^2}{12}$.

**C8. What is valley filling and does it matter here?** Replacing the integrand by its non-increasing envelope; valid
because $P(|e|\ge h/2)$ is non-increasing in $h$. For our pilot (side-lobes ≤ 0.13) it changes the ZZB by < 0.01 %.

**C9. Derive the DLL's noise averaging.** $\varepsilon_{j+1}=(1-K)\varepsilon_j-Kn_j$ ⇒ impulse response $K(1-K)^k$
⇒ variance gain $K/(2-K)$ ⇒ equivalent to averaging $N_{\rm eff}=(2-K)/K$ updates; $B_LT=K/(2(2-K))$ (NB02 §2.3).

**C10. Why is "early minus late" proportional to the error?** $E-L=R(\varepsilon-\frac d2)-R(\varepsilon+\frac d2)\approx-2R'(\frac d2)\varepsilon$ (NB02 §2.1).

---

## D. DLL and the comparison with the MLE

**D1. Why $d = 0.5$ chip?** Narrower spacing approximates the derivative better (lower jitter for band-limited
signals) but shrinks the linear region; 0.5 is a common compromise. Our DLL jitter is ≈1.07× the CRLB; the excess is
the finite-difference approximation (NB05 §B).

**D2. Where does K = 2/17 come from?** We chose $N_{\rm eff}=64/4=16$ updates so the locked DLL averages as much
signal energy as the 64-chip batch: $K=2/(N_{\rm eff}+1)$.
* *Follow-up: "So the DLL ≈ CRLB result is by construction?"* — Yes, the *steady-state* match is by design, and we say
  so. A narrower loop goes below the 64-chip CRLB because it uses more data (NB05 §B).

**D3. Is the DLL–MLE comparison fair?** In steady state yes (equal effective energy). In the threshold region no: the
DLL observes 640 chips and can recover from bad acquisitions. That is why the only *valid* bound for the DLL is the
ZZB for its own 640-chip data, which we also plot (NB05 §C).

**D4. Isn't subtracting the offset $b_j$ cheating?** $b_j$ is computed from the known pilot alone. We recomputed it at
deliberately wrong timing phases and the DLL's MSE did not change (independent audit, archive validation report).
Without it the DLL is biased (MSE ≈ 2×CRLB at 30 dB).

**D5. Why an acquisition stage?** The S-curve restores only within ±1.8 chips; the prior is ±2 chips. Real receivers
acquire before tracking.

**D6. Why first-order and not second-order?** Our delay is constant (no Doppler), so a first-order loop has no
steady-state error. With Doppler we would need second order.

**D7. What is loss of lock and how do you measure it?** The error leaves the S-curve's restoring region and the loop
drifts. We classify each failed run as "acquired outside pull-in" or "lost lock while tracking" (they partition the
total), and we show the DLL's MSE grows with observation time below threshold (NB05 §D).

**D8. Which is better, DLL or MLE?** Above threshold: equal by design. Near threshold: their MSE thresholds are not
significantly different (DLL − MLE = −0.12 dB, CI [−0.87, 0.59]), but the DLL makes *fewer gross errors* (its 1 %
gross-error threshold is 0.68 dB lower, CI [0.40, 0.99]) because it sees 640 chips and recovers from bad
acquisitions. Far below threshold the DLL is much worse (loss of lock, MSE grows
with time). The DLL is cheap and real-time; the MLE needs a search.

---

## E. Experimental method and results

**E1. How many trials, and how do you know the numbers are reliable?** 2000 per point in notebooks 01/05; 20,000
(MLE/MMSE) and 4000 (DLL) in notebook 06; 95 % confidence intervals by bootstrap. Code verified by (i) a test
suite against analytic facts, (ii) an independent re-implementation (CRLB to 0.02 %, ZZB to 0.12 %), (iii) an
adversarial audit.

**E2. Your threshold definition (3 dB) seems arbitrary.** It was fixed before running experiments and is the standard
style of definition; we also report a second, outlier-robust definition (P(gross) > 1 %). Conclusions that hold under
both are the ones we state.

**E3. Why is the threshold so hard to pin down?** The MSE near threshold is dominated by rare large errors; one
2-chip error in 2000 trials adds 0.002 chips², comparable to the CRLB at 10 dB. Hence many trials and bootstrap CIs.

**E4. The ZZB is ≈ 3 dB below your MLE threshold. Is your ZZB wrong?** Both limits are reproduced analytically and
numerically, and an independent implementation agrees to 0.12 %. The gap has two parts: the MMSE estimator (the one the
ZZB really bounds) is 0.85 dB ahead of the MLE (CI 0.25–1.24), so ≈ 2.6 dB of the 3.6 dB gap is bound
looseness, from building it with *pairwise* tests while the estimator must beat *all* competing delays. The size
depends on the pilot (MLE − ZZB 2.7–5.4 dB over five pilots), mainly because the ZZB threshold, defined relative to
each pilot's CRLB, moves with the pilot's bandwidth β².
* *Follow-up: "Mallat et al. found their bound much tighter."* — Different signals (single pulse vs a PN code with
  code-dependent side-lobes) and threshold definitions; we do not claim a contradiction.

**E5. Why does the MLE's MSE exceed the prior variance at very low SNR?** On pure noise the arg-max piles up at the
window edges (~20 % of estimates exactly on an edge), giving ≈3.1 chips² vs $T_a^2/12=1.33$. The MMSE goes to 1.33.

**E6. Why synthetic data only?** The bounds need the exact signal model and ground truth. Real GNSS data adds Doppler,
multipath, data bits and clock effects that change the problem; that is future work.

---

## F. Novelty and contribution (the hardest questions — answer honestly)

**F1. What is your own contribution?**
Up to 70% the *theory* is known (Kay; Ziv–Zakai; Bell et al.; Betz–Kolodziejski). Our contribution is a carefully
verified, self-consistent study for one concrete model, plus three quantitative observations we did not find stated in
the literature for PN pilots: (i) the ZZB's looseness measured against the MMSE estimator (not just the MLE); (ii) its
pilot-to-pilot spread; (iii) a decomposition of DLL failures into acquisition vs loss of lock. The research
contribution we *propose* is the dynamic decomposition in `docs/research_proposal.md`.

**F2. Hasn't the ZZB already been applied to GNSS code tracking?**
Yes — Emmanuele et al. (ICASSP 2012, modified ZZB), Emmanuele's Pisa thesis (2012), Liu et al. (ION ITM 2016). They
bound a single integration period. Our proposal concerns tracking *dynamics* and *loss of lock*.

**F3. Isn't loss of lock already solved by MTLL theory?** MTLL theory gives exit times for a given loop (Holmes &
Biederman 1978; Hasson & Bobrovsky 1995; Groh et al. 2012). We did not find work linking MTLL to an MSE bound or to
the optimal tracker's threshold — that link is what we propose to test (H3).

**F4. How sure are you that it is novel?** Not sure — medium-low confidence, from targeted searches. We missed the
Emmanuele papers in our first survey, which is exactly why we phrase it as "we found no prior work that…".

**F5. Why not use deep learning / RL?** No component of the problem requires it; the optimal tracker is computable
exactly (1-D grid filter). Adaptive-bandwidth RL exists only as proof of concept and would add complexity without a
clear hypothesis.

---

## G. Future work

**G1. What if H2 (lock-conditioning explains the DLL's extra loss) fails?** Then the DLL's near-threshold loss is
partly *jitter inflation inside lock* (nonlinear S-curve region), which is itself a measurable finding; we would
quantify it with the Markov-chain model.

**G2. Why a grid filter, not a particle filter?** The state is one-dimensional, so a fine grid computes the Bayes
recursion essentially exactly (Bucy & Senne 1971) — no Monte Carlo approximation in the "optimal" baseline.

**G3. Connection to robotics/autonomous systems?** The DLL is a feedback estimator; the optimal tracker is the same
Bayes filter used in robot localisation; predicting loss of lock before it happens is an integrity problem for
autonomous navigation (GNSS integrity monitoring). We avoided forcing a robotics dataset because real UWB outliers are
dominated by NLOS, not by the thermal threshold we study.

**G4. How would this extend to a real receiver?** Add Doppler (second-order loop, 2-D state), carrier phase and data
bits; the grid filter becomes 2-D (still feasible) or a particle filter.
