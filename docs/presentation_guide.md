# Presentation guide — 70% evaluation

**Length:** 13 slides, ~15 minutes + questions. **Rule:** every slide says what is *done* vs *proposed*.
Figures are in `figures/`. Threshold numbers come from notebook 06. Question IDs refer to `docs/question_bank.md`.

---

## Slide 1 — Title
**Headline:** Threshold effects in time-delay tracking: CRLB vs Ziv–Zakai bound for a DLL and a batch MLE (BPSK)
**Content:** names, course, mentor, date; one line: *"Status: core study complete (Sprints 1–6); proposal for the
remaining 30%."*
**Say:** "We study when delay estimation suddenly fails at low SNR, and which theory predicts it."

## Slide 2 — Problem and motivation
**Headline:** Delay estimation fails abruptly at low SNR — the classical bound does not see it coming
**Content:** $r(t)=s(t-\tau)+w(t)$; uses (GNSS ranging, radar, synchronisation); the two estimators (batch MLE,
DLL) and two bounds (CRLB, ZZB). Small inset: `s1_objective_one_trial.png` (right panel: wrong peak wins).
**Say:** "At high SNR the estimate sits on the correlation peak. At low SNR a noise bump wins and the estimate jumps
chips away. That jump is the threshold effect."
**Triggers:** A1, A2.

## Slide 3 — Model and why each choice
**Headline:** One model, used identically by every estimator and bound
**Content:** table: 64-chip known pilot, raised-cosine α = 0.5 (why: rectangular ⇒ CRLB = 0), 4 samples/chip
(Nyquist), analytic delay, AWGN, uniform prior $T_a = 4$ chips (why: ZZB needs a prior), SNR = batch $E/N_0$,
metric = MSE. Figure: `s0_signal_model.png` or `s0_autocorrelation.png`.
**Say:** "The autocorrelation is the key picture: its main-lobe *curvature* sets the CRLB, its side-lobes and the
prior width set the ZZB."
**Triggers:** B1–B5, A6.

## Slide 4 — Literature and gap
**Headline:** The theory is known; what is open is the tracking (dynamic) side
**Content:** three columns — *Known:* CRLB/MLE (Kay), ZZB for delay (Ziv–Zakai, Bell et al. 1997), MLE threshold
prediction (Mallat et al. 2014), batch ZZB for GNSS code tracking (Emmanuele 2012, Liu 2016), DLL MTLL (Holmes 1978,
Hasson 1995, Groh 2012). *Partly open:* ZZB looseness vs the MMSE estimator; pilot-to-pilot spread. *Open (as far
as we found):* bounds with tracking dynamics, linked to loss of lock.
**Say:** "We position our 70% as a verified reproduction plus specific observations; the open part is our proposal.
We found no prior work linking a dynamic bound to loss of lock — we don't claim it is first."
**Triggers:** F1–F4. **Be careful:** never say "novel" without "as far as we found".

## Slide 5 — Objectives and hypotheses (original)
**Headline:** What we set out to show
**Content:** (1) MLE attains CRLB above threshold; (2) a threshold exists and the CRLB misses it; (3) the ZZB captures
it and → CRLB at high SNR; (4) DLL and MLE differ for identifiable reasons.
**Say:** "Each was tested; slide 11 says which held and which were challenged."

## Slide 6 — Batch MLE (Sprint 1) — math → code
**Headline:** The MLE is a matched filter
**Content:** Gaussian likelihood → $J(\theta)=\sum r\,s(t-\theta)-\frac12\sum s^2$; arg-max on 1/64-chip grid +
parabola. Noise-free error 0.02 % of CRLB.
**Say:** derive in three lines (C1). "σ² drops out: the MLE needs no SNR knowledge."
**Triggers:** B6, B7, C1.

## Slide 7 — CRLB (Sprint 3)
**Headline:** Precision comes from the waveform's slope
**Content:** score → Fisher information $\frac1{\sigma^2}\sum s'^2$ → $\text{CRLB}=\frac1{2(E/N_0)\beta^2}$; β² checked by
FFT (0.05 %); figure `s3_crlb_vs_mle.png`; pilot dependence (`s3_beta2_pilots.png`, β² from 1.77 to 4.27).
**Say:** "Three reasons it fails below threshold: unbiasedness, locality, and it ignores the prior."
**Triggers:** A4, A5, C2.

## Slide 8 — Ziv–Zakai bound (Sprint 4)
**Headline:** Turning estimation into detection between delays $h$ apart
**Content:** the 6 steps compressed: $\text{MSE}=\int\frac h2P(|e|\ge\frac h2)dh$ → estimator-as-detector ≥ $P_{\min}$ →
$P_{\min}=Q(\sqrt{(E/N_0)(1-\rho(h))})$ → uniform prior $2(1-h/T_a)$ → final formula; limits $T_a^2/12$ and CRLB.
Figure: `s4_zzb_intuition.png` (which $h$ matter at which SNR).
**Say:** the plain-language explanation from NB04 §4.4.
**Triggers:** C3–C8. *Expect to derive one step on the board — practise step 3 and the high-SNR limit.*

## Slide 9 — DLL (Sprint 2)
**Headline:** A feedback tracker: S-curve, loop gain, pull-in
**Content:** $E-L\approx-2R'(d/2)\varepsilon$; S-curve figure `s2_dll_s_curve.png` (pull-in ±1.80 chips); loop
$\theta_{j+1}=\theta_j-KD_j$, $N_{\rm eff}=(2-K)/K=16$; tracking runs `s2_dll_tracking_run.png`.
**Say:** "K is chosen so the locked DLL averages as much signal as the 64-chip batch — so the steady-state match with
the CRLB is by design, not a discovery."
**Triggers:** C9, C10, D1–D6.

## Slide 10 — Results: four curves (Sprint 5)
**Headline:** All methods follow the CRLB at high SNR; below threshold they fail for different reasons
**Content:** `s5_four_curves.png`. Bullet: MLE/CRLB = 0.997 above 14 dB; DLL ≈ 1.07×; MMSE closer to ZZB at low SNR;
DLL far worse below threshold; DLL above the ZZB for its own 640-chip data.
**Say:** "Note the dotted aqua curve: the 64-chip ZZB is *not* a bound for the DLL, because the DLL uses 640 chips.
We plot the bound that actually applies to it."
**Triggers:** D3, E4.

## Slide 11 — Results: the threshold, precisely (Sprint 6)
**Headline:** The threshold, located to ±0.3 dB — and the gap to the bound split in two
**Content:** `s6_threshold_precise.png`; table from NB06 §6.4 (ZZB, MMSE, MLE, DLL with 95 % CIs; the two criteria;
five pilots).
**Say:** the NB06 §6.6 interpretation: what is significant, what is not, how much of the gap is the bound vs the
estimator, and that the gap depends on the pilot.
**Triggers:** E1–E3, D8.
**Key numbers (pilot 2026):** ZZB 5.95 → MMSE 8.54 [8.35, 9.07] → MLE 9.52 [9.00, 9.66] dB; DLL 9.38 [8.72, 9.98].
Bound looseness 2.59 dB, MLE sub-optimality 0.85 dB [0.25, 1.24]. DLL vs MLE: MSE criterion no difference;
gross-error criterion DLL 0.68 dB more robust. Five pilots: MLE − ZZB 2.7–5.4 dB.
**Hypotheses status:** (1) ✓; (2) ✓; (3) ✓ qualitatively, ✗ quantitatively (ZZB ≈ 2.6 dB optimistic even for the
best estimator); (4) ✓ mechanisms identified; DLL vs MLE answer depends on the criterion (more data helps near
threshold, loss of lock hurts below it).

## Slide 12 — Limitations
**Headline:** What our results do *not* show
**Content:** synthetic data; constant delay (no Doppler); coherent, known pilot; one loop design; batch ZZB not a
DLL bound; pilot-to-pilot spread measured but not explained; novelty of the extension uncertain.
**Say:** "We'd rather list these than have you find them."

## Slide 13 — Proposal for the remaining 30%
**Headline:** Where does a delay *tracker* lose information? A dynamic decomposition
**Content:** state $d_k=d_{k-1}+w_k$; optimal tracker = exact grid Bayes filter; PCRB recursion
$J_k=J_D+1/(q+1/J_{k-1})$ and $K_\infty$ (DLL = fixed-gain Kalman in the linear regime); Markov-chain MTLL from the
S-curve; hypotheses H1–H4; 3-week roadmap (docs/research_proposal.md §3–4).
**Say:** "Our static study split the gap into bound looseness and estimator sub-optimality. The proposal does the same
for tracking, and adds the third ingredient only a loop has — loss of lock — and tests whether its exit-time statistics
predict the threshold."
**Triggers:** F2–F5, G1–G4. **Do not** show any proposal experiment as a result.

---

## One-page equation sheet (memorise the logic, not the symbols)
| Quantity | Formula | Where from |
|---|---|---|
| Model | $r[n]=s(t_n-\tau)+w[n]$, $w\sim\mathcal N(0,\sigma^2)$ | NB00 |
| SNR | $E/N_0=\sum s^2/(2\sigma^2)$ | $\sigma^2=N_0f_s/2$, $E=T_s\sum s^2$ |
| MLE | $\arg\max_\theta \sum r\,s(t-\theta)-\frac12\sum s(t-\theta)^2$ | expand Gaussian log-likelihood |
| CRLB | $1/(2(E/N_0)\beta^2)$, $\beta^2=\sum s'^2/\sum s^2$ | Fisher information $\frac1{\sigma^2}\sum s'^2$ |
| ZZB | $\int_0^{T_a}h\,\mathcal V\{(1-\frac h{T_a})Q(\sqrt{(E/N_0)(1-\rho(h))})\}dh$ | tail integral + detection bound |
| ZZB limits | $T_a^2/12$ (no signal), CRLB (strong signal) | $Q(0)=\frac12$; $1-\rho\approx\frac12\beta^2h^2$ |
| DLL | $\theta_{j+1}=\theta_j-K(E-L-b_j)/\bar k$, $N_{\rm eff}=(2-K)/K$ | Taylor of $R$; geometric impulse response |
| Proposal | $J_k=J_D+1/(q+1/J_{k-1})$; $P_\infty=\frac12(-q+\sqrt{q^2+4qR})$ | PCRB for a random walk; scalar Riccati |

## Study order tonight (≈ 3–4 h)
1. NB00 (30 min) → NB01 (40) → NB03 (40) → NB04 (60, derive steps 1–6 yourself on paper) → NB02 (30) → NB06 (20).
2. Question bank A, C, F (45 min) — say answers out loud.
3. Research proposal §0, §3.1–3.4 (30 min).
