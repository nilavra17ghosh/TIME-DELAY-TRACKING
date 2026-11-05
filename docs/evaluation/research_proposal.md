# Research directions for the remaining 30% — comparison and proposal

*Prepared 2026-11-05 for the 70% evaluation. Builds on `docs/surveys/02_literature_survey_2026-10-13.md` (Oct 13) plus a new check
(Nov 5) of the interdisciplinary angles. Tags: **[V]** = source opened or listed with matching abstract this or last
round; **[B]** = standard reference not re-opened.*

---

## 0. What the literature already covers — read this first

Our batch results (MLE, MMSE, CRLB, ZZB, threshold) are **reproductions** of known theory (see the survey). One new
finding from this round tightens the picture further:

* **ZZB-based bounds for GNSS code tracking already exist in batch (per-integration) form.**
  * Emmanuele, Luise, Zanier & Crisci, "Modified Ziv–Zakai bound for TOA estimation of GNSS signal-in-space,"
    ICASSP 2012, pp. 3041–3044 [V, listing];
  * Emmanuele, *Signal design and theoretical bounds for TOA estimation in GNSS applications*, PhD thesis, Univ. of
    Pisa, 2012 — ZZB used as a benchmark for acquisition **and tracking** [V, CORE pdf listing];
  * Liu et al., "The Ziv–Zakai bound for TOA estimation of new generation GNSS signals," ION ITM 2016 [V, listing];
  * a McGill thesis, *GNSS modulation: a unified statistical description with application to tracking bounds*,
    which derives code-tracking thresholds from the ZZB [V, listing].

  Our October survey had missed these — a reminder that "we found nothing" is weak evidence of novelty.

**What remains unaddressed, as far as both rounds of searching show:** these works bound a *single* integration
period (a static delay). None models the **loop's dynamics** (the delay evolving and the tracker carrying
information from epoch to epoch), and none connects a performance bound to **loss of lock / mean time to lose lock
(MTLL)**, which is analysed separately with exit-time theory (Holmes & Biederman 1978; Hasson & Bobrovsky 1995;
Groh et al. 2012). Recursive bounds that *could* do this exist only for generic filtering problems (posterior CRB,
Tichavský et al. 1998; sequential WWB, Rapoport & Oshman 2007; marginal WWB, Fritsche et al. 2015). Confidence that
this is open: **medium-low** — say "we found no prior work that…" and never "this is the first…".

---

## 1. Five candidate directions

Each direction is judged on whether it grows *naturally* out of the problem. The connection to your background
(robotics, state estimation, ML, autonomous systems) is noted honestly, including where it would be forced.

### D1 — Dynamic decomposition of the tracking threshold (estimation theory ↔ state estimation / control)
* **Question.** When the delay *moves*, how much of a DLL's low-SNR failure is (a) fundamental, (b) due to
  sub-optimal estimation, (c) due to the loop's structure (loss of lock)?
* **Why natural.** It is exactly our 70% result (static: ZZB vs MMSE vs MLE) extended to the dynamic setting the
  project title is about ("time-delay *tracking*").
* **Interdisciplinary link.** A DLL *is* a feedback controller around an estimator; the optimal tracker is a Bayes
  filter — the same machinery as robot state estimation (Kalman / grid / particle filters).

### D2 — Physics-informed loss-of-lock risk prediction (safety-critical autonomy, uncertainty estimation)
* **Question.** Can a tracker warn *before* it loses lock, better than classical lock detectors (C/N0 estimators,
  correlator-power tests)?
* **Why natural.** Loss of lock is our observed failure mode; an autonomous vehicle needs to know when its range
  measurement is about to become garbage.
* **Status of literature.** Loss-of-lock *detection* is classical (correlator power / C/N0 tests; patent literature
  notes their weakness under attenuation) [V]; ML is used mostly for integrity/protection levels and NLOS detection
  ("Smart GNSS integrity monitoring for road vehicles," survey 2025 [V]); vector-tracking work defines loss-of-lock
  metrics [V]. We found no work predicting DLL loss of lock from a posterior/physics-based risk measure.

### D3 — Loop bandwidth as a constrained decision problem (control / RL)
* **Question.** Choose the loop gain online to minimise tracking error *subject to* a loss-of-lock risk constraint.
* **Status.** Adaptive bandwidth control is well studied: Fraunhofer's loop-bandwidth control algorithm (LBCA) [V],
  fuzzy and adaptive-gain PLLs [V], adaptive Kalman tracking [V], and an RL framework that tunes the LBCA
  (ION, proof of concept) [V]. **Novelty is low**; an RL agent here would be the "add a network" pattern you asked
  us to avoid unless it is derived from a constrained-optimal-control formulation.

### D4 — Pilot/code selection for threshold robustness (information theory / optimisation)
* **Question.** Our ZZB–MLE gap varies 2.7–5.4 dB across five random 64-chip pilots (notebook 06), mostly because the ZZB threshold moves with β². Can per-code theory (MIE, ZZB
  integrand) predict it, and pick codes that push the threshold down?
* **Status.** ZZB-driven *waveform/pilot* design exists (Xiong & Liu, IEEE SPL 2023; Graff & Humphreys 2024 for OFDM)
  [V]; per-pulse MIE thresholds exist (Mallat et al. 2014) [V]. Discrete PN-code selection by threshold is not
  covered — **incremental**.

### D5 — Threshold-aware ranging likelihoods for robot localization (robotics)
* **Idea.** Use the physics-derived outlier probability $P_{\rm out}(\text{SNR})$ as the mixture weight of a
  heavy-tailed range likelihood in a robot's filter; test on the UTIL UWB dataset (Zhao et al., IJRR 2024) [V].
* **Honest critique.** In real UWB, outliers are dominated by **NLOS/multipath**, not by the thermal-noise threshold
  we study — the physics does not transfer cleanly, and UTIL's localisation benchmark already uses χ² gating [V].
  This connection would be **forced**. (Vision-language models and imitation learning have no natural link to this
  problem; we do not propose them.)

---

## 2. Comparison

| | D1 dynamic decomposition | D2 loss-of-lock prediction | D3 adaptive bandwidth | D4 code selection | D5 robot UWB |
|---|---|---|---|---|---|
| Novelty | **medium** (gap narrowed, see §0) | medium | low | low–medium | low (physics mismatch) |
| Technical depth | **high** | medium | medium | medium | medium |
| Mathematical sophistication | **high** (Bayes filter, PCRB recursion, Markov-chain exit times) | medium | medium–high if done properly | medium | medium |
| Feasible in remaining time | **yes for the core** (1-D state ⇒ exact grid filter) | yes (reuses D1) | partly | **yes (≈1 week)** | risky (real data, NLOS) |
| Compute | low–moderate (CPU) | low | moderate | low | moderate |
| Data | synthetic (ground truth required by theory) | synthetic | synthetic | synthetic | real (UTIL) |
| Validation possible | **strong** (bounds + exact optimum + simulation) | ROC / lead time | MSE vs constraint | threshold vs prediction | localisation error |
| Paper potential | **conference/workshop** (ION GNSS+, ICASSP) if the bound–MTLL link holds | short paper / add-on | low | workshop | low |
| Builds on our 70% | **directly** | via D1 | via S2 | via S6 | weakly |

---

## 3. Recommendation: D1, with D2 as its applied component

**Title:** *Where does a delay tracker lose information? A Bayesian decomposition of the code-tracking threshold.*

### 3.1 Central research question
For a delay that evolves in time, how does the mean-squared error of (i) the optimal Bayesian tracker, (ii) its
linearised version, and (iii) the early–late DLL compare with the posterior CRB as $E_c/N_0$ falls — and can the
DLL's extra loss be predicted quantitatively from its loss-of-lock statistics?

### 3.2 Hypotheses (each testable, each could fail)
* **H1 (linear regime).** Above threshold, a DLL whose gain equals the steady-state Kalman gain attains the
  posterior CRB; so does the optimal tracker.
* **H2 (threshold).** Below some $E_c/N_0$ the optimal tracker leaves the PCRB; the DLL leaves it **earlier**, and the
  extra loss disappears when we condition on the DLL staying in lock.
* **H3 (loss of lock is predictable).** The DLL's MTLL computed from a Markov-chain model built on our measured
  S-curve matches simulation, and the DLL's threshold is where its MTLL falls to the observation horizon.
* **H4 (applied, D2).** The optimal tracker's posterior gives a loss-of-lock risk score that predicts DLL loss of lock
  earlier (better ROC) than a classical correlator-power lock detector.

### 3.3 Intuition
In the static problem we already decomposed the gap to the bound: ZZB → MMSE (bound looseness, ≈ 2.6 dB) → MLE
(estimator sub-optimality, 0.85 dB, CI 0.25–1.24; notebook 06). Tracking adds two things: **memory** (the prior at each step is the previous
posterior, which shrinks the effective search window — so the threshold should move to lower SNR) and **feedback**
(the DLL commits to one estimate and can be pushed outside its pull-in range — loss of lock). D1 measures each effect
separately.

### 3.4 Mathematical formulation
**State and observation.** Per update $k$ (one block of $M$ chips):
$$ d_k = d_{k-1} + w_k,\quad w_k\sim\mathcal N(0,q);\qquad r_k[n] = s_k(t_n-d_k) + v_k[n],\quad v_k[n]\sim\mathcal N(0,\sigma^2). $$
**Likelihood** (same algebra as Sprint 1): $\ln p(\mathbf r_k\mid d) = J_k(d)/\sigma^2 + \text{const}$,
$J_k(d)=\sum_n r_k[n]s_k(t_n-d)-\tfrac12\sum_n s_k(t_n-d)^2$.

**Optimal tracker (MMSE).** Bayes recursion
$$ p(d_k\mid\mathbf r_{1:k}) \propto p(\mathbf r_k\mid d_k)\int \mathcal N(d_k;d_{k-1},q)\,p(d_{k-1}\mid\mathbf r_{1:k-1})\,dd_{k-1}. $$
Because the state is **one-dimensional**, this can be computed essentially exactly on a fine grid (a point-mass
filter, Bucy & Senne 1971 [B]): prediction = convolution with a Gaussian, update = multiplication by
$\exp(J_k(d)/\sigma^2)$. This removes particle-filter approximation error from the comparison.

**Posterior CRB** (Tichavský et al. 1998). For this model the measurement information is constant,
$J_D=\sum_n s_k'(t_n)^2/\sigma^2$ (Sprint 3), and the recursion is
$$ J_k = \big(q + J_{k-1}^{-1}\big)^{-1} + J_D , $$
whose steady state is the scalar Kalman–Riccati solution: with $R=1/J_D$,
$P_\infty = \tfrac12\big(-q+\sqrt{q^2+4qR}\big)$, $\;K_\infty = (P_\infty+q)/(P_\infty+q+R)$.
**Link to the DLL:** our loop $\theta_{k+1}=\theta_k-KD_k$ with $D_k\approx\varepsilon_k+n_k$ *is* a fixed-gain
Kalman filter in the linear regime; choosing $K=K_\infty$ makes it optimal there (H1). This turns our Sprint-2 gain
choice into a principled one.

**Loss of lock as a Markov chain (our own contribution to compute).** Using the noise-free S-curve $D(\varepsilon)$
(Sprint 2) and the discriminator noise variance $\sigma_D^2$ (Sprint-5 theory),
$\varepsilon_{k+1} = \varepsilon_k - K D(\varepsilon_k) - K n_k + w_{k+1}$ is a 1-D Markov chain. Discretising
$\varepsilon$ on a grid inside the lock region $|\varepsilon|<R_{\rm pull}$ ($\approx1.8$ chips), the mean exit times
$m(\varepsilon)$ satisfy the linear system
$$ m(\varepsilon) = 1 + \sum_{\varepsilon'} P(\varepsilon\to\varepsilon')\,m(\varepsilon'),\qquad m=0 \text{ outside the lock region}, $$
so $\text{MTLL}=m(0)$ follows from one matrix solve. (Most published MTLL results are continuous-time, second-order
loops at small-noise asymptotics; this discrete first-order version uses our measured S-curve directly.)

**Risk score (H4).** From the optimal tracker's posterior:
$\rho_k = P\big(|d_k-\hat\theta^{\rm DLL}_k| > R_{\rm pull} \,\big|\, \mathbf r_{1:k}\big)$ — the posterior
probability that the DLL is already outside its pull-in range.

**(Stretch) a sequential large-error bound.** Apply the Bell et al. (1997) ZZB construction to the one-step predictive
prior at each $k$, or compute the marginal WWB of Fritsche et al. (2015) with our grid posterior. This is the most
novel and the most uncertain part; it is optional.

### 3.5 What we would implement ourselves
1. Time-varying delay in the signal model (random walk per block) — small change to `sdet/model.py`.
2. Grid (point-mass) Bayes tracker — new, ~60 lines.
3. Posterior CRB recursion and steady-state Kalman gain — ~15 lines.
4. DLL with $K=K_\infty$ and other gains — existing class, new parameter.
5. Markov-chain MTLL solver from the S-curve — new, ~40 lines.
6. Risk score + classical lock detector baseline — new, ~40 lines.
7. Experiments as notebooks 07–08.

### 3.6 Experiments
| # | Experiment | Tests | Metric |
|---|---|---|---|
| E1 | RMSE vs $E_c/N_0$ for PCRB, grid tracker, linearised Kalman, DLL($K_\infty$) at 2–3 values of $q$ | H1, H2 | RMSE, threshold with bootstrap CI |
| E2 | Same, DLL RMSE *conditioned on lock* | H2 | conditional RMSE |
| E3 | MTLL: Markov chain vs simulation, over $E_c/N_0$ and $K$ | H3 | relative error; MTLL = horizon crossing vs DLL threshold |
| E4 | Loss-of-lock prediction: risk score vs correlator-power detector | H4 | ROC/AUC, detection lead time |
| A1 | Ablation: early–late spacing $d$ | sensitivity | threshold shift |
| A2 | Ablation: pilot (5 codes) | robustness | threshold spread |

**Baselines:** DLL with standard gains; linearised (extended) Kalman tracker; classical lock detector. All share our
verified signal model, so differences are attributable to the estimator.

### 3.7 Expected outcomes and failure cases
* **Likely:** H1 holds (it is close to textbook); the grid tracker's threshold is lower than the DLL's.
* **Uncertain:** H2's "conditioning on lock removes the gap" may only partly hold (a locked DLL near threshold also
  has inflated jitter). H3 may need the continuous-to-discrete approximation refined. H4 could fail if the classical
  detector already reacts quickly at high SNR — then report where each wins.
* **Honest limitation:** synthetic data only; one-dimensional delay (no Doppler); coherent tracking.

### 3.8 How it differs from existing work
| Existing | What it does | What D1 adds |
|---|---|---|
| Emmanuele 2012; Liu 2016 (batch ZZB for GNSS) | bound one integration period | dynamics across epochs; decomposition into bound / estimator / loop |
| MTLL theory (Holmes 1978; Hasson 1995; Groh 2012) | exit times of given loops | link exit time to the MSE threshold and to the optimal tracker |
| PCRB / sequential WWB (Tichavský 1998; Rapoport 2007; Fritsche 2015) | generic recursive bounds | applied to PN code tracking, benchmarked against an exact optimal tracker |
| Optimal-MTLL phase tracking (Fischler & Schuss 2007) | phase (periodic) | code delay (non-periodic, code-dependent sidelobes) |

---

## 4. Roadmap for the remaining 30% (assumes ≈ 3 weeks; adjust to the real calendar)

| Week | Work | Deliverable | Priority | Depends on |
|---|---|---|---|---|
| 1 | Random-walk model; PCRB recursion; grid tracker; DLL with $K_\infty$; **E1** | notebook 07, RMSE-vs-SNR figure | **P0** | — |
| 1 | Verify grid tracker (noise-free, high-SNR = PCRB) | tests added to `tests/` | **P0** | grid tracker |
| 2 | Lock-conditioned RMSE (**E2**); Markov-chain MTLL (**E3**) | notebook 08, MTLL figure | **P0** | E1 |
| 2 | Risk score + classical detector (**E4**) | ROC figure | P1 | grid tracker |
| 2–3 | Ablations A1–A2 | table | P1 | E1 |
| 3 | (stretch) sequential ZZB / marginal WWB | derivation + curve | P2 | E1 |
| 3 | Report and final presentation | report, slides | **P0** | all P0 |

**Tomorrow we can show:** the completed static study (Sprints 1–6), the precise threshold measurement (notebook 06),
the literature position (§0), and this proposal with its formulation. **We must not show** any of E1–E4 as done.
