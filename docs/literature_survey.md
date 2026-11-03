# Literature survey and novelty assessment — Project 13 (2026-10-13)

**Question:** Is what we have done so far *new*, or a reproduction of known work? If something is genuinely open,
how could it be done?

**Method.** (1) A deep-research workflow: 5 search angles, 20 primary sources fetched, 92 claims extracted, 25
claims checked by 3 independent verifiers each (23 confirmed, 2 refuted and excluded). (2) Targeted follow-up searches
on the candidate gap, because the workflow itself flagged its "open" conclusion as resting on absence of evidence.

**Tags:** **[V]** = verified this run (source opened, at least the abstract) · **[B]** = standard background reference,
not re-verified this run.

---

## 1. Bottom line

* **Everything up to Sprint 4, and most of Sprints 5–6, is a *reproduction* of known theory.** This is normal and
  expected for a course project: the value is in deriving and verifying it correctly for a specific model, not in
  novelty.
* **Two small *partial* gaps** exist (incremental, publishable at workshop level at best).
* **One genuine *open* direction** survives the targeted check, but in a **narrower form** than first assumed:
  > a **large-error, recursive Bayesian MSE bound** (sequential/marginal Weiss–Weinstein, or a new *sequential*
  > Ziv–Zakai bound) for **PN code-delay tracking**, and a quantitative link between that bound's threshold and the
  > **DLL loss-of-lock (MTLL)** threshold.
* **What already exists close to it** (so it must *not* be claimed as new): MTLL theory for DLLs (since 1978), and a
  fundamental-limit "optimal-MTLL filtering" theory for **phase** tracking (Schuss, Bobrovsky, Fischler). The open
  part is the *code-delay* setting and the *MSE-bound ↔ MTLL* connection.
* **Confidence that it is open: medium-low.** It is an absence-of-evidence finding (searches listed in §4). Check
  the items in §6 before claiming novelty in a report.

---

## 2. Our project component by component

| Our component | Status | Already done by |
|---|---|---|
| BPSK signal model, matched-filter MLE (Sprint 1) | **Done (textbook)** | Kay 1993 Ch. 7 [B] |
| CRLB for delay, Gabor bandwidth (Sprint 3) | **Done (textbook)** | Kay 1993 §3.11 [B] |
| ZZB derivation for delay; ZZB captures the sidelobe-driven threshold the CRLB misses (Sprint 4) | **Done** | Ziv–Zakai 1969, Chazan–Zakai–Ziv 1975, Bell et al. 1997 [B]; Weiss–Weinstein 1983, Ianniello et al. 1983 [B]; Graff–Humphreys 2024 [V]; Gusi-Amigó et al., ZZB for GNSS direct position estimation [V] |
| MLE threshold location and its prediction (method of interval errors, MIE) | **Done** | Mallat et al. 2014 Part 2 (closed-form thresholds, incl. raised cosine) [V]; Athley 2002 [V]; Richmond 2006 [B] |
| Threshold depends on roll-off and prior width, not on chip rate | **Done** | Mallat et al. 2014 Part 2, Sec. III-C [V] |
| Using the ZZB to rank or design ranging signals | **Done** | Graff–Humphreys 2024 (OFDM pilots) [V]; Xiong–Liu, IEEE SPL 2023 (ZZB-optimal waveform) [V] |
| Linear DLL theory (loop bandwidth, jitter); it fails in the large-error regime | **Done** | Betz–Kolodziejski 2009 [B]; Stephens–Thomas 1995 [B]; Cassel, NAVIGATION 2023 [V] |
| DLL loss of lock, mean time to lose lock (MTLL) | **Done** | Holmes–Biederman 1978 [B]; Welti–Bobrovsky 1989 [V]; Hasson–Bobrovsky 1995 [V]; Groh–Gentner–Selva 2012 [V]; Schuss et al., SIAM 2006 [V]; Olmos–Agustí c. 1992 [V] |
| Bayesian (particle-filter) delay tracking in GNSS | **Done** | Lentmaier et al. 2008 [V]; Closas 2012 [V] |
| Large-error Bayesian bound for *static* GNSS parameters | **Done** | WWB for weak-GNSS frequency estimation, NAVIGATION 71(3) 2024 [V] |
| **Spread of the ZZB–MLE gap across random PN pilots (2.5–5.5 dB)** | **Partially open** | MIE theory is per-autocorrelation (Mallat), so it *can* explain it, but no source studies the code-to-code spread |
| **ZZB looseness measured against the MMSE estimator, not the MLE** | **Partially open (low confidence)** | ZZB papers compare against the MLE or optimise the bound; none found that isolates looseness with the MMSE for TDE |
| **Recursive large-error bound for code tracking, linked to DLL loss of lock** | **Open (medium-low confidence)** | see §3–4 |

**One claim to avoid:** that "Mallat et al. 2014 already show the ZZB mislocates the MLE threshold" was **refuted
0–3** by the verifiers. Do not cite it as precedent for our 3 dB gap.

---

## 3. Annotated survey by theme

### 3.1 Static bounds and the threshold effect for time-delay estimation
* **Ziv & Zakai (1969); Chazan, Zakai & Ziv (1975); Bell, Steinberg, Ephraim & Van Trees (1997)** [B]: the ZZB and
  its improved/valley-filled form. This is exactly what Sprint 4 derives.
* **Weiss & Weinstein (1983); Ianniello (1982); Ianniello, Weinstein & Weiss (ICASSP 1983)** [B]: ZZB and large-error
  analysis for delay estimation; ZZB vs correlator performance.
* **Sadler & Kozick (2006)** [B]: survey of TDE performance bounds.
* **Mallat, Gezici, Dardari, Craeye & Vandendorpe (2014), Parts 1–2** [V] (arXiv:1405.4645, 1405.4647): MLE
  statistics and closed-form ambiguity/asymptotic thresholds for baseband and passband pulses via MIE. Closest prior
  work to our Sprint 1/6 threshold analysis. Limitation relative to us: a single-pulse autocorrelation, not a PN code
  with code-dependent sidelobes.
* **Athley (ICASSP 2002)** [V]: outlier-probability/MSE approximations for single-snapshot ML that remain accurate
  where the CRB fails (DOA; carries over to TDE by analogy).
* **Graff & Humphreys (2024)** [V] (arXiv:2403.12905): ZZB for TOA with OFDM pilots; the ZZB captures
  sidelobe-induced thresholding the CRLB misses; uses the ZZB as a convex design objective.
* **Xiong & Liu (IEEE SPL 2023)** [V] (arXiv:2310.02963): SNR-adaptive ZZB-optimal ranging waveform design.
* **Gusi-Amigó et al., "Ziv–Zakai bound for direct position estimation"** [V]: ZZB in a GNSS setting (static).

### 3.2 DLL code tracking and loss of lock
* **Betz & Kolodziejski (2009); Stephens & Thomas (1995); Kaplan & Hegarty** [B]: linearised code-tracking theory and
  digital loop design. This is the basis of our Sprint 2.
* **Cassel, NAVIGATION 70(4), 2023** [V]: states explicitly that loop bandwidth does not describe the nonlinear,
  large-error regime.
* **MTLL literature** [V]: Welti & Bobrovsky (1989, second-order coherent PN DLL); Hasson & Bobrovsky (1995,
  closed-form MTLL and an *analytical tracking threshold derived from MTLL*); Groh, Gentner & Selva (IEEE T-COM 2012,
  analytic MTLL via OU recasting); Schuss et al. (SIAM 2006, exit problem with 1/f noise); Olmos & Agustí
  (Monte Carlo MTLL). Mostly second-order or continuous-time loops at small-noise asymptotics; ours is a discrete
  first-order loop at low SNR.

### 3.3 Optimal filtering and loss of lock (the closest existing work to our gap)
* **Fischler & Schuss, "Nonlinear Filtering with Optimal MTLL" (arXiv math/0703524, 2007)** [V]: for small noise,
  the causal filter that *maximises* MTLL is Bellman's minimum-noise-energy filter. **Phase estimation only.** No
  MSE bounds.
* **Ezri, Bobrovsky & Schuss, "The exit problem in optimal non-causal estimation" (arXiv 0802.0414)** [V]: optimal
  smoothers keep lock much longer than the EKF.
* **Fischler & Bobrovsky, "Mean time to lose lock of phase tracking by particle filtering", Signal Processing (2006)**
  [V]: a particle filter improves MTLL over the PLL/EKF (≈ 0.94 dB in C/N0). Phase only.
* **Schuss, *Nonlinear Filtering and Optimal Phase Tracking* (Springer)** [V, catalogue entry]: book-length treatment
  of optimal lock time for **phase**.
* **Llorente et al. (arXiv 2606.23925, 2026)** [V]: Kalman-filter tracking loop tuned for *minimum pull-out
  probability*. Carrier phase, not code; no Bayesian bound.

### 3.4 Recursive Bayesian bounds for nonlinear filtering (the tools)
* **Tichavský, Muravchik & Nehorai (IEEE T-SP 1998)** [V, PDF found]: posterior CRB, recursive form. Small-error bound,
  so it cannot show a threshold.
* **Reece & Nicholson (FUSION 2005)** [V]: recursive Bobrovsky–Zakai, Bhattacharyya and WWB bounds for discrete-time
  filtering, tighter than the PCRB.
* **Rapoport & Oshman (IEEE T-SP 2007)** [V]: sequential WWB for Markovian systems. It requires linear dynamics but
  allows a **nonlinear measurement**, which is exactly a random-walk delay seen through a correlator.
* **Fritsche, Özkan, Orguner & Gustafsson (ICASSP 2015)** [V]: marginal WWBs, computable with particles, tighter
  than joint WWBs.
* **Nitzan, Routtenberg & Tabrikian (IEEE SPL 2018)** [V]: recursive Bobrovsky–Zakai-type bound for periodic
  filtering, validated on target tracking.
* **Van Trees & Bell (eds.), *Bayesian Bounds for Parameter Estimation and Nonlinear Filtering/Tracking* (2007)** [B].

---

## 4. The gap, stated precisely, with the evidence of absence

**What exists:**
1. MTLL / exit-time descriptions of DLL loss of lock (for given loops), and
2. optimal-MTLL theory for *phase* tracking, and
3. recursive large-error MSE bounds (WWB, BZB) validated on generic target tracking.

**What we could not find:** any work that
* applies a recursive/marginal large-error Bayesian bound (sequential WWB, BZB, or ZZB) to **PN code-delay
  tracking**, or
* relates such a bound's threshold to the **DLL loss-of-lock** threshold or MTLL, or
* derives a **sequential Ziv–Zakai bound** at all.

**Searches that returned nothing relevant** (2026-10-13): "posterior Cramér-Rao bound GNSS code tracking DLL
Bayesian bound"; "Weiss-Weinstein bound phase-locked loop cycle slip threshold tracking"; "Ziv-Zakai bound
sequential / recursive / nonlinear filtering tracking threshold"; "Bayesian lower bound code delay tracking GNSS
particle filter loss of lock MSE bound"; "Schuss optimal delay/code tracking MTLL". Results returned only static ZZB
papers, generic PCRB tracking papers, phase-only MTLL work, and patents.

**Why it is plausible that the gap is real:** the two communities use different languages. Estimation theorists
bound the MSE (CRB/WWB/ZZB); tracking-loop theorists compute the MTLL (exit times). Code tracking has a feature that
phase tracking lacks — a **non-periodic** measurement function with **code-dependent sidelobes** — so the
phase-tracking optimal-MTLL results do not transfer directly.

---

## 5. How it could be done: a concrete research plan

**Working title:** *A large-error Bayesian bound for PN code tracking and its relation to DLL loss of lock.*

**Model (one extension of our current model):** the delay becomes a random walk,
$d_k = d_{k-1} + w_k$, $w_k\sim\mathcal N(0,q)$; each epoch observes our 64-chip raised-cosine BPSK pilot at
$C/N_0$ with integration time $T$ (the measurement is the nonlinear correlator / likelihood we already have).

**Methods (in order of difficulty):**
1. **PCRB** (Tichavský 1998): the small-error reference; it should *not* show a threshold.
2. **Particle filter** with many particles: the near-MMSE estimator, and our analogue of the "optimal" tracker.
3. **DLL family** (our existing loop) across loop bandwidths $B_L$, plus measured **MTLL** and loss-of-lock
   probability over a horizon.
4. **Marginal WWB** (Fritsche 2015) by particle approximation, and/or the **sequential WWB** (Rapoport–Oshman 2007).
   Our model satisfies its linear-dynamics assumption. Test points must span 0.05–2 chips so the bound can "see"
   sidelobe-scale and pull-in-scale errors.
5. *(stretch, the most novel)* **Sequential ZZB**: apply the Bell et al. (1997) construction to the one-step
   predictive prior at each epoch.

**Hypotheses:**
* **H1:** the large-error recursive bound shows a $C/N_0$ threshold that the PCRB misses, set by $qT$ relative to the
  code's main-lobe width.
* **H2:** the DLL's RMSE leaves the bound several dB before the particle filter does, and the excess is explained by
  loss of lock.
* **H3:** the $C/N_0$ at which the DLL departs from the bound matches the $C/N_0$ at which its MTLL drops to the
  observation horizon.
* **H4 (links to our partial gap):** across 100+ random codes, the threshold correlates with the peak sidelobe
  inside the pull-in / prior region. This also tests the per-code MIE explanation of our 2.5–5.5 dB spread.

**Experiments:** sweep $C/N_0$, $q$ and $B_L$; record RMSE, lock-conditioned RMSE, MTLL and loss-of-lock probability;
repeat over many codes.

**Expected contribution:** a quantitative bridge between Bayesian filtering bounds and the MTLL description of the
code-tracking threshold.

**Realistic scope.** This is well beyond Sprints 5–8 of the course. A **minimum version** (PCRB + particle filter +
DLL with measured MTLL, steps 1–3) is about 2–3 weeks of work and already tests H2. The bounds in steps 4–5 are
research-level. Recommended: finish the course sprints as planned, then propose this as an extension.

**Smaller, safer contributions available now** (from the partial gaps):
* explain the 2.5–5.5 dB pilot-to-pilot spread with per-code MIE (Mallat's formulas applied to each code's
  autocorrelation) — about a week of work;
* report ZZB looseness against the MMSE (already measured: ≈ 0.6 dB MMSE gain vs ≈ 3 dB total gap) — essentially done.

---

## 6. Before claiming novelty: checks still to do
1. IEEE Xplore / Google Scholar: "Weiss-Weinstein bound" + "code tracking"; "Bayesian bound" + "delay lock loop";
   "sequential Ziv-Zakai"; "Ziv-Zakai" + "filtering".
2. Table of contents of Schuss, *Nonlinear Filtering and Optimal Phase Tracking*: is there a delay/code-tracking
   chapter?
3. Van Trees & Bell (2007) book: its nonlinear filtering/tracking chapters and examples.
4. Work by the Closas / Vilà-Valls / Fernández-Prades groups on Bayesian tracking bounds (2015–2026), which this run
   did not cover exhaustively.
5. Ask the course mentor whether they know of such work.

**Caveats of this survey:** several sources were verified from abstracts only (Hasson–Bobrovsky, Schuss SIAM,
Fritsche 2015, Nitzan 2018). The foundational references tagged [B] were not re-opened in this run. Coverage of
2025–2026 is thin.

---

## 7. Sources opened in this survey
* Mallat et al. 2014 Part 2: https://arxiv.org/pdf/1405.4647 · Part 1: https://arxiv.org/abs/1405.4645
* Athley 2002: https://research.chalmers.se/en/publication/15050
* Graff & Humphreys 2024: https://arxiv.org/pdf/2403.12905
* Xiong & Liu 2023: https://arxiv.org/pdf/2310.02963
* Gusi-Amigó et al., ZZB for direct position estimation: https://par.nsf.gov/servlets/purl/10119477
* WWB of frequency estimation for very weak GNSS signals (NAVIGATION 71(3), 2024): https://navi.ion.org/content/71/3/navi.654
* Cassel 2023 (NAVIGATION 70(4)): https://navi.ion.org/content/70/4/navi.598
* Welti & Bobrovsky 1989: https://cris.tau.ac.il/en/publications/optimization-of-a-second-order-pn-code-tracking-loop-using-the-me/
* Hasson & Bobrovsky 1995: https://bvbvf1.bib-bvb.de/vufind/PrimoRecord/cdi_ieee_primary_513842
* Groh, Gentner & Selva 2012: https://elib.dlr.de/79495
* Schuss et al. SIAM 2006: https://epubs.siam.org/doi/epdf/10.1137/050627666
* Olmos & Agustí: https://upcommons.upc.edu/bitstream/2117/100778/1/00245436.pdf
* Fischler & Schuss 2007: https://arxiv.org/abs/math/0703524
* Ezri, Bobrovsky & Schuss: https://arxiv.org/pdf/0802.0414
* Fischler & Bobrovsky 2006: https://www.sciencedirect.com/science/article/abs/pii/S0165168406002076
* Schuss, *Nonlinear Filtering and Optimal Phase Tracking*: https://books.google.com/books/about/Nonlinear_Filtering_and_Optimal_Phase_Tr.html?id=TSueVcSUsvYC
* Llorente et al. 2026: https://arxiv.org/abs/2606.23925
* Lentmaier et al. 2008: https://onlinelibrary.wiley.com/doi/10.1155/2008/372651
* Closas 2012: https://onlinelibrary.wiley.com/doi/10.1155/2012/359128
* Tichavský et al. 1998 (PCRB): https://www.ese.wustl.edu/~nehorai/paper/ieeetsp98-1.pdf
* Reece & Nicholson 2005: https://www.robots.ox.ac.uk/~sjrob/Pubs/fusion05b.pdf
* Rapoport & Oshman 2007: https://oshman.technion.ac.il/files/2016/04/Rapoport_Oshman_TSP_Part_1_V55_N5_May07.pdf
* Fritsche et al. 2015: https://avesis.metu.edu.tr/yayin/8a72eaf0-cecb-4992-9540-0096e6b80169/marginal-weiss-weinstein-bounds-for-discrete-time-filtering
* Nitzan, Routtenberg & Tabrikian 2018: https://cris.bgu.ac.il/en/publications/bobrovsky-zakai-type-bound-for-periodic-stochastic-filtering/
* Bobrovsky–Zakai bound for filtering, prediction and smoothing: https://www.researchgate.net/publication/327483968
