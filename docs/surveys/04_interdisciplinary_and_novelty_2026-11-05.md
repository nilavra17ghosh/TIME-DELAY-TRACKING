# Survey 04 — Interdisciplinary directions and a second novelty check (2026-11-05)

**Questions asked:** (a) which interdisciplinary extensions grow naturally from this project (robotics, safety-critical
autonomy, control/RL, information theory)? (b) has anything appeared, or been missed, that affects our novelty claim?

**Method:** targeted web searches (7 queries) on: GNSS loss-of-lock prediction / integrity; adaptive / RL loop
bandwidth; UWB ranging outliers in robot localisation (UTIL dataset); code/waveform design with the ZZB; 2024–2026
Bayesian bounds for code tracking; plus two searches to identify papers cited as prior work.

## Findings

### The most important finding — prior work we had missed
A 2024 NAVIGATION paper on Weiss–Weinstein bounds for weak GNSS signals notes that **"Denis (2009) and Emmanuele (2012)
investigated bounding GNSS code tracking with the ZZB."** Follow-up identified:
* A. Emmanuele, M. Luise, F. Zanier, M. Crisci, "Modified Ziv–Zakai bound for time-of-arrival estimation of GNSS
  signal-in-space," **ICASSP 2012**, pp. 3041–3044.
* A. Emmanuele, *Signal design and theoretical bounds for time-of-arrival estimation in GNSS applications*, PhD thesis,
  University of Pisa, 2012 — ZZB used for **acquisition and tracking** benchmarks (https://core.ac.uk/download/pdf/14704325.pdf).
* Liu et al., "The Ziv–Zakai bound for time-of-arrival estimation of new generation GNSS signals," **ION ITM 2016**
  (BPSK, BOC, MBOC, AltBOC) (https://www.ion.org/publications/abstract.cfm?articleID=13402).
* McGill thesis, *GNSS modulation: a unified statistical description with application to tracking bounds* — derives
  code-tracking thresholds from the ZZB (https://mcgill.scholaris.ca/items/984f6e2f-28ba-4363-a6ed-055e1105a9c4).
* "Denis (2009)" itself was **not** located.

**Consequence:** "ZZB for code tracking" is *not* open. What still appears unaddressed: bounds that include tracking
**dynamics** across epochs and connect to **loss of lock**.

### Interdisciplinary angles
| Angle | What exists | Verdict |
|---|---|---|
| Loss-of-lock prediction / integrity (safety-critical autonomy) | classical lock detectors (C/N0, correlator power) with known weaknesses; ML mainly for protection levels and NLOS detection ("Smart GNSS integrity monitoring for road vehicles", 2025 survey); loss-of-lock metrics for vector tracking | **natural**; posterior-based risk score not found → proposal component H4 |
| Adaptive loop bandwidth (control / RL) | Fraunhofer LBCA (adaptive bandwidth from discriminator statistics); fuzzy / adaptive-gain PLLs; adaptive Kalman tracking; an RL proof-of-concept tuning LBCA (ION) | **low novelty**; RL would be unjustified complexity |
| Code / waveform design via ZZB (information theory) | ZZB-optimal ranging waveform (Xiong & Liu, SPL 2023); ZZB OFDM pilots (Graff & Humphreys 2024) | incremental for discrete PN codes |
| Robot localisation with UWB (robotics) | UTIL dataset (Zhao et al., IJRR 2024) uses χ² gating; NLOS mitigation with NNs; robust outlier gates | **forced**: real UWB outliers are NLOS/multipath, not the thermal threshold we study |
| Vision-language models, imitation learning | — | no natural link; not proposed |
| Bayesian bounds for GNSS, 2024–2026 | WWB for weak-signal *frequency* (static), arXiv:2401.05606; KF loop with minimum pull-out probability (arXiv:2606.23925, carrier phase, no bound) | nothing on dynamic code-tracking bounds found |

## What we learnt
1. **Our novelty margin is narrower than we thought** — the batch ZZB for GNSS code tracking is established
   (Emmanuele 2012, Liu 2016). The defensible contribution is the *dynamic* decomposition (bound / estimator / loop)
   and the bound–loss-of-lock link.
2. The **natural** interdisciplinary bridge is state estimation and safety (Bayes filters, integrity), not vision or
   RL. Forcing robotics datasets would weaken, not strengthen, the project.
3. Novelty must always be phrased "we found no prior work that…" — this survey proved that a careful first survey can
   still miss directly relevant papers.
