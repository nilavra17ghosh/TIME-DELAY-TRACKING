# Surveys — what we checked, and what each one taught us

Four surveys were run during the project. Read them in order: each one changed what we believed after the previous one.

| # | Date | Survey | Method | Main lesson |
|---|---|---|---|---|
| 01 | 2026-10-07 | [Independent validation](01_independent_validation_2026-10-07.md) | 3 independent tracks: literature check of every formula, from-scratch numerical re-implementation, adversarial code audit | **The core theory and code were correct, but several *claims* were wrong** (unfair DLL/ZZB comparison, thresholds without error bars, wrong pull-in range, over-stated efficiency argument). Led to the fixes and to error bars everywhere. *Historical: its numbers refer to the Oct 7 code.* |
| 02 | 2026-10-13 | [Literature survey](02_literature_survey_2026-10-13.md) | annotated bibliography by theme + targeted gap searches | **Everything up to the batch threshold study is a reproduction**; the candidate gap is a bound for DLL loss of lock. |
| 03 | 2026-10-13 | [Deep-research novelty assessment](03_deep_research_novelty_2026-10-13.md) | automated multi-agent search, 25 claims × 3 verifiers | Confirmed 02; **narrowed the gap**: optimal-MTLL theory exists for *phase* tracking, so the claim must be about *code delay* and the *bound ↔ MTLL link*. |
| 04 | 2026-11-05 | [Interdisciplinary directions & second novelty check](04_interdisciplinary_and_novelty_2026-11-05.md) | targeted searches on robotics, safety, control/RL, code design, recent bounds | **Batch ZZB for GNSS code tracking already exists** (Emmanuele 2012, Liu 2016) — missed by 02/03. The natural bridge is state estimation + integrity; robotics/VLM links would be forced. |

## Cumulative conclusions
1. **What is a reproduction:** MLE, CRLB, ZZB for delay, MLE threshold prediction (MIE), ZZB for GNSS code tracking
   (batch), linear DLL theory, DLL MTLL.
2. **What is our own (observations, not claimed as new theory):** ZZB looseness measured against the MMSE estimator
   (≈ 2.6 dB bound looseness + 0.85 dB MLE sub-optimality), its pilot dependence, the DLL failure decomposition and the
   criterion-dependent DLL-vs-MLE comparison (notebooks 05–06).
3. **What appears open (medium-low confidence):** bounds including tracking dynamics, linked to loss of lock, for PN
   code delay → our proposal (`docs/evaluation/research_proposal.md`).
4. **Method lesson:** every novelty statement is phrased "we found no prior work that…"; two separate surveys each
   missed something the next one found.
