# Survey 03 — Deep-research novelty assessment (2026-10-13)

**Question asked:** which parts of our work reproduce published results, which are open, and does any existing bound
describe DLL loss of lock?

**Method:** an automated deep-research workflow — 5 search angles, 20 primary sources fetched, 92 falsifiable claims
extracted, 25 claims put to 3 independent verifiers each (claim survives only if fewer than 2 of 3 refute it):
23 confirmed, 2 refuted. Raw output was saved in the session scratchpad; the verified findings are reproduced here.

## Verified findings

| # | Finding | Confidence (vote) | Sources |
|---|---|---|---|
| 1 | **Done:** closed-form prediction of the MLE threshold for known-pulse delay estimation already exists (method of interval estimation, MIE), including raised-cosine pulses. Our measured MLE threshold is a reproduction. | high (3–0) | Mallat et al. 2014 Part 2, arXiv:1405.4647; Athley, ICASSP 2002 |
| 2 | **Done (scaling law):** the asymptotic MLE threshold depends only on the normalised autocorrelation shape (roll-off), not on chip rate; the a-priori threshold also depends on the prior width. | high (3–0) | Mallat 2014 Part 2, Sec. III-C |
| 3 | **Done:** the ZZB captures the sidelobe-driven threshold the CRLB misses, and is used as a design objective for ranging signals. **Partly open:** ranking discrete PN codes by the ZZB–MLE gap. | high (3–0) | Graff & Humphreys, arXiv:2403.12905; Xiong & Liu, arXiv:2310.02963 |
| 4 | **Done:** DLL loss of lock is long characterised by mean-time-to-lose-lock (MTLL) / exit-time theory, not by MSE bounds. | high (3–0, one 2–1) | Welti & Bobrovsky 1989; Hasson & Bobrovsky 1995; Groh, Gentner & Selva 2012; Schuss et al., SIAM 2006; Olmos & Agustí |
| 5 | **Done (textbook):** linearised DLL theory does not describe the large-error regime. | medium (3–0, one source) | Cassel, NAVIGATION 70(4), 2023 |
| 6 | **Tools exist:** recursive Bayesian bounds for nonlinear filtering — sequential WWB (needs linear dynamics, allows nonlinear measurement — our case), marginal WWB, recursive Bobrovsky–Zakai-type bound. | high (3–0) | Rapoport & Oshman, IEEE T-SP 2007; Fritsche et al., ICASSP 2015; Nitzan, Routtenberg & Tabrikian, IEEE SPL 2018 |
| 7 | **Open (as far as found):** none of these bounds has been applied to PN code tracking or linked to loss of lock / MTLL. | medium (synthesis; absence of evidence) | — |
| 8 | **Partly open (incremental):** pilot-to-pilot spread of the ZZB–MLE gap; ZZB looseness measured against the MMSE estimator. | low | — |

## Refuted claims (do not cite)
* "Mallat et al. 2014 show the ZZB mislocates the MLE threshold" — refuted 0–3.
* "The established design criterion for DLL parameters is MTLL rather than variance" — refuted 1–2.

## Follow-up targeted searches (same day)
Searches for "posterior CRB GNSS code tracking", "Weiss–Weinstein bound PLL cycle slip", "sequential Ziv–Zakai bound",
"Bayesian bound code delay loss of lock" found no direct prior work, but found the **optimal-MTLL filtering theory for
phase tracking** (Fischler & Schuss, arXiv math/0703524, 2007; Ezri, Bobrovsky & Schuss, arXiv 0802.0414; Fischler &
Bobrovsky, Signal Processing 2006; Schuss, *Nonlinear Filtering and Optimal Phase Tracking*), which narrowed the gap to
*code delay* (non-periodic, code-dependent side-lobes) and to the *MSE-bound ↔ MTLL link*.

## What we learnt
1. Almost all of the core project is a reproduction — fine for a course project, but it must be presented as such.
2. The open space lies on the **tracking (dynamic)** side, not the batch side.
3. Phase-tracking results exist and are close; any claim must be specific to **code delay**.
4. Absence-of-evidence claims need targeted follow-up — which Survey 04 did, and which changed the picture again.
