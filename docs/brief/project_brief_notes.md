# Threshold Effects in Time-Delay Tracking: CRLB vs ZZB for DLL and Batch MLE for BPSK Signals

## 1. Project Overview

This project studies **time-delay estimation of a BPSK-modulated received signal**, with a particular focus on what happens when the signal-to-noise ratio (SNR) becomes low.

The central issue is that the **Cramér-Rao Lower Bound (CRLB)** is a useful predictor of estimator performance in the high-SNR/asymptotic regime, but it can become a poor predictor in the low-SNR regime. In this region, estimators can suddenly experience a large increase in estimation error. This behaviour is commonly referred to as a **threshold effect**.

The project therefore compares:

- **Two practical delay estimators**
  - Delay-Locked Loop (DLL)
  - Batch Maximum Likelihood Estimator (MLE)

- **Two theoretical performance bounds**
  - Cramér-Rao Lower Bound (CRLB)
  - Ziv-Zakai Bound (ZZB)

The final objective is to understand and demonstrate, through theory and simulation, how these quantities behave as SNR changes.

The project brief specifically asks the team to:

1. Estimate the delay of a BPSK-modulated received signal.
2. Implement a DLL for real-time tracking.
3. Implement a batch MLE for delay estimation.
4. Obtain the error variance of both methods over a range of SNR values.
5. Compare their performance with the CRLB.
6. Identify the SNR at which the estimator-performance curves begin to diverge from the CRLB.
7. Formulate, derive, and simulate the ZZB for the given time-delay estimation problem.
8. Compare the ZZB with the CRLB and experimental estimator variances.
9. Explain the observed behaviour with theoretical backing.
10. Explain any performance difference between the DLL and batch MLE.

The official project brief also suggests using a simple **Early-Late discriminator DLL** and generating BPSK signals for different SNRs. It recommends beginning with foundational material on ZZB, DLL, and the application of ZZB to time-delay estimation. [Project Brief, Project 13](project_brief_screenshot.png)

---

# 2. The Core Question of the Project

The entire project can be reduced to one question:

> **How accurately can we estimate the time delay of a noisy BPSK signal, especially in the low-SNR region where the CRLB no longer predicts practical estimator performance well?**

To answer this, we create a known BPSK signal, introduce an unknown time delay and noise, and then try to recover that delay using two different estimators.

We then ask:

- How large is the actual estimation error?
- How does the error change with SNR?
- How closely does it follow the CRLB?
- When does the actual estimator performance depart from the CRLB?
- Can the ZZB explain the low-SNR behaviour better?
- How do the DLL and batch MLE compare?

---

# 3. Problem Statement

Assume a transmitter sends a BPSK-modulated signal.

The receiver does not observe the signal perfectly. Instead, it receives a delayed and noisy version of it.

A conceptual continuous-time model is:

\[
r(t) = s(t-\tau) + n(t)
\]

where:

- \(r(t)\) = received signal
- \(s(t)\) = transmitted BPSK waveform
- \(\tau\) = unknown propagation/time delay
- \(n(t)\) = additive noise

The parameter of interest is therefore:

\[
\boxed{\tau}
\]

The receiver must produce an estimate:

\[
\boxed{\hat{\tau}}
\]

The estimation error is:

\[
e = \hat{\tau}-\tau
\]

and the main practical performance quantity in this project is the **variance of this delay-estimation error**.

---

# 4. Why Time-Delay Estimation Matters

Time-delay estimation is a fundamental problem in communication, radar, sonar, synchronization, localization, navigation, and tracking.

If a known waveform is transmitted and reaches a receiver after an unknown delay, estimating that delay can provide information about:

- synchronization offset,
- propagation time,
- relative position or range,
- tracking state,
- timing alignment.

This project intentionally focuses on the statistical estimation aspect rather than on a complete communication-system implementation.

---

# 5. What Is BPSK?

**BPSK** stands for Binary Phase Shift Keying.

It is a binary digital modulation scheme in which two phase states represent the two possible symbols.

For a simplified baseband representation, the transmitted symbols can be represented as:

\[
a_k \in \{-1,+1\}
\]

For example, a bit sequence

\[
0,\;1,\;1,\;0,\;1
\]

can be mapped to something such as

\[
+1,\;-1,\;-1,\;+1,\;-1
\]

depending on the chosen bit-to-symbol convention.

These symbols are then shaped into a waveform using a pulse-shaping scheme.

For the project, the exact pulse shape, sampling structure, and delay representation must be defined consistently and then used throughout the theoretical derivation and simulations.

---

# 6. What Does "Time Delay" Mean?

Suppose the transmitted waveform is:

\[
s(t)
\]

and the received waveform is:

\[
s(t-\tau)
\]

Then \(\tau\) is the time by which the received signal is shifted relative to the transmitted signal.

For example:

```text
Transmitted:

      /\    /\      /\ 
_____/  \__/  \____/  \________
       ↑
      signal


Received:

            /\    /\      /\ 
___________/  \__/  \____/  \__
            ↑
         delayed
```

The receiver knows the general waveform or signal structure, but it does not know the exact value of \(\tau\).

The purpose of the estimator is to recover that unknown shift.

---

# 7. What Makes the Problem Difficult?

The main difficulty is **noise**.

A noisy observation can be represented conceptually as:

\[
r(t)=s(t-\tau)+n(t)
\]

At high SNR:

\[
\text{signal power} \gg \text{noise power}
\]

so the delayed waveform is relatively easy to identify.

At low SNR:

\[
\text{signal power} \approx \text{noise power}
\]

or even becomes smaller relative to the noise.

The signal can then become difficult to distinguish from competing delay hypotheses.

This is where the project becomes interesting.

---

# 8. What Is SNR?

SNR means **Signal-to-Noise Ratio**.

Conceptually:

\[
\mathrm{SNR}
=
\frac{\text{signal power}}{\text{noise power}}
\]

It is commonly represented in decibels:

\[
\mathrm{SNR}_{dB}
=
10\log_{10}
\left(
\frac{P_s}{P_n}
\right)
\]

A higher SNR means cleaner observations.

A lower SNR means noisier observations.

The project specifically requires experiments over a **range of SNR values** so that both the high-SNR and low-SNR regimes can be studied.

---

# 9. The Two Estimation Approaches

The project compares two ways of estimating the unknown delay.

## 9.1 Delay-Locked Loop (DLL)

A **Delay-Locked Loop** is a feedback-based tracking mechanism.

Its purpose is to continuously adjust an internal delay estimate until it aligns with the incoming signal.

Conceptually:

```text
Received signal
      |
      v
   Correlator /
   discriminator
      |
      v
 timing error estimate
      |
      v
   Loop filter
      |
      v
   Delay control
      |
      +------------------+
                         |
                         v
                  Updated estimate
```

The project specifically suggests a simple **Early-Late discriminator DLL**.

### Early-Late idea

The receiver compares samples/correlations taken slightly:

- earlier than the current timing estimate,
- at the current estimate,
- later than the current timing estimate.

The difference between the early and late branches gives a timing error signal.

If the receiver is correctly aligned, the early and late measurements become balanced.

A conceptual discriminator is:

\[
e_D \propto E-L
\]

where:

- \(E\) = early correlation/energy measurement,
- \(L\) = late correlation/energy measurement.

The exact discriminator form depends on the implementation chosen.

### Why DLL is useful

DLL represents a **real-time tracking approach**.

It updates the delay estimate incrementally rather than processing an entire batch and then producing one estimate.

The project asks us to study the variance of its tracking error as SNR changes.

---

## 9.2 Batch Maximum Likelihood Estimation (MLE)

The second method is a **batch MLE**.

Instead of tracking continuously, the receiver collects a batch of observations and determines which delay value is most likely to have generated the received signal.

Conceptually:

```text
Received batch
      |
      v
Evaluate likelihood
for candidate delays
      |
      v
Choose delay giving
maximum likelihood
      |
      v
Estimated delay
```

The estimate is:

\[
\hat{\tau}_{MLE}
=
\arg\max_{\tau}
L(\tau;r)
\]

where \(L(\tau;r)\) is the likelihood of the received observations under delay hypothesis \(\tau\).

Equivalently, one can often work with the log-likelihood:

\[
\hat{\tau}_{MLE}
=
\arg\max_{\tau}
\log L(\tau;r)
\]

The precise likelihood expression depends on the assumed noise model and signal model used for the project.

### Why MLE is useful

MLE gives a statistical estimation benchmark based on the full batch of observations.

It can therefore be compared against the tracking behaviour of the DLL.

---

# 10. What Is Estimation Error?

Suppose the true delay is:

\[
\tau=0.37
\]

and one run gives:

\[
\hat{\tau}=0.39
\]

Then:

\[
e=0.39-0.37=0.02
\]

One simulation is not enough to understand estimator quality.

We therefore perform many independent trials.

For trial \(i\):

\[
e_i=\hat{\tau}_i-\tau
\]

The error variance can then be estimated from the collection:

\[
\{e_1,e_2,\ldots,e_N\}
\]

This gives an empirical measure of estimator performance.

---

# 11. Why Run Many Monte Carlo Trials?

A single noisy signal realization may make an estimator look unusually good or unusually bad.

Monte Carlo simulation reduces this randomness by repeating the experiment many times.

For each SNR:

```text
Generate signal
      ↓
Apply random delay/noise
      ↓
Estimate delay
      ↓
Store error
      ↓
Repeat many times
      ↓
Compute empirical error variance
```

Thus for every SNR we obtain a statistical estimate of the estimator performance.

---

# 12. CRLB: The First Theoretical Benchmark

The **Cramér-Rao Lower Bound (CRLB)** is a classical theoretical lower bound associated with parameter estimation.

For an unbiased estimator under standard regularity conditions, the variance satisfies a relation of the general form:

\[
\operatorname{Var}(\hat{\tau})
\geq
\frac{1}{I(\tau)}
\]

where:

\[
I(\tau)
\]

is the **Fisher Information** for the delay parameter.

The exact CRLB for this project must be derived from the chosen signal, sampling model, noise model, parameterization, and observation interval.

The important conceptual interpretation is:

> CRLB tells us the local/asymptotic precision suggested by the information in the observations.

At sufficiently high SNR, practical estimators can often operate close to this bound.

---

# 13. Why CRLB Is Not the Whole Story

A major point of this project is:

> **CRLB is not necessarily an accurate description of practical estimator performance in the low-SNR region.**

At low SNR, the likelihood or correlation surface can become difficult to distinguish across nearby and distant delay hypotheses.

Instead of making a small estimation error around the true delay, an estimator may make a much larger error.

This produces a sudden increase in error variance.

That sudden degradation is the **threshold effect**.

---

# 14. Understanding the Threshold Effect Intuitively

Imagine the true delay is:

\[
\tau=0.50
\]

At high SNR, the estimator's candidate objective function may look like:

```text
                 strong peak
                    /\
                   /  \
__________________/    \_________________
                 ↑
              true delay
```

The correct delay produces a clearly dominant peak.

Now reduce SNR:

```text
        noisy / ambiguous objective
           /\       /\
      /\  /  \  /\ /  \    /\
_____/  \/    \/  \/    \__/  \____
          ↑
      true delay
```

Now several delay values can appear plausible.

An estimator may sometimes lock onto an incorrect peak.

Most importantly, this is not merely a small perturbation around the true delay. The estimator can make a **large error**.

That causes the empirical variance to increase much more rapidly than the high-SNR CRLB prediction.

---

# 15. What Is the Threshold SNR?

There is not necessarily a single universal numerical threshold independent of the simulation setup.

For this project, the threshold region is to be **identified experimentally from the curves**.

Operationally, it is the region where the actual estimator variance begins to depart substantially from the CRLB/high-SNR prediction.

Therefore the team should define a clear criterion for identifying the divergence point, for example based on:

- visual separation,
- relative error from CRLB,
- or another explicitly stated quantitative criterion.

The chosen criterion should be reported clearly in the final work.

---

# 16. ZZB: The Second Theoretical Benchmark

The **Ziv-Zakai Bound (ZZB)** is another theoretical bound for estimation error.

Its key relevance here is that it can account for **global ambiguity and threshold behaviour**, rather than only the local information captured by the CRLB.

This makes it particularly useful in the low-SNR regime.

The project brief explicitly states that ZZB is expected to provide a tighter bound in low-SNR conditions and to asymptotically approach the CRLB. [Project Brief, Project 13](project_brief_screenshot.png)

Conceptually:

```text
High SNR:

Actual estimator
       |
       | close to
       v
      ZZB ≈ CRLB


Low SNR:

Actual estimator
       |
       | often better described by
       v
      ZZB

while

CRLB
  |
  | may become overly optimistic
```

---

# 17. Why ZZB Can Behave Better at Low SNR

A useful way to understand the distinction is:

### CRLB

CRLB is strongly tied to **local sensitivity** of the likelihood around the true parameter.

It asks, in effect:

> If the estimator is already operating in the right neighbourhood, how much precision can the data provide locally?

### ZZB

ZZB is designed to account more directly for the probability of confusing different parameter values.

It therefore has a stronger connection to:

> How difficult is it to distinguish one delay hypothesis from another?

This distinction becomes important when noise is large enough that the estimator can jump between widely separated candidate delays.

---

# 18. The Main Comparison

The project therefore combines four quantities:

| Category | Quantity | Purpose |
|---|---|---|
| Practical estimator | DLL variance | Measures real-time tracking performance |
| Practical estimator | Batch MLE variance | Measures batch statistical estimation performance |
| Theoretical bound | CRLB | Classical benchmark, especially useful in high-SNR/asymptotic region |
| Theoretical bound | ZZB | Benchmark intended to better capture low-SNR/threshold behaviour |

The central result will be a graph of these quantities versus SNR.

---

# 19. Expected Plot Structure

A major final figure should contain:

- SNR on the x-axis.
- Delay-estimation error variance on the y-axis.
- DLL variance.
- Batch MLE variance.
- CRLB.
- ZZB.

For example, the conceptual structure is:

```text
Error variance
     ^
     |
     |                         DLL / MLE
     |                       /
     |                    __/
     |                 __/
     |              __/
     |-------------/---------------- CRLB
     |          __/
     |       __/
     |    __/                    ZZB
     |
     +--------------------------------------> SNR
          low                  high
```

The actual curves are expected to be determined by simulation and should not be assumed in advance.

---

# 20. What We Need to Demonstrate

The final project should establish the following experimental/theoretical story.

### High-SNR region

We expect the estimators to behave more regularly.

The practical estimator variance should become increasingly consistent with the asymptotic theoretical prediction.

The ZZB should approach the CRLB in this region.

### Low-SNR region

The practical estimation error can increase sharply because the receiver has difficulty distinguishing the correct delay from incorrect delay hypotheses.

This is the region where the CRLB may no longer accurately predict the practical estimator variance.

The ZZB should be investigated as the bound that better represents the low-SNR regime.

These statements are the conceptual behaviour the project asks the team to investigate and validate; the exact numerical results must come from the simulations.

---

# 21. End-to-End Project Pipeline

The whole project can be viewed as the following pipeline:

```text
                    ┌─────────────────────┐
                    │ Generate BPSK data  │
                    └──────────┬──────────┘
                               │
                               v
                    ┌─────────────────────┐
                    │ Apply true delay τ  │
                    └──────────┬──────────┘
                               │
                               v
                    ┌─────────────────────┐
                    │ Add noise at chosen │
                    │        SNR          │
                    └──────────┬──────────┘
                               │
                               v
                    ┌─────────────────────┐
                    │  Received signal    │
                    └───────┬─────┬───────┘
                            │     │
                   ┌────────┘     └─────────┐
                   v                         v
          ┌─────────────────┐      ┌─────────────────┐
          │      DLL        │      │    Batch MLE    │
          │ Early-Late      │      │                 │
          │ discriminator   │      │ delay estimate  │
          └────────┬────────┘      └────────┬────────┘
                   │                        │
                   v                        v
          estimated delay           estimated delay
                   │                        │
                   └──────────┬─────────────┘
                              v
                     empirical errors
                              │
                              v
                    error variance vs SNR
                              │
              ┌───────────────┼────────────────┐
              v               v                v
            CRLB             ZZB        estimator variances
              └───────────────┼────────────────┘
                              v
                   compare all curves
                              │
                              v
                  identify threshold region
                              │
                              v
                    explain the behaviour
```

---

# 22. Mathematical Components of the Project

There are several mathematical layers.

## 22.1 Signal model

Define the BPSK transmitted waveform:

\[
s(t)
\]

and received model:

\[
r(t)=s(t-\tau)+n(t)
\]

or an equivalent sampled-time version.

---

## 22.2 Parameter to estimate

The unknown parameter is:

\[
\tau
\]

Everything in the project revolves around estimating this parameter accurately.

---

## 22.3 Estimator

We produce:

\[
\hat{\tau}
\]

using either:

- DLL
- batch MLE

---

## 22.4 Estimation error

\[
e=\hat{\tau}-\tau
\]

---

## 22.5 Performance measure

The main empirical quantity is:

\[
\operatorname{Var}(e)
\]

or the corresponding empirical variance estimate over Monte Carlo trials.

---

## 22.6 CRLB

The delay-parameter Fisher information is derived from the selected statistical model, and the CRLB follows from the corresponding information quantity.

The final derivation must be consistent with the exact signal/noise model used in simulation.

---

## 22.7 ZZB

The ZZB must be formulated specifically for the chosen delay-estimation problem.

The project requires both:

- a theoretical derivation,
- and a numerical implementation.

The derivation should make the connection between pairwise delay discrimination and expected estimation error clear.

---

# 23. Why the Project Has Both Theory and Simulation

This is not intended to be a purely coding project.

There are two complementary goals.

## Theoretical goal

Derive and understand:

- the CRLB,
- the ZZB,
- the conditions under which they are meaningful,
- why their predictions differ in different SNR regimes.

## Experimental goal

Implement:

- BPSK signal generation,
- noisy observation generation,
- DLL,
- batch MLE,
- Monte Carlo evaluation,
- CRLB calculation,
- ZZB calculation,
- plotting and comparison.

The value of the project comes from connecting both sides.

---

# 24. What "Validate the Theoretical ZZB" Means

The project asks the team to validate the theoretically derived ZZB using simulations.

That means the work should not stop at:

> "Here is a ZZB equation."

Instead, the process should be:

```text
derive ZZB
   ↓
implement ZZB numerically
   ↓
run estimator simulations
   ↓
calculate empirical estimator variance
   ↓
plot empirical variance and ZZB together
   ↓
check whether the observed low-SNR behaviour
is consistent with the theoretical prediction
```

The discussion should explain where agreement is good and where differences arise.

---

# 25. DLL vs Batch MLE: What We Are Studying

The project is not simply asking which method has the lower variance everywhere.

The comparison should explain **why** their performance differs.

The two methods have fundamentally different operating principles:

### DLL

- feedback-based,
- continuously tracks timing,
- depends on loop/discriminator behaviour,
- designed for real-time operation.

### Batch MLE

- processes a collection of observations together,
- directly optimizes a likelihood criterion,
- can exploit the complete batch,
- produces an estimate after batch processing.

Therefore the project should study both:

\[
\operatorname{Var}(e_{DLL})
\]

and

\[
\operatorname{Var}(e_{MLE})
\]

as functions of SNR.

---

# 26. Important Experimental Variables

The simulation should control and document variables such as:

### Signal-related

- number of BPSK symbols,
- pulse shape,
- samples per symbol,
- observation duration,
- true delay.

### Noise-related

- noise distribution,
- noise variance,
- SNR range.

### DLL-related

- early-late spacing,
- loop filter parameters,
- update mechanism,
- initialization,
- locking/transient handling.

### MLE-related

- candidate delay search grid,
- interpolation strategy if sub-sample delay is used,
- likelihood expression,
- batch length.

### Monte Carlo

- number of trials per SNR,
- random seed policy,
- variance estimation procedure.

These choices are not all specified by the project brief, so they must be explicitly selected, justified, and documented by the team.

---

# 27. A Critical Point: Keep the Model Consistent

A common source of incorrect conclusions in this project would be using one model for theory and another for simulation.

For example:

```text
Simulation:
BPSK + one noise model + discrete delay

Theory:
different waveform + different noise assumptions + continuous delay
```

Then the curves may not be comparable.

Therefore:

> **The signal model, delay model, noise model, and observation assumptions used for the derivations must be clearly matched to the simulation setup.**

This consistency is essential.

---

# 28. Recommended Experimental Progression

The project brief recommends starting with the basics and building progressively.

A sensible progression is:

### Stage 1 — Simplest signal model

Use a controlled BPSK signal and a simple known-delay simulation.

Verify that the received signal is being generated correctly.

### Stage 2 — Batch MLE

Implement delay estimation using likelihood/correlation-based processing.

Test at high SNR first.

### Stage 3 — DLL

Implement the Early-Late DLL.

First test whether the loop can lock correctly in clean/noisy but sufficiently high-SNR conditions.

### Stage 4 — Monte Carlo evaluation

Repeat both estimators over many trials and measure variance across SNR values.

### Stage 5 — CRLB

Derive and implement the CRLB for the exact model.

### Stage 6 — Identify threshold region

Compare empirical estimator variances against CRLB and determine the SNR region where the curves diverge.

### Stage 7 — ZZB derivation

Derive the ZZB specifically for the chosen BPSK delay-estimation problem.

### Stage 8 — ZZB simulation

Implement the bound and compare it with empirical performance.

### Stage 9 — Final analysis

Explain:

- high-SNR behaviour,
- low-SNR behaviour,
- threshold effect,
- CRLB limitations,
- ZZB behaviour,
- DLL vs MLE differences.

---

# 29. What the Final Results Should Look Like

The final report should contain enough evidence to answer the project's questions rather than just showing one graph.

Useful results include:

## Plot 1 — Example BPSK signal

Show:

- clean transmitted waveform,
- delayed version,
- noisy received waveform.

This establishes the signal model visually.

## Plot 2 — DLL tracking behaviour

Show how the DLL's estimated delay evolves over time.

For example:

\[
\hat{\tau}(k)
\]

versus time/sample index.

## Plot 3 — Estimation error distribution

At selected SNRs, show the distribution of:

\[
e=\hat{\tau}-\tau
\]

This can help demonstrate how errors change from small local fluctuations to occasional large errors.

## Plot 4 — Main performance plot

Plot:

- DLL variance,
- MLE variance,
- CRLB,
- ZZB

against SNR.

This is the key project figure.

## Plot 5 — Threshold analysis

Zoom into the SNR region where estimator performance diverges from the CRLB.

---

# 30. What We Should Look for in the Curves

The analysis should focus on the following observations.

### Observation A: High SNR

Does the estimator variance approach the CRLB?

### Observation B: Transition region

At what SNR does the behaviour stop following the high-SNR trend?

### Observation C: Threshold region

Does the empirical error suddenly increase?

### Observation D: ZZB

Does ZZB remain closer to the empirical behaviour in the low-SNR region?

### Observation E: High-SNR convergence

Does ZZB approach CRLB as SNR increases?

### Observation F: DLL versus MLE

Which method performs better under the chosen conditions, and why?

The conclusions must be based on actual simulation results.

---

# 31. What the Threshold Effect Means Physically/Statistically

The threshold effect can be understood as a change from:

### Local estimation regime

The signal is sufficiently clear that the estimator is usually operating near the true delay.

Errors are small.

CRLB is informative.

### Global ambiguity regime

The signal becomes so noisy that several delay hypotheses become difficult to distinguish.

The estimator sometimes makes a much larger jump away from the correct delay.

The error variance increases sharply.

This is exactly why a bound that only captures local behaviour may become insufficient for describing practical performance.

---

# 32. Important Terminology

### Time-delay estimation

Estimating the unknown temporal shift between a transmitted/reference waveform and its received observation.

### BPSK

Binary Phase Shift Keying, a modulation scheme with two phase states.

### SNR

Signal-to-Noise Ratio, describing signal strength relative to noise strength.

### DLL

Delay-Locked Loop, a feedback tracking system used to estimate/track timing or delay.

### Early-Late discriminator

A timing-error discriminator that compares early and late correlation/measurement branches.

### MLE

Maximum Likelihood Estimation, a statistical method that selects the parameter value maximizing the likelihood of the observations.

### CRLB

Cramér-Rao Lower Bound, a classical theoretical lower bound associated with estimator variance under standard assumptions.

### ZZB

Ziv-Zakai Bound, a bound on estimation error designed to remain informative in regimes where global ambiguity and threshold effects become important.

### Threshold effect

A sharp degradation in practical estimator performance as SNR decreases past a critical region.

### Monte Carlo simulation

Repeated randomized experiments used to estimate statistical performance numerically.

---

# 33. What Is Actually Being Asked by the Mentor?

The mentor remarks for Project 13 give a clear intended route:

1. Start by reading about the ZZB foundational material.
2. Study DLLs, with Kaplan's textbook suggested as a source.
3. Study applications of ZZB to time-delay estimation.
4. Use a simple Early-Late discriminator DLL.
5. Derive the ZZB for the specific case.
6. Implement the DLL and batch MLE.
7. Generate BPSK signals for different SNRs.
8. Plot:
   - CRLB,
   - ZZB,
   - DLL tracking-error variance,
   - batch-MLE error variance.
9. Validate the derived ZZB using the simulation.
10. Explain any observed performance difference between DLL and batch MLE.

[Project Brief, Project 13](project_brief_screenshot.png)

---

# 34. What the Project Is NOT

This project should not be interpreted as merely:

- generating BPSK signals,
- implementing one timing recovery algorithm,
- calculating one CRLB number,
- or producing a single performance graph.

The project is fundamentally about the **relationship between estimator performance and theoretical bounds across SNR regimes**.

The threshold phenomenon is the central concept.

---

# 35. Final Research Story

The complete research story should look like this:

```text
We have a delayed BPSK signal in noise.
             |
             v
We want to estimate its unknown delay.
             |
             v
We implement two estimators:
DLL and batch MLE.
             |
             v
We evaluate their error variance across SNR.
             |
             v
At high SNR, CRLB provides a useful benchmark.
             |
             v
As SNR decreases, a threshold region appears.
             |
             v
Practical error can increase sharply.
             |
             v
CRLB may no longer describe this behaviour well.
             |
             v
We derive the ZZB for the specific problem.
             |
             v
We compare ZZB against simulations.
             |
             v
We study whether ZZB better captures
the low-SNR threshold region.
             |
             v
We compare the behaviour of DLL and batch MLE.
```

---

# 36. Final Deliverables Checklist

The final project should be able to provide all of the following:

### Theory

- [ ] Clearly defined BPSK signal and observation model.
- [ ] Clearly defined unknown delay parameter.
- [ ] Derivation of the relevant CRLB.
- [ ] Explanation of why the CRLB is useful in the asymptotic/high-SNR regime.
- [ ] Explanation of the threshold effect.
- [ ] Formulation and derivation of the ZZB for the chosen problem.
- [ ] Explanation of why ZZB is relevant to low-SNR estimation.

### Implementation

- [ ] BPSK signal generator.
- [ ] Delay generation/application.
- [ ] Noise generation for controlled SNR values.
- [ ] Early-Late DLL.
- [ ] Batch MLE.
- [ ] Monte Carlo evaluation.
- [ ] CRLB numerical evaluation.
- [ ] ZZB numerical evaluation.

### Experiments

- [ ] Multiple SNR values.
- [ ] Sufficient Monte Carlo trials.
- [ ] Empirical DLL error variance.
- [ ] Empirical MLE error variance.
- [ ] CRLB curve.
- [ ] ZZB curve.
- [ ] Threshold-region identification.
- [ ] DLL vs MLE comparison.

### Report

- [ ] Problem formulation.
- [ ] Mathematical derivations.
- [ ] Algorithm descriptions.
- [ ] Simulation setup.
- [ ] Results.
- [ ] Plots.
- [ ] Threshold analysis.
- [ ] CRLB vs ZZB discussion.
- [ ] DLL vs MLE discussion.
- [ ] Theoretical explanation for observed results.
- [ ] Final conclusions and limitations.

---

# 37. Suggested Repository Structure

A clean implementation can eventually be organized approximately as:

```text
project13/
│
├── README.md
│
├── docs/
│   ├── theory.md
│   ├── crlb_derivation.md
│   ├── zzb_derivation.md
│   └── experimental_setup.md
│
├── src/
│   ├── bpsk.py
│   ├── channel.py
│   ├── dll.py
│   ├── mle.py
│   ├── crlb.py
│   ├── zzb.py
│   └── metrics.py
│
├── simulations/
│   ├── run_snr_sweep.py
│   ├── monte_carlo.py
│   └── parameters.py
│
├── results/
│   ├── data/
│   ├── figures/
│   └── tables/
│
└── report/
    └── final_report.pdf
```

The exact programming language and repository structure are implementation choices; the project brief itself does not prescribe this structure.

---

# 38. The One-Sentence Understanding

> **This project investigates the low-SNR threshold effect in BPSK time-delay estimation by comparing two practical estimators—an Early-Late DLL and a batch MLE—with two theoretical performance bounds—the CRLB and the ZZB—over a range of SNR values, and then uses theory and simulation to explain where and why their behaviours diverge.**

---

# 39. Recommended Mental Model Before Coding

Before writing any code, make sure the following chain is completely clear:

\[
\boxed{
\text{BPSK signal}
\rightarrow
\text{unknown delay}
\rightarrow
\text{noise}
\rightarrow
\text{received signal}
\rightarrow
\hat{\tau}
\rightarrow
\text{error}
\rightarrow
\text{variance}
\rightarrow
\text{compare with CRLB/ZZB}
}
\]

And remember the central phenomenon:

\[
\boxed{
\text{High SNR: CRLB is useful}
\qquad
\text{Low SNR: threshold effects appear}
}
\]

The role of the ZZB is to investigate whether the theoretical bound remains informative in that difficult low-SNR regime.

---

## Source

This README is based primarily on the official **SDET Projects** document supplied for the course, specifically **Project 13: "Threshold effects in Time-Delay tracking: CRLB vs ZZB for DLL and batch MLE for BPSK signals"**, including its project description and mentor remarks on page 3.

[Official project brief](project_brief_screenshot.png)
