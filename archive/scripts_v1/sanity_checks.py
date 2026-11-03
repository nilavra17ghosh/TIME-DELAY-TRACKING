"""Verify every building block against theory before trusting any results.

Run:  python sanity_checks.py
"""
import sys

import numpy as np

sys.path.insert(0, "src")
from tde import (BatchMLE, BatchScenario, Config, EarlyLateDLL, crlb, crossing_snr,  # noqa: E402
                 prior_variance, rc_pulse, rms_bandwidth_sq, zzb)

failures = []


def check(name, ok, detail):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}")
    if not ok:
        failures.append(name)


cfg = Config()
rng = np.random.default_rng(cfg.seed)
sc = BatchScenario(cfg, rng)

# 1. pulse is a Nyquist pulse: p(0) = 1, p(k) = 0 at other integers
k = np.arange(-cfg.span, cfg.span + 1)
p = rc_pulse(k, cfg.rolloff)
check("RC pulse Nyquist property", abs(p[cfg.span] - 1) < 1e-12 and np.max(np.abs(np.delete(p, cfg.span))) < 1e-12,
      f"p(0)={p[cfg.span]:.6f}, max|p(k!=0)|={np.max(np.abs(np.delete(p, cfg.span))):.1e}")

# 2. energy per chip: int p^2 dt = 1 - alpha/4 for raised cosine
expected = cfg.sps * cfg.n_chips * (1 - cfg.rolloff / 4)
check("Signal energy", abs(sc.energy / expected - 1) < 0.05,
      f"sum s^2 = {sc.energy:.1f}, expected ~{expected:.1f} (random-chip cross terms make it approximate)")

# 3. CRLB by an independent route: continuous-time Gabor bandwidth from the FFT of
#    the waveform sampled at 64 samples/chip, vs the discrete sum used in the code
fs = 64
tt = np.arange(-cfg.span - 2, cfg.n_chips + cfg.span + 2, 1 / fs)
spec = np.abs(np.fft.rfft(sc.wave(tt), n=8 * tt.size)) ** 2
f = np.fft.rfftfreq(8 * tt.size, d=1 / fs)
beta2_fft = np.sum((2 * np.pi * f) ** 2 * spec) / np.sum(spec)
check("CRLB bandwidth, frequency-domain vs discrete", abs(beta2_fft / rms_bandwidth_sq(sc) - 1) < 1e-3,
      f"beta^2 FFT = {beta2_fft:.4f}, discrete = {rms_bandwidth_sq(sc):.4f} 1/chip^2")

# 4. Gabor bandwidth: averaged over random pilots it must equal the pulse's own
tf = np.linspace(-cfg.span, cfg.span, 200001)
pp = rc_pulse(tf, cfg.rolloff)
beta2_pulse = np.trapezoid(np.gradient(pp, tf) ** 2, tf) / np.trapezoid(pp ** 2, tf)
beta2_avg = np.mean([rms_bandwidth_sq(BatchScenario(cfg, np.random.default_rng(s))) for s in range(60)])
check("RMS bandwidth (60-pilot average)", abs(beta2_avg / beta2_pulse - 1) < 0.03,
      f"<beta^2> = {beta2_avg:.3f}, pulse = {beta2_pulse:.3f}; this pilot = {rms_bandwidth_sq(sc):.3f} 1/chip^2")

# 5. MLE is (essentially) exact without noise -> refinement bias floor
mle = BatchMLE(sc)
tau = sc.draw_delays(2000, rng)
r_clean = sc.wave(sc.t[None, :] - tau[:, None])
err = mle.estimate(r_clean) - tau
floor = np.mean(err ** 2)
check("MLE noise-free floor", floor < 0.01 * crlb(sc, 10 ** 3.0),
      f"MSE = {floor:.2e} chips^2 = {floor / crlb(sc, 1e3):.2%} of the CRLB at 30 dB")

# 6. DLL (acquisition + tracking) locks without noise from anywhere in the prior
dll = EarlyLateDLL(cfg, sc.wave.chips)   # same pilot as the MLE, repeated
tau = dll.tau_c + rng.uniform(-cfg.prior_width / 2, cfg.prior_width / 2, 201)
e = dll.run(tau, 0.0, rng)["theta"] - tau
check("DLL noise-free floor", np.mean(e ** 2) < 0.01 * crlb(sc, 1e3),
      f"MSE = {np.mean(e ** 2):.2e} = {np.mean(e ** 2) / crlb(sc, 1e3):.2%} of the CRLB at 30 dB, "
      f"max |err| = {np.max(np.abs(e)):.1e} chips")

# 7. DLL pull-in range must exceed the acquisition grid's worst-case error (0.25 chip)
r_pull = dll.pull_in_range(rng)
check("DLL pull-in range", r_pull > cfg.acq_step / 2 + 0.5,
      f"noise-free loop converges from initial errors up to {r_pull:.2f} chips")

# 8. ZZB limits
lo, hi = zzb(sc, [1e-4, 10 ** 4.0])
check("ZZB -> prior variance at very low SNR", abs(lo / prior_variance(cfg) - 1) < 0.01,
      f"ZZB = {lo:.5f}, Ta^2/12 = {prior_variance(cfg):.5f}")
check("ZZB -> CRLB at very high SNR", abs(hi / crlb(sc, 1e4) - 1) < 0.02,
      f"ZZB/CRLB at 40 dB = {hi / crlb(sc, 1e4):.4f}")

# 9. MLE and MMSE attain the CRLB at high SNR (delays away from the prior edges,
#    where the constrained estimators can beat the unconstrained CRLB)
en0 = 10 ** 2.5
tau = sc.tau_c + rng.uniform(-1.5, 1.5, 4000)
r = sc.observe(tau, en0, rng)
e_ml, e_mm = (x - tau for x in mle.estimate_both(r, sc.sigma2(en0)))
check("MLE efficient at 25 dB", abs(np.mean(e_ml ** 2) / crlb(sc, en0) - 1) < 0.1,
      f"MSE/CRLB = {np.mean(e_ml ** 2) / crlb(sc, en0):.3f}")
check("MMSE efficient at 25 dB", abs(np.mean(e_mm ** 2) / crlb(sc, en0) - 1) < 0.1,
      f"MSE/CRLB = {np.mean(e_mm ** 2) / crlb(sc, en0):.3f}")

# 10. MMSE -> prior mean (MSE -> Ta^2/12) with no signal information
tau = sc.draw_delays(4000, rng)
e_mm = mle.estimate_both(sc.observe(tau, 1e-6, rng), sc.sigma2(1e-6))[1] - tau
check("MMSE -> prior variance at zero SNR", abs(np.mean(e_mm ** 2) / prior_variance(cfg) - 1) < 0.1,
      f"MSE = {np.mean(e_mm ** 2):.3f}, Ta^2/12 = {prior_variance(cfg):.3f}")

# 11. DLL linearised theory matches simulation at high SNR
tau = dll.tau_c + rng.uniform(-cfg.prior_width / 2, cfg.prior_width / 2, 2000)
s2 = sc.sigma2(en0)
err = dll.run(tau, s2, rng)["theta"] - tau
ratio = np.mean(err ** 2) / dll.linear_theory_var(s2)
check("DLL linear theory at 25 dB", abs(ratio - 1) < 0.15, f"MSE/theory = {ratio:.3f}")

# 12. threshold finder on synthetic data with a known crossing at 10 dB
x = np.arange(0.0, 21.0, 2.0)
y = 1 + 2 ** (-(x - 10) / 2)            # = 2 exactly at 10 dB, decreasing in x
y[0] = 1.5                              # non-monotone low end must not matter
check("Threshold finder", abs(crossing_snr(x, y, 2.0) - 10.0) < 1e-9,
      f"found {crossing_snr(x, y, 2.0):.4f} dB for a crossing at 10 dB")

print("\nAll checks passed." if not failures else f"\nFAILED: {failures}")
sys.exit(1 if failures else 0)
