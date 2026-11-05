"""Fast regression tests: the core library still agrees with the theory.

Run from the repository root:   python -m pytest -q tests        (about a minute)

Each test compares the code against an independently derived fact (not against
itself). The full verifications, with plots and explanations, are in the notebooks.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from sdet import (BatchMLE, BatchScenario, Config, EarlyLateDLL, crlb, prior_variance,  # noqa: E402
                  rc_pulse, rms_bandwidth_sq, zzb)


@pytest.fixture(scope="module")
def setup():
    cfg = Config()
    sc = BatchScenario(cfg, np.random.default_rng(cfg.seed))
    return cfg, sc


def test_pulse_is_nyquist(setup):
    cfg, _ = setup
    k = np.arange(-cfg.span, cfg.span + 1)
    p = rc_pulse(k, cfg.rolloff)
    assert abs(p[cfg.span] - 1) < 1e-12
    assert np.max(np.abs(np.delete(p, cfg.span))) < 1e-12


def test_snr_calibration(setup):
    _, sc = setup
    en0 = 10.0
    noise = sc.observe(np.full(3000, sc.tau_c), en0, np.random.default_rng(1)) - sc.wave(sc.t - sc.tau_c)
    assert abs(sc.energy / (2 * noise.var()) / en0 - 1) < 0.01          # E/N0 = sum s^2 / (2 sigma^2)


def test_crlb_bandwidth_frequency_domain(setup):
    cfg, sc = setup
    fs = 64
    tt = np.arange(-cfg.span - 2, cfg.n_chips + cfg.span + 2, 1 / fs)
    spec = np.abs(np.fft.rfft(sc.wave(tt), n=8 * tt.size)) ** 2
    f = np.fft.rfftfreq(8 * tt.size, d=1 / fs)
    beta2_fft = np.sum((2 * np.pi * f) ** 2 * spec) / np.sum(spec)   # Gabor bandwidth via Parseval
    assert abs(beta2_fft / rms_bandwidth_sq(sc) - 1) < 1e-3


def test_mle_noise_free_and_efficient(setup):
    _, sc = setup
    mle, rng = BatchMLE(sc), np.random.default_rng(2)
    tau = sc.draw_delays(500, rng)
    e0 = mle.estimate(sc.wave(sc.t[None, :] - tau[:, None])) - tau
    assert np.mean(e0 ** 2) < 0.01 * crlb(sc, 1e3)                     # interpolation error negligible
    en0 = 10 ** 2.5
    tau = sc.tau_c + rng.uniform(-1.5, 1.5, 3000)                        # away from the prior edges
    e = mle.estimate(sc.observe(tau, en0, rng)) - tau
    assert abs(np.mean(e ** 2) / crlb(sc, en0) - 1) < 0.12               # asymptotic efficiency


def test_zzb_limits(setup):
    cfg, sc = setup
    lo, hi = zzb(sc, [1e-6, 1e4])
    assert abs(lo / prior_variance(cfg) - 1) < 0.01                      # -> Ta^2 / 12
    assert abs(hi / crlb(sc, 1e4) - 1) < 0.01                            # -> CRLB


def test_dll_s_curve_and_lock(setup):
    cfg, sc = setup
    dll, rng = EarlyLateDLL(cfg, sc.wave.chips), np.random.default_rng(3)
    slope = (dll.s_curve(0.01) - dll.s_curve(-0.01))[0] / 0.02
    assert abs(slope - 1) < 0.01                                         # D(eps) ~ eps near zero
    tau = dll.tau_c + rng.uniform(-cfg.prior_width / 2, cfg.prior_width / 2, 50)
    e = dll.run(tau, 0.0, rng)["theta"] - tau
    assert np.max(np.abs(e)) < 1e-4                                      # locks from anywhere in the prior
    assert 1.5 < dll.pull_in_range(rng) < 2.1                            # measured pull-in ~1.8 chips


def test_dll_linear_theory(setup):
    cfg, sc = setup
    dll, rng = EarlyLateDLL(cfg, sc.wave.chips), np.random.default_rng(4)
    en0 = 10 ** 2.5
    tau = sc.draw_delays(1500, rng) - sc.tau_c + dll.tau_c
    e = dll.run(tau, sc.sigma2(en0), rng)["theta"] - tau
    assert abs(np.mean(e ** 2) / dll.linear_theory_var(sc.sigma2(en0)) - 1) < 0.15
