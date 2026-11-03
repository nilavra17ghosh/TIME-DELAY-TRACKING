"""Performance bounds for delay estimation of a known signal in white Gaussian noise.

All bounds are computed from the SAME sampled waveform used in the simulation,
so theory and Monte Carlo refer to exactly one model.
"""
import numpy as np
from scipy.special import erfc


def qfunc(x):
    return 0.5 * erfc(np.asarray(x) / np.sqrt(2.0))


def rms_bandwidth_sq(scenario):
    """beta^2 = sum s'[n]^2 / sum s[n]^2  (in 1/chip^2), the 'Gabor bandwidth'."""
    t = scenario.t - scenario.tau_c
    ds = scenario.wave.derivative(t)
    return float(np.sum(ds ** 2) / np.sum(scenario.wave(t) ** 2))


def crlb(scenario, en0):
    """Classical CRLB on the delay: var >= sigma^2 / sum s'^2 = 1 / (2 (E/N0) beta^2)."""
    return 1.0 / (2.0 * np.asarray(en0) * rms_bandwidth_sq(scenario))


def autocorrelation(scenario, h):
    """Normalised autocorrelation rho(h) = <s_phi, s_{phi+h}> / ||s||^2."""
    base = scenario.tau_c - scenario.cfg.prior_width / 2.0
    s0 = scenario.wave(scenario.t - base)
    h = np.asarray(h, dtype=float)
    out = np.empty(h.size)
    for c in range(0, h.size, 256):                         # chunked: bounds memory
        sh = scenario.wave(scenario.t[None, :] - base - h[c:c + 256, None])
        out[c:c + 256] = sh @ s0
    return out / np.sum(s0 ** 2)


def zzb(scenario, en0, n_h=8001):
    """Ziv-Zakai bound (Bell et al. 1997, scalar, with valley filling) for
    tau ~ U[0, Ta]:

        MSE >= int_0^Ta  h * V{ (1 - h/Ta) * Pmin(h) } dh,
        Pmin(h) = Q( sqrt( (E/N0) * (1 - rho(h)) ) ),

    where Pmin(h) is the minimum error probability of the binary test
    'tau = phi' vs 'tau = phi + h' (equally likely, known signals, AWGN), and
    V{f}(h) = max_{eta >= 0} f(h + eta) is the valley-filling operator.
    rho(h) is treated as independent of phi (verified: varies by ~2e-5).
    With the smoothly truncated pulse, ZZB/CRLB -> 0.997 at 40 dB (converged in
    n_h). Valley filling is numerically inactive for this pilot (sidelobes
    |rho| <= 0.13), changing the ZZB by < 0.01 %.
    """
    ta = scenario.cfg.prior_width
    h = np.linspace(0.0, ta, n_h)
    one_minus_rho = np.clip(1.0 - autocorrelation(scenario, h), 0.0, None)
    en0 = np.atleast_1d(np.asarray(en0, dtype=float))
    f = (1.0 - h / ta)[None, :] * qfunc(np.sqrt(en0[:, None] * one_minus_rho[None, :]))
    f_vf = np.maximum.accumulate(f[:, ::-1], axis=1)[:, ::-1]   # valley filling
    return np.trapezoid(h[None, :] * f_vf, h, axis=1)


def prior_variance(cfg):
    return cfg.prior_width ** 2 / 12.0
