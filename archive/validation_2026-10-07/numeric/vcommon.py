"""Independent re-implementation of the signal model (does NOT import the project's tde package).

Pulse: raised cosine, roll-off b, unit chip period, truncated to |t| <= span (symmetric)
or, optionally, to the project's effective support [-span, span+1) ('asym').
"""
import numpy as np
from scipy.special import erfc


def pilot(seed=2026, n=64):
    return np.random.default_rng(seed).choice([-1.0, 1.0], size=n)


def Q(x):
    return 0.5 * erfc(np.asarray(x) / np.sqrt(2.0))


def rc(t, b=0.5):
    """Raised cosine written via cos/sin directly; singular points handled by limit."""
    t = np.asarray(t, dtype=float)
    x = 2.0 * b * t
    sing = np.isclose(np.abs(x), 1.0, rtol=0.0, atol=1e-10)
    tt = np.where(sing, 0.3, t)  # dummy to avoid 0/0
    with np.errstate(divide="ignore", invalid="ignore"):
        s = np.where(np.abs(tt) < 1e-12, 1.0, np.sin(np.pi * tt) / (np.pi * tt))
        val = s * np.cos(np.pi * b * tt) / (1.0 - (2 * b * tt) ** 2)
    lim = (np.pi / 4.0) * np.sinc(1.0 / (2.0 * b))
    return np.where(sing, lim, val)


def rc_trunc(t, b=0.5, span=6, mode="sym"):
    t = np.asarray(t, dtype=float)
    if mode == "sym":
        m = np.abs(t) <= span
    else:  # project's effective support: chip k contributes for t-k in [-span, span+1)
        m = (t >= -span) & (t < span + 1)
    return np.where(m, rc(t, b), 0.0)


def waveform(t, chips, b=0.5, span=6, mode="sym"):
    """s(t) = sum_k a_k p(t-k), direct sum over all chips (no tap trick)."""
    t = np.asarray(t, dtype=float)
    out = np.zeros(t.shape)
    for k, a in enumerate(chips):
        d = t - k
        sel = (d >= -span - 1) & (d <= span + 1)
        if np.any(sel):
            out[sel] += a * rc_trunc(d[sel], b, span, mode)
    return out


def waveform_deriv(t, chips, b=0.5, span=6, mode="sym", h=1e-6):
    # 4th-order central difference
    f = lambda x: waveform(x, chips, b, span, mode)
    return (-f(t + 2 * h) + 8 * f(t + h) - 8 * f(t - h) + f(t - 2 * h)) / (12 * h)
