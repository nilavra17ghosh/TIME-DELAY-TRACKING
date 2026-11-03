"""BPSK signal model: r[n] = s(t_n - tau) + w[n],  w[n] ~ N(0, sigma^2) i.i.d.

s(t) = sum_k a_k p(t - k), a_k in {-1, +1} (known pilot), p = raised-cosine pulse.
The waveform is evaluated analytically, so any fractional delay is exact (no
interpolation error enters the simulation).
"""
import numpy as np


def rc_pulse(t, alpha):
    """Raised-cosine pulse with unit chip period and roll-off alpha (untruncated)."""
    t = np.asarray(t, dtype=float)
    den = 1.0 - (2.0 * alpha * t) ** 2
    singular = np.abs(den) < 1e-10
    with np.errstate(divide="ignore", invalid="ignore"):
        p = np.sinc(t) * np.cos(np.pi * alpha * t) / den
    # limit at t = +-1/(2 alpha)
    return np.where(singular, np.pi / 4.0 * np.sinc(1.0 / (2.0 * alpha)), p)


def truncated_pulse(t, alpha, span):
    """Raised-cosine pulse truncated symmetrically to |t| <= span chips."""
    t = np.asarray(t, dtype=float)
    return np.where(np.abs(t) <= span, rc_pulse(t, alpha), 0.0)


class BPSKWaveform:
    """Known BPSK pilot waveform s(t) with analytic evaluation at any time."""

    def __init__(self, chips, alpha, span):
        self.chips = np.asarray(chips, dtype=float)
        self.alpha = alpha
        self.span = span
        self._offsets = np.arange(-span - 1, span + 1)   # covers |t - k| <= span

    def __call__(self, t):
        """Evaluate s(t) for an array t of any shape."""
        t = np.asarray(t, dtype=float)
        k = np.floor(t)[..., None] + self._offsets          # contributing chips
        valid = (k >= 0) & (k < self.chips.size)
        a = np.where(valid, self.chips[np.clip(k, 0, self.chips.size - 1).astype(int)], 0.0)
        return np.sum(a * truncated_pulse(t[..., None] - k, self.alpha, self.span), axis=-1)

    def derivative(self, t, h=1e-5):
        """ds/dt via a central difference (pulse is smooth, error ~h^2)."""
        return (self(np.asarray(t) + h) - self(np.asarray(t) - h)) / (2.0 * h)


def random_chips(n, rng):
    return rng.choice([-1.0, 1.0], size=n)


class BatchScenario:
    """Geometry of one MLE batch: observation window, prior window, noise level."""

    def __init__(self, cfg, rng, chips=None):
        """Random pilot of cfg.n_chips drawn from rng, or the given `chips`."""
        self.cfg = cfg
        chips = random_chips(cfg.n_chips, rng) if chips is None else np.asarray(chips, float)
        self.n_chips = chips.size
        self.wave = BPSKWaveform(chips, cfg.rolloff, cfg.span)
        # prior centre placed so the whole pulse train (incl. tails) stays inside
        # the observation window for every tau in the prior -> constant energy.
        self.tau_c = cfg.span + cfg.prior_width / 2.0 + 1.0
        n_win = self.n_chips + 2 * cfg.span + int(np.ceil(cfg.prior_width)) + 3
        self.t = np.arange(n_win * cfg.sps) / cfg.sps
        self.energy = float(np.sum(self.wave(self.t - self.tau_c) ** 2))  # sum s^2

    def sigma2(self, en0):
        """Per-sample noise variance giving total E/N0 = en0 (linear).

        Discrete-time equivalent of E/N0 for a Nyquist-sampled real signal:
        E/N0 = sum_n s[n]^2 / (2 sigma^2).
        """
        return self.energy / (2.0 * en0)

    def draw_delays(self, n, rng):
        w = self.cfg.prior_width
        return self.tau_c + rng.uniform(-w / 2.0, w / 2.0, size=n)

    def observe(self, tau, en0, rng):
        """Received samples for each delay in tau -> array (len(tau), n_samples)."""
        clean = self.wave(self.t[None, :] - np.asarray(tau)[:, None])
        noise = rng.standard_normal(clean.shape) * np.sqrt(self.sigma2(en0))
        return clean + noise
