"""Signal model, parameters and prior (Sprints 0-1).

Observation model:  r[n] = s(t_n - tau) + w[n],   w[n] ~ N(0, sigma^2) i.i.d.
    s(t) = sum_k a_k p(t - k),  a_k in {-1, +1} (known BPSK pilot),
    p    = raised-cosine pulse (roll-off alpha), truncated to |t| <= span chips.

Units used everywhere:
    * time and delay in CHIPS (Tc = 1); variances / MSEs in chips^2,
    * SNR = TOTAL batch energy-to-noise ratio E/N0 of the 64-chip pilot
      (per-chip Ec/N0 = E/N0 / n_chips).
The waveform is evaluated analytically, so any fractional delay is exact.
"""
from dataclasses import dataclass, field

import numpy as np


@dataclass
class Config:
    # ---- signal model -------------------------------------------------------
    n_chips: int = 64            # known BPSK pilot chips in one MLE batch
    sps: int = 4                 # samples per chip (fs = sps / Tc)
    rolloff: float = 0.5         # raised-cosine roll-off (band-limited pulse)
    span: int = 6                # pulse truncated to |t| <= span chips

    # ---- prior on the delay ---------------------------------------------------
    prior_width: float = 4.0     # Ta: tau ~ U[tau_c - Ta/2, tau_c + Ta/2] (chips)

    # ---- batch MLE ------------------------------------------------------------
    mle_grid_step: float = 1 / 64  # template spacing before parabolic refinement

    # ---- DLL ------------------------------------------------------------------
    el_spacing: float = 0.5      # early-late spacing d (chips)
    acq_step: float = 0.5        # acquisition search grid step (chips)
    dll_block_chips: int = 4     # coherent integration per loop update (chips)
    dll_settle_factor: int = 10  # number of updates = factor * N_eff

    # ---- Monte Carlo ----------------------------------------------------------
    snr_db: np.ndarray = field(default_factory=lambda: np.arange(-10.0, 32.0, 2.0))
    n_trials: int = 2000
    seed: int = 2026

    # ---- derived quantities ----------------------------------------------------
    @property
    def dll_n_eff(self) -> float:
        """Number of loop updates the DLL effectively averages.

        Chosen so that the DLL's equivalent integration time 1/(2 B_L) equals the
        MLE batch duration -> both estimators use the same amount of signal energy.
        """
        return self.n_chips / self.dll_block_chips

    @property
    def dll_gain(self) -> float:
        """First-order loop gain K with (2 - K) / K = N_eff."""
        return 2.0 / (self.dll_n_eff + 1.0)

    @property
    def dll_bl(self) -> float:
        """One-sided loop noise bandwidth B_L in 1/chip (first-order loop)."""
        k, t_i = self.dll_gain, self.dll_block_chips
        return k / (2.0 * (2.0 - k) * t_i)

    @property
    def dll_n_updates(self) -> int:
        return int(self.dll_settle_factor * self.dll_n_eff)



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
    """Raised-cosine pulse, smoothly truncated to |t| <= span chips.

    The tail is tapered to zero over the last chip (span-1 < |t| <= span) with a
    raised-cosine window, so the pulse AND its slope are continuous. A hard cut
    leaves a slope jump of ~0.005 at |t| = span; samples landing exactly on it
    made the sampled autocorrelation's curvature disagree with beta^2, so the
    ZZB did not converge cleanly to the CRLB. The taper only touches
    5 < |t| <= 6, where the remaining integer points are zeros of the pulse anyway,
    so the Nyquist property is unchanged.
    """
    t = np.asarray(t, dtype=float)
    a = np.abs(t)
    w = np.where(a <= span - 1, 1.0, 0.5 * (1.0 + np.cos(np.pi * (a - (span - 1)))))
    return np.where(a <= span, rc_pulse(t, alpha) * w, 0.0)


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
