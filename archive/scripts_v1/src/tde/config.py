"""Simulation configuration.

Units convention used everywhere in this project:
    * time and delay are measured in CHIPS (Tc = 1),
    * variances / MSEs are therefore in chips^2,
    * SNR is the TOTAL batch energy-to-noise ratio E/N0 (linear or dB);
      the per-chip value is Ec/N0 = (E/N0) / n_chips.
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
