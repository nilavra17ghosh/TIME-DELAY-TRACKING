"""Delay estimators.

    BatchMLE.estimate        batch maximum-likelihood estimator    (Sprint 1)
    EarlyLateDLL             early-late delay-locked loop          (Sprint 2)
    BatchMLE.estimate_both   + posterior-mean (MMSE) estimator     (beyond Sprint 4)
    EarlyLateDLL.linear_theory_var  DLL variance theory            (beyond Sprint 4, = Sprint 5)
"""
import numpy as np

from .model import BPSKWaveform, truncated_pulse


# ----------------------------------------------------------------------------
# Batch MLE and MMSE
# ----------------------------------------------------------------------------
class BatchMLE:
    """Exact ML (and posterior-mean) delay estimators for a known signal in AWGN.

    log p(r; theta) = const - ||r - s_theta||^2 / (2 sigma^2)
                    = const + ( <r, s_theta> - ||s_theta||^2 / 2 ) / sigma^2
    so  theta_ML = argmax_theta  J(theta) = <r, s_theta> - ||s_theta||^2 / 2,
    searched over the prior window. sigma^2 does not enter the argmax.

    J is evaluated on a fine template grid, then refined by fitting a parabola
    through the best grid point and its two neighbours.

    With the uniform prior, the posterior is p(theta | r) ~ exp(J(theta)/sigma^2)
    on the prior window, and the MMSE estimator is its mean. Unlike the MLE it
    needs sigma^2 (the SNR). It is the estimator the ZZB is closest to, since the
    ZZB bounds the MSE of every estimator, i.e. the MMSE.
    """

    def __init__(self, scenario):
        cfg = scenario.cfg
        lo = scenario.tau_c - cfg.prior_width / 2.0
        hi = scenario.tau_c + cfg.prior_width / 2.0
        n_grid = int(round((hi - lo) / cfg.mle_grid_step)) + 1
        self.grid = np.linspace(lo, hi, n_grid)
        self.step = self.grid[1] - self.grid[0]
        self.templates = scenario.wave(scenario.t[None, :] - self.grid[:, None])
        self.half_energy = 0.5 * np.sum(self.templates ** 2, axis=1)
        self.trap = np.full(n_grid, 1.0)                   # trapezoid weights
        self.trap[[0, -1]] = 0.5

    def objective(self, r):
        return r @ self.templates.T - self.half_energy        # (trials, grid)

    def _ml_from_objective(self, j):
        g = np.argmax(j, axis=1)
        rows = np.arange(j.shape[0])
        interior = (g > 0) & (g < j.shape[1] - 1)
        gm, gp = np.clip(g - 1, 0, None), np.clip(g + 1, None, j.shape[1] - 1)
        jm, j0, jp = j[rows, gm], j[rows, g], j[rows, gp]
        den = jm - 2.0 * j0 + jp
        with np.errstate(divide="ignore", invalid="ignore"):
            delta = np.where(interior & (den < 0), 0.5 * (jm - jp) / den, 0.0)
        theta = self.grid[g] + np.clip(delta, -0.5, 0.5) * self.step
        return np.clip(theta, self.grid[0], self.grid[-1])

    def _mmse_from_objective(self, j, sigma2):
        logp = j / sigma2
        w = np.exp(logp - logp.max(axis=1, keepdims=True)) * self.trap
        return (w @ self.grid) / np.sum(w, axis=1)

    def estimate(self, r):
        return self._ml_from_objective(self.objective(r))

    def estimate_both(self, r, sigma2):
        """(ML estimate, MMSE estimate) from one evaluation of J."""
        j = self.objective(r)
        return self._ml_from_objective(j), self._mmse_from_objective(j, sigma2)


# ----------------------------------------------------------------------------
# Early-late DLL
# ----------------------------------------------------------------------------
class EarlyLateDLL:
    """Coherent early-late DLL with a first-order loop filter.

    For update j (coherent integration over the M chips of block j):
        E_j = sum_n r[n] c_j(t_n - theta_j + d/2)   (early replica)
        L_j = sum_n r[n] c_j(t_n - theta_j - d/2)   (late replica)
        D_j = (E_j - L_j - b_j) / k_bar             (~ theta_j - tau, linear region)
        theta_{j+1} = theta_j - K D_j

    c_j is the local replica built from block j's chips only (as a receiver's
    code generator produces it). It is correlated over a window that covers those
    chips plus the pulse tails, positioned by the DLL's OWN timing theta_j (as a
    receiver's integrate-and-dump follows its local code clock). Summing
    E_j - L_j over all blocks gives exactly the full-stream early-late
    correlation.

    b_j is block j's noise-free discriminator output at zero error. It is not
    zero because the block replica also correlates with neighbouring chips; with
    a periodic pilot this offset pattern repeats every period and, uncorrected,
    biases the loop (MSE/CRLB = 2.07 at 30 dB). The pilot is known, so the
    receiver precomputes b_j and subtracts it. b_j and the slopes k_j depend
    only on the pilot: recomputing them at other sampling phases changes the
    DLL MSE by < 0.1 % (validation/audit/t1_calib.py).

    k_bar is the nominal discriminator gain (mean over blocks of the noise-free
    slope k_j). A constant gain keeps the loop's sum of discriminator outputs
    equal to the full-stream correlation. (In an earlier design with a random,
    non-periodic stream, per-block normalisation by k_j was unstable because
    some k_j were near zero; for this pilot min k_j / k_bar = 0.29.)

    Carrier phase and amplitude are assumed known (coherent DLL), matching the
    MLE's assumptions. The pilot is the MLE's batch pilot repeated periodically
    (like a GNSS PRN code).

    Before tracking, a coarse ACQUISITION stage picks the initial theta: the
    first code period is correlated with the pilot on a grid of step `acq_step`
    over the prior window. Tracking then pulls in from any acquisition error
    inside the S-curve's pull-in range (measured by `pull_in_range`, ~1.9 chips).
    The estimate is theta after the last update.
    """

    GUARD = 32   # chips of slack on both sides so drifting windows stay in the stream

    def __init__(self, cfg, pilot):
        self.cfg = cfg
        m, sps, span, ta = cfg.dll_block_chips, cfg.sps, cfg.span, cfg.prior_width
        n_upd = cfg.dll_n_updates
        period = len(pilot)
        n_stream = int(np.ceil((n_upd * m + 2 * span + self.GUARD) / period)) * period
        self.chips = np.tile(np.asarray(pilot, dtype=float), n_stream // period)
        self.wave = BPSKWaveform(self.chips, cfg.rolloff, span)
        self.tau_c = self.GUARD + span + ta / 2.0 + 1.0
        n_total = int(np.ceil(self.tau_c + ta + n_stream + span + self.GUARD)) * sps
        self.t = np.arange(n_total) / sps
        # block j at zero timing offset: chips j*M ... (j+1)*M - 1, +- (span + 1) chips
        lead = span + 1.0
        start = np.round((self.tau_c + np.arange(n_upd) * m - lead) * sps).astype(int)
        width = int(np.ceil((m + 2 * lead) * sps))
        self.block_idx = start[:, None] + np.arange(width)[None, :]
        self.block_chips = np.arange(n_upd)[:, None] * m + np.arange(m)[None, :]
        self.slope, self.el_diff, self.offset = self._block_constants()
        self.k_bar = float(np.mean(self.slope))
        # acquisition: first code period, grid over the prior window
        n_acq = int(round(ta / cfg.acq_step)) + 1
        self.acq_grid = np.linspace(self.tau_c - ta / 2.0, self.tau_c + ta / 2.0, n_acq)
        lead_acq = span + 1 + ta / 2.0
        lo = int(np.floor((self.tau_c - lead_acq) * sps))
        hi = int(np.ceil((self.tau_c + period + lead_acq) * sps))
        self.acq_idx = np.arange(lo, hi)
        first = BPSKWaveform(np.asarray(pilot, dtype=float), cfg.rolloff, span)
        self.acq_templates = first(self.t[self.acq_idx][None, :] - self.acq_grid[:, None])

    def _replica(self, j, t, theta):
        """Block-j local replica  c_j(t - theta) = sum_{k in block j} a_k p(t - theta - k)."""
        k = self.block_chips[j]
        x = (t - theta)[..., None] - k
        return truncated_pulse(x, self.cfg.rolloff, self.cfg.span) @ self.chips[k]

    def _el(self, j, t, theta):
        d = self.cfg.el_spacing
        return self._replica(j, t, theta - d / 2.0), self._replica(j, t, theta + d / 2.0)

    def acquire(self, r):
        corr = r[:, self.acq_idx] @ self.acq_templates.T
        return self.acq_grid[np.argmax(corr, axis=1)]

    def _block_constants(self, h=1e-4):
        """Noise-free slope k_j, early-minus-late replica g_j and offset b_j per block."""
        slope, diff, offset = [], [], []
        for j, idx in enumerate(self.block_idx):
            t = self.t[idx]
            rx = self.wave(t - self.tau_c)
            e_p, l_p = self._el(j, t, self.tau_c + h)
            e_m, l_m = self._el(j, t, self.tau_c - h)
            slope.append((np.sum(rx * (e_p - l_p)) - np.sum(rx * (e_m - l_m))) / (2 * h))
            e0, l0 = self._el(j, t, self.tau_c)
            diff.append(e0 - l0)
            offset.append(np.sum(rx * (e0 - l0)))
        return np.array(slope), np.array(diff), np.array(offset)

    def received(self, tau, sigma2, rng):
        r = np.empty((tau.size, self.t.size))
        for c in range(0, tau.size, 200):                     # chunked: bounds memory
            r[c:c + 200] = self.wave(self.t[None, :] - tau[c:c + 200, None])
        r += np.sqrt(sigma2) * rng.standard_normal(r.shape)
        return r

    def run(self, tau, sigma2, rng, keep_history=0, checkpoints=(), theta_init=None):
        """Acquire (unless theta_init is given) and track each delay in tau.

        Returns a dict with
            theta   final estimates,
            theta0  initial (acquisition) estimates,
            hist    error trajectories of the first `keep_history` trials,
            cp      errors after each update count in `checkpoints` (len(cp), trials).
        """
        k_loop, sps = self.cfg.dll_gain, self.cfg.sps
        r = self.received(tau, sigma2, rng)
        theta0 = self.acquire(r) if theta_init is None else np.asarray(theta_init, float).copy()
        theta = theta0.copy()
        n_upd = len(self.block_idx)
        hist = np.empty((keep_history, n_upd + 1))
        cp = {}
        if keep_history:
            hist[:, 0] = theta[:keep_history] - tau[:keep_history]
        last = r.shape[1] - 1
        for j in range(n_upd):
            shift = np.round((theta - self.tau_c) * sps).astype(int)   # window follows theta
            idx = np.clip(self.block_idx[j][None, :] + shift[:, None], 0, last)
            early, late = self._el(j, self.t[idx], theta[:, None])
            corr = np.sum(np.take_along_axis(r, idx, axis=1) * (early - late), axis=1)
            theta = theta - k_loop * (corr - self.offset[j]) / self.k_bar
            if keep_history:
                hist[:, j + 1] = theta[:keep_history] - tau[:keep_history]
            if j + 1 in checkpoints:
                cp[j + 1] = theta - tau
        out = dict(theta=theta, theta0=theta0)
        if keep_history:
            out["hist"] = hist
        if checkpoints:
            out["cp"] = np.array([cp[c] for c in checkpoints])
        return out

    def pull_in_range(self, rng, step=0.05, max_err=4.0, tol=0.05):
        """Largest initial error (chips) from which the noise-free loop converges.

        Errors on both sides are tested; the range is the smaller of the two
        contiguous convergent intervals around zero.
        """
        eps = np.arange(step, max_err + step / 2, step)
        rng_out = []
        for sign in (+1.0, -1.0):
            tau = np.full(eps.size, self.tau_c)
            th = self.run(tau, 0.0, rng, theta_init=self.tau_c + sign * eps)["theta"]
            ok = np.abs(th - tau) < tol
            n_ok = np.argmin(ok) if not ok.all() else eps.size
            rng_out.append(eps[n_ok - 1] if n_ok > 0 else 0.0)
        return float(min(rng_out))

    def s_curve(self, eps):
        """Noise-free discriminator characteristic D(eps), averaged over all blocks.

        eps = theta - tau (chips). Near zero D(eps) ~ eps (unit slope by the k_bar
        normalisation); the loop is pulled back towards zero wherever D has the
        same sign as eps, which defines the pull-in region.
        """
        eps = np.atleast_1d(np.asarray(eps, dtype=float))
        rx = self.wave(self.t - self.tau_c)                  # noise-free, tau = tau_c
        theta = self.tau_c + eps
        shift = np.round(eps * self.cfg.sps).astype(int)
        out = np.zeros(eps.size)
        for j in range(len(self.block_idx)):
            idx = np.clip(self.block_idx[j][None, :] + shift[:, None], 0, rx.size - 1)
            early, late = self._el(j, self.t[idx], theta[:, None])
            out += (np.sum(rx[idx] * (early - late), axis=1) - self.offset[j]) / self.k_bar
        return out / len(self.block_idx)

    # ---- beyond Sprint 4 (Sprint 5 material): DLL variance theory --------------
    def equivalent_weighting(self):
        """Weighting waveform G of the linearised loop: eps_final = <w, G> + transient.

        Linearised about zero error, eps_{j+1} = a_j eps_j - (K / k_bar) <w, g_j>
        with a_j = 1 - K k_j / k_bar, so
            eps_final = prod(a) eps_0 - sum_j c_j <w, g_j>,
            c_j = (K / k_bar) prod_{i > j} a_i.
        This is exact for the linear loop and accounts for the noise shared by
        overlapping block windows and for the block-to-block slope variation.
        """
        k = self.cfg.dll_gain
        a = 1.0 - k * self.slope / self.k_bar
        tail = np.append(np.cumprod(a[::-1])[::-1][1:], 1.0)   # prod_{i>j} a_i
        c = k / self.k_bar * tail
        g = np.zeros(self.t.size)
        np.add.at(g, self.block_idx, c[:, None] * self.el_diff)
        return g

    def linear_theory_var(self, sigma2):
        """Noise-induced error variance of the linearised loop: sigma^2 ||G||^2."""
        return np.sum(self.equivalent_weighting() ** 2) * np.asarray(sigma2)
