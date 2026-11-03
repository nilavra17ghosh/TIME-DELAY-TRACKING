"""Threshold estimation with bootstrap confidence intervals."""
import numpy as np

GROSS = 0.5        # gross error: |e| > 0.5 chip (outside the correlation main lobe)


def crossing_snr(snr_db, y, level):
    """Highest SNR at which y (increasing towards low SNR) crosses `level`.

    Scans down from the highest SNR; linear interpolation in (dB, log y).
    Returns nan if never crossed.
    """
    with np.errstate(divide="ignore"):
        ly, lv = np.log10(np.maximum(y, 1e-300)), np.log10(level)
    for i in range(len(snr_db) - 1, 0, -1):
        if ly[i] <= lv < ly[i - 1]:
            f = (lv - ly[i]) / (ly[i - 1] - ly[i])
            return snr_db[i] + f * (snr_db[i - 1] - snr_db[i])
    return np.nan


def thresholds(snr_db, err, crlb, ratio=2.0, gross_level=0.01):
    """Two threshold definitions from an error matrix err (n_snr, n_trials).

    mse:   MSE / CRLB first exceeds `ratio` (3 dB).
    gross: P(|e| > GROSS) first exceeds `gross_level` (1 %). Less sensitive to a
           single outlier than the MSE criterion, because one outlier changes P
           by only 1/n while it can change the MSE by orders of magnitude.
    """
    mse = np.mean(err ** 2, axis=1)
    p_gross = np.mean(np.abs(err) > GROSS, axis=1)
    return crossing_snr(snr_db, mse / crlb, ratio), crossing_snr(snr_db, p_gross, gross_level)


def bootstrap_thresholds(snr_db, err, crlb, n_boot=300, rng=None, **kw):
    """Point estimates and 95 % percentile CIs of both thresholds.

    Trials are resampled independently at each SNR (they are independent runs).
    Returns dict(mse=(est, lo, hi), gross=(est, lo, hi)).
    """
    rng = np.random.default_rng(0) if rng is None else rng
    est = thresholds(snr_db, err, crlb, **kw)
    n = err.shape[1]
    boots = np.empty((n_boot, 2))
    for b in range(n_boot):
        idx = rng.integers(0, n, size=(err.shape[0], n))
        boots[b] = thresholds(snr_db, np.take_along_axis(err, idx, axis=1), crlb, **kw)
    out = {}
    for k, name in enumerate(("mse", "gross")):
        lo, hi = np.nanpercentile(boots[:, k], [2.5, 97.5])
        out[name] = (est[k], lo, hi)
    return out
