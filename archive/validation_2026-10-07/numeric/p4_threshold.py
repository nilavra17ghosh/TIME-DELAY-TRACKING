"""Part 4: threshold robustness with the project's own code (fine SNR grid, bootstrap, 5 pilots, Ta=4/16)."""
import sys
import time
import numpy as np

sys.path.insert(0, r"C:\Users\Nilavra Ghosh\Documents\7THSEM\SDET\project\src")
from tde import BatchMLE, BatchScenario, Config, crlb, zzb, rms_bandwidth_sq  # noqa: E402

SEEDS = [2026, 11, 22, 33, 44]
GRIDS = {4.0: np.arange(2.0, 16.5, 1.0), 16.0: np.arange(6.0, 20.5, 1.0)}


def threshold_snr(snr_db, ratio, level=2.0):
    """Same rule as run_baseline.threshold_snr (highest crossing, log-ratio interpolation)."""
    lr, lv = np.log10(ratio), np.log10(level)
    for i in range(len(snr_db) - 1, 0, -1):
        if lr[i] <= lv < lr[i - 1]:
            f = (lv - lr[i]) / (lr[i - 1] - lr[i])
            return snr_db[i] + f * (snr_db[i - 1] - snr_db[i])
    return np.nan


def job(args):
    ta, seed, n = args
    cfg = Config(prior_width=ta, seed=seed)
    sc = BatchScenario(cfg, np.random.default_rng(cfg.seed))
    mle = BatchMLE(sc)
    rng = np.random.default_rng(seed * 7919 + int(ta))
    grid = GRIDS[ta]
    err = np.empty((grid.size, n))
    t0 = time.time()
    for i, db in enumerate(grid):
        e = 10 ** (db / 10)
        tau = sc.draw_delays(n, rng)
        for c in range(0, n, 500):
            sl = slice(c, c + 500)
            err[i, sl] = mle.estimate(sc.observe(tau[sl], e, rng)) - tau[sl]
    snr_f = np.linspace(grid[0] - 6, grid[-1] + 6, 1201)
    z = zzb(sc, 10 ** (snr_f / 10))
    c = crlb(sc, 10 ** (snr_f / 10))
    return dict(ta=ta, seed=seed, err=err, grid=grid, crlb=crlb(sc, 10 ** (grid / 10)),
                zzb=zzb(sc, 10 ** (grid / 10)), thr_zzb=threshold_snr(snr_f, z / c),
                beta2=rms_bandwidth_sq(sc), dt=time.time() - t0)


if __name__ == "__main__":
    from multiprocessing import Pool
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 4000
    jobs = [(ta, s, n) for ta in (4.0, 16.0) for s in SEEDS]
    with Pool(10) as pool:
        out = pool.map(job, jobs)
    brng = np.random.default_rng(123)
    B = 1000
    save = {}
    print(f"n_trials per SNR = {n}, bootstrap B = {B}")
    print(" Ta  seed  beta2   ZZBthr   MLEthr(1dB)  95%CI           MLEthr(2dB grid)  gap    gap95%CI")
    for d in out:
        ta, seed, err, grid, cr = d["ta"], d["seed"], d["err"], d["grid"], d["crlb"]
        mse = np.mean(err ** 2, 1)
        thr = threshold_snr(grid, mse / cr)
        ev = (grid % 2 == 0)
        thr2 = threshold_snr(grid[ev], (mse / cr)[ev])
        bt = np.empty(B)
        e2 = err ** 2
        for b in range(B):
            idx = brng.integers(0, e2.shape[1], e2.shape)
            bt[b] = threshold_snr(grid, np.take_along_axis(e2, idx, 1).mean(1) / cr)
        lo, hi = np.nanpercentile(bt, [2.5, 97.5])
        print(f"{ta:4.0f} {seed:5d} {d['beta2']:.3f}  {d['thr_zzb']:6.2f}   {thr:6.2f}     [{lo:5.2f},{hi:5.2f}]   {thr2:6.2f}"
              f"          {thr - d['thr_zzb']:5.2f}  [{lo - d['thr_zzb']:4.2f},{hi - d['thr_zzb']:4.2f}]  ({d['dt']:.0f}s)")
        key = f"ta{ta:g}_s{seed}"
        save[key + "_mse"], save[key + "_grid"], save[key + "_crlb"], save[key + "_zzb"] = mse, grid, cr, d["zzb"]
        save[key + "_gross"] = np.mean(np.abs(err) > 0.5, 1)
        save[key + "_boot"] = bt
        save[key + "_thrzzb"] = d["thr_zzb"]
    np.savez("p4_out.npz", **save)
