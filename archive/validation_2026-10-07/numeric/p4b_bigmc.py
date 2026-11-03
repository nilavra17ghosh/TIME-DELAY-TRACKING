"""Part 4b: large-trial MLE runs near threshold for pilot 2026, Ta=4 and Ta=16 (project code)."""
import sys, time, numpy as np
sys.path.insert(0, r"C:\Users\Nilavra Ghosh\Documents\7THSEM\SDET\project\src")
from tde import BatchMLE, BatchScenario, Config, crlb
from multiprocessing import Pool

def job(a):
    ta, db, n = a
    cfg = Config(prior_width=ta); sc = BatchScenario(cfg, np.random.default_rng(cfg.seed)); m = BatchMLE(sc)
    rng = np.random.default_rng(int(db * 10) + int(ta) * 1000 + 99)
    err = np.empty(n)
    for c in range(0, n, 500):
        tau = sc.draw_delays(500, rng)
        err[c:c + 500] = m.estimate(sc.observe(tau, 10 ** (db / 10), rng)) - tau
    return ta, db, err, float(crlb(sc, 10 ** (db / 10)))

def thr(snr_db, ratio, level=2.0):
    lr, lv = np.log10(ratio), np.log10(level)
    for i in range(len(snr_db) - 1, 0, -1):
        if lr[i] <= lv < lr[i - 1]:
            f = (lv - lr[i]) / (lr[i - 1] - lr[i]); return snr_db[i] + f * (snr_db[i - 1] - snr_db[i])
    return np.nan

if __name__ == "__main__":
    jobs = [(4.0, db, 40000) for db in (8.0, 8.5, 9.0, 9.5, 10.0, 11.0)] + \
           [(16.0, db, 60000) for db in (12.0, 12.5, 13.0, 13.5, 14.0, 15.0)]
    with Pool(12) as p: out = p.map(job, jobs)
    rng = np.random.default_rng(4)
    res = {}
    for ta in (4.0, 16.0):
        rows = [o for o in out if o[0] == ta]
        g = np.array([o[1] for o in rows]); cr = np.array([o[3] for o in rows])
        e2 = [o[2] ** 2 for o in rows]
        print(f"Ta = {ta:g}")
        for (t, db, e, c) in rows:
            x = e ** 2
            bt = [x[rng.integers(0, x.size, x.size)].mean() for _ in range(500)]
            lo, hi = np.percentile(bt, [2.5, 97.5])
            print(f"  {db:5.1f} dB n={e.size}  MSE/CRLB {x.mean() / c:6.2f} [{lo / c:5.2f},{hi / c:5.2f}]  P(|e|>0.5) {np.mean(np.abs(e) > 0.5):.5f}  "
                  f"E[e^2 | gross] {x[np.abs(e) > 0.5].mean() if np.any(np.abs(e) > 0.5) else np.nan:.2f}")
        ratio = np.array([x.mean() for x in e2]) / cr
        bts = []
        for b in range(500):
            bts.append(thr(g, np.array([x[rng.integers(0, x.size, x.size)].mean() for x in e2]) / cr))
        print(f"  threshold {thr(g, ratio):.2f} dB, bootstrap 95% CI {np.nanpercentile(bts, [2.5, 97.5]).round(2)}")
