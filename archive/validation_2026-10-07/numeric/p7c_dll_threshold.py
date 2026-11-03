"""Part 7c: DLL threshold re-estimate (project code, default factor 10), 6 x 3000 trials per SNR, bootstrap."""
import sys, numpy as np
sys.path.insert(0, r"C:\Users\Nilavra Ghosh\Documents\7THSEM\SDET\project\src")
from tde import *
from multiprocessing import Pool
SNR = [7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0]
def job(a):
    db, rep = a
    cfg = Config(); sc = BatchScenario(cfg, np.random.default_rng(cfg.seed)); dll = EarlyLateDLL(cfg, sc.wave.chips)
    rng = np.random.default_rng(int(db * 100) + rep * 7 + 31)
    tau = sc.draw_delays(2000, rng)
    th, th0 = dll.run(tau, sc.sigma2(10 ** (db / 10)), rng)
    return db, th - tau, float(crlb(sc, 10 ** (db / 10)))
def thr(s, r, level=2.0):
    lr, lv = np.log10(r), np.log10(level)
    for i in range(len(s) - 1, 0, -1):
        if lr[i] <= lv < lr[i - 1]:
            f = (lv - lr[i]) / (lr[i - 1] - lr[i]); return s[i] + f * (s[i - 1] - s[i])
    return np.nan
if __name__ == "__main__":
    with Pool(12) as p: out = p.map(job, [(db, r) for db in SNR for r in range(3)])
    e2 = {db: np.concatenate([o[1] ** 2 for o in out if o[0] == db]) for db in SNR}
    cr = np.array([[o[2] for o in out if o[0] == db][0] for db in SNR])
    rng = np.random.default_rng(1)
    for db, c in zip(SNR, cr):
        x = e2[db]; bt = [x[rng.integers(0, x.size, x.size)].mean() for _ in range(1000)]
        print(f"DLL {db:5.1f} dB n={x.size}: MSE/CRLB {x.mean() / c:6.2f} 95%CI {np.round(np.percentile(bt, [2.5, 97.5]) / c, 2)}  P(|e|>0.5) {np.mean(x > 0.25):.4f}")
    s = np.array(SNR)
    print("DLL threshold:", round(thr(s, np.array([e2[d].mean() for d in SNR]) / cr), 2))
    bts = [thr(s, np.array([e2[d][rng.integers(0, e2[d].size, e2[d].size)].mean() for d in SNR]) / cr) for _ in range(1000)]
    print("bootstrap 95% CI:", np.nanpercentile(bts, [2.5, 97.5]).round(2))
    # 2 dB grid version for comparison with project
    ev = np.isin(s, [8, 10, 12, 14])
    print("DLL threshold (2 dB grid points only):", round(thr(s[ev], (np.array([e2[d].mean() for d in SNR]) / cr)[ev]), 2))
