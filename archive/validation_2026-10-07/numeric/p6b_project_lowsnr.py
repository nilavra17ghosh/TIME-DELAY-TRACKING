"""Project-code MLE at -10 dB and pure noise (20000 trials) for comparison with own implementation."""
import sys, numpy as np
sys.path.insert(0, r"C:\Users\Nilavra Ghosh\Documents\7THSEM\SDET\project\src")
from tde import *
cfg = Config(); sc = BatchScenario(cfg, np.random.default_rng(cfg.seed)); m = BatchMLE(sc)
rng = np.random.default_rng(777)
for db in [-10.0, -100.0]:
    tau = sc.draw_delays(20000, rng); est = np.empty_like(tau)
    for c in range(0, tau.size, 500):
        est[c:c+500] = m.estimate(sc.observe(tau[c:c+500], 10**(db/10), rng))
    e2 = (est - tau)**2
    bt = [e2[rng.integers(0, e2.size, e2.size)].mean() for _ in range(1000)]
    edge = np.mean((est <= m.grid[0] + 1e-9) | (est >= m.grid[-1] - 1e-9))
    print(f"project code {db:6.0f} dB: MSE {e2.mean():.4f} 95%CI {np.percentile(bt,[2.5,97.5])}, P(edge) {edge:.3f}, Var(est) {est.var():.3f}")
