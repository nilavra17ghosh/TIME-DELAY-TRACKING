"""Part 7b: DLL MSE vs elapsed updates at low SNR (project code, factor 20 full trajectories)."""
import sys, numpy as np
sys.path.insert(0, r"C:\Users\Nilavra Ghosh\Documents\7THSEM\SDET\project\src")
from tde import *
from multiprocessing import Pool
def job(db):
    cfg = Config(dll_settle_factor=20); sc = BatchScenario(cfg, np.random.default_rng(cfg.seed))
    dll = EarlyLateDLL(cfg, sc.wave.chips); rng = np.random.default_rng(int(db) + 500)
    tau = sc.draw_delays(400, rng)
    th, th0, h = dll.run(tau, sc.sigma2(10 ** (db / 10)), rng, keep_history=400)
    return db, h
if __name__ == "__main__":
    with Pool(3) as p: out = p.map(job, [-10.0, 0.0, 2.0])
    cols = [40, 80, 160, 240, 320]
    print("DLL MSE (chips^2) vs number of updates j, 400 trials; j=80 <-> factor 5, 160 <-> factor 10 (project), 320 <-> factor 20")
    for db, h in out:
        print(f"{db:5.0f} dB: " + "  ".join(f"j={c}: {np.mean(h[:, c]**2):7.3f}" for c in cols))
