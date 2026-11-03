"""Part 7: sensitivity of DLL MSE to the stream length (dll_settle_factor), project code."""
import sys
import time
import numpy as np

sys.path.insert(0, r"C:\Users\Nilavra Ghosh\Documents\7THSEM\SDET\project\src")
from tde import BatchScenario, Config, EarlyLateDLL, crlb  # noqa: E402

SNRS = [4.0, 8.0, 10.0, 14.0]
FACTORS = [5, 10, 20]


def job(args):
    db, fac, n = args
    cfg = Config(dll_settle_factor=fac)
    sc = BatchScenario(cfg, np.random.default_rng(cfg.seed))
    dll = EarlyLateDLL(cfg, sc.wave.chips)
    rng = np.random.default_rng(int(1000 * db) + fac)
    tau = sc.draw_delays(n, rng)
    t0 = time.time()
    keep = n if fac == 20 else 0
    out = dll.run(tau, sc.sigma2(10 ** (db / 10)), rng, keep_history=keep)
    theta, theta0 = out[0], out[1]
    hist = out[2] if keep else None
    return dict(db=db, fac=fac, e=theta - tau, e0=theta0 - tau, hist=hist,
                crlb=float(crlb(sc, 10 ** (db / 10))), n_upd=cfg.dll_n_updates, dt=time.time() - t0)


def boot(x2, rng, B=2000):
    m = np.array([x2[rng.integers(0, x2.size, x2.size)].mean() for _ in range(B)])
    return np.percentile(m, [2.5, 97.5])


if __name__ == "__main__":
    from multiprocessing import Pool
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
    jobs = [(db, f, n) for db in SNRS for f in FACTORS]
    with Pool(6) as pool:
        out = pool.map(job, jobs)
    summ = np.genfromtxt(r"C:\Users\Nilavra Ghosh\Documents\7THSEM\SDET\project\results\tables\summary.csv",
                         delimiter=",", names=True)
    rng = np.random.default_rng(9)
    print(f"n = {n} trials per cell")
    print("  dB  factor n_upd   MSE       95%CI               MSE/CRLB  P(|e|>0.5)  P(acq fail)  P(lost lock)  proj(f=10)")
    save = {}
    for d in out:
        e2 = d["e"] ** 2
        lo, hi = boot(e2, rng)
        acq = np.abs(d["e0"]) > 0.5
        i = int(np.where(summ["snr_db"] == d["db"])[0][0])
        print(f"{d['db']:4.0f} {d['fac']:6d} {d['n_upd']:5d}  {e2.mean():.3e} [{lo:.3e},{hi:.3e}] {e2.mean() / d['crlb']:8.2f}"
              f"  {np.mean(np.abs(d['e']) > 0.5):9.3f}  {acq.mean():10.3f}  {np.mean(~acq & (np.abs(d['e']) > 0.5)):10.3f}"
              f"    {summ['mse_dll'][i]:.3e}   ({d['dt']:.0f}s)")
        if d["hist"] is not None:
            save[f"hist_{d['db']:g}"] = d["hist"]
            save[f"crlb_{d['db']:g}"] = d["crlb"]
    # MSE vs number of updates from full trajectories (factor 20 runs)
    print("\nMSE vs elapsed updates (factor-20 runs, same trials):")
    cols = [40, 80, 160, 240, 320]
    print("  dB  " + "  ".join(f"j={c:3d}" for c in cols) + "   [P(|e|>0.5) at same j]")
    for k in sorted(save):
        if not k.startswith("hist"):
            continue
        h = save[k]
        print(f"{k[5:]:>4s}  " + "  ".join(f"{np.mean(h[:, c] ** 2):.3e}" for c in cols) + "   "
              + " ".join(f"{np.mean(np.abs(h[:, c]) > 0.5):.3f}" for c in cols))
    np.savez("p7_out.npz", **save)
