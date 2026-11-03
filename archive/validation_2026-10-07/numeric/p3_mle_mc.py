"""Part 3 (+ data for part 6): independent MLE Monte Carlo.

Own signal generator (tap-sum over chips, symmetric |t|<=6 truncation), own noise calibration,
dense correlation grid (1/32 chip) + vectorised golden-section refinement of the exact
log-likelihood J(theta) = <r, s_theta> - ||s_theta||^2/2 (25 iterations -> ~4e-7 chip).
Prior centre deliberately at a NON-integer sampling phase (different from the project).
"""
import sys
import time
import numpy as np
from vcommon import pilot, rc_trunc

PROJ = r"C:\Users\Nilavra Ghosh\Documents\7THSEM\SDET\project\results\tables\summary.csv"
summ = np.genfromtxt(PROJ, delimiter=",", names=True)

TA, SPS, SPAN = 4.0, 4, 6
chips = pilot()
NCH = chips.size
TAU_C = SPAN + TA / 2 + 1.37          # non-integer -> different sampling phase than project
T = np.arange(int((NCH + 2 * SPAN + TA + 5) * SPS)) / SPS
LO, HI = TAU_C - TA / 2, TAU_C + TA / 2


def sig(theta):
    """s(t_n - theta) for each theta (shape (m,)) -> (m, N)."""
    u = T[None, :] - np.asarray(theta)[:, None]
    k0 = np.round(u)
    out = np.zeros(u.shape)
    for o in range(-SPAN - 1, SPAN + 2):
        k = k0 + o
        d = u - k
        ok = (k >= 0) & (k < NCH) & (np.abs(d) <= SPAN)
        a = np.where(ok, chips[np.clip(k, 0, NCH - 1).astype(int)], 0.0)
        out += a * rc_trunc(d, 0.5, SPAN)
    return out


E_SIG = float(np.sum(sig([TAU_C]) ** 2))
GRID = np.linspace(LO, HI, int(TA * 32) + 1)
STEP = GRID[1] - GRID[0]
TPL = sig(GRID)
HE = 0.5 * np.sum(TPL ** 2, axis=1)


def J(r, theta):
    s = sig(theta)
    return np.sum(r * s, axis=1) - 0.5 * np.sum(s * s, axis=1)


def mle(r, iters=20):
    j = r @ TPL.T - HE
    g = np.argmax(j, axis=1)
    a = np.clip(GRID[g] - STEP, LO, HI)
    b = np.clip(GRID[g] + STEP, LO, HI)
    gr = (np.sqrt(5) - 1) / 2
    c, d = b - gr * (b - a), a + gr * (b - a)
    fc, fd = J(r, c), J(r, d)
    for _ in range(iters):
        left = fc > fd                     # max in [a, d]
        b = np.where(left, d, b)
        a = np.where(left, a, c)
        nc = np.where(left, b - gr * (b - a), d)
        nd = np.where(left, c, a + gr * (b - a))
        fnew = J(r, np.where(left, nc, nd))
        fc, fd = np.where(left, fnew, fd), np.where(left, fc, fnew)
        c, d = nc, nd
    th = 0.5 * (a + b)
    # also consider the window end points (MLE may sit exactly on the boundary)
    cand = np.stack([th, np.full_like(th, LO), np.full_like(th, HI)], 1)
    jc = np.stack([J(r, cand[:, i]) for i in range(3)], 1)
    return cand[np.arange(len(th)), np.argmax(jc, 1)]


def run(db, n, rng, chunk=500, signal=True):
    en0 = 10 ** (db / 10)
    sigma = np.sqrt(E_SIG / (2 * en0))
    tau = rng.uniform(LO, HI, n)
    est = np.empty(n)
    for c in range(0, n, chunk):
        sl = slice(c, c + chunk)
        r = sigma * rng.standard_normal((len(tau[sl]), T.size))
        if signal:
            r += sig(tau[sl])
        est[sl] = mle(r)
    return tau, est


def boot_ci(x2, B=4000, rng=None):
    rng = rng or np.random.default_rng(1)
    n = x2.size
    m = np.array([x2[rng.integers(0, n, n)].mean() for _ in range(B)])
    return np.percentile(m, [2.5, 97.5])


def job(args):
    db, n, seed, signal = args
    rng = np.random.default_rng(seed)
    t0 = time.time()
    tau, est = run(db, n, rng, signal=signal)
    return db, signal, tau, est, time.time() - t0


if __name__ == "__main__":
    from multiprocessing import Pool
    tau0 = np.linspace(LO, HI, 401)
    e0 = mle(sig(tau0)) - tau0
    print(f"noise-free MSE {np.mean(e0 ** 2):.2e} chips^2 (CRLB at 30 dB = 2.04e-4)")
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 5000
    snrs = [4, 6, 8, 10, 12, 20, 30]
    jobs = [(db, n, 1000 + k, True) for k, db in enumerate(snrs)]
    jobs += [(db, 10000, 2000 + k, True) for k, db in enumerate([-10, -20, -40])]
    jobs += [(0, 10000, 3000, False)]
    with Pool(11) as pool:
        out = pool.map(job, jobs)
    res = {}
    brng = np.random.default_rng(5)
    print(" dB    n    myMSE     95%CI(bootstrap)       projMSE   z(diff)  my/CRLB  gross%   n=2000 sampling 95% range")
    for db, signal, tau, est, dt in out:
        key = str(db) if signal else "noise"
        res[f"tau_{key}"], res[f"est_{key}"] = tau, est
        if not signal or db not in snrs:
            continue
        e2 = (est - tau) ** 2
        lo, hi = boot_ci(e2, rng=brng)
        i = int(np.where(summ["snr_db"] == db)[0][0])
        p = summ["mse_mle"][i]
        se_mine = e2.std(ddof=1) / np.sqrt(e2.size)
        se_proj = e2.std(ddof=1) / np.sqrt(2000)
        sub = np.array([e2[brng.integers(0, e2.size, 2000)].mean() for _ in range(4000)])
        plo, phi = np.percentile(sub, [2.5, 97.5])
        z = (p - e2.mean()) / np.hypot(se_mine, se_proj)
        cr = summ["crlb"][i]
        print(f"{db:3d} {e2.size:5d} {e2.mean():.4e} [{lo:.4e},{hi:.4e}] {p:.4e} {z:+6.2f} {e2.mean() / cr:7.3f} "
              f"{100 * np.mean(np.abs(est - tau) > 0.5):6.2f}   [{plo:.4e},{phi:.4e}]  ({dt:.0f}s)")
    np.savez("p3_out.npz", LO=LO, HI=HI, **res)
    print("saved")
