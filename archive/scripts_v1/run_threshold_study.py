"""Threshold study: precise threshold SNRs with bootstrap CIs, across pilots and prior widths.

Run:  python run_threshold_study.py           (~20-30 min)
      python run_threshold_study.py --quick   (smoke test)

Part A (batch estimators): MLE and MMSE, prior widths Ta = 1, 4, 16 chips,
        5 random pilots (Ta = 4, 16), 1 dB SNR grid, 10,000 trials per point.
Part B (DLL): pilot 2026, Ta = 4, 1 dB grid, 4,000 trials per point.

Two threshold definitions, both with 95 % bootstrap CIs:
    mse:   MSE/CRLB first exceeds 2 (3 dB)       - the usual definition
    gross: P(|e| > 0.5 chip) first exceeds 1 %   - robust to single outliers
Outputs: results/tables/thresholds.csv, results/data/threshold_study.npz,
         results/figures/fig6_prior_width.png, fig8_thresholds.png
"""
import argparse
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, "src")
from tde import (BatchMLE, BatchScenario, Config, EarlyLateDLL, bootstrap_thresholds,  # noqa: E402
                 crlb, crossing_snr, rms_bandwidth_sq, zzb)
from tde.plotting import COLORS, INK2, MUTED, TA_RAMP, apply_style, plt  # noqa: E402

OUT = Path("results")
PILOTS = (2026, 11, 22, 33, 44)
GRIDS = {1.0: np.arange(-4.0, 13.0, 1.0), 4.0: np.arange(2.0, 17.0, 1.0), 16.0: np.arange(6.0, 21.0, 1.0)}


def batch_run(ta, seed, snr_db, n, rng):
    cfg = Config(prior_width=ta, seed=seed)
    sc = BatchScenario(cfg, np.random.default_rng(seed))
    mle = BatchMLE(sc)
    en0 = 10 ** (snr_db / 10)
    e_ml, e_mm = np.empty((snr_db.size, n)), np.empty((snr_db.size, n))
    for i, e in enumerate(en0):
        tau = sc.draw_delays(n, rng)
        for c in range(0, n, 1000):
            sl = slice(c, c + 1000)
            ml, mm = mle.estimate_both(sc.observe(tau[sl], e, rng), sc.sigma2(e))
            e_ml[i, sl], e_mm[i, sl] = ml - tau[sl], mm - tau[sl]
    fine = np.linspace(snr_db[0] - 6, snr_db[-1], 600)
    zf, cf = zzb(sc, 10 ** (fine / 10)), crlb(sc, 10 ** (fine / 10))
    return dict(sc=sc, e_ml=e_ml, e_mm=e_mm, cr=crlb(sc, en0), fine=fine, zzb_fine=zf, crlb_fine=cf,
                thr_zzb=crossing_snr(fine, zf / cf, 2.0))


def fmt(t):
    return f"{t[0]:6.2f} [{t[1]:5.2f}, {t[2]:5.2f}]"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    args = ap.parse_args()
    n_batch, n_dll, n_boot = (1000, 400, 50) if args.quick else (10000, 4000, 300)
    pilots = PILOTS[:2] if args.quick else PILOTS
    rng = np.random.default_rng(99)
    for sub in ("data", "figures", "tables"):
        (OUT / sub).mkdir(parents=True, exist_ok=True)
    apply_style()

    rows, store, curves = [], {}, {}
    print("Part A: batch MLE / MMSE")
    for ta in (1.0, 4.0, 16.0):
        for seed in (pilots if ta != 1.0 else (2026,)):
            t0 = time.time()
            r = batch_run(ta, seed, GRIDS[ta], n_batch, rng)
            b_ml = bootstrap_thresholds(GRIDS[ta], r["e_ml"], r["cr"], n_boot=n_boot, rng=rng)
            b_mm = bootstrap_thresholds(GRIDS[ta], r["e_mm"], r["cr"], n_boot=n_boot, rng=rng)
            rows.append(dict(ta=ta, pilot=seed, beta2=rms_bandwidth_sq(r["sc"]),
                             zzb=r["thr_zzb"], mle=b_ml, mmse=b_mm))
            store[f"ta{ta:g}_p{seed}_mle"] = r["e_ml"]
            store[f"ta{ta:g}_p{seed}_mmse"] = r["e_mm"]
            if seed == 2026:
                curves[ta] = r
            print(f"  Ta={ta:4.1f} pilot {seed:4d}: ZZB {r['thr_zzb']:6.2f} | MLE mse {fmt(b_ml['mse'])} "
                  f"gross {fmt(b_ml['gross'])} | MMSE mse {fmt(b_mm['mse'])}  ({time.time() - t0:.0f}s)",
                  flush=True)

    print("Part B: DLL (pilot 2026, Ta = 4)")
    cfg = Config(prior_width=4.0)
    sc = BatchScenario(cfg, np.random.default_rng(cfg.seed))
    dll = EarlyLateDLL(cfg, sc.wave.chips)
    snr_d = GRIDS[4.0]
    e_dll = np.empty((snr_d.size, n_dll))
    for i, db in enumerate(snr_d):
        e = 10 ** (db / 10)
        tau = sc.draw_delays(n_dll, rng) - sc.tau_c + dll.tau_c
        e_dll[i] = dll.run(tau, sc.sigma2(e), rng)["theta"] - tau
    b_dll = bootstrap_thresholds(snr_d, e_dll, crlb(sc, 10 ** (snr_d / 10)), n_boot=n_boot, rng=rng)
    store["ta4_p2026_dll"] = e_dll
    print(f"  DLL: mse {fmt(b_dll['mse'])}  gross {fmt(b_dll['gross'])}")

    # ---- table ----
    lines = ["ta,pilot,beta2,thr_zzb,thr_mle_mse,lo,hi,thr_mle_gross,lo,hi,thr_mmse_mse,lo,hi,"
             "thr_mmse_gross,lo,hi,gap_mle_zzb,gap_mmse_zzb"]
    for r in rows:
        v = [r["ta"], r["pilot"], r["beta2"], r["zzb"], *r["mle"]["mse"], *r["mle"]["gross"],
             *r["mmse"]["mse"], *r["mmse"]["gross"], r["mle"]["mse"][0] - r["zzb"], r["mmse"]["mse"][0] - r["zzb"]]
        lines.append(",".join(f"{x:.4g}" for x in v))
    lines.append(f"4,2026,nan,nan,{','.join(f'{x:.4g}' for x in (*b_dll['mse'], *b_dll['gross']))},"
                 f"nan,nan,nan,nan,nan,nan,nan,nan  # DLL row (columns 5-10 = DLL mse / gross)")
    (OUT / "tables" / "thresholds.csv").write_text("\n".join(lines) + "\n")
    np.savez(OUT / "data" / "threshold_study.npz", **store)

    print("\nGap to ZZB threshold (MSE criterion), dB:")
    for ta in (4.0, 16.0):
        g_ml = [r["mle"]["mse"][0] - r["zzb"] for r in rows if r["ta"] == ta]
        g_mm = [r["mmse"]["mse"][0] - r["zzb"] for r in rows if r["ta"] == ta]
        print(f"  Ta={ta:4.1f}: MLE-ZZB {np.round(g_ml, 2)} mean {np.mean(g_ml):.2f} | "
              f"MMSE-ZZB {np.round(g_mm, 2)} mean {np.mean(g_mm):.2f}")

    # ---- fig 6: prior width (pilot 2026) ----
    fig, ax = plt.subplots(figsize=(9, 6))
    first = curves[4.0]
    ax.semilogy(first["fine"], first["crlb_fine"], color=COLORS["crlb"], ls="--", lw=1.5, label="CRLB")
    for color, ta in zip(TA_RAMP, (1.0, 4.0, 16.0)):
        r = curves[ta]
        ax.semilogy(r["fine"], r["zzb_fine"], color=color, lw=1.8, label=f"ZZB, Ta = {ta:g}")
        ax.semilogy(GRIDS[ta], np.mean(r["e_ml"] ** 2, 1), "o", color=color, markeredgecolor="#fcfcfb",
                    label=f"MLE, Ta = {ta:g}")
        ax.semilogy(GRIDS[ta], np.mean(r["e_mm"] ** 2, 1), "D", color=color, ms=4, mfc="none",
                    label=f"MMSE, Ta = {ta:g}")
    ax.set_xlabel("batch E/N0 (dB)")
    ax.set_ylabel("MSE (chips²)")
    ax.set_title(f"Prior width Ta: wider uncertainty, later threshold ({n_batch} trials/point, pilot 2026)")
    ax.legend(loc="lower left", ncol=2, fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "figures" / "fig6_prior_width.png")
    plt.close(fig)

    # ---- fig 8: thresholds with CIs, per pilot ----
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6), sharey=False)
    for a, ta in zip(ax, (4.0, 16.0)):
        sel = [r for r in rows if r["ta"] == ta]
        y = np.arange(len(sel))
        for off, key, name, mk in ((-0.2, "zzb", "ZZB", "|"), (0.0, "mmse", "MMSE", "D"), (0.2, "mle", "MLE", "o")):
            for k, r in enumerate(sel):
                if key == "zzb":
                    a.plot(r["zzb"], k + off, mk, color=COLORS["zzb"], ms=12, mew=2.5,
                           label=name if k == 0 else None)
                else:
                    est, lo, hi = r[key]["mse"]
                    a.errorbar(est, k + off, xerr=[[est - lo], [hi - est]], fmt=mk, color=COLORS[key], ms=5,
                               capsize=3, label=name if k == 0 else None)
        if ta == 4.0:
            est, lo, hi = b_dll["mse"]
            a.errorbar(est, len(sel) - 0.4, xerr=[[est - lo], [hi - est]], fmt="s", color=COLORS["dll"],
                       capsize=3, label="DLL (pilot 2026)")
        a.set_yticks(y, [f"pilot {r['pilot']}" for r in sel])
        a.set_xlabel("threshold: batch E/N0 where MSE first exceeds 2×CRLB (dB)")
        a.set_title(f"Ta = {ta:g} chips (bars = 95 % bootstrap CI)")
        a.legend(loc="lower right", fontsize=8)
    fig.suptitle("Threshold SNR by pilot: the gap to the ZZB depends on the pilot", color=INK2)
    fig.tight_layout()
    fig.savefig(OUT / "figures" / "fig8_thresholds.png")
    plt.close(fig)
    print(f"Saved {OUT / 'tables' / 'thresholds.csv'}, fig6_prior_width.png, fig8_thresholds.png")


if __name__ == "__main__":
    main()
