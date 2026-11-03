"""Baseline experiment: DLL vs batch MLE / MMSE vs CRLB / ZZB over SNR.

Run:  python run_baseline.py            (full run, ~10 min)
      python run_baseline.py --quick    (fewer trials, for a smoke test)

Outputs: results/data/snr_sweep.npz, results/tables/summary.csv, results/figures/fig1-5,7.
The precise threshold study (fine SNR grid, many trials, several pilots, prior
widths) is in run_threshold_study.py; thresholds printed here (2 dB grid,
2000 trials) are indicative only.
"""
import argparse
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, "src")
from tde import (GROSS, BatchMLE, BatchScenario, Config, EarlyLateDLL, bootstrap_thresholds,  # noqa: E402
                 crlb, crossing_snr, prior_variance, zzb)
from tde.plotting import COLORS, INK2, MUTED, SNR_RAMP, apply_style, plt  # noqa: E402

THRESH_RATIO = 2.0          # MSE threshold criterion: MSE first exceeds 2x CRLB (3 dB)
CHECKPOINTS = (20, 40, 80, 160)   # DLL updates at which the error is recorded
OUT = Path("results")


def snr_sweep(cfg, rng):
    sc = BatchScenario(cfg, np.random.default_rng(cfg.seed))
    mle = BatchMLE(sc)
    dll = EarlyLateDLL(cfg, sc.wave.chips)       # same pilot, repeated periodically
    checkpoints = tuple(c for c in CHECKPOINTS if c <= cfg.dll_n_updates)
    en0 = 10 ** (cfg.snr_db / 10)
    n, m = cfg.n_trials, cfg.snr_db.size
    err = {k: np.empty((m, n)) for k in ("mle", "mmse", "dll", "acq")}
    err_cp = np.empty((m, len(checkpoints), n))
    hist = {}
    for i, (db, e) in enumerate(zip(cfg.snr_db, en0)):
        t0 = time.time()
        tau = sc.draw_delays(n, rng)
        for c in range(0, n, 500):
            sl = slice(c, c + 500)
            ml, mm = mle.estimate_both(sc.observe(tau[sl], e, rng), sc.sigma2(e))
            err["mle"][i, sl], err["mmse"][i, sl] = ml - tau[sl], mm - tau[sl]
        tau = sc.draw_delays(n, rng) - sc.tau_c + dll.tau_c     # same prior, DLL's time origin
        out = dll.run(tau, sc.sigma2(e), rng, keep_history=12, checkpoints=checkpoints)
        err["dll"][i], err["acq"][i] = out["theta"] - tau, out["theta0"] - tau
        err_cp[i] = out["cp"]
        hist[db] = out["hist"]
        print(f"  E/N0 = {db:5.1f} dB   MSE: MLE {np.mean(err['mle'][i] ** 2):.3e}  "
              f"MMSE {np.mean(err['mmse'][i] ** 2):.3e}  DLL {np.mean(err['dll'][i] ** 2):.3e}  "
              f"({time.time() - t0:.1f}s)", flush=True)
    return sc, mle, dll, err, err_cp, checkpoints, hist


def stream_scenario(cfg, sc, dll):
    """The DLL's observation: the pilot repeated over the chips its blocks cover."""
    n_rep = int(np.ceil(cfg.dll_n_updates * cfg.dll_block_chips / cfg.n_chips))
    sc_s = BatchScenario(cfg, None, chips=np.tile(sc.wave.chips, n_rep))
    return sc_s, sc_s.energy / sc.energy     # energy ratio maps batch E/N0 -> stream E/N0


# ----------------------------------------------------------------------------
# figures
# ----------------------------------------------------------------------------
def fig_signal(sc, rng, path):
    tau = sc.tau_c + 0.37
    en0 = 10 ** 2.0
    t = sc.t
    tx = sc.wave(t)                                      # transmitted, chip k centred at t = k
    rx_clean = sc.wave(t - tau)
    rx = sc.observe(np.array([tau]), en0, rng)[0]
    fig, ax = plt.subplots(3, 1, figsize=(9, 6), sharex=True)
    ax[0].plot(t, tx, color=COLORS["mle"], lw=1.6)
    ax[0].stem(np.arange(sc.n_chips), sc.wave.chips, linefmt=MUTED, markerfmt=" ", basefmt=" ")
    ax[0].set_title("Transmitted pilot s(t): raised-cosine BPSK chips (stems = chip values)")
    ax[1].plot(t, rx_clean, color=COLORS["mle"], lw=1.6)
    ax[1].set_title(f"Delayed s(t − τ), τ = {tau:.2f} chips")
    ax[2].plot(t, rx, color=INK2, lw=0.8, label="received r[n]")
    ax[2].plot(t, rx_clean, color=COLORS["mle"], lw=1.6, label="s(t − τ)")
    ec = 10 * np.log10(en0 / sc.n_chips)
    ax[2].set_title(f"Received samples at batch E/N0 = 20 dB (per-chip Ec/N0 = {ec:.1f} dB)")
    ax[2].legend(loc="upper right")
    ax[2].set_xlim(0, 24)
    ax[2].set_xlabel("time (chips)")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def fig_dll_tracking(cfg, hist, crlb_of, path):
    picks = [cfg.snr_db[-6], cfg.snr_db[len(cfg.snr_db) // 2 - 2], cfg.snr_db[3]]
    fig, ax = plt.subplots(1, 3, figsize=(12, 3.6), sharey=True)
    for a, db in zip(ax, picks):
        h = hist[db]
        tchips = np.arange(h.shape[1]) * cfg.dll_block_chips
        for row in h:
            a.plot(tchips, row, color=COLORS["dll"], lw=1.0, alpha=0.7)
        s = np.sqrt(crlb_of(10 ** (db / 10)))
        a.axhline(s, color=COLORS["crlb"], ls="--", lw=1)
        a.axhline(-s, color=COLORS["crlb"], ls="--", lw=1)
        a.set_title(f"E/N0 = {db:.0f} dB")
        a.set_xlabel("time (chips)")
    ax[0].set_ylabel("tracking error θ − τ (chips)")
    ax[0].plot([], [], color=COLORS["dll"], label="DLL error, 12 runs")
    ax[0].plot([], [], color=COLORS["crlb"], ls="--", lw=1, label="± √CRLB (64-chip batch)")
    ax[0].legend(loc="upper right")
    ax[0].set_ylim(-2.5, 2.5)
    fig.suptitle("DLL error after acquisition (t = 0) and during tracking; 12 runs per SNR", color=INK2)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def fig_error_hist(cfg, err, picks, path):
    lim = cfg.prior_width + 0.5
    bins = np.linspace(-lim, lim, 91)
    fig, ax = plt.subplots(3, 3, figsize=(12, 8.5), sharex=True)
    rows = [("mle", "Batch MLE"), ("mmse", "Batch MMSE"), ("dll", "DLL")]
    for col, db in enumerate(picks):
        i = int(np.argmin(np.abs(cfg.snr_db - db)))
        for row, (key, name) in enumerate(rows):
            a = ax[row, col]
            a.hist(np.clip(err[key][i], -lim, lim), bins=bins, color=COLORS[key], edgecolor="#fcfcfb",
                   linewidth=0.5)
            a.set_yscale("log")
            a.set_title(f"{name}, E/N0 = {cfg.snr_db[i]:.0f} dB")
            for x in (-cfg.prior_width, cfg.prior_width):
                a.axvline(x, color=MUTED, ls=":", lw=1)
    for a in ax[-1]:
        a.set_xlabel("error τ̂ − τ (chips)")
    for a in ax[:, 0]:
        a.set_ylabel("count (log)")
    fig.suptitle("Error distributions above, at and below the MLE threshold "
                 f"(dotted = ±Ta, the largest possible batch error; DLL errors beyond ±{lim:g} are "
                 "piled into the edge bins)", color=INK2, fontsize=10)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def fig_main(cfg, res, path):
    fig, ax = plt.subplots(figsize=(9, 6))
    x, xf = cfg.snr_db, res["snr_fine"]
    ax.semilogy(xf, res["crlb_fine"], color=COLORS["crlb"], ls="--", lw=1.5, label="CRLB (64-chip batch)")
    ax.semilogy(xf, res["zzb_fine"], color=COLORS["zzb"], label="ZZB (64-chip batch)")
    ax.semilogy(xf, res["zzb_stream_fine"], color=COLORS["zzb"], ls=":", lw=1.5,
                label="ZZB (DLL's 640-chip data): valid DLL bound")
    lin = xf >= res["thr_dll_gross"]
    ax.semilogy(xf[lin], res["dll_theory_fine"][lin], color=COLORS["dll"], ls="--", lw=1.2,
                label="DLL linear theory (locked region only)")
    ax.semilogy(x, res["mse_mle"], "o-", color=COLORS["mle"], lw=1.5, label="Batch MLE (sim)")
    ax.semilogy(x, res["mse_mmse"], "D-", color=COLORS["mmse"], lw=1.5, ms=5, label="Batch MMSE (sim)")
    ax.semilogy(x, res["mse_dll"], "s-", color=COLORS["dll"], lw=1.5, label="DLL (sim)")
    ax.axhline(prior_variance(cfg), color=COLORS["prior"], ls=":", lw=1.2,
               label="prior variance Ta²/12 (MMSE at zero SNR)")
    ax.axhline(res["mse_mle_noise"], color=COLORS["prior"], ls="-.", lw=1.0,
               label=f"MLE on pure noise ({res['mse_mle_noise']:.2f})")
    ax.set_xlabel("batch E/N0 (dB)   [per-chip Ec/N0 = E/N0 − %.1f dB]" % (10 * np.log10(cfg.n_chips)))
    ax.set_ylabel("MSE of delay estimate (chips²)")
    ax.set_title("Delay-estimation MSE vs SNR: estimators and bounds")
    ax.set_ylim(np.min(res["zzb_stream_fine"]) * 0.5,
                2.0 * max(res["mse_dll"].max(), res["mse_mle"].max()))
    ax.legend(loc="lower left", fontsize=8)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def fig_threshold(cfg, res, path):
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    x = cfg.snr_db
    cr = res["crlb"]
    ax[0].plot(x, 10 * np.log10(res["mse_mle"] / cr), "o-", color=COLORS["mle"], label="Batch MLE")
    ax[0].plot(x, 10 * np.log10(res["mse_mmse"] / cr), "D-", color=COLORS["mmse"], ms=5, label="Batch MMSE")
    ax[0].plot(x, 10 * np.log10(res["mse_dll"] / cr), "s-", color=COLORS["dll"], label="DLL")
    ax[0].plot(res["snr_fine"], 10 * np.log10(res["zzb_fine"] / res["crlb_fine"]),
               color=COLORS["zzb"], label="ZZB")
    ax[0].axhline(10 * np.log10(THRESH_RATIO), color=MUTED, ls=":", lw=1, label="3 dB criterion")
    ax[0].set_xlabel("batch E/N0 (dB)")
    ax[0].set_ylabel("MSE / CRLB (dB)")
    ax[0].set_title("Departure from the 64-chip CRLB")
    ax[0].set_ylim(-6, 20)
    ax[0].legend(loc="upper right")

    floor = 0.5 / cfg.n_trials                     # plot zeros at half a count
    def lp(y):
        return np.maximum(y, floor)
    ax[1].semilogy(x, lp(res["out_mle"]), "o-", color=COLORS["mle"], label="Batch MLE")
    ax[1].semilogy(x, lp(res["out_mmse"]), "D-", color=COLORS["mmse"], ms=5, label="Batch MMSE")
    ax[1].semilogy(x, lp(res["out_dll"]), "s-", color=COLORS["dll"], label="DLL, total")
    ax[1].semilogy(x, lp(res["out_dll_acq"]), "^--", color=COLORS["dll"], lw=1,
                   label=f"DLL, acquired outside pull-in (>{res['pull_in']:.2f} chip)")
    ax[1].semilogy(x, lp(res["out_dll_lock"]), "v:", color=COLORS["dll"], lw=1,
                   label="DLL, acquired inside pull-in, then lost lock")
    ax[1].axhline(0.01, color=MUTED, ls=":", lw=1)
    ax[1].set_xlabel("batch E/N0 (dB)")
    ax[1].set_ylabel("fraction of trials")
    ax[1].set_title(f"Gross-error probability P(|e| > {GROSS} chip); DLL split adds to total",
                    fontsize=10)
    ax[1].set_ylim(floor * 0.8, 1.2)
    ax[1].legend(loc="lower left", fontsize=8)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def fig_dll_time(cfg, checkpoints, err_cp, path):
    picks = [-10.0, 0.0, 4.0, 8.0]
    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    tchips = np.array(checkpoints) * cfg.dll_block_chips
    for color, db in zip(SNR_RAMP, picks):
        i = int(np.argmin(np.abs(cfg.snr_db - db)))
        ax.loglog(tchips, np.mean(err_cp[i] ** 2, axis=1), "s-", color=color,
                  label=f"E/N0 = {cfg.snr_db[i]:.0f} dB")
    ax.set_xlabel("DLL observation time (chips)")
    ax.set_ylabel("DLL MSE (chips²)")
    ax.set_title("Below threshold the DLL MSE grows with observation time (loss of lock)")
    ax.legend(loc="upper left")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    args = ap.parse_args()
    cfg = Config(n_trials=300) if args.quick else Config()
    rng = np.random.default_rng(cfg.seed + 7)
    for sub in ("data", "figures", "tables"):
        (OUT / sub).mkdir(parents=True, exist_ok=True)
    apply_style()

    print(f"Config: {cfg}\n  DLL: K = {cfg.dll_gain:.4f}, B_L = {cfg.dll_bl:.5f}/chip, "
          f"1/(2B_L) = {1 / (2 * cfg.dll_bl):.1f} chips, {cfg.dll_n_updates} updates")
    print("SNR sweep (Ta = %g chips):" % cfg.prior_width)
    sc, mle, dll, err, err_cp, checkpoints, hist = snr_sweep(cfg, rng)
    pull_in = dll.pull_in_range(rng)
    sc_s, e_ratio = stream_scenario(cfg, sc, dll)

    # MLE on pure noise (no signal information): edge pile-up of the argmax
    tau = sc.draw_delays(4000, rng)
    th_noise = mle.estimate(sc.observe(tau, 1e-8, rng))
    mse_mle_noise = float(np.mean((th_noise - tau) ** 2))
    edge_frac = float(np.mean((th_noise <= mle.grid[0]) | (th_noise >= mle.grid[-1])))

    en0 = 10 ** (cfg.snr_db / 10)
    snr_fine = np.linspace(cfg.snr_db[0], cfg.snr_db[-1], 400)
    en0_f = 10 ** (snr_fine / 10)
    acq_out = np.abs(err["acq"]) > pull_in
    gross_dll = np.abs(err["dll"]) > GROSS
    res = dict(
        snr_db=cfg.snr_db, snr_fine=snr_fine, pull_in=pull_in,
        mse_mle_noise=mse_mle_noise, edge_frac_noise=edge_frac, stream_energy_ratio=e_ratio,
        crlb=crlb(sc, en0), zzb=zzb(sc, en0), crlb_fine=crlb(sc, en0_f), zzb_fine=zzb(sc, en0_f),
        crlb_stream=crlb(sc_s, en0 * e_ratio), zzb_stream=zzb(sc_s, en0 * e_ratio),
        zzb_stream_fine=zzb(sc_s, en0_f * e_ratio),
        dll_theory=dll.linear_theory_var(sc.sigma2(en0)),
        dll_theory_fine=dll.linear_theory_var(sc.sigma2(en0_f)),
        out_dll=np.mean(gross_dll, 1),
        out_dll_acq=np.mean(gross_dll & acq_out, 1),         # partition of out_dll
        out_dll_lock=np.mean(gross_dll & ~acq_out, 1),
    )
    for k in ("mle", "mmse", "dll"):
        res[f"mse_{k}"] = np.mean(err[k] ** 2, 1)
        res[f"bias_{k}"] = np.mean(err[k], 1)
        res[f"out_{k}"] = np.mean(np.abs(err[k]) > GROSS, 1)
    thr = {k: bootstrap_thresholds(cfg.snr_db, err[k], res["crlb"], n_boot=200)
           for k in ("mle", "mmse", "dll")}
    for k in thr:
        res[f"thr_{k}_mse"], res[f"thr_{k}_gross"] = thr[k]["mse"][0], thr[k]["gross"][0]
    res["thr_zzb"] = crossing_snr(snr_fine, res["zzb_fine"] / res["crlb_fine"], THRESH_RATIO)
    np.savez(OUT / "data" / "snr_sweep.npz", err_mle=err["mle"], err_mmse=err["mmse"],
             err_dll=err["dll"], err_acq=err["acq"], err_dll_cp=err_cp,
             checkpoints=np.array(checkpoints), **res)

    hdr = ("snr_db,ecn0_db,crlb,zzb,crlb_dll_data,zzb_dll_data,mse_mle,mse_mmse,mse_dll,dll_theory,"
           "bias_mle,bias_mmse,bias_dll,gross_mle,gross_mmse,gross_dll,gross_dll_acq_outside_pullin,"
           "gross_dll_lost_lock")
    rows = np.column_stack([cfg.snr_db, cfg.snr_db - 10 * np.log10(cfg.n_chips), res["crlb"], res["zzb"],
                            res["crlb_stream"], res["zzb_stream"], res["mse_mle"], res["mse_mmse"],
                            res["mse_dll"], res["dll_theory"], res["bias_mle"], res["bias_mmse"],
                            res["bias_dll"], res["out_mle"], res["out_mmse"], res["out_dll"],
                            res["out_dll_acq"], res["out_dll_lock"]])
    np.savetxt(OUT / "tables" / "summary.csv", rows, delimiter=",", header=hdr, comments="", fmt="%.6g")

    print(f"\nDLL pull-in range: {pull_in:.2f} chips; DLL data = {e_ratio:.2f} x batch energy")
    print(f"MLE on pure noise: MSE {mse_mle_noise:.3f} (Ta^2/6 = {2 * prior_variance(cfg):.3f}), "
          f"{edge_frac:.1%} of estimates on a window edge")
    print("Indicative thresholds (2 dB grid, bootstrap 95% CI), batch E/N0:")
    for k in ("mle", "mmse", "dll"):
        m, g = thr[k]["mse"], thr[k]["gross"]
        print(f"  {k.upper():4s}: MSE>2xCRLB {m[0]:5.2f} [{m[1]:5.2f}, {m[2]:5.2f}]   "
              f"P(gross)>1% {g[0]:5.2f} [{g[1]:5.2f}, {g[2]:5.2f}]")
    print(f"  ZZB : {res['thr_zzb']:.2f} dB")
    print(f"  high-SNR MSE/CRLB: MLE {res['mse_mle'][-1] / res['crlb'][-1]:.3f}, "
          f"MMSE {res['mse_mmse'][-1] / res['crlb'][-1]:.3f}, DLL {res['mse_dll'][-1] / res['crlb'][-1]:.3f}, "
          f"DLL theory {res['dll_theory'][-1] / res['crlb'][-1]:.3f}")
    print(f"  DLL MSE / ZZB(DLL data): min {np.min(res['mse_dll'] / res['zzb_stream']):.2f} (must be >= 1)")
    i_noise = int(np.argmin(np.abs(cfg.snr_db + 10)))
    print("  DLL MSE vs time at -10 dB:", ", ".join(
        f"{c * cfg.dll_block_chips} chips: {np.mean(err_cp[i_noise, k] ** 2):.1f}" for k, c in enumerate(checkpoints)))

    figs = OUT / "figures"
    fig_signal(sc, rng, figs / "fig1_signal.png")
    fig_dll_tracking(cfg, hist, lambda e: crlb(sc, e), figs / "fig2_dll_tracking.png")
    t = res["thr_mle_mse"]
    picks = [t + 8, t, t - 8] if np.isfinite(t) else [cfg.snr_db[-3], cfg.snr_db[10], cfg.snr_db[2]]
    fig_error_hist(cfg, err, picks, figs / "fig3_error_hist.png")
    fig_main(cfg, res, figs / "fig4_mse_vs_snr.png")
    fig_threshold(cfg, res, figs / "fig5_threshold.png")
    fig_dll_time(cfg, checkpoints, err_cp, figs / "fig7_dll_mse_vs_time.png")
    print(f"\nSaved figures to {figs}/, data to {OUT / 'data'}/, table to {OUT / 'tables'}/")


if __name__ == "__main__":
    main()
