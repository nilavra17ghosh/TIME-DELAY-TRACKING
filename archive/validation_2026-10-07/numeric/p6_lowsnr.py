"""Part 6: low-SNR MLE MSE vs Ta^2/6; distribution of estimates (uses p3_out.npz)."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

d = np.load("p3_out.npz")
LO, HI = float(d["LO"]), float(d["HI"])
TA = HI - LO
rng = np.random.default_rng(3)
print(f"Ta^2/6 = {TA ** 2 / 6:.4f}, Ta^2/12 = {TA ** 2 / 12:.4f}")
print(" case    n      MSE    95%CI           Var(est)  Cov(est,tau)  bias^2  P(edge)  P(|e|<0.5)  hist of est position (8 bins, fraction)")
fig, ax = plt.subplots(1, 4, figsize=(14, 3.4), sharey=True)
for a, key in zip(ax, ["noise", "-40", "-20", "-10"]):
    tau, est = d[f"tau_{key}"], d[f"est_{key}"]
    e = est - tau
    e2 = e ** 2
    bt = np.array([e2[rng.integers(0, e2.size, e2.size)].mean() for _ in range(2000)])
    lo, hi = np.percentile(bt, [2.5, 97.5])
    edge = np.mean((np.abs(est - LO) < 1e-6) | (np.abs(est - HI) < 1e-6))
    pos = (est - LO) / TA
    hist = np.histogram(pos, bins=8, range=(0, 1))[0] / est.size
    cv = np.cov(est, tau)[0, 1]
    print(f"{key:>6s} {e.size:6d} {e2.mean():.4f} [{lo:.3f},{hi:.3f}]  {est.var():.4f}   {cv:+.4f}      {e.mean() ** 2:.4f}  {edge:.3f}   {np.mean(np.abs(e) < 0.5):.3f}     "
          + " ".join(f"{x:.3f}" for x in hist))
    a.hist(pos * TA - TA / 2, bins=64, range=(-TA / 2, TA / 2), density=True, color="#2a78d6")
    a.axhline(1 / TA, color="k", ls="--", lw=1)
    a.set_title("pure noise" if key == "noise" else f"E/N0 = {key} dB")
    a.set_xlabel("MLE estimate - prior centre (chips)")
ax[0].set_ylabel("density (dashed = uniform)")
fig.tight_layout()
fig.savefig("p6_lowsnr_hist.png", dpi=130)

# compare with project summary at -10 dB
summ = np.genfromtxt(r"C:\Users\Nilavra Ghosh\Documents\7THSEM\SDET\project\results\tables\summary.csv", delimiter=",", names=True)
print("project MLE MSE at -10/-8/-6 dB:", summ["mse_mle"][:3])
