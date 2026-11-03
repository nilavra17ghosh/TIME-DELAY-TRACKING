"""Summary figure: MLE sim (project code, fine grid, 5 pilots) vs ZZB vs CRLB vs MIE predictions."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

p4 = np.load("p4_out.npz")
p5 = np.load("p5_out.npz")
summ = np.genfromtxt(r"C:\Users\Nilavra Ghosh\Documents\7THSEM\SDET\project\results\tables\summary.csv", delimiter=",", names=True)
fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))
for a, ta in zip(ax, (4, 16)):
    snr = p5[f"ta{ta}_snr"]
    b2 = 2.4481
    crl = 1 / (2 * 10 ** (snr / 10) * b2)
    a.semilogy(snr, crl, "k--", lw=1.2, label="CRLB")
    a.semilogy(snr, p5[f"ta{ta}_rice"][:, 0], color="#c43", lw=1.5, label="MIE (Rice/Poisson)")
    a.semilogy(snr, p5[f"ta{ta}_u1"][:, 0], color="#c43", ls=":", lw=1.5, label="MIE (union, Δ=1.335)")
    g = p4[f"ta{ta}_s2026_grid"]
    a.semilogy(g, p4[f"ta{ta}_s2026_zzb"], color="#1baf7a", lw=2, label="ZZB (project code)")
    a.semilogy(g, p4[f"ta{ta}_s2026_mse"], "o", color="#2a78d6", label="MLE sim, pilot 2026 (project code)")
    if ta == 4:
        a.semilogy(summ["snr_db"], summ["mse_mle"], "x", color="#555", label="summary.csv mse_mle")
        try:
            p3 = np.load("p3_out.npz")
            for db in [4, 6, 8, 10, 12, 20, 30]:
                e2 = (p3[f"est_{db}"] - p3[f"tau_{db}"]) ** 2
                a.semilogy([db], [e2.mean()], "s", mfc="none", color="#eb6834",
                           label="independent MLE (this check)" if db == 4 else None)
        except FileNotFoundError:
            pass
    a.set_xlim(snr[0], snr[-1])
    a.set_ylim(1e-3 if ta == 4 else 5e-4, 30 if ta == 16 else 5)
    a.set_title(f"Ta = {ta} chips")
    a.set_xlabel("batch E/N0 (dB)")
    a.grid(alpha=0.3)
ax[0].set_ylabel("MSE (chips²)")
ax[0].legend(fontsize=8, loc="lower left")
fig.tight_layout()
fig.savefig("verify_mse_mie.png", dpi=130)
print("saved verify_mse_mie.png")
