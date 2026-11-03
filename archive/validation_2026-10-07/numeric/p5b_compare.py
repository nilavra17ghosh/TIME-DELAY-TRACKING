"""Compare MIE predictions (p5_out.npz) with simulated MLE (p4_out.npz, pilot 2026)."""
import numpy as np
p5 = np.load("p5_out.npz"); p4 = np.load("p4_out.npz")
for ta in (4, 16):
    snr = p5[f"ta{ta}_snr"]; k = f"ta{ta}_s2026"
    print(f"Ta={ta}:  dB   sim MSE/CRLB  Rice  union1.335 | sim Pgross  Rice Pe  union Pe")
    for db, m, c, g in zip(p4[k + "_grid"], p4[k + "_mse"], p4[k + "_crlb"], p4[k + "_gross"]):
        i = np.argmin(np.abs(snr - db))
        print(f"      {db:5.0f}  {m / c:9.2f}   {p5[f'ta{ta}_rice'][i, 0] / c:6.2f} {p5[f'ta{ta}_u1'][i, 0] / c:8.2f}   | {g:.4f}   {p5[f'ta{ta}_rice'][i, 1]:.4f}  {p5[f'ta{ta}_u1'][i, 1]:.4f}")
