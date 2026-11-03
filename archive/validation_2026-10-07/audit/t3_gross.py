import sys; sys.path.insert(0, r'C:\Users\Nilavra Ghosh\Documents\7THSEM\SDET\project\src')
import numpy as np
from tde import *
d = np.load(r'C:\Users\Nilavra Ghosh\Documents\7THSEM\SDET\project\results\data\snr_sweep.npz')
snr = d['snr_db']; em, ed, ea = d['err_mle'], d['err_dll'], d['err_acq']
print("SNR | MLE gross | MLE |e|>0.75 | acq|e|>0.5 | acq|e|>0.75 | DLL final gross | acq>0.75 & final ok | acq<=0.5 & final gross | acq in(0.5,0.75] final ok")
for i in range(len(snr)):
    if snr[i] < 0 or snr[i] > 14: continue
    a = np.abs(ea[i]); f = np.abs(ed[i]); m = np.abs(em[i])
    print(f"{snr[i]:5.0f} {np.mean(m>0.5):.4f} {np.mean(m>0.75):.4f} {np.mean(a>0.5):.4f} {np.mean(a>0.75):.4f} {np.mean(f>0.5):.4f} "
          f"{np.mean((a>0.75)&(f<=0.5)):.4f} {np.mean((a<=0.5)&(f>0.5)):.4f} {np.mean((a>0.5)&(a<=0.75)&(f<=0.5)):.4f}")
# acquisition error at high SNR: max |acq err| should be <= 0.25
for i in [len(snr)-1, len(snr)-5]:
    print(snr[i], "dB max|acq err| =", np.max(np.abs(ea[i])), " max |dll err|", np.max(np.abs(ed[i])))
# check acq partition: out_dll_acq + out_dll_lock vs out_dll
print("acq+lock - total:", np.round(d['out_dll_acq'] + d['out_dll_lock'] - d['out_dll'],4))
# low SNR DLL final error distribution
i = 0
print("-10 dB: DLL |err| quantiles", np.round(np.quantile(np.abs(ed[i]),[.5,.9,.99,1]),2), " frac |err|>9:", np.mean(np.abs(ed[i])>9))
print("MLE bias at -10..0 dB:", np.round(d['bias_mle'][:6],3))
# variance vs MSE: bias contribution at high SNR
print("bias^2/mse MLE high SNR:", np.round(d['bias_mle'][-8:]**2/d['mse_mle'][-8:],4))
print("bias^2/mse DLL high SNR:", np.round(d['bias_dll'][-8:]**2/d['mse_dll'][-8:],4))
