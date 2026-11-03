import sys; sys.path.insert(0, r'C:\Users\Nilavra Ghosh\Documents\7THSEM\SDET\project\src')
import numpy as np, importlib.util
from tde import *
spec = importlib.util.spec_from_file_location("rb", r'C:\Users\Nilavra Ghosh\Documents\7THSEM\SDET\project\run_baseline.py'); rb = importlib.util.module_from_spec(spec); spec.loader.exec_module(rb)
d = np.load(r'C:\Users\Nilavra Ghosh\Documents\7THSEM\SDET\project\results\data\snr_sweep.npz')
snr, cr = d['snr_db'], d['crlb']; em, ed = d['err_mle'], d['err_dll']
rng = np.random.default_rng(0); tm, td, gap = [], [], []
for b in range(1000):
    i = rng.integers(0, em.shape[1], em.shape[1]); k = rng.integers(0, ed.shape[1], ed.shape[1])
    a = rb.threshold_snr(snr, np.mean(em[:, i]**2, 1)/cr); c = rb.threshold_snr(snr, np.mean(ed[:, k]**2, 1)/cr)
    tm.append(a); td.append(c); gap.append(c-a)
q = lambda v: np.round(np.nanpercentile(v,[2.5,50,97.5]),2)
print("MLE thr 95% CI", q(tm), " DLL thr", q(td), " DLL-MLE gap", q(gap))
# prior sweep repeat with different rng seeds -> spread of gap
cfg0 = Config()
for seed in [11, 12, 13, 14]:
    rng = np.random.default_rng(seed); out = rb.prior_sweep(cfg0, [4.0, 16.0], 1000, rng)
    snr_f = np.linspace(-10, 30, 400); res=[]
    for ta, dd in out.items():
        z = zzb(dd['sc'], 10**(snr_f/10)); c = crlb(dd['sc'], 10**(snr_f/10))
        res.append(rb.threshold_snr(cfg0.snr_db, dd['mse']/crlb(dd['sc'], 10**(cfg0.snr_db/10))) - rb.threshold_snr(snr_f, z/c))
    print(f"seed {seed}: MLE-ZZB gap Ta=4 {res[0]:.2f} dB, Ta=16 {res[1]:.2f} dB")
