import sys; sys.path.insert(0, r'C:\Users\Nilavra Ghosh\Documents\7THSEM\SDET\project\src'); sys.path.insert(0, r'C:\Users\Nilavra Ghosh\Documents\7THSEM\SDET\project')
import numpy as np
from tde import *
from tde.bounds import autocorrelation, qfunc
import importlib.util
spec = importlib.util.spec_from_file_location("rb", r'C:\Users\Nilavra Ghosh\Documents\7THSEM\SDET\project\run_baseline.py'); rb = importlib.util.module_from_spec(spec); spec.loader.exec_module(rb)
cfg = Config(); sc = BatchScenario(cfg, np.random.default_rng(cfg.seed)); dll = EarlyLateDLL(cfg, sc.wave.chips)
# (a) window geometry: noise gain ||e-l||^2 and S-curve vs theta offset for a mid block
j = 80; idx = dll.block_idx[j]; t = dll.t[idx]
print("block window (chips rel. block start):", t[0]-(dll.tau_c+4*j), t[-1]-(dll.tau_c+4*j))
for x in [0, 2, 5, 8, 10, 12, 15, 20]:
    e,l = dll._el(j, t, dll.tau_c + x); print(f"theta-tau_c={x:3d}: ||e-l||^2 in window = {np.sum((e-l)**2):.4f}")
# full S-curve averaged over blocks (noise-free), error eps
eps = np.linspace(-2,2,81); S=[]
for ep in eps:
    v=[]
    for j in range(16, 32):
        idx = dll.block_idx[j]; t = dll.t[idx]; rx = dll.wave(t - dll.tau_c)
        e,l = dll._el(j,t,dll.tau_c+ep); v.append((np.sum(rx*(e-l))-dll.offset[j])/dll.k_bar)
    S.append(np.mean(v))
S=np.array(S); pos = eps[(S>0)]
print("S-curve: D(eps) for eps in", list(zip(eps[40:81:4].round(2), S[40:81:4].round(3))))
print("first eps>0 where D<=0:", eps[(eps>0.05)&(S<=0)][:1], " peak at", eps[np.argmax(S)])
# (b) threshold_snr synthetic tests
x = np.arange(0, 31, 2.0)
r = 1 + 10**((10 - x)/5)    # ratio=2 at x=10 exactly
print("thr synthetic (expect 10):", rb.threshold_snr(x, r))
r2 = r.copy(); r2[0] = 1.5  # non-monotone at low end
print("thr nonmonotone low end (expect 10):", rb.threshold_snr(x, r2))
r3 = np.full_like(x, 2.0); print("ratio exactly 2 everywhere:", rb.threshold_snr(x, r3))
r4 = r.copy(); r4[-1] = 2.5  # noisy high-SNR point above 2
print("high-SNR noisy point >2 (returns crossing at top):", rb.threshold_snr(x, r4))
# (c) ZZB independent recomputation: rho at different base, finer h, brute-force valley filling
h = np.linspace(0, 4, 4001)
rho_base = autocorrelation(sc, h)
s0 = sc.wave(sc.t - sc.tau_c); sh = sc.wave(sc.t[None,:] - sc.tau_c - h[:,None]); rho_c = sh@s0/np.sum(s0**2)
print("rho(base) vs rho(tau_c) max diff:", np.max(np.abs(rho_base-rho_c)))
def zzb_bf(en0):
    f = (1-h/4)*qfunc(np.sqrt(en0*np.clip(1-rho_c,0,None)))
    v = np.array([f[i:].max() for i in range(len(f))])
    return np.trapezoid(h*v, h), np.trapezoid(h*f, h)
for db in [0, 6, 10, 20]:
    zb, znv = zzb_bf(10**(db/10)); print(f"{db} dB ZZB code {zzb(sc,10**(db/10))[0]:.5g}, brute {zb:.5g}, no-valley-fill {znv:.5g}")
# beta^2 range over pilots
b = [rms_bandwidth_sq(BatchScenario(cfg, np.random.default_rng(s))) for s in range(400)]
print("beta^2 over 400 pilots: min %.3f max %.3f mean %.3f" % (min(b), max(b), np.mean(b)))
