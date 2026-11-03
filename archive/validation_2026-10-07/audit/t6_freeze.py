import sys; sys.path.insert(0, r'C:\Users\Nilavra Ghosh\Documents\7THSEM\SDET\project\src')
import numpy as np
from tde import *
cfg = Config(); sc = BatchScenario(cfg, np.random.default_rng(cfg.seed)); dll = EarlyLateDLL(cfg, sc.wave.chips)
rng = np.random.default_rng(5)
tau = dll.tau_c + rng.uniform(-2,2,400)
th, th0, H = dll.run(tau, sc.sigma2(0.1), rng, keep_history=400)
dev = np.abs(th - dll.tau_c)
far = dev > 13
step_last = np.abs(np.diff(H[:, -41:], axis=1)).mean(1)
print("-10 dB: frac |theta-tau_c|>13:", far.mean(), " max dev", dev.max().round(2))
print("mean |step| over last 40 updates: far trials %.2e, near trials %.2e" % (step_last[far].mean(), step_last[~far].mean()))
# MSE growth with time (should grow linearly for a true random walk)
for n in [40, 80, 120, 160]: print(f"MSE after {n} updates: {np.mean(H[:, n]**2):.1f}")
# MLE edge zones at 20 dB
mle = BatchMLE(sc); e = 10**2
for lo_off in [0.0, 1.0]:
    t = sc.tau_c - 2 + lo_off + rng.uniform(0, 0.05, 4000)
    err = mle.estimate(sc.observe(t, e, rng)) - t
    print(f"MLE 20 dB tau within [{lo_off},{lo_off+0.05}] of prior low edge: MSE/CRLB {np.mean(err**2)/crlb(sc,e):.3f}, bias {np.mean(err):+.4f}")
# rho phase invariance at non-integer base
from tde.bounds import autocorrelation
h = np.linspace(0,4,401); s0 = sc.wave(sc.t - sc.tau_c - 0.37); sh = sc.wave(sc.t[None,:]-sc.tau_c-0.37-h[:,None])
print("rho phase variation:", np.max(np.abs(sh@s0/np.sum(s0**2) - autocorrelation(sc,h))))
# naive independent-update DLL theory vs exact (doc: off by 10x)
g2 = np.mean(np.sum(dll.el_diff**2,1)); s2 = sc.sigma2(1e3)
naive = s2*g2/dll.k_bar**2 * cfg.dll_gain/(2-cfg.dll_gain)
print("naive/exact linear theory:", naive/dll.linear_theory_var(s2))
