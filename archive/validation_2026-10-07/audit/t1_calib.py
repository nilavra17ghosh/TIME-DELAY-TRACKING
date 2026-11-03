import sys; sys.path.insert(0, r'C:\Users\Nilavra Ghosh\Documents\7THSEM\SDET\project\src')
import numpy as np, copy
from tde import *
cfg = Config(); sc = BatchScenario(cfg, np.random.default_rng(cfg.seed))
dll = EarlyLateDLL(cfg, sc.wave.chips)
print("tau_c", dll.tau_c, "k_bar", dll.k_bar)
def consts(dll, tref, h=1e-4):
    slope, offset = [], []
    for j, idx in enumerate(dll.block_idx):
        t = dll.t[idx]; rx = dll.wave(t - tref)
        e_p,l_p = dll._el(j,t,tref+h); e_m,l_m = dll._el(j,t,tref-h)
        slope.append((np.sum(rx*(e_p-l_p))-np.sum(rx*(e_m-l_m)))/(2*h))
        e0,l0 = dll._el(j,t,tref); offset.append(np.sum(rx*(e0-l0)))
    return np.array(slope), np.array(offset)
for dt in [0.37, 0.125, -1.63, 1.9]:
    s2, o2 = consts(dll, dll.tau_c+dt)
    print(f"dt={dt}: max|db|={np.max(np.abs(o2-dll.offset)):.3e} (rms b={np.sqrt(np.mean(dll.offset**2)):.3e}, max|b|={np.max(np.abs(dll.offset)):.3f}); "
          f"max|dk|/kbar={np.max(np.abs(s2-dll.slope))/dll.k_bar:.3e}; kbar change {np.mean(s2)/dll.k_bar-1:.2e}")
# sum of b over period
print("sum b_j per period (16 blocks):", [round(float(np.sum(dll.offset[p*16:(p+1)*16])),4) for p in range(10)])
print("slope min/max/kbar", dll.slope.min(), dll.slope.max(), dll.k_bar)
# impact on MSE at high SNR: swap calibration to wrong phase
rng = np.random.default_rng(1)
for db in [25, 30]:
    s2v = sc.sigma2(10**(db/10))
    tau = dll.tau_c + rng.uniform(-2,2,2000)
    seed = 123
    e_ok = dll.run(tau, s2v, np.random.default_rng(seed))[0]-tau
    d2 = copy.copy(dll); sl, of = consts(dll, dll.tau_c+0.37); d2.offset = of; d2.slope=sl; d2.k_bar=float(np.mean(sl))
    e_bad = d2.run(tau, s2v, np.random.default_rng(seed))[0]-tau
    # no calibration at all
    d3 = copy.copy(dll); d3.offset = np.zeros_like(dll.offset)
    e_none = d3.run(tau, s2v, np.random.default_rng(seed))[0]-tau
    cr = crlb(sc, 10**(db/10))
    print(f"{db} dB: MSE/CRLB calib@tau_c {np.mean(e_ok**2)/cr:.4f}, calib@tau_c+0.37 {np.mean(e_bad**2)/cr:.4f}, no b_j {np.mean(e_none**2)/cr:.4f}")
# noise-free floor with wrong-phase calibration
tau0 = dll.tau_c + np.linspace(-2,2,201)
print("noise-free floor wrong-phase calib:", np.mean((d2.run(tau0,0.0,rng)[0]-tau0)**2), "vs crlb30", crlb(sc,1e3))
# noise-free floor at many random taus (not just 201 grid)
tr = dll.tau_c + rng.uniform(-2,2,2000)
e = dll.run(tr,0.0,rng)[0]-tr
print("noise-free floor 2000 random taus:", np.mean(e**2), np.max(np.abs(e)))
