import sys; sys.path.insert(0, r'C:\Users\Nilavra Ghosh\Documents\7THSEM\SDET\project\src')
import numpy as np
from tde import *
cfg = Config(); sc = BatchScenario(cfg, np.random.default_rng(cfg.seed))
dll = EarlyLateDLL(cfg, sc.wave.chips)
# template energy & beta^2 vs fractional delay
mle = BatchMLE(sc)
he = 2*mle.half_energy
print("template energy over grid: min %.6f max %.6f rel spread %.2e" % (he.min(), he.max(), he.max()/he.min()-1))
b2 = []
for dt in np.linspace(-2,2,33):
    t = sc.t - sc.tau_c - dt
    b2.append(np.sum(sc.wave.derivative(t)**2)/np.sum(sc.wave(t)**2))
b2=np.array(b2); print("beta^2 vs tau: min %.5f max %.5f (at tau_c %.5f)" % (b2.min(), b2.max(), rms_bandwidth_sq(sc)))
# derivative accuracy: compare h=1e-5 with h=1e-3 Richardson
t = sc.t - sc.tau_c
d1 = sc.wave.derivative(t, 1e-5); d2 = sc.wave.derivative(t, 1e-3); d3 = sc.wave.derivative(t,5e-4)
rich = (4*d3 - d2)/3
print("deriv rel err (h=1e-5 vs Richardson):", np.linalg.norm(d1-rich)/np.linalg.norm(rich))
# periodic stream: per-period energy and beta^2 in the middle of the stream
tt = dll.t; s = dll.wave(tt - dll.tau_c); ds = dll.wave.derivative(tt - dll.tau_c)
p0 = int((dll.tau_c + 5*64)*4); p1 = p0 + 64*4
print("periodic per-period energy %.3f vs batch energy %.3f; per-period beta^2 %.4f vs batch %.4f" %
      (np.sum(s[p0:p1]**2), sc.energy, np.sum(ds[p0:p1]**2)/np.sum(s[p0:p1]**2), rms_bandwidth_sq(sc)))
# CRLB-equivalent: Fisher info of 64 chips periodic vs batch
print("Fisher ratio periodic-period/batch: %.4f" % (np.sum(ds[p0:p1]**2)/np.sum(sc.wave.derivative(t)**2)))
# Cauchy-Schwarz claim
G = dll.equivalent_weighting()
print("<s',G> (stream) = %.5f" % np.dot(ds, G))
sp2_stream = np.sum(ds[:]**2)
eff_stream = np.sum(G**2)*sp2_stream/np.dot(ds,G)**2
print("||G||^2 ||s'_stream||^2/<s',G>^2 = %.3f   (doc says this equals the 1.07 high-SNR loss)" % eff_stream)
print("sigma^2||G||^2 / CRLB_batch = %.4f" % (np.sum(G**2)*sc.sigma2(1e3)/crlb(sc,1e3)))
print("stream CRLB (all %d samples) / batch CRLB = %.4f" % (tt.size, (sc.sigma2(1e3)/sp2_stream)/crlb(sc,1e3)))
# which part of the stream does the weighting actually use: distribution of G energy per period
gp = [np.sum(G[int((dll.tau_c+64*p)*4):int((dll.tau_c+64*(p+1))*4)]**2) for p in range(12)]
print("fraction of ||G||^2 per period:", np.round(np.array(gp)/np.sum(G**2),3))
# narrower loop -> beats batch CRLB?
for f in [2,4]:
    c2 = Config(dll_block_chips=4); 
    class C(Config): pass
    c2 = Config(); c2.n_chips_eff=None
    import dataclasses
    c3 = dataclasses.replace(cfg, dll_settle_factor=10*f)
    # emulate n_eff = 16*f by subclass property override
    class CF(Config):
        @property
        def dll_n_eff(self): return 16.0*f
        @property
        def dll_n_updates(self): return int(160*f)
    cf = CF()
    d = EarlyLateDLL(cf, sc.wave.chips)
    print(f"N_eff={16*f}: linear theory var / batch CRLB = {d.linear_theory_var(sc.sigma2(1e3))/crlb(sc,1e3):.3f}")
