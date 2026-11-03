"""Part 5: Method-of-interval-errors (MIE) prediction of MLE MSE through the threshold.

Normalised log-likelihood: J(theta)/(sigma ||s||) = d rho(theta - tau) + X(theta),
d = sqrt(2 E/N0), X unit-variance stationary Gaussian process with correlation rho(h).
An outlier occurs when the max of J outside the main lobe (|theta - tau| > h0) beats the
main-lobe peak (~ d + X(tau)).

Model R (Rice/Poisson clumping, sidelobe mean included through a locally-varying level):
  P(no outlier | X(tau)=x) = prod_segments Phi(u(start)) * exp(-int_O nu(u(theta)) dtheta),
  u(theta) = d (1 - rho(theta - tau)) + x,  nu(u) = sqrt(lambda2)/(2 pi) exp(-u^2/2),
  lambda2 = beta^2 (Rice upcrossing rate).  X(theta) taken independent of X(tau) outside the lobe.
Model U (union bound over pairwise ZZB-type error probabilities, candidates spaced Delta):
  P_e(tau) = min(1, sum_k Q( sqrt(E/N0 (1 - rho(h_k))) )).
MSE(tau) = (1 - P_e) CRLB + P_e E[(theta - tau)^2 | outlier], outlier location ~ U(O(tau))
(Model R: weighted by the local crossing density).  Averaged over tau ~ U[0, Ta].
"""
import numpy as np
from scipy.stats import norm
from vcommon import pilot, Q, waveform

chips = pilot()


def make_rho(fs=16, span=6):
    t = np.arange(-span - 2, 64 + span + 2, 1.0 / fs)
    s = waveform(t, chips, span=span)
    n = 1 << 15
    S = np.fft.rfft(s, n)
    f = np.fft.rfftfreq(n, 1.0 / fs)
    w = np.abs(S) ** 2
    w[1:] *= 2
    keep = f < 4
    f, w = f[keep], w[keep]
    b2 = np.sum((2 * np.pi * f) ** 2 * w) / np.sum(w)

    def rho(h):
        h = np.atleast_1d(np.asarray(h, float))
        out = np.concatenate([np.cos(2 * np.pi * np.outer(h.ravel()[i:i + 400], f)) @ w
                              for i in range(0, h.size, 400)]) / np.sum(w)
        return out.reshape(h.shape)
    return rho, b2


rho, B2 = make_rho()
LAM = np.sqrt(B2)


def first_zero():
    h = np.linspace(0, 3, 30001)
    return h[np.argmax(rho(h) < 0)]


H0 = first_zero()


def mie_rice(en0_list, ta, n_tau=81, n_th=None, h0=H0):
    n_th = n_th or int(ta * 64) + 1
    th = np.linspace(0, ta, n_th)
    dth = th[1] - th[0]
    taus = np.linspace(0, ta, n_tau)
    xs = np.linspace(-6, 6, 241)
    wx = norm.pdf(xs) * (xs[1] - xs[0])
    R = rho(th[:, None] - taus[None, :])       # (n_th, n_tau)  depends on |h|
    out = []
    for en0 in np.atleast_1d(en0_list):
        d = np.sqrt(2 * en0)
        crl = 1 / (2 * en0 * B2)
        mse_t = np.empty(n_tau)
        pe_t = np.empty(n_tau)
        for j, tau in enumerate(taus):
            outside = np.abs(th - tau) > h0
            if not outside.any():
                mse_t[j], pe_t[j] = crl, 0.0
                continue
            u = d * (1 - R[outside, j])[None, :] + xs[:, None]      # (nx, n_out)
            nu = LAM / (2 * np.pi) * np.exp(-0.5 * u ** 2)
            lam_tot = nu.sum(1) * dth
            # boundary terms: window end points that lie outside the lobe
            ends = [i for i in (0, n_th - 1) if outside[i]]
            pb = np.ones(xs.size)
            for i in ends:
                pb *= norm.cdf(d * (1 - R[i, j]) + xs)
            p_ok = np.sum(wx * pb * np.exp(-lam_tot))
            pe = 1 - p_ok
            # location of the outlier: crossing density, plus end-point masses
            wloc = np.sum(wx[:, None] * nu, 0) * dth
            e2 = (th[outside] - tau) ** 2
            num = np.sum(wloc * e2)
            den = np.sum(wloc)
            for i in ends:
                pm = np.sum(wx * norm.sf(d * (1 - R[i, j]) + xs))
                num += pm * (th[i] - tau) ** 2
                den += pm
            cond = num / den
            pe_t[j] = pe
            mse_t[j] = (1 - pe) * crl + pe * cond
        out.append((mse_t.mean(), pe_t.mean()))
    return np.array(out)


def mie_union(en0_list, ta, delta=H0, n_tau=81, h0=H0):
    taus = np.linspace(0, ta, n_tau)
    out = []
    for en0 in np.atleast_1d(en0_list):
        crl = 1 / (2 * en0 * B2)
        m, pes = [], []
        for tau in taus:
            cand = np.concatenate([np.arange(tau + h0, ta + 1e-9, delta), np.arange(tau - h0, -1e-9, -delta)])
            if cand.size == 0:
                m.append(crl); pes.append(0); continue
            p = Q(np.sqrt(en0 * np.clip(1 - rho(cand - tau), 0, None)))
            pe = min(1.0, p.sum())
            cond = np.sum(p * (cand - tau) ** 2) / p.sum()
            m.append((1 - pe) * crl + pe * cond); pes.append(pe)
        out.append((np.mean(m), np.mean(pes)))
    return np.array(out)


def threshold_snr(snr_db, ratio, level=2.0):
    lr, lv = np.log10(ratio), np.log10(level)
    for i in range(len(snr_db) - 1, 0, -1):
        if lr[i] <= lv < lr[i - 1]:
            f = (lv - lr[i]) / (lr[i - 1] - lr[i])
            return snr_db[i] + f * (snr_db[i - 1] - snr_db[i])
    return np.nan


if __name__ == "__main__":
    print(f"beta^2 = {B2:.4f}, main-lobe first zero h0 = {H0:.3f} chips")
    summ = np.genfromtxt(r"C:\Users\Nilavra Ghosh\Documents\7THSEM\SDET\project\results\tables\summary.csv",
                         delimiter=",", names=True)
    try:
        p4 = np.load("p4_out.npz")
    except FileNotFoundError:
        p4 = None
    res = {}
    for ta in (4.0, 16.0):
        snr_f = np.arange(-4.0, 24.01, 0.5)
        en = 10 ** (snr_f / 10)
        crl = 1 / (2 * en * B2)
        r_rice = mie_rice(en, ta)
        r_u1 = mie_union(en, ta, delta=H0)
        r_u2 = mie_union(en, ta, delta=1.0)
        r_rice_h = mie_rice(en, ta, h0=1.0)
        res[ta] = (snr_f, r_rice, r_u1, r_u2, r_rice_h)
        print(f"\nTa = {ta:g}: predicted threshold (MSE > 2 CRLB)")
        for name, r in [("MIE-Rice h0=1.335", r_rice), ("MIE-Rice h0=1.0", r_rice_h),
                        ("MIE-union Delta=1.335", r_u1), ("MIE-union Delta=1.0", r_u2)]:
            print(f"   {name:24s}: {threshold_snr(snr_f, r[:, 0] / crl):6.2f} dB")
        # comparison table at the simulated SNRs
        if p4 is not None:
            key = f"ta{ta:g}_s2026"
            g, mse = p4[key + "_grid"], p4[key + "_mse"]
            gross = p4[key + "_gross"]
            print("   dB   sim_MSE   Rice     union1.335  union1.0   sim_Pgross  Rice_Pe  union_Pe")
            for k, db in enumerate(g):
                i = np.argmin(np.abs(snr_f - db))
                print(f"  {db:4.0f} {mse[k]:.3e} {r_rice[i, 0]:.3e} {r_u1[i, 0]:.3e}  {r_u2[i, 0]:.3e}   "
                      f"{gross[k]:.4f}     {r_rice[i, 1]:.4f}   {r_u1[i, 1]:.4f}")
    np.savez("p5_out.npz", **{f"ta{ta:g}_{n}": v for ta, tup in res.items()
                              for n, v in zip(["snr", "rice", "u1", "u2", "riceh"], tup)})
