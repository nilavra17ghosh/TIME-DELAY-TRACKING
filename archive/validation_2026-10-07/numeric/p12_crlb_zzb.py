"""Parts 1 & 2: independent CRLB (continuous + discrete) and ZZB."""
import numpy as np
from vcommon import pilot, Q, waveform, waveform_deriv

PROJ = r"C:\Users\Nilavra Ghosh\Documents\7THSEM\SDET\project\results\tables\summary.csv"
summ = np.genfromtxt(PROJ, delimiter=",", names=True)
snr = summ["snr_db"]
en0 = 10 ** (snr / 10)
chips = pilot()
TA = 4.0

print("pilot first 10 chips:", chips[:10])

# ---------------------------------------------------------------- 1(a) continuous beta^2 via FFT
def beta2_fft(span=6, mode="sym", fs=64, nfft=2 ** 18):
    t = np.arange(-span - 2, 64 + span + 2, 1.0 / fs)
    s = waveform(t, chips, span=span, mode=mode)
    S = np.fft.rfft(s, nfft)
    f = np.fft.rfftfreq(nfft, d=1.0 / fs)
    P = np.abs(S) ** 2
    return np.sum((2 * np.pi * f) ** 2 * P) / np.sum(P), s, t


def beta2_ideal_rc(b=0.5):
    """Untruncated RC: |S(f)|^2 = |P(f)|^2 |A(f)|^2, analytic P(f)."""
    f = np.linspace(-(1 + b) / 2, (1 + b) / 2, 400001)
    af = np.abs(f)
    P = np.where(af <= (1 - b) / 2, 1.0, 0.5 * (1 + np.cos(np.pi / b * (af - (1 - b) / 2))))
    A = np.exp(-2j * np.pi * np.outer(f, np.arange(64))) @ chips if False else None
    # chunked |A(f)|^2
    A2 = np.empty_like(f)
    k = np.arange(64)
    for i in range(0, f.size, 20000):
        A2[i:i + 20000] = np.abs(np.exp(-2j * np.pi * f[i:i + 20000, None] * k) @ chips) ** 2
    w = P ** 2 * A2
    return np.trapezoid((2 * np.pi * f) ** 2 * w, f) / np.trapezoid(w, f)


def beta2_time_cont(span=6, mode="sym", fs=256):
    t = np.arange(-span - 1, 64 + span + 1, 1.0 / fs)
    s = waveform(t, chips, span=span, mode=mode)
    ds = waveform_deriv(t, chips, span=span, mode=mode)
    return np.sum(ds ** 2) / np.sum(s ** 2)


def beta2_discrete(sps=4, span=6, mode="sym", frac=0.0):
    t = np.arange(-span - 3, 64 + span + 3, 1.0 / sps) + frac
    s = waveform(t, chips, span=span, mode=mode)
    ds = waveform_deriv(t, chips, span=span, mode=mode)
    return np.sum(ds ** 2) / np.sum(s ** 2)


b_fft6, s_f, t_f = beta2_fft(6)
b_fft6a, _, _ = beta2_fft(6, "asym")
b_fft12, _, _ = beta2_fft(12)
b_ideal = beta2_ideal_rc()
b_tc6 = beta2_time_cont(6)
b_d4 = beta2_discrete(4, 6)
b_d4a = beta2_discrete(4, 6, "asym")
b_d16 = beta2_discrete(16, 6)
phases = np.linspace(0, 0.25, 11)[:-1]
b_phase = np.array([beta2_discrete(4, 6, "asym", p) for p in phases])
b_proj_from_csv = 1.0 / (2 * en0[0] * summ["crlb"][0])
print("\n=== Part 1: beta^2 (1/chip^2) ===")
print(f"project (from summary.csv crlb)      : {b_proj_from_csv:.5f}")
print(f"continuous FFT, span6 sym (64 sps)   : {b_fft6:.5f}")
print(f"continuous FFT, span6 asym [-6,7)    : {b_fft6a:.5f}")
print(f"continuous FFT, span12 sym           : {b_fft12:.5f}")
print(f"continuous ideal RC (no truncation)  : {b_ideal:.5f}")
print(f"continuous time-domain s' (256 sps)  : {b_tc6:.5f}")
print(f"discrete sps4 sym, phase 0           : {b_d4:.5f}")
print(f"discrete sps4 asym, phase 0          : {b_d4a:.5f}")
print(f"discrete sps16 sym                   : {b_d16:.5f}")
print(f"discrete sps4 asym over sampling phase: min {b_phase.min():.5f} max {b_phase.max():.5f}")
for name, b in [("FFT span6", b_fft6), ("ideal RC", b_ideal), ("disc sps4", b_d4a)]:
    print(f"  rel diff project vs {name}: {(b_proj_from_csv / b - 1) * 100:+.3f} %  (CRLB rel diff {(b / b_proj_from_csv - 1) * 100:+.3f} %)")
crlb_mine = 1 / (2 * en0 * b_fft6)
print("max |CRLB_proj/CRLB_mine(FFT) - 1| over grid:", np.max(np.abs(summ["crlb"] / crlb_mine - 1)))

# ---------------------------------------------------------------- 2. ZZB
def rho_cont(h, span=6, mode="sym", fs=16):
    """Continuous-time normalised autocorrelation, evaluated EXACTLY at each lag h via
    R(h) = int |S(f)|^2 cos(2 pi f h) df (no lag interpolation)."""
    t = np.arange(-span - 2, 64 + span + 2, 1.0 / fs)
    s = waveform(t, chips, span=span, mode=mode)
    n = 1 << 15
    S = np.fft.rfft(s, n)
    f = np.fft.rfftfreq(n, 1.0 / fs)
    w = np.abs(S) ** 2
    w[1:] *= 2.0
    keep = f < 4.0
    f, w = f[keep], w[keep]
    h = np.atleast_1d(h)
    R = np.concatenate([np.cos(2 * np.pi * np.outer(h[i:i + 500], f)) @ w for i in range(0, h.size, 500)])
    return R / np.sum(w)


def rho_disc(h, sps=4, span=6, mode="sym", phi=0.0):
    t = np.arange(-span - 3, 64 + span + TA + 3, 1.0 / sps)
    s0 = waveform(t - phi, chips, span=span, mode=mode)
    out = np.empty(len(h))
    for i, hh in enumerate(h):
        out[i] = np.dot(s0, waveform(t - phi - hh, chips, span=span, mode=mode))
    return out / np.dot(s0, s0)


def zzb_from_rho(h, rho, en0, vf=True):
    f = (1 - h / TA)[None, :] * Q(np.sqrt(np.atleast_1d(en0)[:, None] * np.clip(1 - rho, 0, None)[None, :]))
    if vf:
        f = np.maximum.accumulate(f[:, ::-1], axis=1)[:, ::-1]
    return np.trapezoid(h[None, :] * f, h, axis=1)


def thr(snr_f, ratio, level=2.0):
    lr, lv = np.log10(ratio), np.log10(level)
    for i in range(len(snr_f) - 1, 0, -1):
        if lr[i] <= lv < lr[i - 1]:
            f = (lv - lr[i]) / (lr[i - 1] - lr[i])
            return snr_f[i] + f * (snr_f[i - 1] - snr_f[i])
    return np.nan


h = np.linspace(0, TA, 4001)
rc6 = rho_cont(h, 6)
z_vf = zzb_from_rho(h, rc6, en0)
z_novf = zzb_from_rho(h, rc6, en0, vf=False)
print("\n=== Part 2: ZZB (continuous rho, n_h=4001, span6 sym) vs project ===")
print(" snr   proj_zzb    my_zzb     rel%    my_noVF   VF-gain%  my/CRLB")
for i in range(len(snr)):
    print(f"{snr[i]:5.0f} {summ['zzb'][i]:.5e} {z_vf[i]:.5e} {100 * (z_vf[i] / summ['zzb'][i] - 1):+7.3f} "
          f"{z_novf[i]:.5e} {100 * (z_vf[i] / z_novf[i] - 1):7.2f} {z_vf[i] / crlb_mine[i]:7.3f}")
print("max |rel diff| ZZB:", np.max(np.abs(z_vf / summ["zzb"] - 1)))

# sensitivity
snr_f = np.linspace(-10, 30, 801)
en_f = 10 ** (snr_f / 10)
variants = {}
hh = np.linspace(0, TA, 2001)
variants["cont rho, span6, n_h=2001"] = (zzb_from_rho(hh, rho_cont(hh, 6), en_f), b_fft6)
hh = np.linspace(0, TA, 501)
variants["cont rho, span6, n_h=501"] = (zzb_from_rho(hh, rho_cont(hh, 6), en_f), b_fft6)
hh = np.linspace(0, TA, 20001)
variants["cont rho, span6, n_h=20001"] = (zzb_from_rho(hh, rho_cont(hh, 6), en_f), b_fft6)
hh = np.linspace(0, TA, 2001)
variants["cont rho, span12"] = (zzb_from_rho(hh, rho_cont(hh, 12), en_f), b_fft12)
variants["cont rho, span6 asym"] = (zzb_from_rho(hh, rho_cont(hh, 6, "asym"), en_f), b_fft6a)
variants["disc rho sps4, span6 asym (=project)"] = (zzb_from_rho(hh, rho_disc(hh, 4, 6, "asym"), en_f), b_d4a)
variants["disc rho sps16, span6 sym"] = (zzb_from_rho(hh, rho_disc(hh, 16, 6), en_f), b_d16)
variants["disc rho sps4, span12 sym"] = (zzb_from_rho(hh, rho_disc(hh, 4, 12), en_f), beta2_discrete(4, 12))
variants["cont rho, NO valley filling"] = (zzb_from_rho(hh, rho_cont(hh, 6), en_f, vf=False), b_fft6)
ref = variants["cont rho, span6, n_h=2001"][0]
print("\n=== ZZB / CRLB sensitivity ===")
print(f"{'variant':40s} {'beta2':>8s} {'ZZBthr dB':>9s} {'max|dZZB|%':>10s} {'ZZB@6dB':>9s} {'ZZB@10dB':>9s}")
i6, i10 = np.argmin(np.abs(snr_f - 6)), np.argmin(np.abs(snr_f - 10))
for k, (z, b2) in variants.items():
    c = 1 / (2 * en_f * b2)
    print(f"{k:40s} {b2:8.4f} {thr(snr_f, z / c):9.3f} {100 * np.max(np.abs(z / ref - 1)):10.3f} {z[i6]:9.5f} {z[i10]:9.5f}")
print("ZZB low-SNR (-40 dB):", zzb_from_rho(hh, rho_cont(hh, 6), 1e-4), "Ta^2/12 =", TA ** 2 / 12)
# first zero of rho / main-lobe width
r = rho_cont(np.linspace(0, TA, 4001), 6)
hz = np.linspace(0, TA, 4001)[np.argmax(r < 0)]
print("first zero crossing of rho(h):", hz, " max |rho| for h>1.5:", np.max(np.abs(r[np.linspace(0, TA, 4001) > 1.5])))
np.savez("p12_out.npz", h=np.linspace(0, TA, 4001), rho=r, b2=b_fft6, snr_f=snr_f, zzb_f=ref,
         zzb_novf=variants["cont rho, NO valley filling"][0])
