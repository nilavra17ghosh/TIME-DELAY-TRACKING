"""
SDET project - Sprint 1 preliminary simulation
Threshold effects in time-delay tracking: batch MLE vs CRLB for BPSK

WHAT THIS SCRIPT SHOWS (read this first - no signal-processing background needed):

  1. We build a BPSK signal. "BPSK" just means a waveform made of +1 / -1
     pieces called chips. The specific +/- pattern here is a PRN code (a
     pseudo-random m-sequence), the same kind GPS uses.

  2. We look at its autocorrelation: slide a copy of the code over itself and
     measure the match at every shift. A good code gives one tall spike at
     "perfectly aligned" and almost nothing everywhere else (a "thumbtack").
     That tall spike is what an estimator hunts for.

  3. The estimation experiment. We take the signal, delay it by an unknown
     amount tau, add noise, and try to recover tau with the BATCH MLE:
        -> correlate the noisy signal with a clean replica,
        -> the delay is wherever the correlation peaks.
     We repeat this many times at each SNR (noise level) and measure how much
     the estimate scatters -> that scatter is the error variance.

  4. The point of the whole project - the THRESHOLD EFFECT:
        - High SNR: the estimate sits right on the true peak. Its variance
          follows the CRLB (the classical "best possible" limit).
        - Low SNR: noise occasionally makes a WRONG location look taller than
          the true peak, so the estimate jumps far away. The variance shoots
          up and the CRLB no longer describes it. The SNR where this happens
          is the threshold.

  The CRLB is drawn as the reference line. Watch the measured curve hug it at
  high SNR and break away below the threshold - that break-away is exactly the
  phenomenon the ZZB (a later sprint) is designed to predict.

Next sprints add: the early-late DLL, the CRLB derivation, and the ZZB.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.signal import max_len_seq

# ----------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------
RNG = np.random.default_rng(42)     # fixed seed -> reproducible numbers

NBITS   = 7                         # m-sequence degree -> code length 2^7-1 = 127 chips
SPS     = 16                        # samples per chip (oversampling for fractional delay)
TRUE_TAU_CHIPS = 50.37             # the delay we hide and try to recover (in chips)

SNR_DB  = np.arange(-30, 21, 2)     # noise levels to sweep (per-sample SNR, in dB)
TRIALS  = 600                       # Monte-Carlo trials per SNR point

CODE_LEN = (2 ** NBITS) - 1         # 127 chips
N        = CODE_LEN * SPS           # total samples per code period (2032)
TRUE_TAU_SAMP = TRUE_TAU_CHIPS * SPS


# ----------------------------------------------------------------------
# Build the BPSK / PRN waveform
# ----------------------------------------------------------------------
def make_bpsk_waveform():
    """m-sequence code mapped to +/-1, each chip repeated SPS times, then
    band-limited to a finite bandwidth (real BPSK is never infinite-bandwidth).
    Band-limiting also makes the signal's derivative well-defined, so the CRLB
    below is accurate rather than relying on instantaneous square edges."""
    seq, _ = max_len_seq(NBITS)          # 0/1 m-sequence of length 127
    code = 2 * seq.astype(float) - 1     # map 0->-1, 1->+1
    wave = np.repeat(code, SPS).astype(float)   # hold each chip for SPS samples

    freqs = np.fft.fftfreq(len(wave))    # cycles per sample
    fc = 1.5 * (1.0 / SPS)               # cutoff ~1.5x the chip rate
    H = np.exp(-(freqs / fc) ** 2)       # smooth (Gaussian) low-pass
    wave = np.real(np.fft.ifft(np.fft.fft(wave) * H))
    wave /= np.sqrt(np.mean(wave ** 2))  # normalize to unit power
    return code, wave


def frac_delay(x, d):
    """Delay a periodic signal x by d samples (d may be fractional), via FFT."""
    n = len(x)
    k = np.fft.fftfreq(n)
    X = np.fft.fft(x)
    return np.real(np.fft.ifft(X * np.exp(-1j * 2 * np.pi * k * d)))


def circ_xcorr(received, replica):
    """Circular cross-correlation (all delays at once) using the FFT."""
    return np.real(np.fft.ifft(np.fft.fft(received) * np.conj(np.fft.fft(replica))))


def estimate_delay(received, replica):
    """Batch MLE: delay = location of the correlation peak, refined by a
    parabolic fit to the three samples around the peak (sub-sample accuracy)."""
    c = circ_xcorr(received, replica)
    p = int(np.argmax(c))
    ym1, y0, yp1 = c[(p - 1) % N], c[p], c[(p + 1) % N]
    denom = (ym1 - 2 * y0 + yp1)
    offset = 0.5 * (ym1 - yp1) / denom if denom != 0 else 0.0
    return (p + offset) % N


# ----------------------------------------------------------------------
# CRLB for time delay:   var(tau_hat) >= sigma^2 / sum( (ds/dt)^2 )
# ----------------------------------------------------------------------
def crlb_samples(wave, sigma2):
    deriv = np.gradient(wave)            # discrete time-derivative of the signal
    return sigma2 / np.sum(deriv ** 2)   # CRLB in samples^2


# ----------------------------------------------------------------------
# Monte-Carlo sweep
# ----------------------------------------------------------------------
def run():
    code, wave = make_bpsk_waveform()
    replica = wave.copy()                        # clean reference (tau = 0)
    clean = frac_delay(wave, TRUE_TAU_SAMP)      # the delayed signal (noise added later)
    p_sig = np.mean(wave ** 2)                   # signal power (= 1 here)

    emp_var_chip = np.zeros(len(SNR_DB))         # measured variance (chips^2)
    crlb_chip    = np.zeros(len(SNR_DB))         # CRLB reference   (chips^2)

    for i, snr_db in enumerate(SNR_DB):
        sigma2 = p_sig / (10 ** (snr_db / 10.0))     # noise variance for this SNR
        sigma = np.sqrt(sigma2)

        errs = np.empty(TRIALS)
        for t in range(TRIALS):
            received = clean + sigma * RNG.standard_normal(N)
            tau_hat = estimate_delay(received, replica)
            e = tau_hat - TRUE_TAU_SAMP
            e = (e + N / 2) % N - N / 2               # wrap into [-N/2, N/2)
            errs[t] = e
        emp_var_chip[i] = np.var(errs) / (SPS ** 2)   # samples^2 -> chips^2
        crlb_chip[i]    = crlb_samples(wave, sigma2) / (SPS ** 2)

    ceiling_chip = (CODE_LEN ** 2) / 12.0            # a-priori ceiling: uniform over the code
    return code, wave, emp_var_chip, crlb_chip, ceiling_chip


# ----------------------------------------------------------------------
# Figure (one consolidated 3-panel figure)
# ----------------------------------------------------------------------
def plot(code, wave, emp_var_chip, crlb_chip, ceiling_chip):
    C_SIG, C_CRLB, C_EMP, C_CEIL = "#16233F", "#3B7DD8", "#D9622B", "#9AA6B8"
    fig = plt.figure(figsize=(12, 8))
    gs = GridSpec(2, 2, figure=fig, height_ratios=[1, 1.25], hspace=0.38, wspace=0.26)

    # Panel A: a slice of the BPSK waveform
    axA = fig.add_subplot(gs[0, 0])
    nshow = 10 * SPS
    axA.step(np.arange(nshow) / SPS, wave[:nshow], where="post", color=C_SIG, lw=1.8)
    axA.set_title("BPSK waveform  (first 10 chips)", fontsize=11, weight="bold")
    axA.set_xlabel("chips"); axA.set_ylabel("amplitude")
    axA.set_ylim(-1.5, 1.5); axA.grid(alpha=0.25)

    # Panel B: autocorrelation (the thumbtack)
    axB = fig.add_subplot(gs[0, 1])
    full = np.repeat(code, 1).astype(float)
    ac = np.real(np.fft.ifft(np.abs(np.fft.fft(full)) ** 2)) / CODE_LEN
    lags = np.arange(CODE_LEN); lags[lags > CODE_LEN // 2] -= CODE_LEN
    order = np.argsort(lags)
    axB.plot(lags[order], ac[order], color=C_SIG, lw=1.6)
    axB.set_title("Code autocorrelation  (thumbtack)", fontsize=11, weight="bold")
    axB.set_xlabel("lag (chips)"); axB.set_ylabel("normalized")
    axB.grid(alpha=0.25)

    # Panel C: the threshold effect - variance vs SNR
    axC = fig.add_subplot(gs[1, :])
    axC.semilogy(SNR_DB, crlb_chip, "--", color=C_CRLB, lw=2.2, label="CRLB (classical bound)")
    axC.semilogy(SNR_DB, emp_var_chip, "o-", color=C_EMP, lw=2.2, ms=5,
                 label="Batch MLE variance (measured)")
    axC.axhline(ceiling_chip, color=C_CEIL, ls=":", lw=1.8,
                label="a-priori ceiling (code length)")
    axC.set_title("The threshold effect: batch MLE vs CRLB",
                  fontsize=12, weight="bold")
    axC.set_xlabel("SNR (dB)  ->  higher"); axC.set_ylabel(r"delay-error variance (chips$^2$)")
    axC.grid(which="both", alpha=0.25); axC.legend(loc="upper right", fontsize=10)

    fig.suptitle("SDET Sprint 1 - preliminary results (BPSK, 127-chip m-sequence)",
                 fontsize=13, weight="bold")
    fig.savefig("sdet_sprint1_results.png", dpi=150,
                bbox_inches="tight")
    print("saved figure -> sdet_sprint1_results.png")


if __name__ == "__main__":
    code, wave, emp_var_chip, crlb_chip, ceiling_chip = run()
    plot(code, wave, emp_var_chip, crlb_chip, ceiling_chip)

    # quick text summary for the meeting
    print("\n  SNR(dB)   MLE var (chips^2)     CRLB (chips^2)")
    for s, e, c in zip(SNR_DB, emp_var_chip, crlb_chip):
        print(f"   {s:5.0f}     {e:14.4g}     {c:14.4g}")
    # find roughly where the measured curve leaves the CRLB (x5)
    ratio = emp_var_chip / crlb_chip
    above = np.where(ratio > 5)[0]
    if len(above):
        thr = SNR_DB[above[-1]]
        print(f"\n  Approx threshold SNR (MLE var > 5x CRLB up to): ~{thr:.0f} dB")
