"""SDET Project 13: threshold effects in BPSK time-delay estimation.

    model.py       signal, noise, prior            (Sprints 0-1)
    estimators.py  batch MLE (S1), DLL (S2), DLL theory (S5)   [+ MMSE: extension]
    bounds.py      CRLB (S3), ZZB (S4; numerics S5)
    style.py       plot colours / style
"""
from .model import Config, BPSKWaveform, BatchScenario, rc_pulse, truncated_pulse, random_chips
from .estimators import BatchMLE, EarlyLateDLL
from .bounds import crlb, zzb, autocorrelation, prior_variance, rms_bandwidth_sq, qfunc
