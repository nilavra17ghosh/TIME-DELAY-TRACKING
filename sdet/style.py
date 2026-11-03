"""Shared figure style: validated categorical palette, recessive chrome."""
import matplotlib.pyplot as plt

# fixed entity -> colour mapping, used identically in every figure
COLORS = {
    "mle": "#2a78d6",      # categorical slot 1 (blue)
    "dll": "#eb6834",      # categorical slot 2 (orange)
    "zzb": "#1baf7a",      # categorical slot 3 (aqua)
    "mmse": "#eda100",     # categorical slot 4 (yellow; low contrast -> always legended)
    "crlb": "#0b0b0b",     # reference line: primary ink, dashed
    "prior": "#898781",    # reference line: muted ink, dotted
}
TA_RAMP = ["#86b6ef", "#2a78d6", "#104281"]   # ordinal blue ramp (steps 250/450/650)
SNR_RAMP = ["#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]   # ordinal blue ramp (250/400/550/700)
INK, INK2, MUTED, GRID, AXIS = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
SURFACE = "#fcfcfb"


def apply_style():
    plt.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
        "axes.edgecolor": AXIS, "axes.labelcolor": INK2, "axes.titlecolor": INK,
        "axes.titlesize": 11, "axes.titleweight": "semibold", "axes.labelsize": 10,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
        "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelcolor": INK2, "ytick.labelcolor": INK2,
        "legend.frameon": False, "legend.fontsize": 9, "legend.labelcolor": INK2,
        "lines.linewidth": 2.0, "lines.markersize": 6, "font.size": 10,
        "savefig.dpi": 160, "savefig.bbox": "tight",
    })
