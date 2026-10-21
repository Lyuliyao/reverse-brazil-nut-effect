"""
Figure 9 — force-balance diagnostics for the canonical CFD-DEM case
(d_L = 600 um, d_S = 200 um, rho_p = 2470, phi = 0.10, alpha = 30 deg,
cor = 0.1, mu = 0.5).

Layout: 1 row x 3 cols sharing x-axis (time)
   (A) species mean heights z_L(t), z_S(t)             -- physical heights
   (B) <fz_contact> / W' per particle, per species     -- contact force
   (C) <fz_fluid excess> / W' per particle, per species -- fluid force

Phase boundaries t1, t2 are taken from zero crossings of (z_L - z_S)
(same convention as Figures 6-8). Interval labels above panel A name the
states directly: RBNE phase = [t1, t2], BNE phase = after t2 (the narrow
mixed crossover at t2 is not labelled).

Reads the canonical case's time_series.tsv.
"""
import os
from pathlib import Path
from bundle_io import REPO, CASES, OUTPUT
import numpy as np
import pandas as pd
import json
import matplotlib as mpl
import matplotlib.pyplot as plt

# Font sizes chosen so the printed size matches Figure 6:
# this figure renders ~12.5 in wide and is inserted at 0.9\linewidth
# (jfm.cls textwidth = 32pc = 5.31 in), so printed = pt x 0.38.
mpl.rcParams.update({
    "font.family":      "STIXGeneral",
    "mathtext.fontset": "stix",
    "axes.titlesize":   19,
    "axes.labelsize":   19,
    "xtick.labelsize":  17,
    "ytick.labelsize":  17,
    "legend.fontsize":  16,
})

ROOT = str(REPO)
CASE = os.path.join(str(CASES),
                    "bne_vf0p1_rs100um_rl300um_rho2470_cor0p1_mu0p5_ang30")
OUT = str(OUTPUT)
os.makedirs(OUT, exist_ok=True)

with open(REPO / "cases" / Path(CASE).name / "case_params.json") as f:
    params = json.load(f)
RL    = params["large_radius"]      # 3e-4 m
RS    = params["small_radius"]      # 1e-4 m
RHOP  = params["rhop"]
RHOF  = 971.0                       # force-normalization parameter
GZ    = abs(params["gz"])           # z-component of gravity magnitude

# Submerged weight (z-component) W' = (rho_p - rho_f) V |g_z|, per particle.
# Natural scale for the z-force balance: in pure terminal settling the
# drag balances W', and contact forces of order W' are what flip the
# net vertical force during the BNE recovery.
V_S = (4.0/3.0) * np.pi * RS**3
V_L = (4.0/3.0) * np.pi * RL**3
W_S = (RHOP - RHOF) * V_S * GZ
W_L = (RHOP - RHOF) * V_L * GZ

# colors — IDENTICAL to the OVITO-style snapshot rasterizer used in
# Figures 6-8, so the line colour for "small" (blue) and "large" (orange)
# matches what the reader sees in the snapshots.
C_S = (0.27, 0.41, 0.78)            # small  royal blue
C_L = (0.93, 0.55, 0.40)            # large  warm coral / peach

# ---------------------------- load ----------------------------------
ts = pd.read_csv(os.path.join(CASE, "time_series.tsv"), sep="\t")
ts = ts.sort_values("time").reset_index(drop=True)

t  = ts["time"].to_numpy()
zS_mm = ts["z_mean_small"].to_numpy() * 1e3        # mm
zL_mm = ts["z_mean_large"].to_numpy() * 1e3
gap   = ts["z_gap_large_minus_small"].to_numpy() * 1e3   # mm
nS    = ts["count_small"].to_numpy()
nL    = ts["count_large"].to_numpy()

# columns in time_series.tsv are already per-particle MEAN forces
# (analyze_bne_case.py uses np.mean). Normalize by submerged weight W'.
# fz_fluid in the CFDEM output includes the pressure-gradient force,
# which contains the hydrostatic (buoyancy) component +rho_f V g_z.
# Subtract it so panel (C) shows the genuine hydrodynamic excess
# (drag + viscous + dynamic-pressure-gradient) -- in pure terminal
# settling this excess equals +W' (drag balances submerged weight).
buoy_S = RHOF * V_S * GZ           # > 0, upward (z-component of buoyancy)
buoy_L = RHOF * V_L * GZ
fzc_S = ts["fz_contact_small"].to_numpy() / W_S
fzc_L = ts["fz_contact_large"].to_numpy() / W_L
fzf_S = (ts["fz_fluid_small"].to_numpy() - buoy_S) / W_S
fzf_L = (ts["fz_fluid_large"].to_numpy() - buoy_L) / W_L

# ---------------- moving-time-average smoothing -------------------
# Time series is sampled at dt = 0.005 s. Use a 1.0 s centered window
# (~200 samples) to suppress per-frame stochastic fluctuations while
# preserving the slow regime-change trend.
DT_SAMPLE = float(np.median(np.diff(t)))           # ~0.005 s
WIN_S     = 1.0                                    # smoothing window in seconds
WIN_N     = max(3, int(round(WIN_S / DT_SAMPLE)))

def smooth(y):
    return (pd.Series(y)
            .rolling(window=WIN_N, center=True, min_periods=1)
            .mean()
            .to_numpy())

zS_mm_s, zL_mm_s = smooth(zS_mm), smooth(zL_mm)
fzc_S_s, fzc_L_s = smooth(fzc_S), smooth(fzc_L)
fzf_S_s, fzf_L_s = smooth(fzf_S), smooth(fzf_L)

# ---------------------------- t1, t2 --------------------------------
sign = np.sign(gap)
idx  = np.where(np.diff(sign) != 0)[0]
zc   = []
for i in idx:
    x0, x1 = t[i], t[i+1]
    y0, y1 = gap[i], gap[i+1]
    if y1 != y0:
        zc.append(x0 - y0 * (x1 - x0) / (y1 - y0))
zc    = np.array(zc)
i_min = int(np.argmin(gap))
t_min = t[i_min]
t1    = float(zc[zc < t_min][-1])
t2    = float(zc[zc > t_min][0])

# ---------------------------- figure --------------------------------
# Render width chosen so the fonts look proportionally right when the
# figure is inserted at 0.9\linewidth in the paper.
fig = plt.figure(figsize=(12.5, 3.9))
gs  = fig.add_gridspec(
    nrows=1, ncols=3,
    wspace=0.30,
    left=0.060, right=0.995, top=0.87, bottom=0.16,
)

panel_letters = ["A", "B", "C"]
axes = []

def add_phase_lines(ax, label_phases=False):
    """Vertical t1/t2 dashed lines on every panel; interval labels above
    panel A name the states directly (RBNE phase = [t1, t2]; BNE phase =
    after t2; the narrow mixed crossover at t2 is not labelled as an
    interval)."""
    ax.axvline(t1, ls="--", color="0.4", lw=1.0)
    ax.axvline(t2, ls="--", color="0.4", lw=1.0)
    if label_phases:
        kw = dict(transform=ax.transAxes, ha="center",
                  va="bottom", fontsize=15.5)
        xmin, xmax = ax.get_xlim()
        def xf(tt): return (tt - xmin) / (xmax - xmin)
        ax.text(xf(0.5*(t1+t2)),    1.02, "RBNE phase",
                color="#7a1f15", **kw)
        ax.text(xf(0.5*(t2+t[-1])), 1.02, "BNE phase",
                color="#1c4e74", **kw)
        ax.text(xf(t1)+0.006, 0.965, r"$t_1$", transform=ax.transAxes,
                fontsize=15.5, color="0.25", ha="left", va="top")
        ax.text(xf(t2)+0.006, 0.965, r"$t_2$", transform=ax.transAxes,
                fontsize=15.5, color="0.25", ha="left", va="top")

def plot_pair(ax, t, raw_L, raw_S, smooth_L, smooth_S,
              label_L, label_S):
    """Draw a (small, large) species pair: raw as thin transparent line,
    smoothed bold on top."""
    ax.plot(t, raw_L,    color=C_L, lw=0.6, alpha=0.30)
    ax.plot(t, raw_S,    color=C_S, lw=0.6, alpha=0.30)
    ax.plot(t, smooth_L, color=C_L, lw=2.2, label=label_L)
    ax.plot(t, smooth_S, color=C_S, lw=2.2, label=label_S)

# (A) species mean heights
axA = fig.add_subplot(gs[0, 0]); axes.append(axA)
plot_pair(axA, t, zL_mm, zS_mm, zL_mm_s, zS_mm_s,
          r"$\bar z_L$ (large)", r"$\bar z_S$ (small)")
axA.set_xlim(0, t[-1])
axA.set_xlabel("time (s)")
axA.set_ylabel("species mean height (mm)")
axA.grid(alpha=0.18)
# lower LEFT is the only reliably empty region: both curves descend into
# the lower right; keep the box clear of the t1 dashed line as well
axA.legend(loc="lower left", bbox_to_anchor=(0.07, 0.01), frameon=False)
add_phase_lines(axA, label_phases=True)

# (B) <fz_contact>/W'
axB = fig.add_subplot(gs[0, 1]); axes.append(axB)
plot_pair(axB, t, fzc_L, fzc_S, fzc_L_s, fzc_S_s,
          "large", "small")
axB.axhline(0, color="k", lw=0.6)
axB.set_xlim(0, t[-1])
axB.set_xlabel("time (s)")
axB.set_ylabel(r"$\langle f_z^{\rm contact}\rangle / W'$")
axB.grid(alpha=0.18)
axB.legend(loc="lower left", frameon=False)
add_phase_lines(axB)

# (C) <fz_fluid excess>/W'
axC = fig.add_subplot(gs[0, 2]); axes.append(axC)
plot_pair(axC, t, fzf_L, fzf_S, fzf_L_s, fzf_S_s,
          "large", "small")
axC.axhline(0, color="k", lw=0.6)
axC.set_xlim(0, t[-1])
axC.set_xlabel("time (s)")
axC.set_ylabel(r"$\langle f_z^{\rm drag}\rangle / W'$")
axC.grid(alpha=0.18)
axC.legend(loc="upper left", frameon=False)
add_phase_lines(axC)

# panel letters outside upper-left of each axes
for letter, ax in zip(panel_letters, axes):
    ax.text(-0.015, 1.06, letter, transform=ax.transAxes,
            fontsize=19, fontweight="bold", ha="right", va="bottom",
            color="0.10")

out_pdf = os.path.join(OUT, "figure09.pdf")
out_png = os.path.join(OUT, "figure09.png")
fig.savefig(out_pdf, bbox_inches="tight")
fig.savefig(out_png, bbox_inches="tight", dpi=180)
print(f"saved {out_pdf}")
print(f"saved {out_png}")
print(f"t1 = {t1:.3f} s   t_min(gap) = {t_min:.3f} s   t2 = {t2:.3f} s")
print(f"final fzc_L/W'_L = {fzc_L[-1]:+.3f}, fzc_S/W'_S = {fzc_S[-1]:+.3f}")
print(f"final (fzf-buoy)_L/W'_L = {fzf_L[-1]:+.3f}, (fzf-buoy)_S/W'_S = {fzf_S[-1]:+.3f}")
print(f"W'_L = {W_L:.3e} N    W'_S = {W_S:.3e} N    W'/mg = {(RHOP-RHOF)/RHOP:.3f}")
