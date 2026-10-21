"""
Figure 10 — mechanism-isolation controls.

Overlays Delta z_bar(t) = z_L_bar - z_S_bar for four runs of the canonical
parameters (d_L = 600 um, d_S = 200 um, phi = 0.10, alpha = 30 deg; a separate
realisation from the Figure 6 case, insertion seed 11000081):

    baseline (full incline)   g = (g sin a, 0, -g cos a)
    settling-only             g_par = 0 (no downslope shear)
    shear-only                g_perp = 0 (no wall-normal settling drive)
    equal-settling (Stokes)   full gravity + constant upward body force on
                              each large particle
                              F_add = W'_L [1-(d_S/d_L)^2] = (8/9) W'_L
                              = 1.278e-6 N, so the isolated-sphere Stokes
                              settling velocities of the two species coincide

Reads data/mechanisms/*_meanz.csv and writes PNG and PDF to the output directory.
"""
import csv
import os
from bundle_io import REPO, OUTPUT
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams.update({
    "font.family":      "STIXGeneral",
    "mathtext.fontset": "stix",
    "axes.labelsize":   17,
    "xtick.labelsize":  15,
    "ytick.labelsize":  15,
    "legend.fontsize":  14.5,
})

RESULTS = str(REPO / "data" / "mechanisms")
OUT = str(OUTPUT)
os.makedirs(OUT, exist_ok=True)

# palette shared with Figures 6-9: royal blue / warm coral
C_S = (0.27, 0.41, 0.78)
C_L = (0.93, 0.55, 0.40)

SERIES = [
    ("mech_settling_only_dS200_dL600_phi0p1_ang30_seed1",
     r"settling-only ($g_\parallel = 0$)", "#b0413e", "-", 1.9, 2),
    ("mech_shear_only_dS200_dL600_phi0p1_ang30_seed1",
     r"shear-only ($g_\perp = 0$)", "0.40", "--", 1.7, 3),
    ("mech_baseline_dS200_dL600_phi0p1_ang30_seed1",
     "baseline (full incline)", "black", "-", 2.3, 4),
    ("mech_equal_settling_stokes_dS200_dL600_phi0p1_ang30_seed1",
     r"equal settling (Stokes $F_{\mathrm{add}} = \frac{8}{9}W_L'$)",
     C_S, "-", 1.9, 5),
]


def read_dz_mm(stem):
    path = os.path.join(RESULTS, f"{stem}_meanz.csv")
    t, dz = [], []
    with open(path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            try:
                t.append(float(row["time"]))
                dz.append(1.0e3 * float(row["dz"]))
            except (KeyError, ValueError):
                continue
    if not t:
        raise SystemExit(f"empty CSV: {path}")
    return np.array(t), np.array(dz)


fig, ax = plt.subplots(figsize=(8.2, 5.0))

curves = [(lbl, c, ls, lw, z, *read_dz_mm(stem))
          for stem, lbl, c, ls, lw, z in SERIES]
ymin = min(dz.min() for *_, dz in curves)
ymax = max(dz.max() for *_, dz in curves)
pad = 0.08 * (ymax - ymin)

# RBNE / BNE bands, same convention as Figure 6E
ax.axhspan(ymin - pad, 0.0, color=C_L, alpha=0.14, lw=0, zorder=0)
ax.axhspan(0.0, ymax + pad, color=C_S, alpha=0.10, lw=0, zorder=0)
ax.axhline(0.0, color="k", lw=0.6, zorder=1)

for lbl, c, ls, lw, z, t, dz in curves:
    ax.plot(t, dz, color=c, ls=ls, lw=lw, label=lbl, zorder=z)

ax.set_xlim(0, max(t[-1] for *_, t, _dz in curves))
ax.set_ylim(ymin - pad, ymax + pad)
ax.set_xlabel("time (s)")
ax.set_ylabel(r"$\Delta \bar z = \bar z_{L} - \bar z_{S}$  (mm)")
ax.text(0.025, 0.045, r"RBNE  ($\Delta\bar z<0$)", transform=ax.transAxes,
        ha="left", va="bottom", color="#7a1f15", fontsize=15)
ax.text(0.977, 0.955, r"BNE  ($\Delta\bar z>0$)", transform=ax.transAxes,
        ha="right", va="top", color="#1c4e74", fontsize=15)
# middle-right band (t > 20, -1 < dz < -0.35 mm) is empty of data
ax.legend(loc="center right", bbox_to_anchor=(0.985, 0.33),
          frameon=False)
ax.grid(alpha=0.20)

for ext in ("pdf", "png"):
    fig.savefig(os.path.join(OUT, f"figure10.{ext}"),
                bbox_inches="tight",
                **({"dpi": 200} if ext == "png" else {}))
plt.close(fig)

for lbl, *_rest, t, dz in curves:
    print(f"  {lbl[:46]:46s} t:[{t[0]:.2f},{t[-1]:.2f}]s "
          f"min={dz.min():+.3f}mm  end={dz[-1]:+.3f}mm")
print("saved figure10.{png,pdf}")
