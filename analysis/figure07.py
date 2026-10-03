"""
Paper simulation figure 2 — effect of total volume fraction phi and
incline angle alpha on the transient RBNE -> BNE evolution.

Fixed parameters: d_L = 600 um, d_S = 200 um, rho_p = 2470 kg/m^3,
                  cor = 0.1, mu = 0.5
Scanned         : phi in {0.05, 0.10, 0.20, 0.30, 0.40}
                  alpha in {20, 30} deg
                  (alpha = 10 deg runs were terminated early and are
                  excluded by request)

Layout: 3 rows x 5 columns
        row 1: peak-RBNE snapshot for alpha = 30 deg (canonical)
        row 2: end-state snapshot for alpha = 30 deg
        row 3: Delta zbar(t) overlay; alpha = 30 (solid), alpha = 20 (dashed)

Each trace uses archived per-frame species means; the selected particle
snapshots were extracted from packed forcebalance.h5 segments.
"""
import os, glob
import numpy as np
import h5py
from pathlib import Path
from bundle_io import REPO, CASES, OUTPUT, load_snapshot, read_height_series, read_fluid_profile
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

mpl.rcParams.update({
    "font.family":      "STIXGeneral",
    "mathtext.fontset": "stix",
    "axes.titlesize":   18,
    "axes.labelsize":   18,
    "xtick.labelsize":  15,
    "ytick.labelsize":  15,
    "legend.fontsize":  15,
})

ROOT = str(REPO)
SBATCH_C = str(CASES)
SCR_C = str(CASES)
OUT = str(OUTPUT)
os.makedirs(OUT, exist_ok=True)

DEM_DT  = 1.0e-7
C_S     = np.array([0.27, 0.41, 0.78])    # small  royal blue
C_L     = np.array([0.93, 0.55, 0.40])    # large  warm coral
EPS_MM  = 0.10                             # |Delta z_bar| threshold for mixed (mm)

PHIS = [0.05, 0.10, 0.20, 0.30, 0.40]
ALPHAS = [20, 30]

# ---------- locate case directories across SBATCH and SCR ----------
def case_dir(phi, alpha):
    """Find the case directory for a given (phi, alpha) — searches
    /mnt/home SBATCH first, then /mnt/scratch."""
    name = f"bne_vf0p{int(round(phi*100)):02d}".replace("vf0p05", "vf0p05")
    name = name.replace("vf0p10", "vf0p1")
    name = name.replace("vf0p20", "vf0p2")
    name = name.replace("vf0p30", "vf0p3")
    name = name.replace("vf0p40", "vf0p4")
    name = f"{name}_rs100um_rl300um_rho2470_cor0p1_mu0p5_ang{alpha}"
    for root in (SBATCH_C, SCR_C):
        p = os.path.join(root, name)
        if os.path.isdir(p):
            return p
    raise FileNotFoundError(f"case not found: {name}")


def case_segments(case_dir):
    seg_root = os.path.join(case_dir, "DEM", "post", "forcebalance_segments")
    return sorted(glob.glob(os.path.join(seg_root, "t*_to_*", "forcebalance.h5")))


def compute_dzbar(case_dir):
    return read_height_series(case_dir)


def load_frame(case_dir, t_target):
    return load_snapshot(case_dir, t_target)


# --------------------- OVITO-style sphere rasterizer -----------------
def render_snapshot(xyz_m, types, radii_m, box,
                    px_per_mm=80, light=(-0.40, 0.40, 0.85),
                    ambient=0.40, diff_w=0.50, spec_w=0.40, spec_pow=18,
                    bg=(1.0, 1.0, 1.0)):
    x_mm = xyz_m[:, 0] * 1e3
    y_mm = xyz_m[:, 1] * 1e3
    z_mm = xyz_m[:, 2] * 1e3
    r_mm = radii_m * 1e3
    Lx_mm = (box[0, 1] - box[0, 0]) * 1e3
    H_mm  = (box[2, 1] - box[2, 0]) * 1e3
    W = int(round(Lx_mm * px_per_mm))
    H = int(round(H_mm  * px_per_mm))
    rgb  = np.tile(np.array(bg, dtype=np.float32), (H, W, 1))
    zbuf = np.full((H, W), -np.inf, dtype=np.float32)
    lx, ly, lz = light
    n = (lx*lx + ly*ly + lz*lz) ** 0.5
    lx, ly, lz = lx / n, ly / n, lz / n
    order = np.argsort(y_mm)
    for i in order:
        rp = r_mm[i] * px_per_mm
        cx = (x_mm[i] - box[0, 0]*1e3) * px_per_mm
        cy = (z_mm[i] - box[2, 0]*1e3) * px_per_mm
        x0 = max(0, int(np.floor(cx - rp)) - 1)
        x1 = min(W, int(np.ceil (cx + rp)) + 1)
        y0 = max(0, int(np.floor(cy - rp)) - 1)
        y1 = min(H, int(np.ceil (cy + rp)) + 1)
        if x0 >= x1 or y0 >= y1:
            continue
        yy, xx = np.mgrid[y0:y1, x0:x1].astype(np.float32)
        nx = (xx - cx) / rp
        ny = (yy - cy) / rp
        d2 = nx*nx + ny*ny
        mask = d2 <= 1.0
        if not mask.any():
            continue
        nz = np.sqrt(np.clip(1.0 - d2, 0.0, 1.0))
        depth = y_mm[i] + r_mm[i] * nz
        diff = np.clip(nx*lx + ny*ly + nz*lz, 0.0, 1.0)
        hx, hy, hz = lx, ly, lz + 1.0
        nh = (hx*hx + hy*hy + hz*hz) ** 0.5
        hx, hy, hz = hx/nh, hy/nh, hz/nh
        spec = np.clip(nx*hx + ny*hy + nz*hz, 0.0, 1.0) ** spec_pow
        base = (C_S if types[i] == 1 else C_L).astype(np.float32)
        col  = (ambient + diff_w * diff)[..., None] * base \
               + spec_w * spec[..., None]
        col  = np.clip(col, 0.0, 1.0)
        rim  = np.clip((1.0 - nz)**4, 0.0, 1.0)[..., None]
        col  = col * (1.0 - 0.22 * rim)
        sub_z = zbuf[y0:y1, x0:x1]
        sub_c = rgb [y0:y1, x0:x1]
        m2 = mask & (depth > sub_z)
        sub_z[m2] = depth[m2]
        sub_c[m2] = col[m2]
    return rgb, (0.0, Lx_mm, 0.0, H_mm)


# ---------------------- gather data ---------------------------------
print("computing dz_bar(t) and locating snapshot frames per case...")
data = {}                                    # data[(phi, alpha)] = dict
for phi in PHIS:
    for ang in ALPHAS:
        cdir = case_dir(phi, ang)
        t, dz = compute_dzbar(cdir)
        i_min = int(np.argmin(dz))
        rec = dict(case_dir=cdir, t=t, dz=dz,
                   t_min=float(t[i_min]),
                   dz_min=float(dz.min()),
                   t_end=float(t[-1]))
        data[(phi, ang)] = rec
        print(f"  phi={phi}, a={ang}deg:  t<=[{t[0]:.2f},{t[-1]:.2f}]  "
              f"dz_min={rec['dz_min']:+.3f} mm @ t={rec['t_min']:.2f} s  "
              f"end dz={dz[-1]:+.3f} mm")

# ---------------------- figure layout -------------------------------
NCOL = len(PHIS)                              # 5
# match Fig 1's 13-inch render width so font size on the printed page is
# the same as Fig 1; height scales with the 3-row content
fig = plt.figure(figsize=(13.0, 5.6))
gs  = fig.add_gridspec(
    nrows=3, ncols=NCOL,
    height_ratios=[1.0, 1.0, 1.55],
    hspace=0.32, wspace=0.05,
    left=0.055, right=0.99, top=0.935, bottom=0.085,
)

# common axis ranges
t_xmax = max(d["t_end"] for d in data.values())
dz_max_mag = max(max(abs(d["dz"].min()), abs(d["dz"].max()))
                 for d in data.values())
dz_lim = 1.10 * dz_max_mag

# --- snapshots (rows 1, 2) use alpha = 30 deg ---
SNAP_ALPHA = 30
letters_row1 = ["A", "B", "C", "D", "E"]
letters_row2 = ["F", "G", "H", "I", "J"]
letters_row3 = ["K", "L", "M", "N", "O"]

for j, phi in enumerate(PHIS):
    rec = data[(phi, SNAP_ALPHA)]

    # ---- row 1: peak RBNE snapshot ----
    ax1 = fig.add_subplot(gs[0, j])
    xyz, typ, rad, t_act, bb = load_frame(rec["case_dir"], rec["t_min"])
    img, ext = render_snapshot(xyz, typ, rad, bb, px_per_mm=85)
    ax1.imshow(img, extent=ext, origin="lower", interpolation="bilinear")
    ax1.add_patch(Rectangle((0, 0), ext[1], ext[3],
                            fill=False, edgecolor="0.15", linewidth=1.0))
    ax1.plot([0, ext[1]], [0, 0], color="k", lw=1.4)
    ax1.set_xlim(0, ext[1]); ax1.set_ylim(0, ext[3])
    ax1.set_aspect("equal"); ax1.set_xticks([]); ax1.set_yticks([])
    for s in ax1.spines.values(): s.set_visible(False)
    ax1.text(0.985, 0.93, "peak RBNE", transform=ax1.transAxes,
             fontsize=15.5, color="0.10", ha="right", va="top",
             bbox=dict(boxstyle="round,pad=0.22", fc="white",
                       ec="0.55", lw=0.5))
    ax1.text(0.5, -0.10, rf"$t = {t_act:.2f}$ s",
             transform=ax1.transAxes, fontsize=16.5, color="0.10",
             ha="center", va="top")
    if j == 0:
        ax1.plot([0.6, 1.6], [0.32, 0.32], color="k", lw=2.4,
                 solid_capstyle="butt")
        ax1.text(1.1, 0.10, "1 mm", ha="center", va="top",
                 fontsize=14, color="k")
    ax1.set_title(rf"$\phi = {phi:.2f}$", fontsize=18, pad=8)
    ax1.text(-0.015, 1.06, letters_row1[j], transform=ax1.transAxes,
             fontsize=19, fontweight="bold", ha="right", va="bottom",
             color="0.10")

    # ---- row 2: end-state snapshot ----
    ax2 = fig.add_subplot(gs[1, j])
    xyz, typ, rad, t_act, bb = load_frame(rec["case_dir"], rec["t_end"])
    img, ext = render_snapshot(xyz, typ, rad, bb, px_per_mm=85)
    ax2.imshow(img, extent=ext, origin="lower", interpolation="bilinear")
    ax2.add_patch(Rectangle((0, 0), ext[1], ext[3],
                            fill=False, edgecolor="0.15", linewidth=1.0))
    ax2.plot([0, ext[1]], [0, 0], color="k", lw=1.4)
    ax2.set_xlim(0, ext[1]); ax2.set_ylim(0, ext[3])
    ax2.set_aspect("equal"); ax2.set_xticks([]); ax2.set_yticks([])
    for s in ax2.spines.values(): s.set_visible(False)
    end_dz = rec["dz"][-1]
    end_lbl = ("BNE"   if end_dz > +EPS_MM else
               "RBNE"  if end_dz < -EPS_MM else "mixed")
    ax2.text(0.985, 0.93, f"end ({end_lbl})", transform=ax2.transAxes,
             fontsize=15.5, color="0.10", ha="right", va="top",
             bbox=dict(boxstyle="round,pad=0.22", fc="white",
                       ec="0.55", lw=0.5))
    ax2.text(0.5, -0.10, rf"$t = {t_act:.2f}$ s",
             transform=ax2.transAxes, fontsize=16.5, color="0.10",
             ha="center", va="top")
    ax2.text(-0.015, 1.02, letters_row2[j], transform=ax2.transAxes,
             fontsize=19, fontweight="bold", ha="right", va="bottom",
             color="0.10")

    # ---- row 3: Delta zbar(t) overlay (both alphas) ----
    ax3 = fig.add_subplot(gs[2, j])
    rec30 = data[(phi, 30)]
    rec20 = data[(phi, 20)]
    # symmetric shading
    ax3.fill_between(rec30["t"], rec30["dz"], 0,
                     where=(rec30["dz"] <  0), color=C_L, alpha=0.20,
                     linewidth=0)
    ax3.fill_between(rec30["t"], rec30["dz"], 0,
                     where=(rec30["dz"] >= 0), color=C_S, alpha=0.20,
                     linewidth=0)
    ax3.plot(rec30["t"], rec30["dz"], color="k", lw=1.6,
             label=r"$\alpha = 30^\circ$")
    ax3.plot(rec20["t"], rec20["dz"], color="0.30", lw=1.2,
             ls="--", label=r"$\alpha = 20^\circ$")
    ax3.axhline(0, color="k", lw=0.6)
    ax3.axhline(+EPS_MM, color="0.55", ls=":", lw=0.7)
    ax3.axhline(-EPS_MM, color="0.55", ls=":", lw=0.7)
    # mark the snapshot times for alpha=30
    for t_mark in (rec30["t_min"], rec30["t_end"]):
        ax3.axvline(t_mark, color="k", lw=0.6, alpha=0.30)
    ax3.scatter([rec30["t_min"]],
                [np.interp(rec30["t_min"], rec30["t"], rec30["dz"])],
                s=36, marker="o", facecolor="white",
                edgecolor="k", linewidth=1.0, zorder=5)
    ax3.scatter([rec30["t_end"]],
                [np.interp(rec30["t_end"], rec30["t"], rec30["dz"])],
                s=36, marker="o", facecolor="white",
                edgecolor="k", linewidth=1.0, zorder=5)
    ax3.set_xlim(0, t_xmax)
    ax3.set_ylim(-dz_lim, dz_lim)
    ax3.grid(alpha=0.18)
    ax3.set_xlabel("time (s)", fontsize=18)
    if j == 0:
        ax3.set_ylabel(r"$\Delta\bar{z} = \bar{z}_L - \bar{z}_S$  (mm)",
                       fontsize=18)
    else:
        ax3.set_yticklabels([])
    ax3.text(-0.015, 1.02, letters_row3[j], transform=ax3.transAxes,
             fontsize=19, fontweight="bold", ha="right", va="bottom",
             color="0.10")
    if j == 0:
        ax3.legend(loc="upper left", frameon=False, fontsize=15,
                   bbox_to_anchor=(0.0, 0.92))

out_pdf = os.path.join(OUT, "figure07.pdf")
out_png = os.path.join(OUT, "figure07.png")
fig.savefig(out_pdf, bbox_inches="tight")
fig.savefig(out_png, bbox_inches="tight", dpi=180)
print(f"saved {out_pdf}")
print(f"saved {out_png}")
