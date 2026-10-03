"""
Same as fig_sim3_radius_phi_ang20.py but at alpha = 30 deg.

   row 1: d_L = 600  um peak-RBNE snapshot
   row 2: d_L = 1200 um peak-RBNE snapshot
   row 3: Delta zbar(t) overlay for both d_L values

5 columns: phi in {0.05, 0.10, 0.20, 0.30, 0.40}.
Fixed: d_S = 200 um, rho_p = 2470 kg/m^3, cor = 0.1, mu = 0.5,
       alpha = 30 deg.
"""
import os, glob
import numpy as np
import h5py
from pathlib import Path
from bundle_io import REPO, CASES, OUTPUT, load_snapshot, read_height_series, read_fluid_profile
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

mpl.rcParams.update({
    "font.family":      "STIXGeneral",
    "mathtext.fontset": "stix",
    "axes.titlesize":   18,
    "axes.labelsize":   18,
    "xtick.labelsize":  15,
    "ytick.labelsize":  15,
})

ROOT = str(REPO)
SBATCH_C = str(CASES)
SCR_C = str(CASES)
OUT = str(OUTPUT)
os.makedirs(OUT, exist_ok=True)

DEM_DT  = 1.0e-7
C_S     = np.array([0.27, 0.41, 0.78])
C_L     = np.array([0.93, 0.55, 0.40])

PHIS    = [0.05, 0.10, 0.20, 0.30, 0.40]
ALPHA   = 30
RL_TAGS = [("rl300um", 600), ("rl600um", 1200)]    # (folder tag, d_L in um)


def case_dir(phi, rl_tag):
    name = f"bne_vf0p{int(round(phi*100)):02d}"
    name = name.replace("vf0p10", "vf0p1") \
               .replace("vf0p20", "vf0p2") \
               .replace("vf0p30", "vf0p3") \
               .replace("vf0p40", "vf0p4")
    name = f"{name}_rs100um_{rl_tag}_rho2470_cor0p1_mu0p5_ang{ALPHA}"
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


def render_snapshot(xyz_m, types, radii_m, box,
                    px_per_mm=90, light=(-0.40, 0.40, 0.85),
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


# ---------------- gather snapshots per (rl, phi) -----------------
print("computing peak-RBNE snapshots at alpha=30 deg, two large diameters...")
records = {}                                          # records[(rl_tag, phi)] = dict
for rl_tag, dL_um in RL_TAGS:
    for phi in PHIS:
        cdir = case_dir(phi, rl_tag)
        t, dz = compute_dzbar(cdir)
        i_min = int(np.argmin(dz))
        t_min, dz_min = float(t[i_min]), float(dz.min())
        xyz, typ, rad, t_act, bb = load_frame(cdir, t_min)
        img, ext = render_snapshot(xyz, typ, rad, bb, px_per_mm=90)
        records[(rl_tag, phi)] = dict(
            phi=phi, dL_um=dL_um, t_act=t_act, dz_min=dz_min,
            t=t, dz=dz, t_end=float(t[-1]),
            img=img, ext=ext)
        print(f"  d_L={dL_um}um  phi={phi}  "
              f"t_min={t_min:.2f}s  dz_min={dz_min:+.3f} mm")

# ---------------- figure ------------------------------------------
NCOL = len(PHIS)
NROW = len(RL_TAGS) + 1                  # 2 snapshot rows + 1 curve row
# match Fig 1's 13-inch render width so font size on the printed page is
# the same as Fig 1
fig = plt.figure(figsize=(13.0, 6.4))
gs  = fig.add_gridspec(
    nrows=NROW, ncols=NCOL,
    height_ratios=[1.0, 1.0, 1.55],
    hspace=0.34, wspace=0.05,
    left=0.07, right=0.99, top=0.94, bottom=0.075,
)

letters = [["A","B","C","D","E"],
           ["F","G","H","I","J"],
           ["K","L","M","N","O"]]
EPS_MM = 0.10                             # mixed-band threshold

for i, (rl_tag, dL_um) in enumerate(RL_TAGS):
    for j, phi in enumerate(PHIS):
        rec = records[(rl_tag, phi)]
        ax = fig.add_subplot(gs[i, j])
        ax.imshow(rec["img"], extent=rec["ext"], origin="lower",
                  interpolation="bilinear")
        Lx_mm, H_mm = rec["ext"][1], rec["ext"][3]
        ax.add_patch(Rectangle((0, 0), Lx_mm, H_mm,
                               fill=False, edgecolor="0.15", linewidth=1.0))
        ax.plot([0, Lx_mm], [0, 0], color="k", lw=1.4)
        ax.set_xlim(0, Lx_mm); ax.set_ylim(0, H_mm)
        ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values(): s.set_visible(False)
        # column header (only on top row)
        if i == 0:
            ax.set_title(rf"$\phi = {phi:.2f}$", fontsize=18, pad=6)
        # row label on the left of the first column
        if j == 0:
            ax.text(-0.06, 0.5,
                    rf"$d_L = {dL_um}\,\mu$m",
                    transform=ax.transAxes, fontsize=17,
                    rotation=90, ha="right", va="center", color="0.10")
            # 1 mm scale bar (only on first cell of each row)
            ax.plot([0.6, 1.6], [0.32, 0.32], color="k", lw=2.4,
                    solid_capstyle="butt")
            ax.text(1.1, 0.10, "1 mm", ha="center", va="top",
                    fontsize=14, color="k")
        ax.text(0.985, 0.93,
                rf"$\Delta\bar{{z}}_{{\min}}={rec['dz_min']:+.2f}$ mm",
                transform=ax.transAxes, fontsize=15.5, color="0.10",
                ha="right", va="top",
                bbox=dict(boxstyle="round,pad=0.22", fc="white",
                          ec="0.55", lw=0.5))
        ax.text(0.5, -0.10, rf"$t = {rec['t_act']:.2f}$ s",
                transform=ax.transAxes, fontsize=16.5, color="0.10",
                ha="center", va="top")
        # panel letter outside upper-left; row 0 has a title above so use a
        # slightly larger y-offset to clear it
        y_off = 1.06 if i == 0 else 1.02
        ax.text(-0.015, y_off, letters[i][j], transform=ax.transAxes,
                fontsize=19, fontweight="bold", ha="right", va="bottom",
                color="0.10")

# -------- row 3: Delta zbar(t) overlay (both d_L) ----------
t_xmax = max(rec["t_end"] for rec in records.values())
dz_max_mag = max(max(abs(rec["dz"].min()), abs(rec["dz"].max()))
                 for rec in records.values())
dz_lim = 1.10 * dz_max_mag

for j, phi in enumerate(PHIS):
    ax3 = fig.add_subplot(gs[2, j])
    rec600  = records[("rl300um", phi)]      # d_L = 600 um
    rec1200 = records[("rl600um", phi)]      # d_L = 1200 um
    # background shading from the d_L=1200um curve (deeper RBNE so it bounds the band)
    ax3.fill_between(rec1200["t"], rec1200["dz"], 0,
                     where=(rec1200["dz"] <  0), color=C_L, alpha=0.18,
                     linewidth=0)
    ax3.fill_between(rec1200["t"], rec1200["dz"], 0,
                     where=(rec1200["dz"] >= 0), color=C_S, alpha=0.18,
                     linewidth=0)
    ax3.plot(rec1200["t"], rec1200["dz"], color="k", lw=1.6,
             label=r"$d_L = 1200\,\mu$m")
    ax3.plot(rec600["t"],  rec600["dz"],  color="0.30", lw=1.2,
             ls="--", label=r"$d_L = 600\,\mu$m")
    ax3.axhline(0,        color="k",    lw=0.6)
    ax3.axhline(+EPS_MM,  color="0.55", ls=":", lw=0.7)
    ax3.axhline(-EPS_MM,  color="0.55", ls=":", lw=0.7)
    # mark each row's snapshot time on its own curve
    for rec, color in [(rec1200, "k"), (rec600, "0.30")]:
        t_snap = float(rec["t"][int(np.argmin(rec["dz"]))])
        ax3.axvline(t_snap, color=color, lw=0.6, alpha=0.30)
        ax3.scatter([t_snap], [rec["dz"].min()],
                    s=36, marker="o", facecolor="white",
                    edgecolor=color, linewidth=1.0, zorder=5)
    ax3.set_xlim(0, t_xmax)
    ax3.set_ylim(-dz_lim, dz_lim)
    ax3.grid(alpha=0.18)
    ax3.set_xlabel("time (s)", fontsize=18)
    if j == 0:
        ax3.set_ylabel(r"$\Delta\bar{z} = \bar{z}_L - \bar{z}_S$  (mm)",
                       fontsize=18)
    else:
        ax3.set_yticklabels([])
    ax3.text(-0.015, 1.02, letters[2][j], transform=ax3.transAxes,
             fontsize=19, fontweight="bold", ha="right", va="bottom",
             color="0.10")
    if j == 0:
        ax3.legend(loc="upper left", frameon=False, fontsize=15,
                   bbox_to_anchor=(0.0, 0.92))

out_pdf = os.path.join(OUT, "figure08.pdf")
out_png = os.path.join(OUT, "figure08.png")
fig.savefig(out_pdf, bbox_inches="tight")
fig.savefig(out_png, bbox_inches="tight", dpi=180)
print(f"saved {out_pdf}")
print(f"saved {out_png}")
