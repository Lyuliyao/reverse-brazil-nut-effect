"""
Figure 6 — canonical three-phase evolution (RBNE -> mixed -> BNE),
with horizontally averaged fluid-velocity profiles near t1 and t2.

Layout
  Row 1 : A-D  x-z side-view OVITO-style snapshots (well mixed, RBNE phase,
               mixed, BNE phase).
  Row 2 : E    Delta z(t) = z_L_bar - z_S_bar with phase shading, snapshot
               markers, and dashed lines at the labelled zero-crossing
               times t1 and t2;
          F,G  horizontally averaged fluid streamwise velocity profiles
               <u_x>(z) from the CFD side, at the CFD output times nearest
               to the two zero crossings t1 and t2 (fluid fields are written
               every 0.5 s: t = 1.5 s for t1 = 1.30 s, t = 28.5 s for
               t2 = 28.69 s). Both panels share the same velocity axis so
               the flow acceleration between t1 and t2 is directly visible.

Case: canonical baseline cor0p1 mu0p5 (phi=0.1, alpha=30 deg,
      d_L=600 um, d_S=200 um) — the case shown in the paper.
"""
import os, json, re, datetime
import numpy as np
import pandas as pd
import h5py
from pathlib import Path
from bundle_io import REPO, CASES, OUTPUT, load_snapshot, read_height_series, read_fluid_profile
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.patches import Rectangle

# JFM-style typography
mpl.rcParams.update({
    "font.family":      "STIXGeneral",
    "mathtext.fontset": "stix",
    "axes.titlesize":   20,
    "axes.labelsize":   20,
    "xtick.labelsize":  18,
    "ytick.labelsize":  18,
    "legend.fontsize":  18,
})

CASES_ROOT = str(CASES)
OUT = str(OUTPUT)
os.makedirs(OUT, exist_ok=True)

COR_STR, MU_STR = "cor0p1", "mu0p5"      # canonical case used in the paper

DEM_DT = 1.0e-7
NX, NY, NZ = 32, 32, 16                  # CFD mesh (blockMeshDict)
LZ_M = 0.005

# species colours shared by Figures 6-10: royal blue (small), warm coral (large)
C_S = np.array([0.27, 0.41, 0.78])   # small  royal blue
C_L = np.array([0.93, 0.55, 0.40])   # large  warm coral / peach


def load_frame(seg_dir, seg_label, t_target):
    return load_snapshot(Path(seg_dir).parents[2], t_target)


# ---------------- fluid velocity profiles from the CFD side ----------------
def read_openfoam_U(path):
    """Parse an ascii OpenFOAM volVectorField internalField."""
    txt = open(path).read()
    m = re.search(r"internalField\s+nonuniform\s+List<vector>\s*\n(\d+)\s*\n\(",
                  txt)
    n = int(m.group(1))
    body = txt[m.end():]
    vals = re.findall(r"\(([-0-9.eE+ ]+)\)", body[:body.index("\n)\n")])
    U = np.array([[float(x) for x in v.split()] for v in vals[:n]])
    assert U.shape == (n, 3), U.shape
    return U


def fluid_profile(cfd_dir, t_target):
    return read_fluid_profile(Path(cfd_dir).parent, t_target)


# ------------------- OVITO-style sphere rasterizer ------------------
def render_snapshot(xyz_m, types, radii_m, box,
                    px_per_mm=90, light=(-0.40, 0.40, 0.85),
                    ambient=0.40, diff_w=0.50, spec_w=0.40, spec_pow=18,
                    bg=(1.0, 1.0, 1.0)):
    """Render an x-z view (camera looking along -y) using Phong shading
    with depth buffering. Returns an RGB image."""
    x_mm = xyz_m[:, 0] * 1e3
    y_mm = xyz_m[:, 1] * 1e3        # depth
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

    order = np.argsort(y_mm)        # ascending: back to front
    for i in order:
        rp_x = r_mm[i] * px_per_mm
        rp_y = r_mm[i] * px_per_mm
        cx = (x_mm[i] - box[0, 0]*1e3) * px_per_mm
        cy = (z_mm[i] - box[2, 0]*1e3) * px_per_mm
        x0 = max(0, int(np.floor(cx - rp_x)) - 1)
        x1 = min(W, int(np.ceil (cx + rp_x)) + 1)
        y0 = max(0, int(np.floor(cy - rp_y)) - 1)
        y1 = min(H, int(np.ceil (cy + rp_y)) + 1)
        if x0 >= x1 or y0 >= y1:
            continue
        yy, xx = np.mgrid[y0:y1, x0:x1].astype(np.float32)
        nx = (xx - cx) / rp_x
        ny = (yy - cy) / rp_y
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


def find_zero_crossings(t, gap):
    """Linear-interpolated zero crossings of gap(t)."""
    sign = np.sign(gap)
    idx  = np.where(np.diff(sign) != 0)[0]
    zc   = []
    for i in idx:
        x0, x1 = t[i], t[i+1]
        y0, y1 = gap[i], gap[i+1]
        if y1 != y0:
            zc.append(x0 - y0 * (x1 - x0) / (y1 - y0))
    return np.array(zc)


def make_figure():
    case_name = (
        f"bne_vf0p1_rs100um_rl300um_rho2470_{COR_STR}_{MU_STR}_ang30"
    )
    case_dir = os.path.join(CASES_ROOT, case_name)
    seg_dir  = os.path.join(case_dir, "DEM", "post", "forcebalance_segments")
    cfd_dir  = os.path.join(case_dir, "CFD")

    # ---- delta z(t) ----
    ts = pd.read_csv(os.path.join(case_dir, "time_series.tsv"),
                     sep="\t")
    ts = ts.sort_values("time").reset_index(drop=True)
    t   = ts["time"].to_numpy()
    gap = ts["z_gap_large_minus_small"].to_numpy() * 1e3   # mm

    # phase boundaries from gap zero crossings
    zc    = find_zero_crossings(t, gap)
    i_min = int(np.argmin(gap))
    t_min = t[i_min]
    zc_lo = zc[zc < t_min]
    zc_hi = zc[zc > t_min]
    t1 = float(zc_lo[-1]) if len(zc_lo) else None
    t2 = float(zc_hi[0])  if len(zc_hi) else None
    if t1 is None or t2 is None:
        raise RuntimeError("could not locate both zero crossings t1/t2")

    # ---- fluid velocity profiles near t1 and t2 ----
    prof = {}
    for key, tc in (("F", t1), ("G", t2)):
        t_used, zcen, ux, uz = fluid_profile(cfd_dir, tc)
        prof[key] = dict(t_target=tc, t_used=t_used, z=zcen, ux=ux, uz=uz)
        print(f"  fluid profile {key}: CFD t = {t_used:g} s "
              f"(nearest to {tc:.3f} s), max <u_x> = {ux.max():.2f} mm/s")
    u_max = max(p["ux"].max() for p in prof.values())
    u_lim = float(np.ceil(u_max / 5.0) * 5.0)

    # ---- snapshot specs ----
    # labels name the states directly (RBNE phase / mixed / BNE phase);
    # the "mixed" snapshot is taken just before the second zero
    # crossing t2 = 28.69 s, not after it
    snap_specs = [
        ("t0_to_1",   0.005,  "well mixed"),
        ("t9_to_10",  10.0,   "RBNE phase"),
        ("t27_to_28", 28.0,   "mixed"),
        ("t39_to_40", t[-1],  "BNE phase"),
    ]

    # ---- figure ----
    fig = plt.figure(figsize=(13.0, 6.3))
    gs  = fig.add_gridspec(
        nrows=2, ncols=1,
        height_ratios=[1.0, 1.25],
        hspace=0.34,
        left=0.07, right=0.99, top=0.95, bottom=0.115,
    )
    gs1 = gs[0].subgridspec(1, 4, wspace=0.06)
    gs2 = gs[1].subgridspec(1, 3, width_ratios=[2.35, 0.80, 0.80],
                            wspace=0.26)

    snap_times = []
    for j, (seg, t_target, phase_lbl) in enumerate(snap_specs):
        ax = fig.add_subplot(gs1[0, j])
        xyz, typ, rad, t_act, bb = load_frame(seg_dir, seg, t_target)
        snap_times.append(t_act)

        img, extent = render_snapshot(xyz, typ, rad, bb, px_per_mm=90)
        ax.imshow(img, extent=extent, origin="lower", interpolation="bilinear")

        Lx_mm = extent[1]
        H_mm  = extent[3]
        ax.add_patch(Rectangle((0, 0), Lx_mm, H_mm,
                               fill=False, edgecolor="0.15", linewidth=1.0))
        ax.plot([0, Lx_mm], [0, 0], color="k", lw=1.6)

        ax.set_xlim(0, Lx_mm)
        ax.set_ylim(0, H_mm)
        ax.set_aspect("equal")
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)

        ax.text(0.5, -0.10, f"$t = {t_act:.2f}$ s", transform=ax.transAxes,
                fontsize=20, color="0.10", ha="center", va="top")

        if j == 0:
            sb_x0, sb_y0, sb_len = 0.35, 0.10, 1.0
            ax.plot([sb_x0, sb_x0 + sb_len],
                    [sb_y0 * H_mm, sb_y0 * H_mm],
                    color="k", lw=2.6, solid_capstyle="butt",
                    path_effects=[pe.Stroke(linewidth=5.0,
                                            foreground="white"),
                                  pe.Normal()], zorder=6)
            ax.text(sb_x0 + sb_len/2, sb_y0 * H_mm + 0.14,
                    "1 mm", ha="center", va="bottom", fontsize=18,
                    color="k", zorder=6,
                    path_effects=[pe.Stroke(linewidth=3.0,
                                            foreground="white"),
                                  pe.Normal()])

        ax.text(-0.015, 1.04, chr(ord('A')+j), transform=ax.transAxes,
                fontsize=20, fontweight="bold", fontfamily="STIXGeneral",
                ha="left", va="baseline", color="0.10")
        ax.text(0.07, 1.04, phase_lbl, transform=ax.transAxes,
                fontsize=20, fontweight="normal", fontfamily="STIXGeneral",
                ha="left", va="baseline", color="0.10")

    # ---- row 2, panel E: delta z(t) ----
    ax_b = fig.add_subplot(gs2[0, 0])
    ax_b.fill_between(t, gap, 0, where=(gap < 0),  color=C_L, alpha=0.28,
                      label="RBNE (large below small)", linewidth=0)
    ax_b.fill_between(t, gap, 0, where=(gap >= 0), color=C_S, alpha=0.28,
                      label="BNE (large above small)", linewidth=0)
    ax_b.plot(t, gap, color="k", lw=1.6)
    ax_b.axhline(0, color="k", lw=0.6)

    for tx in (t1, t2):
        ax_b.axvline(tx, ls="--", color="0.35", lw=1.0)
    ymin, ymax = ax_b.get_ylim()
    ax_b.text(t1+0.5, ymax*0.84, rf"$t_1 = {t1:.2f}$ s",
              ha="left", fontsize=18)
    ax_b.text(t2+0.5, ymax*0.84, rf"$t_2 = {t2:.2f}$ s",
              ha="left", fontsize=18)
    # region labels name the states directly (RBNE phase / BNE phase);
    # the narrow mixed crossover around t2 is identified by snapshot C
    ax_b.text(0.5*(t1+t2), 0.55*np.min(gap), "RBNE phase",
              ha="center", color="#7a1f15", fontsize=19)
    ax_b.text(0.5*(t2+t[-1])+1.4, ymax*0.26, "BNE phase",
              ha="center", color="#1c4e74", fontsize=19)

    for k, ts_snap in enumerate(snap_times):
        ax_b.scatter([ts_snap], [np.interp(ts_snap, t, gap)],
                     s=55, marker="o", clip_on=False,
                     facecolor="white", edgecolor="k", linewidth=1.2, zorder=5)
        is_last = (k == len(snap_times) - 1)
        dx = -8 if is_last else 6
        dy = -14 if is_last else 8
        ha = "right" if is_last else "left"
        va = "top" if is_last else "baseline"
        ax_b.annotate(chr(ord('A')+k),
                      xy=(ts_snap, np.interp(ts_snap, t, gap)),
                      xytext=(dx, dy), textcoords="offset points",
                      fontsize=19, fontweight="bold", ha=ha, va=va)

    ax_b.set_xlim(0, t[-1])
    ax_b.set_xlabel("time (s)", fontsize=20)
    ax_b.set_ylabel(
        r"$\Delta \bar z = \bar z_{\mathrm{L}} - \bar z_{\mathrm{S}}$  (mm)",
        fontsize=20)
    ax_b.legend(loc="upper center", bbox_to_anchor=(0.46, 1.005),
                frameon=False, fontsize=17,
                handlelength=1.4, handletextpad=0.6)
    ax_b.grid(alpha=0.20)
    ax_b.text(-0.015, 1.03, "E", transform=ax_b.transAxes,
              ha="left", va="bottom", fontsize=20, fontweight="bold",
              color="0.10")

    # ---- row 2, panels F/G: fluid velocity profiles near t1 and t2 ----
    axes_fg = []
    for col, key in ((1, "F"), (2, "G")):
        p  = prof[key]
        ax = fig.add_subplot(gs2[0, col])
        ax.plot(p["ux"], p["z"], color="k", lw=2.0, marker="o",
                markersize=4.5, markerfacecolor="white",
                markeredgewidth=1.0, clip_on=False, zorder=4)
        ax.fill_betweenx(p["z"], 0, p["ux"], color="0.60", alpha=0.18,
                         linewidth=0)
        ax.set_xlim(0, u_lim)
        ax.set_xticks(np.arange(0, u_lim + 1, 10))
        ax.set_ylim(0, LZ_M * 1e3)
        ax.set_xlabel(r"$\langle u_x \rangle$ (mm/s)", fontsize=20)
        if col == 1:
            ax.set_ylabel(r"$z$ (mm)", fontsize=20)
        else:
            ax.set_yticklabels([])
        ax.tick_params(labelsize=18)
        ax.grid(alpha=0.20)

        ax.text(-0.03, 1.03, key, transform=ax.transAxes,
                fontsize=20, fontweight="bold", fontfamily="STIXGeneral",
                ha="left", va="bottom", color="0.10")
        sub = "1" if key == "F" else "2"
        ax.text(0.13, 1.03, rf"fluid, $t \approx t_{sub}$",
                transform=ax.transAxes, fontsize=19, ha="left",
                va="bottom", color="0.10")
        ax.text(0.94, 0.06, rf"$t = {p['t_used']:g}$ s",
                transform=ax.transAxes, fontsize=17, ha="right",
                va="bottom", color="0.25")
        axes_fg.append(ax)

    out_pdf = os.path.join(OUT, "figure06.pdf")
    out_png = os.path.join(OUT, "figure06.png")
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_png, bbox_inches="tight", dpi=200)
    plt.close(fig)

    meta = {
        "generated": datetime.datetime.now().isoformat(timespec="seconds"),
        "script": os.path.abspath(__file__),
        "case": case_name,
        "t1_s": t1, "t2_s": t2,
        "panels_FG": "horizontally averaged fluid streamwise velocity "
                     "<u_x>(z) from OpenFOAM U fields (average over the "
                     "periodic x and y directions), at the CFD output "
                     "times nearest to t1 and t2",
        "fluid_times_s": {k: prof[k]["t_used"] for k in ("F", "G")},
        "max_ux_mm_s": {k: float(prof[k]["ux"].max()) for k in ("F", "G")},
        "velocity_axis_mm_s": u_lim,
        "outputs": [out_pdf, out_png],
    }
    with open(os.path.join(OUT, "figure06_metadata.json"),
              "w") as fjs:
        json.dump(meta, fjs, indent=2)

    print(f"saved {out_png}")
    print(f"  t1 = {t1:.3f} s   t_min(gap) = {t_min:.3f} s   t2 = {t2:.3f} s")
    print(f"  snapshot times: {[f'{x:.3f}' for x in snap_times]}")


if __name__ == "__main__":
    make_figure()
