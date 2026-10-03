#!/usr/bin/env python3
"""Compute the analytic Stokes-flow force for the equal-settling control.

For an isolated sphere of diameter d in an unbounded Stokes flow,

    V_St = (rho_p-rho_f) g_perp d^2 / (18 mu_f).

The small particle is left unforced. A constant upward force is applied to the
large particle so that its theoretical Stokes terminal speed equals that of the
small particle. The required force is

    F_add,L = W'_L - 3 pi mu_f d_L V_St,S
            = W'_L [1 - (d_S/d_L)^2],

where

    W'_L = (rho_p-rho_f) (pi d_L^3/6) g_perp.

Thus the force does not require a dilute or single-particle simulation and, for
equal particle densities, is independent of viscosity. Positive force is upward
(+z). This is a theoretical free-space Stokes control; residual differential
settling in the finite-concentration CFD-DEM mixture is an output, not something
to tune away.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


DEFAULT_RHO_P = 2470.0
DEFAULT_RHO_F = 971.0  # density actually used by the coupled CFD field
DEFAULT_NU = 1.2e-3    # 1200 cSt = 1.2e-3 m^2/s
DEFAULT_G = 9.796
DEFAULT_ANGLE_DEG = 30.0
DEFAULT_D_SMALL_UM = 200.0
DEFAULT_D_LARGE_UM = 600.0


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rho-p", type=float, default=DEFAULT_RHO_P, help="particle density [kg/m^3]")
    parser.add_argument("--rho-f", type=float, default=DEFAULT_RHO_F, help="fluid density [kg/m^3]")
    parser.add_argument("--nu", type=float, default=DEFAULT_NU, help="fluid kinematic viscosity [m^2/s]")
    parser.add_argument("--g", type=float, default=DEFAULT_G, help="gravity magnitude [m/s^2]")
    parser.add_argument("--angle", type=float, default=DEFAULT_ANGLE_DEG, help="incline angle [deg]")
    parser.add_argument("--g-perp", type=float, help="override positive wall-normal gravity [m/s^2]")
    parser.add_argument("--d-small-um", type=float, default=DEFAULT_D_SMALL_UM)
    parser.add_argument("--d-large-um", type=float, default=DEFAULT_D_LARGE_UM)
    parser.add_argument("--value-only", action="store_true", help="print only F_add,L [N]")
    parser.add_argument("--write-json", type=Path, help="optional output JSON path")
    args = parser.parse_args()

    if args.rho_p <= args.rho_f:
        raise SystemExit("rho_p must exceed rho_f for downward gravitational settling")
    for name, value in (("nu", args.nu), ("g", args.g), ("d_small_um", args.d_small_um),
                        ("d_large_um", args.d_large_um)):
        if value <= 0:
            raise SystemExit(f"{name} must be positive, got {value}")
    if args.d_large_um <= args.d_small_um:
        raise SystemExit("d_large_um must exceed d_small_um for this equal-size-material control")

    g_perp = args.g_perp if args.g_perp is not None else args.g * math.cos(math.radians(args.angle))
    if g_perp <= 0:
        raise SystemExit(f"g_perp must be positive, got {g_perp}")

    d_small = args.d_small_um * 1e-6
    d_large = args.d_large_um * 1e-6
    delta_rho = args.rho_p - args.rho_f
    mu = args.rho_f * args.nu

    v_small_stokes = delta_rho * g_perp * d_small**2 / (18.0 * mu)
    v_large_stokes = delta_rho * g_perp * d_large**2 / (18.0 * mu)
    volume_large = math.pi * d_large**3 / 6.0
    buoyant_weight_large = delta_rho * volume_large * g_perp

    size_ratio_squared = (d_small / d_large) ** 2
    addforce_large_z = buoyant_weight_large * (1.0 - size_ratio_squared)
    remaining_downward_force = buoyant_weight_large - addforce_large_z
    stokes_drag_large_at_target = 3.0 * math.pi * mu * d_large * v_small_stokes

    record = {
        "definition": "analytic free-space Stokes equal-settling control",
        "rho_p_kg_per_m3": args.rho_p,
        "rho_f_kg_per_m3": args.rho_f,
        "nu_m2_per_s": args.nu,
        "mu_Pa_s": mu,
        "g_m_per_s2": args.g,
        "angle_deg": args.angle,
        "g_perp_m_per_s2": g_perp,
        "d_small_m": d_small,
        "d_large_m": d_large,
        "V_small_Stokes_m_per_s": v_small_stokes,
        "V_large_Stokes_unforced_m_per_s": v_large_stokes,
        "unforced_speed_ratio": v_large_stokes / v_small_stokes,
        "buoyant_weight_large_N": buoyant_weight_large,
        "addforce_large_z_N": addforce_large_z,
        "addforce_small_z_N": 0.0,
        "remaining_downward_force_large_N": remaining_downward_force,
        "stokes_drag_large_at_target_N": stokes_drag_large_at_target,
        "force_fraction_of_buoyant_weight": addforce_large_z / buoyant_weight_large,
        "closed_form": "F_add,L = Wprime_L * (1 - (d_small/d_large)^2)",
        "interpretation": (
            "Do not fit or retune this force using finite-concentration mixture velocities, "
            "mean-height drift, or segregation flux."
        ),
    }

    # Algebraic consistency check of the two equivalent expressions.
    scale = max(abs(remaining_downward_force), abs(stokes_drag_large_at_target), 1e-300)
    if abs(remaining_downward_force - stokes_drag_large_at_target) / scale > 1e-12:
        raise RuntimeError("internal Stokes force-balance check failed")

    if args.write_json is not None:
        args.write_json.parent.mkdir(parents=True, exist_ok=True)
        args.write_json.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    if args.value_only:
        print(f"{addforce_large_z:.15g}")
        return

    print(json.dumps(record, indent=2, sort_keys=True))
    print("\nLIGGGHTS setting:")
    print("fix                 addsmall gsmall addforce 0.0 0.0 0.0")
    print(f"fix                 addlarge glarge addforce 0.0 0.0 {addforce_large_z:.15g}")
    print("\nGenerator argument:")
    print(f"--addforce-large {addforce_large_z:.15g}")


if __name__ == "__main__":
    main()
