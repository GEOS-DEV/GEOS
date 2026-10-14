#!/usr/bin/env python3
"""Generate the prescribed cooling field T(x, y, z, t) for the 2D lab-block wellbore deck.

The 2D model is plane strain, so there is no out-of-plane direction for the fluid to lose heat
into; the cooling has to be prescribed rather than solved.

PHYSICS. Before a fracture exists the fluid moves through POROUS rock, where fluid and grains are
in local thermal equilibrium. The thermal front is then set by retardation, not by conduction into
fracture walls:

    R_T  = (rho c)_bulk / (phi (rho c)_f)          thermal retardation, ~11 here
    r_f  = sqrt( Q t / (phi pi H) )                fluid front (volume balance)
    r_T  = r_f / sqrt( R_T )                       thermal front, lagging the fluid front

and the front is smeared by bulk conduction over 2 sqrt(alpha_b t), giving

    T = T_rock + (T_inj - T_rock) * 0.5 * erfc( (r - r_T) / (2 sqrt(alpha_b t)) )

NOTE the smear EXCEEDS the front position at every time here (29.5 mm vs 13.9 mm at t = 200 s), so
this is a smooth diffusive profile, NOT a sharp front. The 1D column work showed that a smooth
driver spreads precipitation almost uniformly while only a sharp front localizes it, so expect weak
localization at these conditions. That is a physical result of the lab scale, not a modelling
choice - do not sharpen the profile to get a nicer answer.

An earlier version of this script used a wall-conduction balance, rc = sqrt(rho_f c_f Q / 2 pi h).
That applies to flow in an OPEN FRACTURE, not to porous flow, and is wrong for the pre-fracture
stage this field has to cover.

SHIFTED FRONT (--front-offset L_c, optional --front-width delta). Places the front at a fixed
distance L_c from the wall with width delta, so precipitation sits downstream of the inlet rather than
at it. This is a sensitivity input for the plug position, not a prediction. delta defaults to the
conduction smear 2 sqrt(alpha_b t).

    python3 make_cooling_field.py --front-offset 0.025 --front-width 0.003 --out ../coolingField_Lc25_w3

The tables are read back and plotted to cooling_field.png in the output directory.
"""
import argparse
import math
import numpy as np
import os

# scipy is not available in this environment; erfc is in the stdlib via math.erfc.
_erfc = np.vectorize(math.erfc)

parser = argparse.ArgumentParser()
parser.add_argument("--front-offset", type=float, default=None, help="front distance from the wall [m]")
parser.add_argument("--front-width", type=float, default=None, help="front width delta [m]")
parser.add_argument("--out", default=os.path.dirname(os.path.abspath(__file__)), help="output directory")
args = parser.parse_args()

OUT = os.path.abspath(args.out)
os.makedirs(OUT, exist_ok=True)

RC_F, RC_R = 4.2e6, 2.2e6     # J/m3/K, fluid and rock volumetric heat capacity
LAM, PHI = 2.5, 0.05          # W/m/K rock conductivity, rock porosity
H = 0.001                     # m, modelled slice thickness

# Injection: 5 mL/min over a 50 mm interval, so the 1 mm slice takes 1/50 of it.
Q_TOTAL = 5.0e-6 / 60.0
INTERVAL = 0.050
Q = Q_TOTAL * (H / INTERVAL)

T_INJ, T_ROCK = 333.15, 298.15    # 60 C prepared solution into 25 C rock
R_BORE = 0.00635                  # m, meshed borehole radius (Source/Rock interface, = 12.7 mm dia)

RC_B = PHI * RC_F + (1.0 - PHI) * RC_R
R_T = RC_B / (PHI * RC_F)
AL_B = LAM / RC_B

# The profile is steep near the wall at early times, and a cell there interpolates between table
# nodes in both x and y, so every coordinate out to 16 mm gets 0.2 mm spacing (cells are 0.17 mm at the
# wall); beyond that the front is smeared over tens of mm and coarser spacing is enough. The early
# times resolve the conduction front, which moves most in the first second.
# A shifted front gets 0.5 mm spacing out to three widths past it.
if args.front_offset is None:
    _half = np.concatenate([np.arange(0.0, 0.016, 0.0002), np.arange(0.016, 0.040, 0.001), np.arange(0.040, 0.0781, 0.002)])
else:
    _width = args.front_width if args.front_width is not None else 2.0 * math.sqrt(LAM / (PHI * RC_F + (1.0 - PHI) * RC_R) * 200.0)
    _mid = min(max(R_BORE + args.front_offset + 3.0 * _width, 0.016), 0.0781)
    _half = np.concatenate([np.arange(0.0, 0.016, 0.0002), np.arange(0.016, _mid, 0.0005), np.arange(_mid, 0.0781, 0.002)])
X = np.concatenate([-_half[:0:-1], _half])
Y = X.copy()
Z = np.array([0.0])
T = np.array([0.0, 0.1, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 40.0, 70.0, 110.0, 160.0, 200.0])


def front_position(t):
    """Front distance from the wall: the retardation front, or the fixed --front-offset."""
    if args.front_offset is not None:
        return args.front_offset
    return np.sqrt(Q * t / (PHI * np.pi * H)) / np.sqrt(R_T)


def profile(r, t):
    """Temperature at radius r and time t.

    The borehole itself is held at the preparation temperature: it is continuously refilled with
    fresh 60 C solution, so it never cools. The front is therefore measured from the borehole WALL
    into the rock, not from the axis. Getting this wrong puts 60 C-saturated fluid into cells the
    temperature field calls 25 C, i.e. supersaturated at the injection point by construction.
    """
    tt = max(t, 1.0e-2)
    r_T = front_position(tt)
    smear = args.front_width if args.front_width is not None else 2.0 * np.sqrt(AL_B * tt)
    s = np.asarray(r) - R_BORE                      # distance into the rock
    temp = T_ROCK + (T_INJ - T_ROCK) * 0.5 * _erfc((s - r_T) / smear)
    return np.where(np.asarray(r) <= R_BORE, T_INJ, temp)


# Fortran order: x fastest, then y, z, t
vals = []
for ti in T:
    for _ in Z:
        for yi in Y:
            vals.append(profile(np.hypot(X, yi), ti))
vals = np.concatenate(vals)
assert vals.size == X.size * Y.size * Z.size * T.size

# NOTE: macOS is case-insensitive, so the time axis must not be called t.geos -- T.geos clobbers it.
for name, arr in [("x", X), ("y", Y), ("z", Z), ("time", T)]:
    np.savetxt(os.path.join(OUT, name + ".geos"), arr, fmt="%.7f")
np.savetxt(os.path.join(OUT, "T.geos"), vals, fmt="%.4f")

print(f"Q per slice = {Q:.4e} m3/s   R_T = {R_T:.2f}   alpha_b = {AL_B:.3e} m2/s")
print(f"grid {X.size} x {Y.size} x {Z.size} x {T.size} = {vals.size} values\n")
print(f"{'t[s]':>6} {'r_T':>9} {'smear':>9} {'T at wall':>11} {'T at 30mm':>11}")
for ti in T[1:]:
    r_T = front_position(max(ti, 1e-2))
    smear = args.front_width if args.front_width is not None else 2*np.sqrt(AL_B*ti)
    print(f"{ti:>6.0f} {r_T*1e3:>7.1f}mm {smear*1e3:>7.1f}mm "
          f"{profile(np.array([R_BORE]), ti)[0]:>10.1f}K {profile(np.array([0.030]), ti)[0]:>10.1f}K")


def plot_tables():
    """Plot the field as written, so the figure shows what GEOS reads, not the analytic profile."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    xs, ys, ts = (np.loadtxt(os.path.join(OUT, n + ".geos"), ndmin=1) for n in ("x", "y", "time"))
    field = np.loadtxt(os.path.join(OUT, "T.geos")).reshape(ts.size, ys.size, xs.size) - 273.15

    # Supersaturation of the 60 C-saturated feed, with ln K linear in 1/T between 25 and 60 C
    def omega(temp_c):
        inv = 1.0 / (temp_c + 273.15)
        slope = (7.0575 - 6.6495) / (1.0 / 333.15 - 1.0 / 298.15)
        return np.exp(slope * (1.0 / 333.15 - inv))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.8))

    # Along the crack: x = 0, y > 0
    i0 = np.argmin(np.abs(xs))
    pos = ys >= R_BORE
    s_mm = (ys[pos] - R_BORE) * 1e3
    ax1b = ax1.twinx()
    for ti, color in zip(ts[1:], plt.cm.viridis(np.linspace(0, 1, ts.size - 1))):
        it = np.where(ts == ti)[0][0]
        temp = field[it, pos, i0]
        ax1.plot(s_mm, temp, color=color, label=f"{ti:g} s")
        ax1b.plot(s_mm, omega(temp), color=color, ls=":", lw=1)
    ax1.set_xlabel("distance from borehole wall along the crack (mm)")
    ax1.set_ylabel("temperature (C), solid")
    ax1b.set_ylabel(r"$\Omega$ of the 60 C-saturated feed, dotted")
    ax1.legend(fontsize=7, ncol=2, title="time")
    ax1.grid(alpha=0.3)

    # Map at the latest time
    im = ax2.pcolormesh(xs * 1e3, ys * 1e3, field[-1], shading="auto", cmap="coolwarm")
    ax2.add_patch(plt.Circle((0, 0), R_BORE * 1e3, fill=False, color="k", lw=0.8))
    if args.front_offset is not None:
        ax2.add_patch(plt.Circle((0, 0), (R_BORE + args.front_offset) * 1e3, fill=False, color="k", ls="--", lw=0.8))
    ax2.set_aspect("equal")
    ax2.set_xlabel("x (mm)")
    ax2.set_ylabel("y (mm)")
    ax2.set_title(f"t = {ts[-1]:g} s")
    fig.colorbar(im, ax=ax2, label="temperature (C)")

    if args.front_offset is None:
        fig.suptitle("Retardation front")
    else:
        width = "2 sqrt(alpha t)" if args.front_width is None else f"{args.front_width*1e3:g} mm"
        fig.suptitle(f"Shifted front: L_c = {args.front_offset*1e3:g} mm, width = {width}")
    fig.tight_layout()
    path = os.path.join(OUT, "cooling_field.png")
    fig.savefig(path, dpi=150)
    print(f"\nplot: {path}")


plot_tables()
