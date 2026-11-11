#!/usr/bin/env python3
"""Heterogeneous Gc field for GEOS: correlated background + local weak spots.

Background (same recipe as the porosity/surface-area notebook in
playground/ChemoMechanics/Uno/surface_area_porosity_field_generator.ipynb):
1. Gaussian random field with Gaussian covariance (FFT spectral simulation)
2. Normalize to [0,1] via the Gaussian CDF (preserves spatial structure)
3. Power transform f^alpha so the mean lands at bg_mean
4. Rescale linearly to [bg_min, bg_max]

Weak spots: placed ON the bottom (inlet) boundary at y = 0, at random x
positions, imprinted on the background by blending log10(Gc) to log10(gcmin)
inside a Gaussian envelope of std sigma (x) and sigma_y (y); sigma_y > sigma
elongates the spots into weak channels along y. Centers are snapped to table
nodes so they hit gcmin exactly. With --spot-ymax > 0 the spots are instead
placed at random (x, y) inside the strip, spot_ymin <= y <= spot_ymax, with no
copies on the boundaries. --no-top-spots drops the copies on the top boundary
(y = Y_MAX), so each inlet spot counts once toward --nspots.

Optional inlet floor: the background is raised to at least inlet_floor within
inlet_depth of the inlet (y = 0) and tapered (cosine in log10 Gc) to the
background over a further inlet_taper, before the weak spots are imprinted.

Table grid spacing is 2 mm — no finer than the 1 mm uniform XML mesh.
Writes a GEOS TableFunction voxel table (Fortran order: x fastest, then y,
then z) into ./gcField/ plus a preview PNG and the spot-center list.

Usage: python3 generate_gc_field.py [--gcmin 1e-3] [--seed 42]
                                    [--bg-min 1.0] [--bg-max 3.0] [--bg-mean 1.9]
                                    [--len-x 0.02] [--len-y 0.02] [--angle 0]
                                    [--nspots 6] [--sigma 0.008] [--sigma-y 0.02]
                                    [--inlet-floor 10] [--inlet-depth 0.01] [--inlet-taper 0.01]
                                    [--spot-ymin 0.01] [--spot-ymax 0.05] [--no-top-spots]
                                    [--outdir gcField]
"""
import argparse
import math
import pathlib

import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("--gcmin", type=float, default=1e-2, help="Gc at weak-spot centers")
parser.add_argument("--bg-min", type=float, default=1.0, help="background Gc minimum")
parser.add_argument("--bg-max", type=float, default=20.0, help="background Gc maximum")
parser.add_argument("--bg-mean", type=float, default=10.0, help="target background Gc mean")
parser.add_argument("--len-x", type=float, default=0.005, help="background correlation length in x [m]")
parser.add_argument("--len-y", type=float, default=0.005, help="background correlation length in y [m]")
parser.add_argument("--angle", type=float, default=0.0, help="background anisotropy rotation [deg]")
parser.add_argument("--nspots", type=int, default=7, help="number of weak spots on the bottom boundary")
parser.add_argument("--sigma", type=float, default=0.001, help="std dev of the weak-spot Gaussian falloff [m]")
parser.add_argument("--sigma-y", type=float, default=0.003, help="std dev of the weak-spot falloff in y [m] (default: sigma)")
parser.add_argument("--inlet-floor", type=float, default=20.0, help="minimum background Gc near the inlet (0 = off)")
parser.add_argument("--inlet-depth", type=float, default=0.005, help="depth of the inlet floor from y = 0 [m]")
parser.add_argument("--inlet-taper", type=float, default=0.005, help="taper length below the inlet floor [m]")
parser.add_argument("--spot-ymin", type=float, default=0.0, help="min depth of interior spots [m] (with --spot-ymax)")
parser.add_argument("--spot-ymax", type=float, default=0.0, help="max depth of interior spots [m] (0 = spots on y = 0 and y = Y_MAX)")
parser.add_argument("--no-top-spots", action="store_true", help="no weak-spot copies on the top boundary y = Y_MAX")
parser.add_argument("--seed", type=int, default=30, help="RNG seed for reproducibility")
parser.add_argument("--outdir", type=str, default=None, help="output directory (default: ./gcField next to this script)")
args = parser.parse_args()
if args.sigma_y is None:
    args.sigma_y = args.sigma

assert args.bg_min < args.bg_mean < args.bg_max

# Domain of the InternalMesh in the main XMLs
X_MIN, X_MAX = 0.0, 1.0
Y_MIN, Y_MAX = -0.0001, 0.2
Z_MIN, Z_MAX = 0.0, 0.005

DX = 0.001  # table spacing: equal to the mesh size (the allowed minimum),
            # needed to resolve the small weak spots

# Spot centers sit on the bottom boundary (y = 0), spread along x
MARGIN_X = 0.05
MIN_SPACING = 0.08  # min center distance along the inlet [m]

rng = np.random.default_rng(args.seed)

# ---------------------------------------------------------------- grid
x = np.arange(0.0, X_MAX + DX / 2, DX)
y = np.arange(0.0, Y_MAX + DX / 2, DX)
nx, ny = x.size, y.size
XX, YY = np.meshgrid(x, y, indexing="ij")

# ------------------- background GRF (FFT, Gaussian covariance) -> [bg_min, bg_max]
cos_t, sin_t = np.cos(np.radians(args.angle)), np.sin(np.radians(args.angle))
mx = int(2 ** np.ceil(np.log2(nx + int(6 * args.len_x / DX))))
my = int(2 ** np.ceil(np.log2(ny + int(6 * args.len_y / DX))))
fx = (np.arange(mx) - mx // 2) * DX
fy = (np.arange(my) - my // 2) * DX
FX, FY = np.meshgrid(fx, fy, indexing="ij")
u = FX * cos_t + FY * sin_t
v = -FX * sin_t + FY * cos_t
kernel = np.exp(-((u / args.len_x) ** 2 + (v / args.len_y) ** 2))
spectrum = np.fft.fft2(np.fft.ifftshift(kernel)).real
spectrum[spectrum < 0] = 0.0
noise = rng.standard_normal((mx, my))
raw = np.fft.ifft2(np.sqrt(spectrum) * np.fft.fft2(noise)).real[:nx, :ny]

erf = np.vectorize(math.erf)
f01 = 0.5 * (1.0 + erf((raw - raw.mean()) / (raw.std() * np.sqrt(2.0))))
mean_rel = (args.bg_mean - args.bg_min) / (args.bg_max - args.bg_min)
alpha = 1.0 / mean_rel - 1.0
base = args.bg_min + f01 ** alpha * (args.bg_max - args.bg_min)

# ------------------------------- optional tougher background near the inlet
if args.inlet_floor > 0.0:
    s_taper = np.clip((YY - args.inlet_depth) / max(args.inlet_taper, DX), 0.0, 1.0)
    w = 0.5 * (1.0 + np.cos(np.pi * s_taper))  # 1 within inlet_depth, 0 beyond the taper
    floor = 10.0 ** (w * np.log10(args.inlet_floor) + (1.0 - w) * np.log10(np.minimum(base, args.inlet_floor)))
    base = np.maximum(base, floor)

# ---------------------- weak spots on the bottom (inlet) boundary, y = 0
centers = []
attempts = 0
while len(centers) < args.nspots:
    attempts += 1
    if attempts > 200000:
        raise RuntimeError("could not place all spots; lower --nspots or MIN_SPACING")
    cx = rng.uniform(X_MIN + MARGIN_X, X_MAX - MARGIN_X)
    if args.spot_ymax > 0.0:
        # interior spots: random depth, spacing measured in 2D
        cy = rng.uniform(args.spot_ymin, args.spot_ymax)
        if all(math.hypot(cx - px, cy - py) >= MIN_SPACING for px, py in centers):
            centers.append((cx, cy))
    elif all(abs(cx - px) >= MIN_SPACING for px, _ in centers):
        centers.append((cx, 0.0))
        if not args.no_top_spots:
            centers.append((cx, Y_MAX))
# snap to table nodes so spot centers hit gcmin exactly
centers = np.round(np.array(centers) / DX) * DX

# circular weak spots: blend log10(Gc) to log10(gcmin) inside each envelope
envelope = np.zeros_like(base)
for cx, cy in centers:
    r2 = ((XX - cx) / args.sigma) ** 2 + ((YY - cy) / args.sigma_y) ** 2
    envelope = np.maximum(envelope, np.exp(-r2 / 2.0))
gc2d = 10.0 ** ((1.0 - envelope) * np.log10(base) + envelope * np.log10(args.gcmin))

# --------------------------------------- write GEOS table (Fortran order)
y_tab = np.concatenate(([Y_MIN], y))               # extra inlet-reservoir row
gc_tab = np.concatenate((gc2d[:, :1], gc2d), axis=1)
z_tab = np.array([Z_MIN, Z_MAX])
gc3d = np.repeat(gc_tab[:, :, None], z_tab.size, axis=2)
values = gc3d.ravel(order="F")

outdir = pathlib.Path(args.outdir) if args.outdir else pathlib.Path(__file__).parent / "gcField"
outdir.mkdir(exist_ok=True)
np.savetxt(outdir / "x.csv", x, fmt="%.6g")
np.savetxt(outdir / "y.csv", y_tab, fmt="%.6g")
np.savetxt(outdir / "z.csv", z_tab, fmt="%.6g")
np.savetxt(outdir / "gc.csv", values, fmt="%.6e")
np.savetxt(outdir / "spot_centers.txt", centers, fmt="%.6f",
           header=f"x y  (gcmin={args.gcmin}, bg=[{args.bg_min}, {args.bg_max}] "
                  f"mean {args.bg_mean}, len=({args.len_x}, {args.len_y}), "
                  f"angle={args.angle}, sigma=({args.sigma}, {args.sigma_y}), "
                  f"inlet_floor={args.inlet_floor} to y={args.inlet_depth}+{args.inlet_taper}, "
                  f"spot_y=[{args.spot_ymin}, {args.spot_ymax}], seed={args.seed})")

print(f"wrote {values.size} voxel values to {outdir}")
print(f"alpha = {alpha:.3f}")
print(f"background: min {base.min():.4f}, max {base.max():.4f}, mean {base.mean():.4f}")
print(f"Gc range in table: [{values.min():.3e}, {values.max():.3e}]")
print("spot centers (x, y):")
for cx, cy in centers:
    print(f"  ({cx:.4f}, {cy:.4f})")

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8, 4.5))
    pc = ax.pcolormesh(XX, YY, gc2d, vmin=0.0, vmax=args.bg_max,
                       cmap="Blues", shading="auto")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_aspect("equal")
    ax.set_ylim(0, 0.2)
    ax.set_title(f"Gc: bg [{args.bg_min}, {args.bg_max}] mean {args.bg_mean}, "
                 f"len {args.len_x} m, {args.nspots} weak spots @ {args.gcmin:g}",
                 fontsize=10)
    # hollow red boxes at the weak-spot centers, +-2 sigma wide
    from matplotlib.patches import Rectangle
    for cx, cy in centers:
        ax.add_patch(Rectangle((cx - 2 * args.sigma, cy - 2 * args.sigma_y), 4 * args.sigma, 4 * args.sigma_y,
                               fill=False, ec="red", lw=1.0))
    fig.colorbar(pc, ax=ax, label="Gc [J/m$^2$]")
    fig.tight_layout()
    fig.savefig(outdir / "gc_field_preview.png", dpi=150)
    print(f"preview: {outdir / 'gc_field_preview.png'}")
except ImportError:
    print("matplotlib not available; skipped preview plot")