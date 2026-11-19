#!/usr/bin/env python3
"""
Analytical solution for the dislocated reservoir with a vertical fault.

Reference: Frigo, Castelletto, Cusini et al., "Augmented Lagrangian Method for
thin faults in poromechanics", J. Comput. Phys. 561 (2026) 114988, Section 5.1.2.

The domain has a vertical fault at x=0 separating two offset reservoirs:
  - Left reservoir:  x < 0, z in [-b, a]  (pressurized)
  - Right reservoir:  x > 0, z in [-a, b]  (pressurized)

Parameters:
  a = inner reservoir half-extent (overlap region: |z| < a)
  b = outer reservoir half-extent (fault extends from -b to b)

Tangential traction (Eq. 48):
  t_T(z) = C/2 * ln( (z-a)^2*(z+a)^2 / ((z-b)^2*(z+b)^2) )

Tangential slip (Eq. 50):
              | 0              if z <= -b
              | -(z+b)         if -b < z <= -a
  g_T(z) = C/A * | (a-b)          if -a < z < a
              | (z-b)          if  a <= z < b
              | 0              if z >= b

where:
  C = (1-2*nu) * alpha * dp / (2*pi*(1-nu))     (Eq. 49)
  A = G / (2*pi*(1-nu))                          (Eq. 51)
  G = E / (2*(1+nu))  (shear modulus)
  alpha = Biot coefficient = 1 - K_drained / K_grain
  dp = pressure change in the reservoir (negative for depletion)

Usage:
  python generate_analytical.py
  python generate_analytical.py --npoints 5000 --zmin -500 --zmax 500
  python generate_analytical.py --dp -25e6 --a 75 --b 150  # explicit parameters
"""

import argparse
import numpy as np
import os


def compute_parameters(E, nu, K_grain, dp, a, b):
    """Compute derived parameters from material properties."""
    G = E / (2 * (1 + nu))
    K_drained = E / (3 * (1 - 2*nu))
    alpha = 1.0 - K_drained / K_grain

    C = (1 - 2*nu) * alpha * dp / (2 * np.pi * (1 - nu))
    A = G / (2 * np.pi * (1 - nu))

    params = {
        'E': E, 'nu': nu, 'G': G,
        'K_drained': K_drained, 'K_grain': K_grain,
        'alpha': alpha,
        'dp': dp, 'a': a, 'b': b,
        'C': C, 'A': A,
    }
    return params


def analytical_traction(z, params):
    """
    Tangential traction on the fault (Eq. 48).

    t_T(z) = C/2 * ln( (z-a)^2*(z+a)^2 / ((z-b)^2*(z+b)^2) )

    Has logarithmic singularities at z = +/-a and z = +/-b.
    """
    a = params['a']
    b = params['b']
    C = params['C']

    z = np.asarray(z, dtype=float)
    t_T = np.zeros_like(z)

    for i, zi in enumerate(z):
        # Avoid singularities at z = +/-a, +/-b
        if abs(abs(zi) - a) < 1e-12 or abs(abs(zi) - b) < 1e-12:
            t_T[i] = np.nan
        else:
            num = (zi - a)**2 * (zi + a)**2
            den = (zi - b)**2 * (zi + b)**2
            t_T[i] = C / 2 * np.log(num / den)

    return t_T


def analytical_slip(z, params):
    """
    Tangential slip on the fault (Eq. 50).

    Piecewise linear profile; returns absolute slip |g_T|.
    """
    a = params['a']
    b = params['b']
    C = params['C']
    A = params['A']

    z = np.asarray(z, dtype=float)
    g_T = np.zeros_like(z)

    for i, zi in enumerate(z):
        if zi <= -b:
            g_T[i] = 0.0
        elif zi <= -a:
            g_T[i] = -(zi + b)
        elif zi < a:
            g_T[i] = (a - b)
        elif zi < b:
            g_T[i] = (zi - b)
        else:
            g_T[i] = 0.0

    g_T *= C / A

    return np.abs(g_T)


def default_params():
    """Return the default parameter dict for the dislocated reservoir benchmark."""
    return compute_parameters(
        E=14.95e9,
        nu=0.15,
        K_grain=7.119e10,
        dp=-25e6,
        a=75.0,
        b=150.0,
    )


def write_analytical_file(filename, z, values, header="Var1\tVar2"):
    """Write 2-column analytical data file (tab-separated)."""
    with open(filename, 'w') as f:
        f.write(header + '\n')
        for zi, vi in zip(z, values):
            f.write(f'{zi}\t{vi}\n')
    print(f"  Written: {filename} ({len(z)} points)")


def main():
    parser = argparse.ArgumentParser(
        description="Generate analytical solutions for dislocated reservoir benchmark",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    # Material parameters
    parser.add_argument('--E', type=float, default=14.95e9,
                        help="Young's modulus [Pa]")
    parser.add_argument('--nu', type=float, default=0.15,
                        help="Poisson ratio")
    parser.add_argument('--K_grain', type=float, default=7.1190476190476195e10,
                        help="Grain bulk modulus [Pa]")

    # Loading
    parser.add_argument('--dp', type=float, default=-25e6,
                        help="Pressure change in reservoir [Pa] (negative=depletion)")

    # Geometry
    parser.add_argument('--a', type=float, default=75.0,
                        help="Inner reservoir half-extent [m]")
    parser.add_argument('--b', type=float, default=150.0,
                        help="Outer reservoir half-extent [m]")

    # Output grid
    parser.add_argument('--npoints', type=int, default=2000,
                        help="Number of output points")
    parser.add_argument('--zmin', type=float, default=-250.0,
                        help="Minimum z coordinate [m]")
    parser.add_argument('--zmax', type=float, default=250.0,
                        help="Maximum z coordinate [m]")

    # Output files
    parser.add_argument('--slip-file', type=str, default='slip_Analytical.txt',
                        help="Output file for slip profile")
    parser.add_argument('--traction-file', type=str, default='Sxy_Analytical.txt',
                        help="Output file for traction profile")
    parser.add_argument('--outdir', type=str, default='.',
                        help="Output directory")

    args = parser.parse_args()

    # Compute parameters
    params = compute_parameters(args.E, args.nu, args.K_grain, args.dp, args.a, args.b)

    print("=" * 60)
    print("Analytical Solution Generator")
    print("  Dislocated Reservoir with Vertical Fault")
    print("  (Frigo et al., JCP 2026, Eqs. 48-51)")
    print("=" * 60)
    print(f"  E = {params['E']:.4e} Pa")
    print(f"  nu = {params['nu']}")
    print(f"  G = {params['G']:.4e} Pa")
    print(f"  K_drained = {params['K_drained']:.4e} Pa")
    print(f"  K_grain = {params['K_grain']:.4e} Pa")
    print(f"  alpha (Biot) = {params['alpha']:.6f}")
    print(f"  dp = {params['dp']:.4e} Pa")
    print(f"  a = {params['a']} m  (inner half-extent)")
    print(f"  b = {params['b']} m  (outer half-extent)")
    print(f"  C = {params['C']:.6e} Pa")
    print(f"  A = {params['A']:.6e} Pa")
    print(f"  C/A = {params['C']/params['A']:.10e}")
    print(f"  Max slip |C/A*(a-b)| = {abs(params['C']/params['A']*(args.a-args.b)):.15f} m")

    # Generate z-grid
    z = np.linspace(args.zmin, args.zmax, args.npoints)

    # Compute analytical solutions
    slip = analytical_slip(z, params)
    traction = analytical_traction(z, params)

    # Write output files
    os.makedirs(args.outdir, exist_ok=True)
    slip_path = os.path.join(args.outdir, args.slip_file)
    trac_path = os.path.join(args.outdir, args.traction_file)

    write_analytical_file(slip_path, z, slip)
    write_analytical_file(trac_path, z, traction)

    print(f"\n  Slip range: [{slip.min():.6e}, {slip.max():.6e}] m")
    print(f"  Traction range: [{np.nanmin(traction):.6e}, {np.nanmax(traction):.6e}] Pa")
    print("=" * 60)


if __name__ == "__main__":
    main()
