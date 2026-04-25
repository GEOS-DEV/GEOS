"""
wellboreTemperatureSolutions.py

Shared analytic and numerical reference solutions for wellbore thermal
diffusion validation tests.

Covers:
  - Linear analytic transient (Wang & Papamichos 1994)
  - Log steady-state
  - Non-linear FDM reference
"""

import numpy as np
from scipy import special
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve


# ---------------------------------------------------------------------------
# Linear analytic solutions
# ---------------------------------------------------------------------------

def steadyState(Tin, Tout, Rin, Rout, r):
    """Log-linear steady-state radial temperature profile.

    Parameters
    ----------
    Tin, Tout : float
        Temperature imposed at inner (r=Rin) and outer (r=Rout) boundaries (°C or K).
    Rin, Rout : float
        Inner and outer radii (m).
    r : array_like
        Radial coordinates at which to evaluate the solution (m).

    Returns
    -------
    numpy.ndarray
        Temperature at each r.
    """
    r = np.asarray(r)
    return Tin + (Tout - Tin) * (np.log(r) - np.log(Rin)) / (np.log(Rout) - np.log(Rin))


def transientTemperature(Tin, Tout, Rin, r, thermalDiffusivity, t):
    """Approximate transient radial temperature (Wang & Papamichos 1994).

    Valid for an infinite domain (Rout → ∞) in the early-time regime.

    Reference
    ---------
    Wang, Y. and Papamichos, E. (1994). Conductive heat flow and thermally
    induced fluid flow around a wellbore in a poroelastic medium.
    Water Resources Research, 30(12), 3375–3384.

    Parameters
    ----------
    Tin : float
        Temperature at the wellbore wall (r = Rin).
    Tout : float
        Initial / far-field temperature.
    Rin : float
        Wellbore radius (m).
    r : array_like
        Radial coordinates (m).
    thermalDiffusivity : float
        κ = λ / (ρ c_eff)  (m²/s).
    t : float
        Elapsed time (s).

    Returns
    -------
    numpy.ndarray
        Temperature at each r.
    """
    r  = np.asarray(r)
    xi = (r - Rin) / (2.0 * np.sqrt(thermalDiffusivity * t))
    return Tout + (Tin - Tout) * np.sqrt(Rin / r) * special.erfc(xi)


# ---------------------------------------------------------------------------
# Temperature-dependent material properties
# ---------------------------------------------------------------------------

def thermalConductivity(lambda0, lambda_gradient, T, Tref):
    """Linear T-dependent thermal conductivity: λ(T) = λ0 + dλ/dT * (T - Tref).

    Parameters
    ----------
    lambda0 : float or array
        Reference thermal conductivity at Tref (W/m/K).
    lambda_gradient : float or array
        dλ/dT (W/m/K²).
    T : array_like
        Current temperature field (°C or K).
    Tref : float
        Reference temperature (°C or K).

    Returns
    -------
    numpy.ndarray
        λ(T) at each point.
    """
    return np.asarray(lambda0) + np.asarray(lambda_gradient) * (np.asarray(T) - Tref)


def effectiveVolumetricHeatCapacity(c0, c_gradient, T, Tref, porosity):
    """Porosity-weighted solid volumetric heat capacity: c_eff(T) = (1-φ)(c0 + dc/dT*(T-Tref)).

    Consistent with the GEOS energy-balance assembly, which stores the
    solid-only heat capacity and multiplies by the solid volume fraction
    (1 - φ) when computing the effective capacity.

    Parameters
    ----------
    c0 : float
        Reference solid volumetric heat capacity at Tref (J/m³/K).
    c_gradient : float
        dc/dT (J/m³/K²).
    T : array_like
        Current temperature field (°C or K).
    Tref : float
        Reference temperature (°C or K).
    porosity : float
        Rock porosity φ ∈ [0, 1].

    Returns
    -------
    numpy.ndarray
        Effective volumetric heat capacity at each point.
    """
    return (1.0 - porosity) * (c0 + c_gradient * (np.asarray(T) - Tref))


# ---------------------------------------------------------------------------
# Non-linear FD solver  (non-uniform radial stencil, implicit, lagged coefficients)
# ---------------------------------------------------------------------------

def _buildMatrix(lam, cap, r, dt, N):
    """Assemble the implicit backward-Euler system matrix for one time step.

    Grid layout (N+1 nodes):
        r[0]    = Rin  — Dirichlet boundary face (inner)
        r[1:-1] = GEOS cell centres (N-1 interior nodes)
        r[N]    = Rout — Dirichlet boundary face (outer)

    Interior node i (1 ≤ i ≤ N-1):
        dr_L = r[i] - r[i-1],  dr_R = r[i+1] - r[i],  S = dr_L + dr_R

        Discretised PDE (cylindrical):
            c_eff * dT/dt = d/dr[λ dT/dr] + λ/r * dT/dr

        d²T/dr²  ≈ (2/S) * [(T[i+1]-T[i])/dr_R − (T[i]-T[i-1])/dr_L]
        dT/dr    ≈ (T[i+1] - T[i-1]) / S
        dλ/dr    ≈ (λ[i+1] - λ[i-1]) / S   (lagged)

    Material coefficients are lagged at T^n (explicit treatment of non-linearity).
    """
    A = np.zeros((N + 1, N + 1))

    for i in range(1, N):
        dr_L = r[i]     - r[i - 1]
        dr_R = r[i + 1] - r[i]
        S    = dr_L + dr_R
        r_i  = r[i]

        alpha_i   = lam[i] / cap[i]                          # thermal diffusivity
        dlambda_i = (lam[i + 1] - lam[i - 1]) / S            # dλ/dr (lagged)
        beta_i    = dlambda_i / cap[i]

        fd1 = (alpha_i / r_i + beta_i) / S                   # first-derivative coefficient

        A[i, i - 1] = -dt * (2.0 * alpha_i / (S * dr_L) - fd1)
        A[i, i    ] =  1.0 + dt * 2.0 * alpha_i / (dr_L * dr_R)
        A[i, i + 1] = -dt * (2.0 * alpha_i / (S * dr_R) + fd1)

    # Dirichlet BCs at boundary faces
    A[0, 0] = 1.0
    A[N, N] = 1.0
    return A


def solveRadialDiffusion(r_cells, Rin, Rout, tmax, dt,
                         Tin, Tout,
                         lambda0, lambda_gradient,
                         c0, c_gradient, Tref, porosity):
    """Solve the non-linear radial heat equation with an implicit FD scheme.

    The computational grid augments the GEOS cell centres with the actual
    boundary faces so that Dirichlet conditions are imposed at the correct
    radial positions (not at the first/last cell centres).

    Parameters
    ----------
    r_cells : array_like
        GEOS cell-centre radial coordinates (m), sorted ascending.
    Rin, Rout : float
        Inner and outer boundary radii (m).
    tmax : float
        Total simulation time (s).
    dt : float
        Time-step size (s). A smaller dt improves accuracy.
    Tin, Tout : float
        Dirichlet temperatures at Rin and Rout (°C or K).
    lambda0, lambda_gradient : float
        Thermal conductivity reference value and gradient (W/m/K and W/m/K²).
    c0, c_gradient : float
        Solid heat capacity reference value and gradient (J/m³/K and J/m³/K²).
    Tref : float
        Reference temperature for material property linearisation (°C or K).
    porosity : float
        Rock porosity φ.

    Returns
    -------
    numpy.ndarray
        Temperature at each GEOS cell centre after tmax seconds.
    """
    r_cells = np.asarray(r_cells, dtype=float)
    r       = np.concatenate([[Rin], r_cells, [Rout]])
    N       = len(r) - 1
    n_steps = max(1, int(round(tmax / dt)))

    T      = np.full(N + 1, Tout, dtype=float)
    T[0]   = Tin

    for _ in range(n_steps):
        lam = thermalConductivity(lambda0, lambda_gradient, T, Tref)
        cap = effectiveVolumetricHeatCapacity(c0, c_gradient, T, Tref, porosity)
        A   = csr_matrix(_buildMatrix(lam, cap, r, dt, N))
        T   = spsolve(A, T)

    return T[1:-1]   # strip the two boundary-face nodes
