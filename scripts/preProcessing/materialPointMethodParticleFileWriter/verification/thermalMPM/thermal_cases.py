"""Definitions for the GEOS-MPM thermal verification suite (SI units).

No GEOS modules are imported, allowing reference calculations and post-processing
independently of the solver and HPC environment.
"""
from __future__ import annotations

import math

L = 0.010                   # m
Y = 0.00125                 # m
Z = L / 32.0               # m, one cell in plane strain
RHO = 1000.0               # kg/m^3
CP = 10.0                  # J/(kg K)
K = 5.0                    # W/(m K)
ALPHA = K / (RHO * CP)    # m^2/s
T0 = 300.0                 # K
AMP = 100.0               # K
SIGMA = 0.001             # m
RC_GLOBAL = 0.025         # m^2 K/W
RC_MATERIAL = 0.05        # m^2 K/W
RC_GROUP = 0.0125         # m^2 K/W
Q_SPECIFIC = 20000.0      # W/kg
THERMAL_EXPANSION = 1.0e-4 # 1/K
CAPACITY_SLOPE = 0.05      # J/(kg K^2)
CONDUCTIVITY_SLOPE = 0.05  # W/(m K^2)
RC_COMBINED = 0.020        # m^2 K/W

# Base cases and variants. Each ID has its own pfw_input and postProcess scripts.
# Tests use the <Thermal> child schema from geos_mpm_thermal_extension.patch.
CASES = {
    "uniform": dict(kind="uniform", description="Spatially uniform temperature stays constant", stop=0.008, requires=["transport"]),
    "cosine": dict(kind="cosine", description="Exact periodic Fourier-mode decay", stop=0.012, periodic=True, requires=["transport"]),
    "cosine_coarse": dict(kind="cosine", description="Fourier convergence: coarse mesh", stop=0.012, periodic=True, nx=18, requires=["transport"]),
    "cosine_fine": dict(kind="cosine", description="Fourier convergence: fine mesh", stop=0.012, periodic=True, nx=66, requires=["transport"]),
    "gaussian": dict(kind="gaussian", description="Infinite-line Gaussian diffusion at early time", stop=0.0007, requires=["transport"]),
    "step": dict(kind="step", description="Early-time error-function diffusion from hot/cold step", stop=0.0007, requires=["transport"]),
    "source_local": dict(kind="source_local", description="Localized particle selection: adiabatic source inside x-band only", stop=0.012, conductivity=0.0, requires=["transport", "source_local"]),
    "source": dict(kind="source", description="Uniform particle source, exact adiabatic dT/dt", stop=0.012, conductivity=0.0, requires=["transport", "source"]),
    "capacity_linear": dict(kind="capacity_linear", description="Exact nonlinear caloric inversion under adiabatic uniform heating", stop=0.012, conductivity=0., requires=["transport", "source", "linear_capacity"]),
    "conductivity_linear": dict(kind="conductivity_linear", description="Exact Kirchhoff-transformed steady profile for k(T)=k0+k1(T-Tr)", stop=0.008, requires=["transport", "linear_conductivity", "isothermal"]),
    "sink": dict(kind="sink", description="Uniform particle cooling, exact adiabatic dT/dt", stop=0.012, conductivity=0.0, requires=["transport", "source"]),
    "dirichlet_steady": dict(kind="dirichlet_steady", description="Steady linear temperature with isothermal ends", stop=0.010, requires=["transport", "isothermal"]),
    "dirichlet_transient": dict(kind="dirichlet_transient", description="1D transient Dirichlet slab Fourier series", stop=0.010, requires=["transport", "isothermal"]),
    "contact_perfect": dict(kind="contact_perfect", description="Perfect contact limit equals erfc diffusion through a continuous solid", stop=0.0007, requires=["transport", "contact_global"]),
    "contact_gap": dict(kind="contact", description="Positive initial interface gap remains thermally insulating", stop=0.03, requires=["transport", "contact_global"]),
    "contact_global": dict(kind="contact", description="Two-lump contact exchange; global resistance", stop=0.20, requires=["transport", "contact_global"]),
    "contact_material": dict(kind="contact", description="Material-pair override of global resistance", stop=0.20, requires=["transport", "contact_global", "contact_material"]),
    "contact_combined": dict(kind="contact", description="Combined group-and-material pair precedence", stop=0.20, requires=["transport", "contact_global", "contact_material", "contact_group"]),
    "contact_same_material": dict(kind="contact", description="Identical material different contact groups with group-specific resistance", stop=0.20, requires=["transport", "contact_group"]),
    "contact_same_group": dict(kind="contact", description="Different materials in same mechanical contact group have material-pair resistance", stop=0.20, requires=["transport", "contact_material"]),
    "contact_group": dict(kind="contact", description="Contact-group override of material and global resistances", stop=0.20, requires=["transport", "contact_global", "contact_material", "contact_group"]),
    "contact_insulating": dict(kind="contact", description="No thermal transfer across separated contact groups", stop=0.03, requires=["transport", "contact_global"]),
    "contact_dfg": dict(kind="contact", description="Same-material DFG crack thermal resistance", stop=0.20, requires=["transport", "contact_global", "contact_dfg"]),
    "contact_three_field": dict(kind="three_field", description="Three-field junction with only A-C and B-C exchange", stop=0.20, requires=["transport", "contact_global"]),
    "energy_closed": dict(kind="energy_closed", description="Conservation for closed, nonuniform thermally conducting solid", stop=0.010, requires=["transport"]),
    "thermoelastic_free": dict(kind="thermoelastic", description="Free thermal strain alpha DeltaT with near-zero stress", stop=0.015, requires=["transport", "expansion"]),
    "plastic_adiabatic": dict(kind="plastic", description="Adiabatic plastic work-to-heat, diagnostic beta*sigma*epsp", stop=0.02, conductivity=0.0, requires=["transport", "plastic"]),
    "timestep": dict(kind="cosine", description="Thermal time-step stability stress test", stop=0.004, conductivity=25.0, periodic=True, requires=["transport", "timestep"]),
    "mpi_reference": dict(kind="cosine", description="Serial reference for MPI parity", stop=0.008, periodic=True, requires=["transport"]),
    "mpi_2rank": dict(kind="cosine", description="2 MPI ranks, same resolution as serial", stop=0.008, periodic=True, xpar=2, ranks=2, requires=["transport"]),
    "gpu_reference": dict(kind="cosine", description="CPU reference for GPU parity", stop=0.008, periodic=True, requires=["transport"]),
    "gpu_device": dict(kind="cosine", description="HIP execution comparison; requires GPU build", stop=0.008, periodic=True, requires=["transport"]),
}


def diffusivity(case: dict) -> float:
    return float(case.get("conductivity", K)) / (RHO * CP)


def temperature(kind: str, x: float, t: float, *, case: dict | None = None) -> float:
    """Exact or explicitly documented idealized analytical references."""
    case = case or {}
    a = diffusivity(case)
    z = x - L / 2
    if kind == "uniform" or kind == "energy_closed":
        return T0 if kind == "uniform" else (T0 + 0.5 * AMP if x < L / 2 else T0 - 0.5 * AMP)
    if kind == "cosine":
        return T0 + 0.5 * AMP * math.cos(2 * math.pi * x / L) * math.exp(-a * (2 * math.pi / L) ** 2 * t)
    if kind == "gaussian":
        var = SIGMA ** 2 + 2 * a * max(t, 0.0)
        return T0 + AMP * SIGMA / math.sqrt(var) * math.exp(-z * z / (2 * var))
    if kind == "step":
        if t <= 0:
            return T0 + (AMP if x < L / 2 else 0.0)
        return T0 + 0.5 * AMP * math.erfc(z / (2 * math.sqrt(a * t)))
    if kind == "capacity_linear":
        # u=c0*(T-Tr)+c1/2*(T-Tr)^2 = q*t, T(0)=Tr.
        return T0 + 2*Q_SPECIFIC*t/(CP+math.sqrt(CP*CP+2*CAPACITY_SLOPE*Q_SPECIFIC*t))
    if kind == "source":
        return T0 + Q_SPECIFIC / CP * t
    if kind == "sink":
        return T0 - Q_SPECIFIC / CP * t
    if kind == "source_local":
        return T0 + (Q_SPECIFIC / CP * t if 0.003 <= x <= 0.004 else 0.0)
    if kind == "contact_perfect":
        return temperature("step", x, t, case=case)
    if kind in ("dirichlet_steady", "dirichlet_transient", "conductivity_linear"):
        nx=int(case.get('nx',34))
        boundary=L/(nx-2)
        span=L-2*boundary
        xp=min(span,max(0.,x-boundary))
        if kind == "conductivity_linear":
            # Kirchhoff transformation U(T)=int k(T)dT, exact fixed-flux profile.
            total=K*AMP+0.5*CONDUCTIVITY_SLOPE*AMP*AMP
            potential=(1-xp/span)*total
            theta=2*potential/(K+math.sqrt(K*K+2*CONDUCTIVITY_SLOPE*potential))
            return T0+theta
        if kind == "dirichlet_steady":
            return T0 + AMP * (1 - xp / span)
        if t <= 0:
            return T0
        if xp <= 0: return T0+AMP
        if xp >= span: return T0
        # Exact Fourier series for the initial uniform interior with fixed
        # hot/cold reservoirs at the interior strip boundaries.
        s=1-xp/span
        for n in range(1,301):
            s-=2/(n*math.pi)*math.sin(n*math.pi*xp/span)*math.exp(-a*(n*math.pi/span)**2*t)
        return T0+AMP*s
    if kind in ("contact", "three_field", "thermoelastic", "plastic"):
        return T0
    raise KeyError(kind)


def contact_resistance(case_id: str) -> float:
    if case_id == "contact_material":
        return RC_MATERIAL
    if case_id == "contact_combined":
        return RC_COMBINED
    if case_id in ("contact_group", "contact_same_material"):
        return RC_GROUP
    if case_id == "contact_same_group":
        return RC_MATERIAL
    if case_id in ("contact_insulating", "contact_gap"):
        return math.inf
    if case_id == "contact_perfect":
        return 0.0
    return RC_GLOBAL


def contact_lumped_temperatures(t: float, resistance: float) -> tuple[float, float]:
    """Exact ODE for two initially uniform, insulated, equal-capacity bodies.

    This approximates spatially resolved MPM when Bi = (L/2)/(K*R'') << 1.
    """
    if math.isinf(resistance):
        return T0 + AMP, T0
    vol = (L / 2) * Y * Z
    cap = RHO * CP * vol
    g = (Y * Z) / resistance
    decay = math.exp(-2 * g * t / cap)
    return T0 + 0.5 * AMP * (1 + decay), T0 + 0.5 * AMP * (1 - decay)


def three_field_lumped_temperatures(t: float) -> tuple[float, float, float]:
    """Exact star-graph ODE A--C--B (no A--B), equal face conductances."""
    # A and B are each half the height of C; same x thickness.
    ca = RHO * CP * (L / 2) * (Y / 2) * Z
    cc = 2 * ca
    g = (Y / 2 * Z) / RC_GLOBAL
    a0, b0, c0 = T0 + AMP, T0, T0 + 0.25 * AMP
    s0, d0 = 0.5 * (a0 + b0), 0.5 * (a0 - b0)
    mean = (ca * a0 + ca * b0 + cc * c0) / (2 * ca + cc)
    u = (s0 - c0) * math.exp(-(g / ca + 2 * g / cc) * t)
    d = d0 * math.exp(-(g / ca) * t)
    s = mean + cc / (2 * ca + cc) * u
    c = mean - 2 * ca / (2 * ca + cc) * u
    return s + d, s - d, c
