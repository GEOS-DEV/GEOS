"""Spherical inclusion in the positive octant, with three symmetry planes.

Default: spherical specimen under the same virtual-fluid pressure event used by
weakTraceCircularInclusionRadial. The static reference is genuinely 3D; do NOT
reuse the plane-strain circular-inclusion solution.

WEAK_TRACE_LOADING=ftableCube is a separate kinematic smoke test: a sphere cuts
three Moving/FTable faces of a cube. It is NOT the radial-pressure problem.
Geometry flags and particle positions are identical between singleField and
traceContactGroups; only the contact-group labels and solver options change.
"""
from __future__ import annotations

import math
import os
import numpy as np
import pfw_geometryObjects as geom
import pfw_tracerPoints as tracers


def positive_int(name, default):
    value = int(os.environ.get(name, str(default)))
    if value < 1:
        raise ValueError(name + " must be positive")
    return value


def positive_float(name, default):
    value = float(os.environ.get(name, str(default)))
    if not math.isfinite(value) or value <= 0:
        raise ValueError(name + " must be finite and positive")
    return value


def bool_env(name, default):
    raw = os.environ.get(name)
    if raw is None:
        return bool(default)
    value = raw.strip().lower()
    if value in ("1", "true", "yes", "on"):
        return True
    if value in ("0", "false", "no", "off"):
        return False
    raise ValueError(name + " must be a boolean")


def marker_search_mode_env():
    """Select legacy (0), indexed (1), or checked indexed (2) marker search.

    Preserve the existing true/false spellings for modes 1/0. Mode 2 runs both
    searches and verifies their neighbor sequences; use it for short correctness
    runs, not for performance measurements. This changes solver XML only.
    """
    name = "WEAK_TRACE_FAST_MARKER_SEARCH"
    if os.environ.get(name, "").strip() == "2":
        return 2
    return int(bool_env(name, True))


loading = os.environ.get("WEAK_TRACE_LOADING", "pressure")
variant = os.environ.get("WEAK_TRACE_VARIANT", "traceContactGroups")
if loading not in ("pressure", "ftableCube"):
    raise ValueError("WEAK_TRACE_LOADING must be pressure or ftableCube")
if variant not in ("singleField", "traceContactGroups"):
    raise ValueError("WEAK_TRACE_VARIANT must be singleField or traceContactGroups")

cpp = positive_int("WEAK_TRACE_CPP", 24)  # total grid cells, including one buffer per end
ppc = positive_int("WEAK_TRACE_PPC", 2)  # uncut material keeps the standard two particles per cell
if cpp < 8:
    raise ValueError("Use at least eight total cells per direction for this test")
partitions = [positive_int("WEAK_TRACE_" + d + "PAR", 1) for d in "XYZ"]
if any(cpp % n for n in partitions):
    raise ValueError("WEAK_TRACE_CPP must be divisible by each solver partition count")
outer_radius = positive_float("WEAK_TRACE_OUTER_RADIUS", 0.72)
gap_cells = positive_int("WEAK_TRACE_GAP_CELLS", 3)
if gap_cells >= cpp - 2:
    raise ValueError("Positive-side empty gap must be smaller than the physical grid")
# L - b = gap_cells*h, h=L/(cpp-2). Negative faces remain exactly at zero.
length = outer_radius / (1.0 - gap_cells / float(cpp - 2)) if loading == "pressure" else 0.72
radius = positive_float("WEAK_TRACE_INCLUSION_RADIUS", 0.22 if loading == "pressure" else 1.2 * length)
if loading == "pressure" and radius >= outer_radius:
    raise ValueError("Inclusion radius must be below the specimen radius")
if loading == "ftableCube" and not length < radius < math.sqrt(3.0) * length:
    raise ValueError("FTable smoke requires L < interface radius < sqrt(3)*L so it cuts the moving faces")
h = length / (cpp - 2)
surface_skin_cells = positive_float("WEAK_TRACE_SURFACE_SKIN_CELLS", 0.5)
ramp_time = positive_float("WEAK_TRACE_RAMP_TIME", 8.0)
stop_time = positive_float("WEAK_TRACE_STOP_TIME", 12.0)
if stop_time < ramp_time:
    raise ValueError("Stop time must include the complete loading ramp")
pressure = positive_float("WEAK_TRACE_EXTERNAL_PRESSURE", 0.004)
young_matrix = positive_float("WEAK_TRACE_MATRIX_E", 1.0)
young_inclusion = positive_float("WEAK_TRACE_INCLUSION_E", 4.0)
poisson = float(os.environ.get("WEAK_TRACE_POISSON", "0.25"))
if not -1.0 < poisson < 0.5:
    raise ValueError("Poisson ratio must be between -1 and 0.5")
density = positive_float("WEAK_TRACE_DENSITY", 1.0)
case = os.environ.get("VV_CASE_NAME", "weakTraceSphericalInclusionOctant_" + loading + "_" + variant)
run_debug = bool_env("WEAK_TRACE_RUN_DEBUG", True)

pfw = dict(
    caseName=case, runDebug=run_debug, mBatch=True, mSubmitJobs=False,
    mWallTime=os.environ.get("WEAK_TRACE_WALLTIME", "00:45:00"),
    mPartition=os.environ.get("WEAK_TRACE_PARTITION", "pdebug" if run_debug else "pbatch"),
    logLevel=int(os.environ.get("WEAK_TRACE_LOG_LEVEL", "2" if variant == "traceContactGroups" else "0")),
    planeStrain=0, periodic=[False, False, False],
    xpar=partitions[0], ypar=partitions[1], zpar=partitions[2],
    nI=cpp, nJ=cpp, nK=cpp, ppc=ppc,
    xmin=0.0, ymin=0.0, zmin=0.0, xmax=length, ymax=length, zmax=length,
    surfaceRefinementLevel=positive_int("WEAK_TRACE_SURFACE_REFINEMENT_LEVEL", 2),
    surfaceRefinementMethod=os.environ.get("WEAK_TRACE_SURFACE_REFINEMENT_METHOD", "simplex"),
    surfaceRefinementFallback="cartesian",
    particleReferenceGeometry="fromRefinement", particleDomainOutput="native",
    endTime=stop_time, plotInterval=positive_float("WEAK_TRACE_PLOT_INTERVAL", stop_time / 40.0),
    restartInterval=10.0 * stop_time, timeIntegrationOption="ExplicitDynamic",
    cflFactor=positive_float("WEAK_TRACE_CFL", 0.15), initialDt=1.0e-16,
    writeStatistics="all", reactionHistory=1,
    reactionWriteInterval=positive_float("WEAK_TRACE_HISTORY_INTERVAL", stop_time / 200),
    boxAverageHistory=1,
    boxAverageWriteInterval=positive_float("WEAK_TRACE_HISTORY_INTERVAL", stop_time / 200),
    updateMethod=os.environ.get("WEAK_TRACE_UPDATE_METHOD", "PIC"),
    cpdiDomainScaling=0, outputType="silo", plotGridFields=1,
    frictionCoefficient=0.0, contactNormalType="Aligned", useSurfacePositionForContact=1,
    damageFieldPartitioning=0, weakDiscontinuityInsertion="off",
    weakInterfaceTraceGeometry="surface", enableWeakInterfacePrimalEnrichment=0,
    prescribedFTable=0, prescribedBoundaryFTable=0,
    # Ordering: x-, x+, y-, y+, z-, z+. Zero is outflow, NOT a wall.
    boundaryConditionTypes=[1, 0, 1, 0, 1, 0],
    particleFileFields=["Velocity", "MaterialType", "ContactGroup", "SurfaceFlag",
                        "RVector", "SurfaceNormal", "SurfacePosition"],
    plottableFields=["particleStress", "particleVelocity", "particleDisplacement"],
)
if pfw["updateMethod"].upper() == "XPIC":
    pfw["updateOrder"] = positive_int("WEAK_TRACE_UPDATE_ORDER", 2)

pfw["gridFieldNames"] = [
    "gridMass", "gridActive", "gridVelocity", "gridInternalForce", "gridExternalForce",
    "gridContactForce", "gridWeakInterfaceTraceActive", "gridWeakInterfaceTraceSkipReason",
    "gridWeakInterfaceTraceSupportWeight", "gridWeakInterfaceTraceVelocityJumpPost",
    "gridWeakInterfaceTraceForce", "gridWeakInterfaceTraceBoundaryForce",
]
if variant == "traceContactGroups":
    qr_backend = os.environ.get("WEAK_TRACE_QR_BACKEND", "auto")
    if qr_backend not in ("auto", "reference", "native", "lapack"):
        raise ValueError("WEAK_TRACE_QR_BACKEND must be auto, reference, native, or lapack")
    pfw.update(
        enableContact=1, enableWeakInterfaceTraceProjection=1,
        weakInterfaceTraceBoundaryCoupling=1, weakInterfaceTraceProjectionScale=1.0,
        weakInterfaceTraceQRBackend=qr_backend,
        weakInterfaceTraceComponentRanks=positive_int("WEAK_TRACE_COMPONENT_RANKS", 3),
        weakInterfaceTraceFastMarkerSearch=marker_search_mode_env(),
        weakInterfaceTracePerformanceLog=1 if bool_env("WEAK_TRACE_PERFORMANCE_LOG", False) else 0,
        weakInterfaceTraceSuppressNodalContact=1, weakInterfaceTracePairs=[[0, 1]],
        weakInterfaceTraceSurfaceFitNeighbors=12, weakInterfaceTraceSurfaceFitRadius=2.5,
        weakInterfaceTraceMinNormalCoherence=0.0, weakInterfaceTraceMaxFitResidual=1.0e6,
        weakInterfaceTraceMortarQuadratureOrder=2, weakInterfaceTraceMassCutoffRatio=1.0e-12,
        weakInterfaceTraceMinWeight=1.0e-12, weakInterfaceTraceCutTolerance=1.0e-8,
        weakInterfaceTraceMinSegmentFraction=1.0e-4,
        weakInterfaceTraceMinFacetAreaFraction=0.0, weakInterfaceTraceMinTriangleQuality=0.0,
        weakInterfaceTraceFallbackToSingleField=1, weakInterfaceTraceGapStabilization=0.0,
        weakInterfaceTraceMaxGapCorrectionVelocity=0.0,
        weakInterfaceTraceLinearSolverTolerance=1.0e-10,
        weakInterfaceTraceAbsoluteTolerance=positive_float("WEAK_TRACE_ABSOLUTE_TOLERANCE", 1.0e-12),
        weakInterfaceTraceQRRankTolerance=positive_float("WEAK_TRACE_QR_RANK_TOLERANCE", 1.0e-12),
        weakInterfaceTraceQRMaxEntries=positive_int("WEAK_TRACE_QR_MAX_ENTRIES", 16000000),
    )
else:
    pfw.update(enableContact=0, enableWeakInterfaceTraceProjection=0)


def elastic_xml(name, young):
    return ('<ElasticIsotropic name="%s" defaultDensity="%.17g" '
            'defaultYoungModulus="%.17g" defaultPoissonRatio="%.17g"/>'
            % (name, density, young, poisson))


pfw["materials"] = ["octantMatrix", "octantInclusion"]
pfw["materialPropertyString"] = elastic_xml("octantMatrix", young_matrix) + "\n" + elastic_xml("octantInclusion", young_inclusion)
if loading == "pressure":
    pressure_radius = positive_float("WEAK_TRACE_CONFINING_SPHERE_RADIUS",
                                     0.5 * (radius + outer_radius))
    if not radius < pressure_radius < outer_radius:
        raise ValueError("WEAK_TRACE_CONFINING_SPHERE_RADIUS must lie between the inclusion and specimen radii")
    # Apply the virtual-fluid background pressure outside an origin-centered
    # spherical surface, not outside an axis-aligned cube. Keeping this transition
    # between the inclusion and outer specimen surface avoids putting pressure
    # background nodes inside the support of the weak material interface.
    pfw["useEvents"] = 1
    pfw["mpmEventsString"] = f'''
<ConfiningPressure name="octantSphericalOuterPressure"
  startTime="0.0" endTime="{ramp_time}"
  confiningPressureSphereRadius="{pressure_radius}"
  startPressure="0.0" endPressure="{pressure}" interpType="1"/>
'''
else:
    # Different geometry/loading: cube cut by a spherical material interface.
    # This tests arbitrary intersections with Moving faces as well as symmetry.
    # It has no concentric-sphere analytic stress reference.
    pfw["prescribedBoundaryFTable"] = 1
    pfw["boundaryConditionTypes"] = [1, 2, 1, 2, 1, 2]
    pfw["fTableInterpType"] = "Cosine"
    pfw["fTable"] = [[0, 1, 1, 1], [0.5*ramp_time, .998, .998, .998],
                     [ramp_time, 1, 1, 1], [stop_time, 1, 1, 1]]
    if stop_time == ramp_time:
        pfw["fTable"].pop()  # Avoid duplicate table times.


class RadialRegion(geom.Geometry):
    """Analytic region; reuses the existing membership/normal/offset interface.

    The PFW generation box clips to the octant. The three symmetry cuts are NOT
    flagged as material free surfaces. This also makes halo queries symmetric.
    """
    def __init__(self, inclusion):
        self.inclusion = inclusion
        material = 1 if inclusion else 0
        super().__init__("octant_inclusion" if inclusion else "octant_matrix",
                         vel=[0., 0., 0.], mat=material,
                         group=material if variant == "traceContactGroups" else 0,
                         particleType=1)

    def _surface(self, pt):
        r = float(np.linalg.norm(pt))
        outer = (not self.inclusion and loading == "pressure" and
                 abs(r - outer_radius) < abs(r - radius))
        return r, outer_radius if outer else radius, outer

    def isInterior(self, pt, skinDepth):
        r, nearest, outer = self._surface(pt)
        if self.inclusion:
            if r >= radius:
                return -1
        elif r < radius or (loading == "pressure" and r >= outer_radius):
            return -1
        flag_depth = min(float(skinDepth), surface_skin_cells * h)
        if abs(r - nearest) <= flag_depth:
            return int(geom.SurfaceFlag.Surface if outer else geom.SurfaceFlag.WeakDiscontinuity)
        return 0

    def getSurfaceNormal(self, pt):
        r, _, outer = self._surface(pt)
        direction = np.asarray(pt, dtype=float) / r if r > 0 else np.array([1., 0., 0.])
        return direction if self.inclusion or outer else -direction

    def getSurfacePosition(self, pt):
        r, surface_radius, _ = self._surface(pt)
        direction = np.asarray(pt, dtype=float) / r if r > 0 else np.array([1., 0., 0.])
        return (surface_radius - r) * direction


pfw["objects"] = [RadialRegion(False), RadialRegion(True)]
# The stock PFW can run rank-sliced without a factory; it owns lattice rows and
# queries read-only halos against these complete analytic objects.
tracer_points = []
if loading == "pressure":
    directions = np.array([[1., .05, .05], [.05, 1., .05], [.05, .05, 1.],
                           [1., 1., .05], [1., .05, 1.], [.05, 1., 1.], [1., 1., 1.]])
    directions /= np.linalg.norm(directions, axis=1)[:, None]
    for direction in directions:
        for r in (0.5*radius, .90*radius, 1.10*radius, .65*outer_radius, .90*outer_radius):
            tracer_points.append(tuple(r * direction))
    tracers.set_tracers(pfw, tracer_points,
        ["velocityX", "velocityY", "velocityZ", "stressXX", "stressYY", "stressZZ",
         "stressXY", "stressXZ", "stressYZ", "materialType"],
        write_interval=stop_time/200, output_prefix="sphericalOctantTracer")

pfw_expected = dict(case_id="weakTraceSphericalInclusionOctant", loading=loading,
    variant=variant, inclusion_radius=radius, outer_radius=outer_radius,
    domain_length=length, grid_spacing=h, positive_gap_cells=gap_cells if loading == "pressure" else 0,
    cells_per_direction=cpp, particles_per_cell_direction=ppc,
    surface_refinement_level=pfw["surfaceRefinementLevel"],
    surface_refinement_method=pfw["surfaceRefinementMethod"],
    surface_skin_cells=surface_skin_cells,
    matrix_E=young_matrix, inclusion_E=young_inclusion, poisson=poisson, density=density,
    external_pressure=pressure if loading == "pressure" else None,
    confining_pressure_sphere_radius=pressure_radius if loading == "pressure" else None,
    pressure_ramp_time=ramp_time, stop_time=stop_time,
    reference="3D small-strain static pressure solution after equilibration" if loading == "pressure" else "none: FTable kinematic smoke",
    exact_octant_volume=math.pi*outer_radius**3/6 if loading == "pressure" else length**3)
