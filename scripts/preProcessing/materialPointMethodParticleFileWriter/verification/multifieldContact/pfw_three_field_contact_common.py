"""Shared builders for the three-field projected-contact verification decks."""

from __future__ import annotations

from typing import Iterable, Sequence

import pfw_geometryObjects as geom


CELL_SIZE = 0.125
XMIN = -1.0
XMAX = 1.0
YMIN = -0.5
YMAX = 0.5
SOLID_YMIN = -0.375
SOLID_YMAX = 0.375


def _material_block(densities: Sequence[float]) -> tuple[list[str], str]:
    """Return density-scaled elastic materials with a common wave speed."""
    names = []
    blocks = []
    for index, density in enumerate(densities):
        name = f"threeFieldElastic{index}"
        names.append(name)
        # Scaling E with density keeps the elastic wave speed approximately
        # unchanged, so the mass-ratio case isolates contact conditioning.
        young_modulus = 100.0 * float(density)
        blocks.append(
            f"""<ElasticIsotropic
    name="{name}"
    defaultDensity="{float(density):.16g}"
    defaultYoungModulus="{young_modulus:.16g}"
    defaultPoissonRatio="0.25"/>"""
        )
    return names, "\n".join(blocks)


def _base_pfw(
    case_name: str,
    *,
    densities: Sequence[float] = (1.0, 1.0, 1.0),
    friction_coefficient: float = 0.0,
    gap_correction: str = "Implicit",
    use_surface_positions: bool = True,
    explicit_normals: bool = True,
    end_time: float = 0.02,
) -> dict:
    if len(densities) != 3:
        raise ValueError("Three-field tests require exactly three material densities")

    material_names, material_xml = _material_block(densities)
    pfw = {
        "runDebug": True,
        "caseName": case_name,
        "xmin": XMIN,
        "xmax": XMAX,
        "ymin": YMIN,
        "ymax": YMAX,
        "zmin": -0.5 * CELL_SIZE,
        "zmax": 0.5 * CELL_SIZE,
        "planeStrain": 1,
        "periodic": [False, False, False],
        # (nI-2) and (nJ-2) give 16 x 8 physical cells, both with h=0.125.
        "xpar": 1,
        "ypar": 1,
        "zpar": 1,
        "nI": 18,
        "nJ": 10,
        "nK": 3,
        "ppc": 2,
        "mBatch": True,
        "mCores": 1,
        "mWallTime": "00:05:00",
        "mSubmitJobs": False,
        "autoRestart": False,
        "endTime": end_time,
        "plotInterval": end_time / 10.0,
        "restartInterval": 2.0 * end_time,
        "outputType": "silo",
        "timeIntegrationOption": "ExplicitDynamic",
        "updateMethod": "PIC",
        "cflFactor": 0.20,
        "initialDt": 1.0e-8,
        "cpdiDomainScaling": 0,
        "subdivideParticles": 0,
        "damageFieldPartitioning": 0,
        "needsNeighborList": 0,
        "enableContact": 1,
        "boundaryConditionTypes": [0, 0, 0, 0, 1, 1],
        "frictionCoefficient": float(friction_coefficient),
        "contactGapActivationRelativeTolerance": 1.0e-10,
        "contactGapCorrection": gap_correction,
        "contactNormalType": "Difference",
        "contactNormalExponent": 1.0,
        "useSurfacePositionForContact": 1 if use_surface_positions else 0,
        "explicitSurfaceNormalInfluence": 1000.0 if explicit_normals else 0.0,
        "disableSurfaceNormalsAndPositionsOnCPDIScaling": 0,

        "contactSolver": "ProjectedGaussSeidel",
        "contactPGSMaximumIterations": 200,
        "contactPGSVelocityTolerance": 1.0e-10,
        "contactPGSRelaxation": 1.0,


        # # Newton-Raphson
        # "contactSolver": "NewtonRaphson",
        # "contactNRMaximumIterations": 50,
        # "contactNRVelocityTolerance": 1.0e-10,
        # "contactNRFiniteDifferenceRelativeStep": 1.0e-6,
        # "contactNRLineSearchMinimumScale": 1.0e-4,
        # "contactNRRegularization": 1.0e-12,

        "contactPGSRequireConvergence": 1,
        "contactPGSUseLogisticRegressionForMultifield": 0,

        "contactSolverDiagnostics": 1,
        "contactSolverDiagnosticMaxNodes": 20,
        "contactSolverFailureDiagnosticMaxNodes": 20,
        "contactSolverFailureDiagnostics": 1,

        "maxParticleVelocity": 2.0,
        "minParticleJacobian": 0.01,
        "maxParticleJacobian": 10.0,
        "solverProfiling": 0,
        "reactionHistory": 0,
        "boxAverageHistory": 0,
        "writeStatistics": "all",
        "plotGridFields": 1,
        "gridFieldNames": [
            "gridActive",
            "gridMass",
            "gridUncontactedVelocity",
            "gridVelocity",
            "gridContactForce",
            "gridCenterOfMass",
            "gridSurfaceNormal",
            "gridSurfacePosition",
        ],
        "plottableFields": [
            "particleID",
            "particleMass",
            "particleMaterialType",
            "particleGroup",
            "particleSurfaceFlag",
            "particleCenter",
            "particleVelocity",
            "particleSurfaceNormal",
            "particleSurfacePosition",
            "gridActive",
            "gridMass",
            "gridUncontactedVelocity",
            "gridVelocity",
            "gridContactForce",
            "gridCenterOfMass",
            "gridSurfaceNormal",
            "gridSurfacePosition",
        ],
        "materials": material_names,
        "materialPropertyString": material_xml,
    }

    particle_fields = [
        "Velocity",
        "MaterialType",
        "ContactGroup",
        "SurfaceFlag",
        "RVector",
    ]
    if explicit_normals or use_surface_positions:
        particle_fields.extend(["SurfaceNormal", "SurfacePosition"])
    pfw["particleFileFields"] = particle_fields
    return pfw


def _fixed_normal_box(
    name: str,
    x0: Sequence[float],
    x1: Sequence[float],
    *,
    velocity: Sequence[float],
    material: int,
    group: int,
    flagged_surfaces: Sequence[bool],
    normal: Sequence[float] | None = None,
):
    obj = geom.box(
        name,
        list(x0),
        list(x1),
        vel=list(velocity),
        mat=int(material),
        group=int(group),
        particleType=2,
        dim=2,
        flaggedSurfaces=list(flagged_surfaces),
    )
    if normal is not None:
        obj = geom.surfaceNormalWrapper(f"{name}FixedNormal", obj, list(normal))
    return obj


def make_shared_interface_case(
    case_name: str,
    *,
    gap: float,
    velocities: Sequence[Sequence[float]],
    groups: Sequence[int] = (0, 1, 2),
    densities: Sequence[float] = (1.0, 1.0, 1.0),
    friction_coefficient: float = 0.0,
    gap_correction: str = "Implicit",
    end_time: float = 0.02,
) -> dict:
    """Two left fields contact one right field on a common planar interface.

    ``gap > 0`` separates the nominal surfaces, ``gap == 0`` makes them
    coincident, and ``gap < 0`` gives a small prescribed geometric overlap.
    The lower-left and upper-left bodies do not contact each other; they create
    two coupled constraints through the right-hand field at the junction node.
    """
    if len(velocities) != 3 or len(groups) != 3:
        raise ValueError("velocities and groups must each contain three entries")
    if sorted(int(group) for group in groups) != [0, 1, 2]:
        raise ValueError("groups must be a permutation of (0, 1, 2)")

    pfw = _base_pfw(
        case_name,
        densities=densities,
        friction_coefficient=friction_coefficient,
        gap_correction=gap_correction,
        use_surface_positions=True,
        explicit_normals=True,
        end_time=end_time,
    )

    left_surface = -0.5 * float(gap)
    right_surface = 0.5 * float(gap)
    z0 = pfw["zmin"]
    z1 = pfw["zmax"]

    lower_left = _fixed_normal_box(
        "lowerLeftField",
        (-0.80, SOLID_YMIN, z0),
        (left_surface, 0.0, z1),
        velocity=velocities[0],
        material=0,
        group=groups[0],
        flagged_surfaces=(False, False, True, False),
        normal=(1.0, 0.0, 0.0),
    )
    upper_left = _fixed_normal_box(
        "upperLeftField",
        (-0.80, 0.0, z0),
        (left_surface, SOLID_YMAX, z1),
        velocity=velocities[1],
        material=1,
        group=groups[1],
        flagged_surfaces=(False, False, True, False),
        normal=(1.0, 0.0, 0.0),
    )
    right = _fixed_normal_box(
        "rightField",
        (right_surface, SOLID_YMIN, z0),
        (0.80, SOLID_YMAX, z1),
        velocity=velocities[2],
        material=2,
        group=groups[2],
        flagged_surfaces=(True, False, False, False),
        normal=(-1.0, 0.0, 0.0),
    )

    # For gap < 0 the nominal geometries overlap.  Priority is intentional:
    # the left fields own the overlap particles while the surface positions
    # still describe interpenetrating nominal surfaces.
    pfw["objects"] = [lower_left, upper_left, right]
    return pfw


def make_collinear_chain_case(
    case_name: str,
    *,
    gap_correction: str = "Implicit",
    use_surface_positions: bool = False,
    groups: Sequence[int] = (0, 1, 2),
    end_time: float = 0.02,
) -> dict:
    """Create the A|B|C one-node chain used by the PGS three-field unit test.

    The middle strip is one cell wide and centered on x=0, so the common node
    receives nonzero CPDI support from all fields.  Its fixed +x explicit normal
    prevents the two opposing B faces from cancelling.  With
    ``use_surface_positions=False``, A-B and B-C are active from mapped centers
    while the farther A-C pair should remain inactive under Implicit correction.
    """
    if sorted(int(group) for group in groups) != [0, 1, 2]:
        raise ValueError("groups must be a permutation of (0, 1, 2)")

    pfw = _base_pfw(
        case_name,
        densities=(1.0, 1.0, 1.0),
        friction_coefficient=0.0,
        gap_correction=gap_correction,
        use_surface_positions=use_surface_positions,
        explicit_normals=True,
        end_time=end_time,
    )

    left_interface = -0.50 * CELL_SIZE
    right_interface = 0.50 * CELL_SIZE
    z0 = pfw["zmin"]
    z1 = pfw["zmax"]
    velocities = ((0.10, 0.0, 0.0), (0.0, 0.0, 0.0), (-0.10, 0.0, 0.0))

    left = _fixed_normal_box(
        "leftChainField",
        (-0.80, -0.25, z0),
        (left_interface, 0.25, z1),
        velocity=velocities[0],
        material=0,
        group=groups[0],
        flagged_surfaces=(False, False, True, False),
        normal=(1.0, 0.0, 0.0),
    )
    middle = _fixed_normal_box(
        "middleChainField",
        (left_interface, -0.25, z0),
        (right_interface, 0.25, z1),
        velocity=velocities[1],
        material=1,
        group=groups[1],
        flagged_surfaces=(True, False, True, False),
        # One velocity field stores one mapped normal at a node.  This fixed
        # orientation keeps the chain active and makes the integration test
        # match the abstract unit-test constraints n_AB=n_BC=+x.
        normal=(1.0, 0.0, 0.0),
    )
    right = _fixed_normal_box(
        "rightChainField",
        (right_interface, -0.25, z0),
        (0.80, 0.25, z1),
        velocity=velocities[2],
        material=2,
        group=groups[2],
        flagged_surfaces=(True, False, False, False),
        normal=(-1.0, 0.0, 0.0),
    )
    pfw["objects"] = [left, middle, right]
    return pfw


def make_corner_case(case_name: str, *, penetration: float = 0.05 * CELL_SIZE) -> dict:
    """Create three fields meeting at a corner with nonparallel pair normals."""
    if penetration <= 0.0:
        raise ValueError("penetration must be positive")

    pfw = _base_pfw(
        case_name,
        densities=(1.0, 1.0, 1.0),
        friction_coefficient=0.0,
        gap_correction="Implicit",
        use_surface_positions=True,
        explicit_normals=True,
        end_time=0.02,
    )
    z0 = pfw["zmin"]
    z1 = pfw["zmax"]
    half_overlap = 0.5 * penetration

    lower_left = geom.box(
        "cornerDrivingField",
        [-0.75, -0.375, z0],
        [half_overlap, half_overlap, z1],
        vel=[0.10, 0.10, 0.0],
        mat=0,
        group=0,
        particleType=2,
        dim=2,
        flaggedSurfaces=[False, False, True, True],
    )
    lower_right = _fixed_normal_box(
        "cornerRightField",
        (-half_overlap, -0.375, z0),
        (0.75, half_overlap, z1),
        velocity=(0.0, 0.0, 0.0),
        material=1,
        group=1,
        flagged_surfaces=(True, False, False, False),
        normal=(-1.0, 0.0, 0.0),
    )
    upper_left = _fixed_normal_box(
        "cornerUpperField",
        (-0.75, -half_overlap, z0),
        (half_overlap, 0.375, z1),
        velocity=(0.0, 0.0, 0.0),
        material=2,
        group=2,
        flagged_surfaces=(False, True, False, False),
        normal=(0.0, -1.0, 0.0),
    )
    pfw["objects"] = [lower_left, lower_right, upper_left]
    return pfw
