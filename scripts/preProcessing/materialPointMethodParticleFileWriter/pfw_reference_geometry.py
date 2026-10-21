"""Canonical affine particle geometry, independent of interpolation.

The two-pass refiner already owns the cut vertices. Encoding them adds no object
queries and does not change particle positions, physical volumes, or rank ownership.
Three row vectors encode a tetrahedron (full edges), a triangular wedge (two full
in-plane edges and a full thickness vector), or a box (full edges).
"""
from __future__ import annotations

import numpy as np

PARALLELEPIPED, TRIANGLE_WEDGE, TETRAHEDRON = 0, 1, 2
COLUMN_NAMES = ("GeometryVersion", "GeometryShape", "ParticleVolume") + tuple(
    "Geometry" + ij for ij in ("XX", "XY", "XZ", "YX", "YY", "YZ", "ZX", "ZY", "ZZ"))


def encode_geometry(point, widths, vertices=None, plane_strain=False):
    """Return (shape, generators) without changing any mechanical input fields.

    widths are the EXISTING equal-volume box full widths. Their product is an
    independent check, not a target used to rescale or repair the actual simplex.
    """
    point = np.asarray(point, dtype=float)
    widths = np.asarray(widths, dtype=float)
    if point.shape != (3,) or widths.shape != (3,):
        raise ValueError("Native reference geometry requires 3D point and widths")
    if not np.all(np.isfinite(point)) or not np.all(np.isfinite(widths)) or np.any(widths <= 0):
        raise ValueError("Nonfinite point or nonpositive particle widths")
    if vertices is None:
        return PARALLELEPIPED, np.diag(widths)

    v = np.asarray(vertices, dtype=float)
    expected = 3 if plane_strain else 4
    if v.shape != (expected, 3) or not np.all(np.isfinite(v)):
        raise ValueError("Expected three triangle or four tetrahedron vertices")
    tol = 256 * np.finfo(float).eps * max(float(np.max(np.abs(v))), float(np.max(widths)), 1.0)
    if np.max(np.abs(v.mean(axis=0) - point)) > tol:
        raise ValueError("Native simplex centroid disagrees with particle center")
    if plane_strain:
        if np.max(np.abs(v[:, 2] - point[2])) > tol:
            raise ValueError("Triangle reference geometry must be in the plane-strain midplane")
        shape = TRIANGLE_WEDGE
        g = np.array([v[1]-v[0], v[2]-v[0], [0.0, 0.0, widths[2]]])
        factor = 0.5
    else:
        shape = TETRAHEDRON
        g = v[1:] - v[0]
        factor = 1.0 / 6.0
    det = float(np.linalg.det(g))
    if det < 0:
        g[[0, 1]] = g[[1, 0]]  # reorder only; no geometric change
        det = -det
    volume = factor * det
    expected_volume = float(np.prod(widths))
    if not np.isfinite(volume) or volume <= 0 or not np.isclose(volume, expected_volume, rtol=1e-9, atol=0):
        raise ValueError("Native simplex volume disagrees with equivalent-box particle volume")
    return shape, g


def geometry_columns(geometry, volume, delimiter="\t"):
    """Schema 2: one shape, physical volume, and nine full-edge components.

    Initial files carry one geometry. GEOS initializes G0=G, then advances G
    independently of constitutive resets. Restarts carry both matrices.
    """
    shape, g = geometry
    if not np.isscalar(shape) or not np.isfinite(shape) or shape not in (0, 1, 2):
        raise ValueError("Unsupported geometry shape")
    g = np.asarray(g, dtype=float)
    if g.shape != (3, 3) or not np.all(np.isfinite(g)):
        raise ValueError("Geometry requires finite 3x3 full-edge generators")
    measure = np.linalg.det(g) * (1.0, 0.5, 1.0/6.0)[int(shape)]
    if not np.isfinite(volume) or volume <= 0 or measure <= 0 or not np.isclose(measure, volume, rtol=1e-9, atol=0):
        raise ValueError("Initial geometry measure must agree with positive ParticleVolume")
    if shape == TRIANGLE_WEDGE:
        tolerance = 1e-12 * np.max(np.abs(g))
        if np.any(np.abs(g[[0,1,2,2], [2,2,0,1]]) > tolerance):
            raise ValueError("Wedge must be an in-plane triangle extruded in z")
    return delimiter.join(["2", str(int(shape)), format(float(volume), ".17g")]
                          + [format(float(x), ".17g") for x in g.flat])
