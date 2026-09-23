"""Optional native particle geometry, independent of the legacy solver R-vectors.

The two-pass refiner already owns the cut vertices. Encoding them adds no object
queries and does not change particle positions, physical volumes, or rank ownership.
Three row vectors encode a tetrahedron (full edges), a triangular wedge (two full
in-plane edges and a full thickness vector), or a box (half edges).
"""
from __future__ import annotations

import numpy as np

PARALLELEPIPED, TRIANGLE_WEDGE, TETRAHEDRON = 0, 1, 2
COLUMN_NAMES = ("ParticleReferenceGeometry",) + tuple(
    "GeometryRVector" + ij for ij in ("XX", "XY", "XZ", "YX", "YY", "YZ", "ZX", "ZY", "ZZ"))


def encode_reference_geometry(point, widths, vertices=None, plane_strain=False):
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
        return PARALLELEPIPED, np.diag(widths * 0.5)

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


def geometry_columns(geometry, delimiter="\t"):
    """Round-trip-safe precision for the optional columns; legacy columns unchanged."""
    shape, generators = geometry
    return delimiter.join([str(int(shape))] + [format(float(x), ".17g") for x in np.asarray(generators).flat])
