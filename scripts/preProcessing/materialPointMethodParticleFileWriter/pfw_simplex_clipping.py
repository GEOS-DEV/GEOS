"""Small, dependency-light plane cutter used by PFW surface refinement.

The temporary geometry lives in a unit square/cube.  This avoids tolerances
that depend on the units, aspect ratio, or global position of the specimen.
Each half of the same cut is triangulated/tetrahedralized independently, but
uses the same plane and the same fan origin on shared polygonal faces.

This is *not* a general mesh generator.  It cuts one convex box with one plane.
Only NumPy is needed, and no connectivity is retained by the MPM solver.
"""
from __future__ import annotations

import math
import numpy as np

_EPS = 64.0 * np.finfo(float).eps


def _unique(points):
    """Remove repeated vertices, without a geometry-scale-dependent threshold."""
    result = []
    for point in points:
        if not any(np.linalg.norm(point - old, ord=np.inf) <= _EPS for old in result):
            result.append(np.asarray(point, dtype=float))
    return result


def _canonical_polygon(points):
    """Use the lexicographically first vertex as a deterministic fan origin."""
    points = _unique(points)
    if not points:
        return points
    first = min(range(len(points)), key=lambda i: tuple(points[i]))
    return points[first:] + points[:first]


def _clip_polygon(points, normal, offset):
    """Sutherland-Hodgman clipping to normal.dot(x) <= offset."""
    if not points:
        return []
    result = []
    previous = points[-1]
    previous_distance = float(np.dot(normal, previous) - offset)
    if abs(previous_distance) <= _EPS:
        previous_distance = 0.0
    for current in points:
        distance = float(np.dot(normal, current) - offset)
        if abs(distance) <= _EPS:
            distance = 0.0
        if (distance <= 0.0) != (previous_distance <= 0.0):
            fraction = previous_distance / (previous_distance - distance)
            result.append(previous + fraction * (current - previous))
        if distance <= 0.0:
            result.append(current)
        previous, previous_distance = current, distance
    return _canonical_polygon(result)


def _cube_faces():
    """Cyclic vertex order on the six faces of the unit cube."""
    return [np.asarray(face, dtype=float) for face in (
        ((0, 0, 0), (0, 1, 0), (1, 1, 0), (1, 0, 0)),
        ((0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)),
        ((0, 0, 0), (1, 0, 0), (1, 0, 1), (0, 0, 1)),
        ((0, 1, 0), (0, 1, 1), (1, 1, 1), (1, 1, 0)),
        ((0, 0, 0), (0, 0, 1), (0, 1, 1), (0, 1, 0)),
        ((1, 0, 0), (1, 1, 0), (1, 1, 1), (1, 0, 1)),
    )]


def _clip_cube(normal, offset):
    faces = [_clip_polygon(list(face), normal, offset) for face in _cube_faces()]
    faces = [face for face in faces if len(face) >= 3]
    cap = _unique([point for face in faces for point in face
                   if abs(np.dot(normal, point) - offset) <= 4.0 * _EPS])
    if len(cap) >= 3:
        center = np.mean(cap, axis=0)
        axis = np.eye(3)[int(np.argmin(np.abs(normal)))]
        u = np.cross(normal, axis)
        u /= np.linalg.norm(u)
        v = np.cross(normal, u)
        cap.sort(key=lambda point: math.atan2(np.dot(point - center, v),
                                             np.dot(point - center, u)))
        cap = _canonical_polygon(cap)
        # A cut exactly on a cube face has already supplied that face.
        if not any(len(face) == len(cap) and all(
                any(np.linalg.norm(p - q, ord=np.inf) <= _EPS for q in face)
                for p in cap) for face in faces):
            faces.append(cap)
    return faces


def _unit_simplices(normal, offset, dimension):
    if dimension == 2:
        square = [np.array(p, dtype=float) for p in ((0, 0), (1, 0), (1, 1), (0, 1))]
        polygon = _clip_polygon(square, normal, offset)
        for j in range(1, len(polygon) - 1):
            vertices = np.array((polygon[0], polygon[j], polygon[j + 1]))
            edges = vertices[1:] - vertices[0]
            area = abs(float(np.linalg.det(edges))) / 2.0
            if area > 0.0:
                yield vertices, area
    else:
        faces = _clip_cube(normal, offset)
        points = _unique([point for face in faces for point in face])
        if len(points) < 4:
            return
        anchor = min(points, key=lambda point: tuple(point))
        # Coning the boundary triangles to a boundary vertex partitions a convex
        # polyhedron.  Faces through the anchor contribute zero-volume tets.
        for face in faces:
            if any(np.linalg.norm(anchor - vertex, ord=np.inf) <= _EPS for vertex in face):
                continue  # Geometrically zero; do not emit roundoff-sized tets.
            for j in range(1, len(face) - 1):
                vertices = np.array((anchor, face[0], face[j], face[j + 1]))
                volume = abs(float(np.linalg.det(vertices[1:] - anchor))) / 6.0
                if volume > 0.0:
                    yield vertices, volume


def split_box(lower, widths, plane_point, normal, plane_strain=False):
    """Return (negative_simplices, positive_simplices) for one shared plane.

    Each item is (physical_vertices, physical_volume).  Plane-strain items are
    triangles at the midplane and their volume is area times widths[2].  This
    routine does not discard void or assign materials; the caller does that.
    It verifies closure and never rescales weights to hide a failed partition.
    """
    lower = np.asarray(lower, dtype=float)
    widths = np.asarray(widths, dtype=float)
    normal = np.asarray(normal, dtype=float)
    point = np.asarray(plane_point, dtype=float)
    dimension = 2 if plane_strain else 3
    if not np.all(np.isfinite(widths)) or np.any(widths <= 0.0):
        raise ValueError("Box widths must be positive and finite")
    scaled_normal = normal[:dimension] * widths[:dimension]
    scale = np.linalg.norm(scaled_normal)
    if not np.isfinite(scale) or scale <= 0.0:
        raise ValueError("A cutting plane requires a finite, nonzero normal")
    unit_normal = scaled_normal / scale
    offset = float(np.dot(normal[:dimension], (point - lower)[:dimension]) / scale)
    jacobian = float(np.prod(widths))
    sides = []
    unit_volumes = []
    for sign in (1.0, -1.0):
        simplices = []
        for unit_vertices, measure in _unit_simplices(sign * unit_normal, sign * offset, dimension):
            vertices = np.tile(lower, (dimension + 1, 1))
            vertices[:, :dimension] += unit_vertices * widths[:dimension]
            if plane_strain:
                vertices[:, 2] += 0.5 * widths[2]
            volume = measure * jacobian
            if not np.isfinite(volume) or volume <= 0.0:
                raise ValueError("A cut simplex has a nonpositive or nonfinite volume")
            simplices.append((vertices, volume))
            unit_volumes.append(measure)
        sides.append(simplices)
    if abs(math.fsum(unit_volumes) - 1.0) > 1.0e-11:
        raise ValueError("Plane-cut simplex partition does not close to the parent volume")
    return tuple(sides)
