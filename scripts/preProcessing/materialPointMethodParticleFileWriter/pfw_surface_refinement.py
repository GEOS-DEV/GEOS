"""Two-pass, narrow-band surface refinement for particleFileWriter.

Pass 1 queries ordinary coarse particle centers.  Occupancy, object/material/
contact-group labels, and surface flags select the domains needing more work.
Pass 2 replaces only those domains, including empty-center domains next to
material.  The passes are pipelined through three x-slices, so memory does not
scale with the rank's entire particle count and interior geometry is not
queried again merely to find the surface band.

Ranks own coarse domains, not refined particle centers.  A read-only one-row
halo is *requeried* on each rank; it is never emitted.  Thus no particle lists,
Python geometry objects, or global occupancy arrays need MPI communication.
The caller must supply objects covering the owned domains plus this halo.
Objects and property functions must be deterministic functions of position.

Geometry objects are unchanged.  The simplex path uses their existing
isInterior / getSurfaceNormal / getSurfacePosition methods.  It is a local
single-plane approximation, not a globally conforming mesh of the object.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from itertools import product
import math
from typing import Optional

import numpy as np

from pfw_simplex_clipping import split_box


def owned_x_range(number_of_rows, rank, size):
    """Zero-based, half-open ownership; matches the writer's legacy ceil rule."""
    if number_of_rows < 1 or size < 1 or not 0 <= rank < size:
        raise ValueError("Invalid lattice size or MPI rank")
    return ((number_of_rows * rank + size - 1) // size,
            (number_of_rows * (rank + 1) + size - 1) // size)


def object_slice_bounds(xmin, xmax, number_of_rows, rank, size, halo=False, periodic=False):
    """Bounds passed to make_objects(xmin, xmax), including any read-only halo.

    A periodic x-halo is disjoint at the seam.  The existing two-scalar factory
    cannot express disjoint intervals, so use the full x extent in that case.
    This costs object memory, but preserves order/priority without new APIs.
    """
    if halo and periodic:
        return float(xmin), float(xmax)
    start, stop = owned_x_range(number_of_rows, rank, size)
    padding = 1 if halo else 0
    dx = (xmax - xmin) / number_of_rows
    return (xmin + max(0, start - padding) * dx,
            xmin + min(number_of_rows, stop + padding) * dx)


@dataclass
class Sample:
    """Geometry result with optional staged ordinary-pass particle text."""
    point: np.ndarray
    object_index: int
    material: int
    group: int
    surface_flag: int
    surface_checked: bool = False
    surface: Optional[tuple] = None  # (outward unit normal, absolute surface point)
    coarse_record: Optional[str] = None  # Optional ordinary-pass output, without ID.

    @property
    def region(self):
        return self.object_index, self.material, self.group


def region(sample):
    return None if sample is None else sample.region


def _integer(value, name):
    number = float(value)
    if not math.isfinite(number) or number != int(number):
        raise ValueError(name + " must be a finite integer")
    return int(number)


def _property(obj, method, attribute, point, default):
    if hasattr(obj, method):
        return getattr(obj, method)(point)
    value = getattr(obj, attribute, default)
    return value(obj, point) if callable(value) else value


class GeometryQueries:
    """Adapter around existing geometry methods, preserving first-object priority."""

    def __init__(self, objects, surface_depth, plane_strain=False, sort_objects=False):
        self.objects = objects
        self.surface_depth = float(surface_depth)
        self.plane_strain = bool(plane_strain)
        self.sort_objects = bool(sort_objects)
        self.stats = Counter()

    def objects_in_interval(self, lower, upper):
        if not self.sort_objects:
            return tuple(range(len(self.objects)))
        # Refined queries use the full parent interval, NOT its center x.  An
        # object wholly inside an empty-center parent's edge must remain visible.
        return tuple(i for i, obj in enumerate(self.objects)
                     if obj.xMin() <= upper and obj.xMax() >= lower)

    def query(self, point, indices, coarse=False):
        point = np.asarray(point, dtype=float)
        self.stats["coarse_center_queries" if coarse else "refinement_point_queries"] += 1
        for index in indices:
            obj = self.objects[index]
            self.stats["isInterior_calls"] += 1
            flag = _integer(obj.isInterior(point, self.surface_depth), "Surface flag")
            if flag < 0:
                continue
            # Query these immediately: some compound objects cache their last
            # membership result for getMat/getGroup.  Do not sample full particle
            # state here (in particular, do not consume random strength draws).
            material = _integer(_property(obj, "getMat", "mat", point, 0), "Material")
            group = _integer(_property(obj, "getGroup", "group", point, 0), "Contact group")
            if material < 0:
                return None  # Preserve the legacy first-object exclusion rule.
            return Sample(point.copy(), index, material, group, flag)
        return None

    def get_surface(self, sample):
        """Lazily cache only the surface information actually used in pass 2."""
        if sample.surface_checked:
            return sample.surface
        sample.surface_checked = True
        if sample.surface_flag == 0:
            return None
        obj = self.objects[sample.object_index]
        if not hasattr(obj, "getSurfaceNormal") or not hasattr(obj, "getSurfacePosition"):
            return None
        # Restore any object's point-local cache before its surface getters.
        self.stats["isInterior_calls"] += 1
        if obj.isInterior(sample.point, self.surface_depth) < 0:
            return None
        self.stats["surface_queries"] += 1
        normal = np.asarray(obj.getSurfaceNormal(sample.point), dtype=float).copy()
        offset = np.asarray(obj.getSurfacePosition(sample.point), dtype=float).copy()
        if normal.shape != (3,) or offset.shape != (3,):
            return None
        if not np.all(np.isfinite(normal)) or not np.all(np.isfinite(offset)):
            return None
        if self.plane_strain:
            normal[2] = 0.0
            offset[2] = 0.0
        length = float(np.linalg.norm(normal))
        if length <= 1.0e-14:
            return None
        sample.surface = (normal / length, sample.point + offset)
        return sample.surface


@dataclass
class ParticlePiece:
    """One final material point, ready for the writer's existing field evaluator."""
    point: np.ndarray
    widths: np.ndarray
    volume: float
    sample: Sample
    parent_index: tuple
    # A cut simplex carries its adopted plane for consistent particle-relative
    # surface offsets.  Ordinary Cartesian/fallback pieces have no plane.
    normal: Optional[np.ndarray] = None
    surface_point: Optional[np.ndarray] = None
    surface_flag: Optional[int] = None
    vertices: Optional[np.ndarray] = None
    coarse_record: Optional[str] = None


class TwoPassRefiner:
    """Generate replacement particles in a deterministic global parent order."""

    def __init__(self, queries, lower, counts, widths, periodic,
                 rank=0, size=1, plane_strain=False, method="cartesian", level=1,
                 material_refinement=None, fallback="cartesian", plane_tolerance=0.05,
                 capture_coarse_record=None):
        self.queries = queries
        self.lower = np.asarray(lower, dtype=float)
        self.counts = tuple(int(n) for n in counts)
        self.widths = np.asarray(widths, dtype=float)
        self.periodic = tuple(bool(p) for p in periodic)
        self.plane_strain = bool(plane_strain)
        self.dimension = 2 if self.plane_strain else 3
        self.method = method.lower().replace("-conforming", "")
        self.level = _integer(level, "Surface refinement level")
        self.fallback = fallback.lower()
        self.plane_tolerance = float(plane_tolerance)
        self.material_refinement = material_refinement or (lambda material: 1)
        self.capture_coarse_record = capture_coarse_record
        if self.method not in ("cartesian", "simplex"):
            raise ValueError("surfaceRefinementMethod must be 'cartesian' or 'simplex'")
        if self.level < 1 or any(n < 1 for n in self.counts):
            raise ValueError("Refinement level and coarse counts must be positive")
        if self.fallback not in ("cartesian", "error"):
            raise ValueError("surfaceRefinementFallback must be 'cartesian' or 'error'")
        if not 0.0 <= self.plane_tolerance < 0.5:
            raise ValueError("surfaceRefinementPlaneTolerance must be in [0, 0.5)")
        if np.any(self.widths <= 0.0) or not np.all(np.isfinite(self.widths)):
            raise ValueError("Coarse particle widths must be positive and finite")
        if self.plane_strain and self.counts[2] != 1:
            raise ValueError("Plane-strain initialization requires one z particle layer")
        self.start, self.stop = owned_x_range(self.counts[0], rank, size)
        self.slices = {}  # x index -> sparse {(j,k): Sample}; missing entries are void
        self.stats = Counter()
        self.neighbor_offsets = tuple(product((-1, 0, 1), (-1, 0, 1),
                                              (0,) if self.plane_strain else (-1, 0, 1)))

    def _canonical_index(self, index):
        result = []
        for d, value in enumerate(index):
            if 0 <= value < self.counts[d]:
                result.append(value)
            elif self.periodic[d] and not (self.plane_strain and d == 2):
                result.append(value % self.counts[d])
            else:
                return None
        return tuple(result)

    def _prepare_slices(self, i):
        needed = set()
        for di in (-1, 0, 1):
            key = self._canonical_index((i + di, 0, 0))
            if key is not None:
                needed.add(key[0])
        for old in list(self.slices):
            if old not in needed:
                del self.slices[old]
        for row in sorted(needed):
            if row in self.slices:
                continue
            x = self.lower[0] + (row + 0.5) * self.widths[0]
            objects = self.queries.objects_in_interval(x, x)
            samples = {}
            for j in range(self.counts[1]):
                for k in range(self.counts[2]):
                    point = self.lower + (np.array((row, j, k)) + 0.5) * self.widths
                    sample = self.queries.query(point, objects, coarse=True)
                    if sample is not None:
                        # The real writer can stage its ordinary first-pass text
                        # here.  Reuse it for untouched parents: no second
                        # membership/property query in the bulk.  Halo samples
                        # never evaluate particle state or create output records.
                        if self.capture_coarse_record is not None and self.start <= row < self.stop:
                            sample.coarse_record = self.capture_coarse_record(sample)
                        samples[j, k] = sample
            self.slices[row] = samples
        self.stats["peak_cached_slices"] = max(self.stats["peak_cached_slices"], len(self.slices))
        self.stats["peak_cached_occupied_centers"] = max(
            self.stats["peak_cached_occupied_centers"], sum(len(s) for s in self.slices.values()))

    def _neighborhood(self, index):
        """Return samples plus periodic image shifts, including the center."""
        neighbors = []
        for delta in self.neighbor_offsets:
            raw = tuple(index[d] + delta[d] for d in range(3))
            key = self._canonical_index(raw)
            if key is None:
                continue  # Outside generation bounds is not an empty target cell.
            sample = self.slices[key[0]].get(key[1:])
            shift = (np.array(raw) - np.array(key)) * self.widths
            neighbors.append((sample, shift))
        return neighbors

    def _is_candidate(self, parent, neighbors):
        for sample, _ in neighbors:
            if sample is not None and sample.surface_flag != 0:
                return True
            if region(sample) != region(parent):
                return True  # Includes BOTH sides of occupied/void and A/B pairs.
        return False

    def _factor(self, parent, neighbors, candidate):
        samples = [s for s, _ in neighbors if s is not None] if candidate else [parent]
        base = max((_integer(self.material_refinement(s.material), "Particle refinement")
                    for s in samples if s is not None), default=1)
        if base < 1:
            raise ValueError("particleRefinement values must be at least 1")
        # One common subdivision for a candidate domain: independently using
        # A's and B's factors could duplicate or omit part of the same domain.
        return base * self.level if candidate else base

    def _cartesian_piece(self, lower, widths, objects, index, sample=None):
        point = lower + 0.5 * widths
        if sample is None:
            sample = self.queries.query(point, objects)
        if sample is None:
            return []
        return [ParticlePiece(point, widths.copy(), float(np.prod(widths)), sample, index)]

    def _fallback(self, reason, lower, widths, objects, index, center):
        self.stats["fallback_" + reason] += 1
        if self.fallback == "error":
            raise ValueError("Simplex surface refinement: " + reason +
                             " at parent " + str(index) + ", subdomain center " +
                             str((lower + 0.5 * widths).tolist()))
        return self._cartesian_piece(lower, widths, objects, index, center)

    def _simplex_pieces(self, lower, widths, objects, index, neighbors):
        """Cut one refined subdomain; fall back as a unit, never drop one side."""
        center_point = lower + 0.5 * widths
        center = self.queries.query(center_point, objects)
        points = [center_point]
        samples = [center]
        for corner in product((0.0, 1.0), repeat=self.dimension):
            point = lower.copy()
            point[:self.dimension] += np.array(corner) * widths[:self.dimension]
            if self.plane_strain:
                point[2] += 0.5 * widths[2]
            points.append(point)
            # Geometry objects commonly use half-open bounds.  Query one ULP
            # toward this subdomain's interior so a corner exactly on xmax is
            # not mistaken for an extra void region.  The cutter still uses the
            # exact corner coordinate stored above.
            samples.append(self.queries.query(np.nextafter(point, center_point), objects))
        labels = {region(s) for s in samples}
        if len(labels) > 2:
            return self._fallback("multiple_regions", lower, widths, objects, index, center)

        # Surface data is fetched only here, never for every object/domain.
        # Prefer a local query, then a nearby coarse surface marker.  Periodic
        # markers keep their canonical query point but translate the plane.
        seeds = [(s, np.zeros(3)) for s in samples if s is not None and s.surface_flag]
        seeds += [(s, shift) for s, shift in neighbors if s is not None and s.surface_flag]
        seeds.sort(key=lambda item: float(np.dot(item[0].point + item[1] - center_point,
                                                item[0].point + item[1] - center_point)))
        plane = None
        distances = None
        scale = float(np.linalg.norm(widths[:self.dimension]))
        epsilon = 128.0 * np.finfo(float).eps * scale
        for source, shift in seeds:
            surface = self.queries.get_surface(source)
            if surface is None:
                continue
            normal, surface_point = surface
            surface_point = surface_point + shift
            trial = np.array([np.dot(normal, p - surface_point) for p in points])
            if trial.min() < -epsilon and trial.max() > epsilon:
                plane = source, normal, surface_point
                distances = trial
                break
            if len(labels) == 1:
                # A usable nearest-surface plane does not cut this uniform
                # subdomain.  Do not query every other flagged marker just to
                # reconfirm that fact.  The single-local-surface assumption is
                # the same one used by the actual cutter.
                self.stats["uncut_subdomains"] += 1
                return self._cartesian_piece(lower, widths, objects, index, center)
        if plane is None:
            if len(labels) == 1:
                self.stats["uncut_subdomains"] += 1
                return self._cartesian_piece(lower, widths, objects, index, center)
            return self._fallback("no_usable_plane", lower, widths, objects, index, center)

        source, normal, surface_point = plane
        negative = samples[int(np.argmin(distances))]
        positive = samples[int(np.argmax(distances))]
        if region(negative) == region(positive):
            if len(labels) == 1:
                return self._cartesian_piece(lower, widths, objects, index, center)
            return self._fallback("unresolved_sides", lower, widths, objects, index, center)
        if region(negative) != source.region:
            return self._fallback("inconsistent_normal", lower, widths, objects, index, center)
        for sample, distance in zip(samples, distances):
            expected = negative if distance < 0.0 else positive
            # Small disagreement is expected for a tangent to a curved surface.
            # Large disagreement means one plane cannot represent this subdomain.
            if abs(distance) > max(epsilon, self.plane_tolerance * scale):
                if region(sample) != region(expected):
                    return self._fallback("inconsistent_plane", lower, widths, objects, index, center)
        try:
            sides = split_box(lower, widths, surface_point, normal, self.plane_strain)
        except ValueError:
            return self._fallback("invalid_partition", lower, widths, objects, index, center)

        pieces = []
        for sign, side, owner in zip((1.0, -1.0), sides, (negative, positive)):
            if owner is None:
                continue  # Void remains void; do not renormalize material volume.
            matching = [s for s in samples if region(s) == owner.region and s.surface_flag]
            flag = source.surface_flag if owner.region == source.region else (
                matching[0].surface_flag if matching else 0)
            for vertices, volume in side:
                point = vertices.mean(axis=0)
                actual = self.queries.query(point, objects)
                # Property getters can be undefined outside an object's actual
                # region.  Reject the WHOLE cut if one centroid has the wrong
                # region; never filter individual simplices and leave holes.
                if region(actual) != owner.region:
                    return self._fallback("centroid_region", lower, widths, objects, index, center)
                equivalent_widths = widths.copy()
                factor = (volume / float(np.prod(widths))) ** (1.0 / self.dimension)
                equivalent_widths[:self.dimension] *= factor
                pieces.append(ParticlePiece(point, equivalent_widths, volume, actual, index,
                                            sign * normal, surface_point.copy(), flag, vertices))
        self.stats["cut_subdomains"] += 1
        self.stats["simplex_particles"] += len(pieces)
        return pieces

    def particles(self, progress=None):
        """Yield final particles; owned parent rows alone may produce output."""
        for i in range(self.start, self.stop):
            self._prepare_slices(i)
            if progress is not None:
                progress(i)
            for j in range(self.counts[1]):
                for k in range(self.counts[2]):
                    index = (i, j, k)
                    parent = self.slices[i].get((j, k))
                    neighbors = self._neighborhood(index)
                    candidate = self._is_candidate(parent, neighbors)
                    self.stats["owned_coarse_domains"] += 1
                    if candidate:
                        self.stats["candidate_domains"] += 1
                    if parent is None and not candidate:
                        continue
                    lower = self.lower + np.array(index) * self.widths
                    factor = self._factor(parent, neighbors, candidate)
                    child_widths = self.widths.copy()
                    child_widths[:self.dimension] /= factor
                    objects = self.queries.objects_in_interval(lower[0], lower[0] + self.widths[0])
                    emitted = False
                    for child in product(range(factor), range(factor),
                                         range(1 if self.plane_strain else factor)):
                        child_lower = lower + np.array(child) * child_widths
                        if candidate and self.method == "simplex":
                            pieces = self._simplex_pieces(child_lower, child_widths, objects, index, neighbors)
                        elif not candidate and factor == 1:
                            # Exactly the original particle; no duplicate geometry query.
                            pieces = self._cartesian_piece(child_lower, child_widths, objects, index, parent)
                            pieces[0].coarse_record = parent.coarse_record
                        else:
                            pieces = self._cartesian_piece(child_lower, child_widths, objects, index)
                        for piece in pieces:
                            emitted = True
                            self.stats["emitted_particles"] += 1
                            yield piece
                    if emitted and parent is None:
                        self.stats["recovered_empty_parents"] += 1
                    if candidate and parent is not None:
                        self.stats["replaced_occupied_parents"] += 1

    def report(self):
        return dict(self.stats + self.queries.stats)
