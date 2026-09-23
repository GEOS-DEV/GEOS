"""Run with: python -m unittest discover -s tests -p 'test_surface_refinement.py' -v

These tests need NumPy only.  Rank tests instantiate independent writer slabs;
they test the actual ownership/halo/geometry algorithm, not an MPI transport.
"""
import math
import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pfw_simplex_clipping import split_box
from pfw_surface_refinement import (GeometryQueries, TwoPassRefiner,
                                    owned_x_range, object_slice_bounds)


class HalfSpace:
    def __init__(self, normal, offset, sign=1, material=0, flagged=True):
        self.normal = np.asarray(normal, dtype=float)
        self.normal /= np.linalg.norm(self.normal)
        self.offset = float(offset)
        self.sign = sign
        self.mat = material
        self.group = material
        self.flagged = flagged

    def isInterior(self, point, depth):
        distance = self.sign * (np.dot(self.normal, point) - self.offset)
        if distance > 0.0:
            return -1
        return 5 if self.flagged and abs(distance) <= depth else 0

    def getSurfaceNormal(self, point):
        return self.sign * self.normal

    def getSurfacePosition(self, point):
        return -(np.dot(self.normal, point) - self.offset) * self.normal

    def xMin(self):
        return -np.inf

    def xMax(self):
        return np.inf


class Ball:
    def __init__(self, center, radius, dimension=3, material=0, periodic_x=False):
        self.center = np.array(center, dtype=float)
        self.radius = radius
        self.dimension = dimension
        self.mat = material
        self.group = material
        self.periodic_x = periodic_x

    def delta(self, point):
        delta = np.asarray(point) - self.center
        if self.periodic_x:
            delta[0] = (delta[0] + 0.5) % 1.0 - 0.5
        if self.dimension == 2:
            delta[2] = 0.0
        return delta

    def isInterior(self, point, depth):
        distance = np.linalg.norm(self.delta(point)) - self.radius
        return -1 if distance > 0 else (2 if distance >= -depth else 0)

    def getSurfaceNormal(self, point):
        delta = self.delta(point)
        length = np.linalg.norm(delta)
        return delta / length if length > 0 else np.zeros(3)

    def getSurfacePosition(self, point):
        delta = self.delta(point)
        return (self.radius - np.linalg.norm(delta)) * self.getSurfaceNormal(point)

    def xMin(self):
        return -np.inf if self.periodic_x else self.center[0] - self.radius

    def xMax(self):
        return np.inf if self.periodic_x else self.center[0] + self.radius


class Everywhere:
    def __init__(self, material=1):
        self.mat = material
        self.group = material

    def isInterior(self, point, depth):
        return 0

    def xMin(self):
        return -np.inf

    def xMax(self):
        return np.inf


def generate(objects, counts=(6, 5, 1), method="simplex", level=2,
             rank=0, size=1, periodic=(False, False, False), **kwargs):
    widths = 1.0 / np.asarray(counts)
    plane_strain = counts[2] == 1
    depth = np.linalg.norm(widths[:2 if plane_strain else 3])
    queries = GeometryQueries(objects, depth, plane_strain, sort_objects=True)
    refiner = TwoPassRefiner(queries, [0, 0, 0], counts, widths, periodic,
                            rank=rank, size=size, plane_strain=plane_strain,
                            method=method, level=level, **kwargs)
    pieces = list(refiner.particles())
    return pieces, refiner.report()


def signature(pieces):
    """Include all geometric output and flags, but not rank-local object indices."""
    return [(p.parent_index, tuple(p.point), tuple(p.widths), p.volume,
             p.sample.material, p.sample.group, p.sample.surface_flag,
             p.surface_flag, None if p.normal is None else tuple(p.normal),
             None if p.surface_point is None else tuple(p.surface_point)) for p in pieces]


class ClippingTests(unittest.TestCase):
    def test_random_plane_volume_and_first_moment(self):
        rng = np.random.default_rng(20260916)
        for dimension in (2, 3):
            for _ in range(100):
                normal = rng.normal(size=3)
                if dimension == 2:
                    normal[2] = 0.0
                sides = split_box(np.zeros(3), np.ones(3), rng.random(3), normal,
                                  plane_strain=dimension == 2)
                volume = math.fsum(v for side in sides for _, v in side)
                moment = sum((v * t.mean(axis=0) for side in sides for t, v in side), np.zeros(3))
                self.assertAlmostEqual(volume, 1.0, places=12)
                np.testing.assert_allclose(moment, [0.5, 0.5, 0.5], atol=2e-13, rtol=0)
                self.assertTrue(all(v > 0 for side in sides for _, v in side))

    def test_planes_through_corners_edges_and_faces(self):
        for dimension in (2, 3):
            for normal, point in [([1, 0, 0], [0, 0, 0]), ([1, 0, 0], [1, 0, 0]),
                                  ([1, -1, 0], [0, 0, 0]), ([1, 1, 0], [1, 0, 0]),
                                  ([1, 0, 0], [1e-8, 0, 0])]:
                sides = split_box(np.zeros(3), np.ones(3), point, normal, dimension == 2)
                self.assertAlmostEqual(math.fsum(v for side in sides for _, v in side), 1.0, places=12)
            self.assertAlmostEqual(sum(v for _, v in sides[0]), 1e-8, places=14)

    def test_translation_anisotropy_and_plane_strain_thickness(self):
        lower = np.array([1e4, -5e3, 7.0])
        widths = np.array([0.125, 0.005, 0.003])
        for strain in (True, False):
            sides = split_box(lower, widths, lower + 0.5 * widths, [1, 0, 0], strain)
            volume = sum(v for side in sides for _, v in side)
            self.assertAlmostEqual(volume / np.prod(widths), 1.0, places=12)
            self.assertAlmostEqual(sum(v for _, v in sides[0]) / volume, 0.5, places=12)
            if strain:
                self.assertTrue(all(np.all(t[:, 2] == lower[2] + 0.5 * widths[2])
                                    for side in sides for t, _ in side))


class TwoPassTests(unittest.TestCase):
    def test_recovery_of_empty_center_parents_both_modes(self):
        for method in ("cartesian", "simplex"):
            pieces, report = generate([HalfSpace([1, 0, 0], 0.30)], counts=(4, 4, 1), method=method, level=4)
            self.assertGreater(report["recovered_empty_parents"], 0)
            self.assertTrue(any(p.parent_index[0] == 1 for p in pieces))
            if method == "simplex":
                self.assertAlmostEqual(sum(p.volume for p in pieces), 0.30, places=12)
                self.assertEqual(sum(v for k, v in report.items() if k.startswith("fallback_")), 0)

    def test_slice_filter_keeps_objects_intersecting_empty_parent_extent(self):
        # The x=.375 parent center lies outside the object's bounding interval,
        # but its [.25,.50] reference domain contains part of the circle.
        for method in ("cartesian", "simplex"):
            pieces, report = generate([Ball([0.60, 0.50, 0.5], 0.19, 2)],
                                      counts=(4, 4, 1), method=method, level=4)
            self.assertGreater(report["recovered_empty_parents"], 0)
            self.assertTrue(any(p.parent_index[0] == 1 for p in pieces))

    def test_staged_bulk_records_are_owned_only_and_reused(self):
        calls = []
        def capture(sample):
            calls.append(tuple(sample.point))
            return "staged-at-" + str(tuple(sample.point))
        pieces, report = generate([Everywhere()], counts=(7, 3, 1), method="cartesian",
                                  rank=1, size=3, capture_coarse_record=capture)
        start, stop = owned_x_range(7, 1, 3)
        self.assertEqual(len(calls), (stop - start) * 3)
        self.assertEqual(len(calls), len(pieces))
        self.assertTrue(all(p.coarse_record.startswith("staged-at-") for p in pieces))
        self.assertEqual(report.get("refinement_point_queries", 0), 0)

    def test_unflagged_interfaces_still_refine(self):
        objects = [HalfSpace([1, 0, 0], 0.47, flagged=False), Everywhere()]
        pieces, report = generate(objects, method="cartesian")
        self.assertGreater(report["candidate_domains"], 0)
        self.assertEqual({p.sample.material for p in pieces}, {0, 1})
        self.assertAlmostEqual(sum(p.volume for p in pieces), 1.0, places=12)

    def test_no_geometry_resampling_away_from_boundary(self):
        pieces, report = generate([Everywhere()], method="cartesian", counts=(12, 8, 3))
        self.assertEqual(report.get("candidate_domains", 0), 0)
        self.assertEqual(report.get("refinement_point_queries", 0), 0)
        self.assertEqual(report["coarse_center_queries"], 12 * 8 * 3)
        self.assertEqual(len(pieces), 12 * 8 * 3)
        self.assertLessEqual(report["peak_cached_slices"], 3)

    def test_planar_two_material_partition_2d_and_3d(self):
        for counts, normal in [((5, 4, 1), [1, 0.31, 0]), ((4, 3, 3), [1, 0.31, 0.19])]:
            objects = [HalfSpace(normal, 0.49, 1, 0), HalfSpace(normal, 0.49, -1, 1)]
            pieces, report = generate(objects, counts=counts)
            self.assertAlmostEqual(math.fsum(p.volume for p in pieces), 1.0, places=12)
            self.assertEqual(sum(v for k, v in report.items() if k.startswith("fallback_")), 0)
            for p in pieces:
                self.assertGreater(p.volume, 0)
                self.assertAlmostEqual(np.prod(p.widths) / p.volume, 1.0, places=12)
                if counts[2] == 1:
                    self.assertEqual(p.widths[2], 1.0)
                if p.normal is not None:
                    offset = -np.dot(p.normal, p.point - p.surface_point) * p.normal
                    self.assertAlmostEqual(np.dot(p.normal, p.point + offset - p.surface_point), 0.0, places=13)

    def test_material_refinement_uses_one_complementary_partition(self):
        objects = [HalfSpace([1, 0.2, 0], 0.5), HalfSpace([1, 0.2, 0], 0.5, -1, 1)]
        pieces, _ = generate(objects, counts=(4, 3, 1), material_refinement=lambda mat: mat + 1)
        self.assertAlmostEqual(sum(p.volume for p in pieces), 1.0, places=12)

    def test_curved_free_surface_and_inclusion_are_positive(self):
        for dimension in (2, 3):
            counts = (6, 6, 1 if dimension == 2 else 5)
            for background in (False, True):
                objects = [Ball([0.47, 0.53, 0.51], 0.29, dimension)]
                if background:
                    objects.append(Everywhere())
                pieces, report = generate(objects, counts=counts)
                self.assertGreater(report.get("cut_subdomains", 0), 0)
                self.assertTrue(all(p.volume > 0 for p in pieces))
                if background:
                    self.assertAlmostEqual(sum(p.volume for p in pieces), 1.0, places=11)
                else:
                    exact = math.pi * 0.29**2 if dimension == 2 else 4/3 * math.pi * 0.29**3
                    self.assertLess(abs(sum(p.volume for p in pieces) / exact - 1), 0.08)

    def test_strict_fallback_reports_missing_surface_information(self):
        with self.assertRaisesRegex(ValueError, "no_usable_plane"):
            generate([HalfSpace([1, 0, 0], 0.47, flagged=False)], fallback="error")

    def test_undetected_island_is_not_claimed_to_be_recovered(self):
        # No occupied coarse center or flagged neighbor near this small island.
        pieces, report = generate([Ball([0.5, 0.5, 0.5], 0.01, 2)], counts=(2, 2, 1))
        self.assertEqual(pieces, [])
        self.assertEqual(report.get("candidate_domains", 0), 0)

    def test_invalid_options_fail(self):
        for values in ({"method": "other"}, {"level": 0}, {"level": 1.5},
                       {"fallback": "ignore"}, {"plane_tolerance": 0.8}):
            with self.assertRaises(ValueError):
                generate([Everywhere()], **values)


class RankSlicingTests(unittest.TestCase):
    def test_integer_ownership_has_no_gaps_or_duplicates(self):
        for rows in (1, 5, 11):
            for ranks in (1, 2, 3, 7, 17):
                union = [i for r in range(ranks) for i in range(*owned_x_range(rows, r, ranks))]
                self.assertEqual(union, list(range(rows)))

    def test_factory_bounds_include_actual_owned_rows_and_halo(self):
        for ranks in (2, 3, 8):
            for rank in range(ranks):
                start, stop = owned_x_range(5, rank, ranks)
                lo, hi = object_slice_bounds(0, 1, 5, rank, ranks, halo=True)
                self.assertLessEqual(lo, max(0, start - 1) / 5 + 1e-15)
                self.assertGreaterEqual(hi, min(5, stop + 1) / 5 - 1e-15)
        self.assertEqual(object_slice_bounds(0, 1, 9, 0, 3, halo=True, periodic=True), (0, 1))

    def test_serial_equals_independent_rank_slabs_and_empty_ranks(self):
        for counts in ((5, 4, 1), (5, 3, 3)):
            dim = 2 if counts[2] == 1 else 3
            for method in ("cartesian", "simplex"):
                objects = [Ball([0.57, 0.51, 0.49], 0.30, dim), Everywhere()]
                serial, _ = generate(objects, counts=counts, method=method)
                for ranks in (2, 3, 7):
                    combined = []
                    for rank in range(ranks):
                        lo, hi = object_slice_bounds(0, 1, counts[0], rank, ranks, halo=True)
                        local_objects = [o for o in objects if o.xMin() <= hi and o.xMax() >= lo]
                        pieces, _ = generate(local_objects, counts=counts, method=method, rank=rank, size=ranks)
                        combined.extend(pieces)
                    self.assertEqual(signature(serial), signature(combined))

    def test_periodic_seam_and_translated_surface_markers(self):
        for method in ("cartesian", "simplex"):
            objects = [Ball([0.96, 0.51, 0.5], 0.22, 2, periodic_x=True)]
            serial, _ = generate(objects, method=method, counts=(7, 6, 1), periodic=(True, False, False))
            combined = []
            for rank in range(3):
                pieces, _ = generate(objects, method=method, counts=(7, 6, 1), rank=rank, size=3,
                                     periodic=(True, False, False))
                combined.extend(pieces)
            self.assertEqual(signature(serial), signature(combined))
            self.assertTrue(any(p.point[0] < 0.1 for p in serial))
            self.assertTrue(any(p.point[0] > 0.9 for p in serial))


if __name__ == "__main__":
    unittest.main()
