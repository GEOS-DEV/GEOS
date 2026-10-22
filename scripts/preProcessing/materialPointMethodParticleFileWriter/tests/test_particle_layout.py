"""Configuration and reduction contracts for the hex-bulk/tet-surface preset."""
import math
from pathlib import Path
import sys
import unittest
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pfw_particle_layout import (resolve_surface_layout, validate_hex_tet_bulk_type,
                                 ParticleLayoutStatistics, combine_layout_statistics)


class LayoutTests(unittest.TestCase):
    def resolve(self, **kwargs):
        args = dict(method='hexTet', geometry=None, reference_geometry='equivalentBox',
                    simplex_type=None, plane_strain=False)
        args.update(kwargs)
        return resolve_surface_layout(**args)

    def test_preset_and_explicit_point_cut_types(self):
        for requested in (None, 0, 1, 3):
            result = self.resolve(simplex_type=requested)
            self.assertEqual((result.method, result.geometry, result.simplex_particle_type,
                              result.hex_tet), ('simplex', 'fromRefinement',
                                               3 if requested is None else requested, True))

    def test_old_inputs_and_preferred_alias_precedence(self):
        for method in ('cartesian', 'simplex', 'simplex-conforming'):
            result = self.resolve(method=method)
            self.assertEqual(result.geometry, 'box')
            self.assertIsNone(result.simplex_particle_type)
            self.assertFalse(result.hex_tet)
        self.assertEqual(self.resolve(method='simplex', reference_geometry='fromRefinement').geometry,
                         'fromRefinement')
        self.assertEqual(self.resolve(geometry='fromRefinement', reference_geometry_explicit=True).geometry,
                         'fromRefinement')

    def test_conflicting_geometry_and_dimension(self):
        for settings in (dict(geometry='box'), dict(plane_strain=True),
                         dict(reference_geometry_explicit=True)):
            with self.assertRaises(ValueError):
                self.resolve(**settings)
        self.resolve(reference_geometry='fromRefinement', reference_geometry_explicit=True)

    def test_bulk_interpolation(self):
        for kind in (0, 1, 2):
            validate_hex_tet_bulk_type(kind)
        for kind in (-1, 3, 4):
            with self.assertRaises(ValueError):
                validate_hex_tet_bulk_type(kind)

    def test_counts_volumes_and_empty_ranks(self):
        first, second = ParticleLayoutStatistics(), ParticleLayoutStatistics()
        first.record(2, 0, .5)
        for _ in range(100):
            first.record(3, 2, .001)
        second.record(3, 2, .05)
        report = combine_layout_statistics([[], first.report(), [], second.report()])
        self.assertEqual(report['total_count'], 102)
        self.assertAlmostEqual(report['total_volume'], .65)
        box, tet = report['by_type_and_shape']
        self.assertEqual((box['shape_name'], box['count']), ('hex', 1))
        self.assertEqual((tet['shape_name'], tet['count']), ('tet', 101))
        self.assertEqual((tet['min_particle_volume'], tet['max_particle_volume']), (.001, .05))
        self.assertEqual(combine_layout_statistics([[], []])['total_count'], 0)

    def test_invalid_volume(self):
        for volume in (0, -1, math.inf, math.nan):
            with self.assertRaises(ValueError):
                ParticleLayoutStatistics().record(2, 0, volume)


if __name__ == '__main__':
    unittest.main()
