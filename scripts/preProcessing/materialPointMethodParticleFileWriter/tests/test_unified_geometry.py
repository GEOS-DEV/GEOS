"""Single-schema geometry contracts; production C++ math has its own runner."""
import sys
from pathlib import Path
import unittest
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pfw_reference_geometry import encode_geometry, geometry_columns, COLUMN_NAMES


class UnifiedGeometryTests(unittest.TestCase):
    def test_box_uses_full_edges_and_explicit_volume(self):
        widths = np.array([.13, .24, .37])
        geometry = encode_geometry([0, 0, 0], widths)
        row = np.array(geometry_columns(geometry, np.prod(widths)).split(), float)
        self.assertEqual(len(row), len(COLUMN_NAMES))
        np.testing.assert_array_equal(row[:3], [2, 0, np.prod(widths)])
        np.testing.assert_array_equal(row[3:].reshape(3, 3), np.diag(widths))

    def test_simplex_shapes_have_no_second_integration_geometry(self):
        for shape, volume in [(1, 12), (2, 4)]:
            g = np.diag([2., 3., 4.])
            row = np.array(geometry_columns((shape, g), volume).split(), float)
            self.assertEqual(len(row), 12)
            np.testing.assert_array_equal(row[:3], [2, shape, volume])
            np.testing.assert_array_equal(row[3:].reshape(3, 3), g)

    def test_rejects_conflicting_or_invalid_data(self):
        invalid = [((0, np.eye(3)), 8), ((2, np.eye(3)), 1),
                   ((0, np.eye(3)), float('nan')), ((2.5, np.eye(3)), 1),
                   ((3, np.eye(3)), 1), ((0, np.full((3, 3), np.nan)), 1),
                   ((0, np.diag([-1, 1, 1])), 1)]
        for geometry, volume in invalid:
            with self.subTest(geometry=geometry, volume=volume):
                with self.assertRaises(ValueError):
                    geometry_columns(geometry, volume)

    def test_sheared_box_has_one_definition(self):
        g = np.array([[1., .2, 0], [0, 1, 0], [0, 0, 1]])
        row = np.array(geometry_columns((0, g), 1).split(), float)
        np.testing.assert_array_equal(row[3:].reshape(3, 3), g)

    def test_wedge_rejects_out_of_plane_edges(self):
        g = np.eye(3); g[0, 2] = .1
        with self.assertRaisesRegex(ValueError, 'in-plane'):
            geometry_columns((1, g), .5)


if __name__ == '__main__':
    unittest.main()
