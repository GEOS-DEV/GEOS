"""PFW CPTI box tessellation checks against independent point ownership."""
import sys
from pathlib import Path
import unittest
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pfw_reference_geometry import box_tetrahedra, encode_geometry, geometry_columns


class CPTIGeometryTests(unittest.TestCase):
    def test_volume_centroid_and_partition(self):
        rng = np.random.default_rng(81023)
        for _ in range(30):
            center = rng.uniform(-10, 10, 3)
            widths = rng.uniform(.1, 2, 3)
            tets = list(box_tetrahedra(center, widths))
            self.assertEqual(len(tets), 6)
            volume = 0.
            first = np.zeros(3)
            samples = center + rng.uniform(-.5, .5, (200, 3))*widths
            owners = np.zeros(len(samples), dtype=int)
            for point, equal_widths, vertices in tets:
                shape, edges = encode_geometry(point, equal_widths, vertices)
                self.assertEqual(shape, 2)
                v = np.linalg.det(edges)/6
                self.assertGreater(v, 0)
                self.assertAlmostEqual(v, np.prod(widths)/6)
                volume += v
                first += v*point
                geometry_columns((shape, edges), np.prod(equal_widths))
                coords = np.linalg.solve((vertices[1:]-vertices[0]).T, (samples-vertices[0]).T).T
                owners += (coords.min(axis=1) >= 0) & (coords.sum(axis=1) <= 1)
            np.testing.assert_allclose(volume, np.prod(widths), rtol=2e-13)
            np.testing.assert_allclose(first/volume, center, atol=2e-13)
            np.testing.assert_array_equal(owners, np.ones(len(samples)))

    def test_shared_box_faces(self):
        def faces(center):
            from itertools import combinations
            answer = []
            for _, _, vertices in box_tetrahedra(center, [1, 1, 1]):
                for face in combinations(vertices, 3):
                    if all(abs(v[0]-.5) < 1e-14 for v in face):
                        answer.append(tuple(sorted(tuple(v) for v in face)))
            return sorted(answer)
        self.assertEqual(faces([0, 0, 0]), faces([1, 0, 0]))
        self.assertEqual(len(faces([0, 0, 0])), 2)

    def test_invalid_box(self):
        for point, widths in [([0, 0], [1, 1, 1]), ([0, 0, 0], [1, 0, 1]),
                              ([0, np.nan, 0], [1, 1, 1])]:
            with self.assertRaises(ValueError):
                list(box_tetrahedra(point, widths))


if __name__ == '__main__':
    unittest.main()
