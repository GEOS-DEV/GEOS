"""PFW native geometry encoding tests; no GEOS/VisIt runtime is required."""
import sys
from pathlib import Path
import unittest
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pfw_reference_geometry import encode_reference_geometry, geometry_columns, COLUMN_NAMES

class ReferenceGeometryTests(unittest.TestCase):
    def test_box_half_edges(self):
        tag, g = encode_reference_geometry([0,0,0], [2,4,6])
        self.assertEqual(tag, 0)
        np.testing.assert_array_equal(g, np.diag([1,2,3]))

    def test_tetrahedron_round_trip_and_orientation(self):
        v = np.array([[0.,0,0],[2,0,0],[0,3,0],[0,0,4]])
        for order in ([0,1,2,3], [0,2,1,3], [2,3,0,1]):
            q = v[order]
            tag, g = encode_reference_geometry(v.mean(0), [1,1,4], q)
            self.assertEqual(tag, 2)
            self.assertAlmostEqual(np.linalg.det(g)/6, 4)
            a = v.mean(0) - g.sum(0)/4
            reconstructed = np.vstack([a, a+g])
            self.assertTrue(all(np.min(np.linalg.norm(reconstructed-x, axis=1)) < 1e-14 for x in v))

    def test_wedge_centroid_and_thickness(self):
        v = np.array([[0.,0,5],[2,0,5],[0,3,5]])
        tag, g = encode_reference_geometry(v.mean(0), [1,3,0.7], v, plane_strain=True)
        self.assertEqual(tag, 1)
        self.assertAlmostEqual(np.linalg.det(g)/2, 2.1)
        a = v.mean(0) - (g[0]+g[1])/3 - g[2]/2
        bottom = np.array([a,a+g[0],a+g[1]])
        wedge = np.vstack([bottom,bottom+g[2]])
        np.testing.assert_allclose(wedge.mean(0), v.mean(0))
        self.assertAlmostEqual(wedge[:,2].max()-wedge[:,2].min(),0.7)

    def test_geometry_columns_round_trip(self):
        tag,g = encode_reference_geometry([1,2,3],[0.11,0.22,0.33])
        row=geometry_columns((tag,g)).split('\t')
        self.assertEqual(len(row), len(COLUMN_NAMES))
        np.testing.assert_array_equal(np.array(row[1:],float).reshape(3,3),g)

    def test_volume_mismatch_is_not_repaired(self):
        v = np.array([[0.,0,0],[1,0,0],[0,1,0],[0,0,1]])
        with self.assertRaisesRegex(ValueError,'volume'):
            encode_reference_geometry(v.mean(0),[1,1,1],v)

    def test_invalid_centroid_and_nonfinite(self):
        v = np.array([[0.,0,0],[1,0,0],[0,1,0],[0,0,1]])
        with self.assertRaisesRegex(ValueError,'centroid'):
            encode_reference_geometry([0,0,0],[1,1,1/6],v)
        with self.assertRaises(ValueError):
            encode_reference_geometry([float('nan'),0,0],[1,1,1])

    def test_reject_out_of_plane_triangle(self):
        v = np.array([[0.,0,0],[1,0,.1],[0,1,0]])
        with self.assertRaisesRegex(ValueError,'midplane'):
            encode_reference_geometry(v.mean(0),[1,.5,1],v,plane_strain=True)

    def test_tiny_but_resolved_tetrahedron(self):
        v = 1e-7*np.array([[0.,0,0],[1,0,0],[0,1,0],[0,0,1]])
        tag,g = encode_reference_geometry(v.mean(0),[1e-7,1e-7,1e-7/6],v)
        self.assertEqual(tag,2)
        self.assertGreater(np.linalg.det(g),0)

if __name__ == '__main__':
    unittest.main()
