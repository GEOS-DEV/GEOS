import unittest
import os
import argparse
import numpy as np
import random
import string
import hdf5_wrapper


def random_string(N):
    return ''.join(random.choices(string.ascii_uppercase + string.ascii_lowercase + string.digits, k=N))


def build_test_dict(depth=0, max_depth=3):
    r = [np.random.randint(2, 20) for x in range(5)]
    test = {'int': np.random.randint(-1000000, 1000000),
            'float': np.random.random(),
            '1d_array': np.random.randn(r[0]),
            '3d_array': np.random.randn(r[1], r[2], r[3]),
            'string': random_string(10),
            'string_array': np.array([random_string(x + 10) for x in range(r[4])])}
    if (depth < max_depth):
        test['child_a'] = build_test_dict(depth + 1, max_depth)
        test['child_b'] = build_test_dict(depth + 1, max_depth)
        test['child_c'] = build_test_dict(depth + 1, max_depth)

    return test


# Test the unit manager definitions
class TestHDF5Wrapper(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.test_dir = 'wrapper_tests'
        os.makedirs(cls.test_dir, exist_ok=True)
        cls.test_dict = build_test_dict()

    def compare_wrapper_dict(self, x, y):
        kx = x.keys()
        ky = y.keys()

        for k in kx:
            if k not in ky:
                raise Exception('y key not in x object (%s)' % (k))

        for k in ky:
            if k not in kx:
                raise Exception('x key not in y object (%s)' % (k))

            vx, vy = x[k], y[k]
            tx, ty = type(vx), type(vy)
            if ((tx != ty) and not (isinstance(vx, (dict, hdf5_wrapper.hdf5_wrapper)) and isinstance(vy, (dict, hdf5_wrapper.hdf5_wrapper)))):
                self.assertTrue(np.issubdtype(tx, ty))

            if isinstance(vx, (dict, hdf5_wrapper.hdf5_wrapper)):
                self.compare_wrapper_dict(vx, vy)
            else:
                if isinstance(vx, np.ndarray):
                    self.assertTrue(np.shape(vx) == np.shape(vy))
                    self.assertTrue((vx == vy).all())
                else:
                    self.assertTrue(vx == vy)

    def test_a_insert_write(self):
        data = hdf5_wrapper.hdf5_wrapper(os.path.join(self.test_dir, 'test_insert.hdf5'), mode='w')
        data.insert(self.test_dict)

    def test_b_manual_write(self):
        data = hdf5_wrapper.hdf5_wrapper(os.path.join(self.test_dir, 'test_manual.hdf5'), mode='w')
        for k, v in self.test_dict.items():
            data[k] = v

    def test_c_link_write(self):
        data = hdf5_wrapper.hdf5_wrapper(os.path.join(self.test_dir, 'test_linked.hdf5'), mode='w')
        for k, v in self.test_dict.items():
            if ('child' in k):
                child_path = os.path.join(self.test_dir, 'test_%s.hdf5' % (k))
                data_child = hdf5_wrapper.hdf5_wrapper(child_path, mode='w')
                data_child.insert(v)
                data.link(k, child_path)
            else:
                data[k] = v

    def test_d_compare_wrapper(self):
        data = hdf5_wrapper.hdf5_wrapper(os.path.join(self.test_dir, 'test_insert.hdf5'))
        self.compare_wrapper_dict(self.test_dict, data)

    def test_e_compare_wrapper_copy(self):
        data = hdf5_wrapper.hdf5_wrapper(os.path.join(self.test_dir, 'test_insert.hdf5'))
        tmp = data.copy()
        self.compare_wrapper_dict(self.test_dict, tmp)

    def test_f_compare_wrapper(self):
        data = hdf5_wrapper.hdf5_wrapper(os.path.join(self.test_dir, 'test_manual.hdf5'))
        self.compare_wrapper_dict(self.test_dict, data)

    def test_g_compare_wrapper(self):
        data = hdf5_wrapper.hdf5_wrapper(os.path.join(self.test_dir, 'test_linked.hdf5'))
        self.compare_wrapper_dict(self.test_dict, data)


def main():
    """Entry point for the geosx_xml_tools unit tests

    Args:
        -v/--verbose (int): Output verbosity
    """

    # Parse the user arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbose', type=int, help='Verbosity level', default=2)
    args = parser.parse_args()

    # Unit manager tests
    suite = unittest.TestLoader().loadTestsFromTestCase(TestHDF5Wrapper)
    unittest.TextTestRunner(verbosity=args.verbose).run(suite)


if __name__ == "__main__":
    main()
