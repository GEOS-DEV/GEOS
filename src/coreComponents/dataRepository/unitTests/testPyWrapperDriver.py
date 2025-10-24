"""
------------------------------------------------------------------------------------------------------------
SPDX-License-Identifier: LGPL-2.1-only

Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
Copyright (c) 2018-2024 TotalEnergies
Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
Copyright (c) 2023-2024 Chevron
Copyright (c) 2019-     GEOS/GEOSX Contributors
All rights reserved

See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
------------------------------------------------------------------------------------------------------------
"""

import unittest
import testPyWrapper


class TestPyWrapperSetValue(unittest.TestCase):
    
    def test_set_string_value(self):
        wrapper = testPyWrapper.get_string_wrapper()
        wrapper.set_value("hello")
        self.assertEqual(wrapper.value(), "hello")
    
    def test_set_int_value(self):
        wrapper = testPyWrapper.get_int_wrapper()
        wrapper.set_value(42)
        self.assertEqual(wrapper.value(), 42)
    
    def test_set_long_value(self):
        wrapper = testPyWrapper.get_long_wrapper()
        wrapper.set_value(123456789)
        self.assertEqual(wrapper.value(), 123456789)
    
    def test_set_double_value(self):
        wrapper = testPyWrapper.get_double_wrapper()
        wrapper.set_value(3.14)
        self.assertAlmostEqual(wrapper.value(), 3.14, places=5)
    
    def test_type_mismatch(self):
        wrapper = testPyWrapper.get_int_wrapper()
        with self.assertRaises(TypeError):
            wrapper.set_value("not an int")
    
    def test_unsupported_type(self):
        wrapper = testPyWrapper.get_int_wrapper()
        with self.assertRaises(TypeError):
            wrapper.set_value({"key": "value"})

if __name__ == '__main__':
    unittest.main()