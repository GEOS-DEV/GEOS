import unittest
import re
import os
import filecmp
from geosx_xml_tools import regex_tools, unit_manager, xml_processor
from geosx_xml_tools.tests import generate_test_xml
import argparse
from parameterized import parameterized

# Create an instance of the unit manager
unitManager = unit_manager.UnitManager()


# Test the unit manager definitions
class TestUnitManager(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.tol = 1e-6

    def test_unit_dict(self):
        unitManager.buildUnits()
        self.assertTrue(bool(unitManager.units))

    # Scale value tests
    @parameterized.expand([['meter', '2', 2.0], ['meter', '1.234', 1.234],
                           ['meter', '1.234e1', 12.34],
                           ['meter', '1.234E1', 12.34],
                           ['meter', '1.234e+1', 12.34],
                           ['meter', '1.234e-1', 0.1234],
                           ['mumeter', '1',
                            1.0e-6], ['micrometer', '1', 1.0e-6],
                           ['kilometer', '1', 1.0e3], ['ms', '1', 1.0e-3],
                           ['millisecond', '1', 1.0e-3], ['Ms', '1', 1.0e6],
                           ['m/s', '1', 1.0], ['micrometer/s', '1', 1.0e-6],
                           ['micrometer/ms', '1', 1.0e-3],
                           ['micrometer/microsecond', '1', 1.0],
                           ['m**2', '1', 1.0], ['km**2', '1', 1.0e6],
                           ['kilometer**2', '1', 1.0e6], ['(km*mm)', '1', 1.0],
                           ['(km*mm)**2', '1', 1.0],
                           ['km^2', '1', 1.0e6, True],
                           ['bbl/day', '1', 0.000001840130728333],
                           ['cP', '1', 0.001]])
    def test_units(self, unit, scale, expected_value, expect_fail=False):
        try:
            val = float(unitManager([scale, unit]))
            self.assertTrue(
                (abs(val - expected_value) < self.tol) != expect_fail)
        except TypeError:
            self.assertTrue(expect_fail)


# Test the behavior of the parameter regex
class TestParameterRegex(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.regexHandler = regex_tools.DictRegexHandler()
        cls.regexHandler.target['foo'] = '1.23'
        cls.regexHandler.target['bar'] = '4.56e7'

    @parameterized.expand(
        [['$:foo*1.234', '1.23*1.234'],
         ['$:foo*1.234/$:bar', '1.23*1.234/4.56e7'],
         ['$:foo*1.234/($:bar*$:foo)', '1.23*1.234/(4.56e7*1.23)'],
         ['$foo*1.234/$bar', '1.23*1.234/4.56e7'],
         ['$foo$*1.234/$bar', '1.23*1.234/4.56e7'],
         ['$foo$*1.234/$bar$', '1.23*1.234/4.56e7'],
         ['$blah$*1.234/$bar$', '1.23*1.234/4.56e7', True],
         ['$foo$*1.234/$bar$', '4.56e7*1.234/4.56e7', True]])
    def test_parameter_regex(self,
                             parameterInput,
                             expectedValue,
                             expect_fail=False):
        try:
            result = re.sub(regex_tools.patterns['parameters'],
                            self.regexHandler, parameterInput)
            self.assertTrue((result == expectedValue) != expect_fail)
        except Exception:
            self.assertTrue(expect_fail)


# Test the behavior of the unit regex
class TestUnitsRegex(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.tol = 1e-6

    @parameterized.expand(
        [['1.234[m**2/s]', '1.234'], ['1.234 [m**2/s]', '1.234'],
         ['1.234[m**2/s]*3.4', '1.234*3.4'],
         ['1.234[m**2/s] + 5.678[mm/s]', '1.234 + 5.678e-3'],
         ['1.234 [m**2/s] + 5.678 [mm/s]', '1.234 + 5.678e-3'],
         ['(1.234[m**2/s])*5.678', '(1.234)*5.678']])
    def test_units_regex(self, unitInput, expectedValue, expect_fail=False):
        try:
            result = re.sub(regex_tools.patterns['units'],
                            unitManager.regexHandler, unitInput)
            self.assertTrue((result == expectedValue) != expect_fail)
        except Exception:
            self.assertTrue(expect_fail)


# Test the symbolic math regex
class TestSymbolicRegex(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.tol = 1e-6

    @parameterized.expand([['`1.234`', '1.234'], ['`1.234*2.0`', '2.468'],
                           ['`10`', '1e1'], ['`10*2`', '2e1'],
                           ['`1.0/2.0`', '5e-1'], ['`2.0**2`', '4'],
                           ['`1.0 + 2.0**2`', '5'], ['`(1.0 + 2.0)**2`', '9'],
                           ['`((1.0 + 2.0)**2)**(0.5)`', '3'],
                           ['`(1.2e3)*2`', '2.4e3'], ['`1.2e3*2`', '2.4e3'],
                           ['`2.0^2`', '4', True], ['`sqrt(4.0)`', '2', True]])
    def test_symbolic_regex(self,
                            symbolicInput,
                            expectedValue,
                            expect_fail=False):
        try:
            result = re.sub(regex_tools.patterns['symbolic'],
                            regex_tools.SymbolicMathRegexHandler,
                            symbolicInput)
            self.assertTrue((result == expectedValue) != expect_fail)
        except Exception:
            self.assertTrue(expect_fail)


# Test the complete xml processor
class TestXMLProcessor(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        generate_test_xml.generate_test_xml_files('.')

    @parameterized.expand(
        [['no_advanced_features_input.xml', 'no_advanced_features_target.xml'],
         ['parameters_input.xml', 'parameters_target.xml'],
         ['included_input.xml', 'included_target.xml'],
         ['symbolic_parameters_input.xml', 'symbolic_parameters_target.xml']])
    def test_xml_processor(self, input_file, target_file, expect_fail=False):
        try:
            tmp = xml_processor.process(input_file,
                                        outputFile=input_file + '.processed',
                                        verbose=0,
                                        keep_parameters=False,
                                        keep_includes=False)
            self.assertTrue(filecmp.cmp(tmp, target_file) != expect_fail)
        except Exception:
            self.assertTrue(expect_fail)


# Main entry point for the unit tests
def run_unit_tests(test_dir, verbose):
    # Create and move to the test directory
    pwd = os.getcwd()
    os.makedirs(test_dir, exist_ok=True)
    os.chdir(test_dir)

    # Unit manager tests
    suite = unittest.TestLoader().loadTestsFromTestCase(TestUnitManager)
    unittest.TextTestRunner(verbosity=verbose).run(suite)

    # Parameter regex handler tests
    suite = unittest.TestLoader().loadTestsFromTestCase(TestParameterRegex)
    unittest.TextTestRunner(verbosity=verbose).run(suite)

    # Regex handler tests
    suite = unittest.TestLoader().loadTestsFromTestCase(TestUnitsRegex)
    unittest.TextTestRunner(verbosity=verbose).run(suite)

    # Symbolic regex handler tests
    suite = unittest.TestLoader().loadTestsFromTestCase(TestSymbolicRegex)
    unittest.TextTestRunner(verbosity=verbose).run(suite)

    # xml processor tests
    suite = unittest.TestLoader().loadTestsFromTestCase(TestXMLProcessor)
    unittest.TextTestRunner(verbosity=verbose).run(suite)

    os.chdir(pwd)


def main():
    """Entry point for the geosx_xml_tools unit tests

    @arg -o/--output Output directory (default = ./test_results)
    """

    # Parse the user arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-t',
                        '--test_dir',
                        type=str,
                        help='Test output directory',
                        default='./test_results')
    parser.add_argument('-v',
                        '--verbose',
                        type=int,
                        help='Verbosity level',
                        default=2)
    args = parser.parse_args()

    # Process the xml file
    run_unit_tests(args.test_dir, args.verbose)


if __name__ == "__main__":
    main()
