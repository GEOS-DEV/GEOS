import unittest
import re
import os
import filecmp
from geosx_xml_tools import regex_tools, unit_manager, xml_processor
from geosx_xml_tools.tests import generate_test_xml


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

  def checkUnits(self, unitStruct, expectedValue):
    val = float(unitManager(unitStruct))
    self.assertTrue(abs(val - expectedValue) < self.tol)

  # Scale value tests
  def test_scale_a(self):
    self.checkUnits(['2', 'meter'], 2.0)

  def test_scale_b(self):
    self.checkUnits(['1.234', 'meter'], 1.234)

  def test_scale_c(self):
    self.checkUnits(['1.234e1', 'meter'], 12.34)

  def test_scale_d(self):
    self.checkUnits(['1.234E1', 'meter'], 12.34)

  def test_scale_e(self):
    self.checkUnits(['1.234e+1', 'meter'], 12.34)

  def test_scale_f(self):
    self.checkUnits(['1.234e-1', 'meter'], 0.1234)

  # Unit prefix tests
  def test_prefix_a(self):
    self.checkUnits(['1.0', 'mumeter'], 1.0e-6)

  def test_prefix_b(self):
    self.checkUnits(['1.0', 'micrometer'], 1.0e-6)

  def test_prefix_c(self):
    self.checkUnits(['1.0', 'kilometer'], 1.0e3)

  def test_prefix_d(self):
    self.checkUnits(['1.0', 'meter'], 1.0)

  def test_prefix_e(self):
    self.checkUnits(['1.0', 'ms'], 1.0e-3)

  def test_prefix_f(self):
    self.checkUnits(['1.0', 'millisecond'], 1.0e-3)

  def test_prefix_g(self):
    self.checkUnits(['1.0', 'Ms'], 1.0e6)

  # Compound unit tests
  def test_compound_a(self):
    self.checkUnits(['1.0', 'm/s'], 1.0)

  def test_compound_b(self):
    self.checkUnits(['1.0', 'micrometer/s'], 1.0e-6)

  def test_compound_c(self):
    self.checkUnits(['1.0', 'micrometer/ms'], 1.0e-3)

  def test_compound_d(self):
    self.checkUnits(['1.0', 'micrometer/microsecond'], 1.0)

  def test_compound_e(self):
    self.checkUnits(['1.0', 'm**2'], 1.0)

  def test_compound_f(self):
    self.checkUnits(['1.0', 'km**2'], 1.0e6)

  def test_compound_g(self):
    self.checkUnits(['1.0', 'kilometer**2'], 1.0e6)

  def test_compound_h(self):
    self.checkUnits(['1.0', '(km*mm)'], 1.0)

  def test_compound_i(self):
    self.checkUnits(['1.0', '(km*mm)**2'], 1.0)

  @unittest.expectedFailure
  def test_compound_j(self):
    self.checkUnits(['1.0', 'km^2'], 1.0e6)

  # Etc unit tests
  def test_etcUnits_a(self):
    self.checkUnits(['1.0', 'bbl/day'], 0.000001840130728333)

  def test_etcUnits_b(self):
    self.checkUnits(['1.0', 'cP'], 0.001)


# Test the behavior of the parameter regex
class TestParameterRegex(unittest.TestCase):

  @classmethod
  def setUpClass(cls):
    cls.regexHandler = regex_tools.DictRegexHandler()
    cls.regexHandler.target['foo'] = '1.23'
    cls.regexHandler.target['bar'] = '4.56e7'

  def evaluateRegex(self, parameterInput, expectedValue):
    result = re.sub(regex_tools.patterns['parameters'], self.regexHandler, parameterInput)
    self.assertEqual(result, expectedValue)

  def test_parameter_regex_a(self):
    self.evaluateRegex('$:foo*1.234', '1.23*1.234')

  def test_parameter_regex_b(self):
    self.evaluateRegex('$:foo*1.234/$:bar', '1.23*1.234/4.56e7')

  def test_parameter_regex_c(self):
    self.evaluateRegex('$:foo*1.234/($:bar*$:foo)', '1.23*1.234/(4.56e7*1.23)')

  def test_parameter_regex_d(self):
    self.evaluateRegex('$foo*1.234/$bar', '1.23*1.234/4.56e7')

  def test_parameter_regex_e(self):
    self.evaluateRegex('$foo$*1.234/$bar', '1.23*1.234/4.56e7')

  def test_parameter_regex_f(self):
    self.evaluateRegex('$foo$*1.234/$bar$', '1.23*1.234/4.56e7')

  @unittest.expectedFailure
  def test_parameter_regex_g(self):
    self.evaluateRegex('$blah$*1.234/$bar$', '1.23*1.234/4.56e7')

  @unittest.expectedFailure
  def test_parameter_regex_h(self):
    self.evaluateRegex('$foo$*1.234/$bar$', '4.56e7*1.234/4.56e7')


# Test the behavior of the unit regex
class TestUnitsRegex(unittest.TestCase):

  @classmethod
  def setUpClass(cls):
    cls.tol = 1e-6

  def evaluateRegex(self, unitInput, expectedValue):
    result = re.sub(regex_tools.patterns['units'], unitManager.regexHandler, unitInput)
    self.assertEqual(result, expectedValue)

  def test_units_regex_a(self):
    self.evaluateRegex('1.234[m**2/s]', '1.234')

  def test_units_regex_b(self):
    self.evaluateRegex('1.234 [m**2/s]', '1.234')

  def test_units_regex_c(self):
    self.evaluateRegex('1.234[m**2/s]*3.4', '1.234*3.4')

  def test_units_regex_d(self):
    self.evaluateRegex('1.234[m**2/s] + 5.678[mm/s]', '1.234 + 5.678e-3')

  def test_units_regex_e(self):
    self.evaluateRegex('1.234 [m**2/s] + 5.678 [mm/s]', '1.234 + 5.678e-3')

  def test_units_regex_f(self):
    self.evaluateRegex('(1.234[m**2/s])*5.678', '(1.234)*5.678')


# Test the symbolic math regex
class TestSymbolicRegex(unittest.TestCase):

  @classmethod
  def setUpClass(cls):
    cls.tol = 1e-6

  def evaluateRegex(self, symbolicInput, expectedValue):
    result = re.sub(regex_tools.patterns['symbolic'], regex_tools.SymbolicMathRegexHandler, symbolicInput)
    self.assertEqual(result, expectedValue)

  def test_symbolic_regex_a(self):
    self.evaluateRegex('`1.234`', '1.234')

  def test_symbolic_regex_b(self):
    self.evaluateRegex('`1.234*2.0`', '2.468')

  def test_symbolic_regex_c(self):
    self.evaluateRegex('`10`', '1e1')

  def test_symbolic_regex_d(self):
    self.evaluateRegex('`10*2`', '2e1')

  def test_symbolic_regex_e(self):
    self.evaluateRegex('`1.0/2.0`', '5e-1')

  def test_symbolic_regex_f(self):
    self.evaluateRegex('`2.0**2`', '4')

  def test_symbolic_regex_g(self):
    self.evaluateRegex('`1.0 + 2.0**2`', '5')

  def test_symbolic_regex_h(self):
    self.evaluateRegex('`(1.0 + 2.0)**2`', '9')

  def test_symbolic_regex_i(self):
    self.evaluateRegex('`((1.0 + 2.0)**2)**(0.5)`', '3')

  def test_symbolic_regex_j(self):
    self.evaluateRegex('`(1.2e3)*2`', '2.4e3')

  def test_symbolic_regex_k(self):
    self.evaluateRegex('`1.2e3*2`', '2.4e3')

  @unittest.expectedFailure
  def test_symbolic_regex_l(self):
    self.evaluateRegex('`2.0^2`', '4')

  @unittest.expectedFailure
  def test_symbolic_regex_m(self):
    self.evaluateRegex('`sqrt(4.0)`', '2')


# Test the complete xml processor
class TestXMLProcessor(unittest.TestCase):
  @classmethod
  def setUpClass(cls):
    # Create a directory to hold tests
    cls.pwd = os.getcwd()
    rand_name = xml_processor.generate_random_name(suffix='')[:10]
    # os.makedirs('./geosx_xml_tools_tests_%s/included' % (rand_name), exist_ok=True)
    os.makedirs('./geosx_xml_tools_tests_%s/included' % (rand_name))
    os.chdir('./geosx_xml_tools_tests_%s' % (rand_name))

    # Generate test xml files to process
    generate_test_xml.generate_test_xml_files('.')

  @classmethod
  def tearDownClass(cls):
    os.chdir(cls.pwd)

  def diff_xml(self, source, target):
    self.assertTrue(filecmp.cmp(source, target))
    # os.remove(source)

  def test_no_advanced_xml(self):
    tmp = xml_processor.process('./no_advanced_features_input.xml',
                                outputFile='./no_advanced_features_processed.xml',
                                verbose=0,
                                keep_parameters=False,
                                keep_includes=False)
    self.diff_xml(tmp, './no_advanced_features_target.xml')

  def test_parameters_xml(self):
    tmp = xml_processor.process('./parameters_input.xml',
                                outputFile='./parameters_processed.xml',
                                verbose=0,
                                keep_parameters=False,
                                keep_includes=False)
    self.diff_xml(tmp, './parameters_target.xml')

  def test_includes_xml(self):
    tmp = xml_processor.process('./included_input.xml',
                                outputFile='./included_processed.xml',
                                verbose=0,
                                keep_parameters=False,
                                keep_includes=False)
    self.diff_xml(tmp, './included_target.xml')

  def test_symbolic_xml(self):
    tmp = xml_processor.process('./symbolic_parameters_input.xml',
                                outputFile='./symbolic_parameters_processed.xml',
                                verbose=0,
                                keep_parameters=False,
                                keep_includes=False)
    self.diff_xml(tmp, './symbolic_parameters_target.xml')


# Main entry point for the unit tests
def run_unit_tests():
  # Unit manager tests
  suite = unittest.TestLoader().loadTestsFromTestCase(TestUnitManager)
  unittest.TextTestRunner(verbosity=2).run(suite)

  # Parameter regex handler tests
  suite = unittest.TestLoader().loadTestsFromTestCase(TestParameterRegex)
  unittest.TextTestRunner(verbosity=2).run(suite)

  # Regex handler tests
  suite = unittest.TestLoader().loadTestsFromTestCase(TestUnitsRegex)
  unittest.TextTestRunner(verbosity=2).run(suite)

  # Symbolic regex handler tests
  suite = unittest.TestLoader().loadTestsFromTestCase(TestSymbolicRegex)
  unittest.TextTestRunner(verbosity=2).run(suite)

  # xml processor tests
  suite = unittest.TestLoader().loadTestsFromTestCase(TestXMLProcessor)
  unittest.TextTestRunner(verbosity=2).run(suite)
