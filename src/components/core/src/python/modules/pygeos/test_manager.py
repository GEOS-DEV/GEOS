import unittest
import re
import os
import filecmp
from . import preprocessGEOSXML, DictRegexHandler, unitManager, symbolicMathRegexHandler, regexConfig
from . import __file__ as modPath


class TestUnitManager(unittest.TestCase):

  def setUp(self):
    self.tol = 1e-6

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


class TestParameterRegex(unittest.TestCase):

  def setUp(self):
    self.regexHandler = DictRegexHandler()
    self.regexHandler.target['foo'] = '1.23'
    self.regexHandler.target['bar'] = '4.56e7'

  def tearDown(self):
    del self.regexHandler

  def evaluateRegex(self, parameterInput, expectedValue):
    result = re.sub(regexConfig.parameters, self.regexHandler, parameterInput)
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


class TestUnitsRegex(unittest.TestCase):

  def setUp(self):
    self.tol = 1e-6

  def evaluateRegex(self, unitInput, expectedValue):
    result = re.sub(regexConfig.units, unitManager.regexHandler, unitInput)
    self.assertEqual(result, expectedValue)

  def test_units_regex_a(self):
    self.evaluateRegex('1.234[m**2/s]', '1.234')

  def test_units_regex_b(self):
    self.evaluateRegex('1.234 [m**2/s]', '1.234')

  def test_units_regex_c(self):
    self.evaluateRegex('1.234[m**2/s]*3.4', '1.234*3.4')

  def test_units_regex_d(self):
    self.evaluateRegex('1.234[m**2/s] + 5.678[mm/s]', '1.234 + 0.005678')

  def test_units_regex_e(self):
    self.evaluateRegex('1.234 [m**2/s] + 5.678 [mm/s]', '1.234 + 0.005678')

  def test_units_regex_f(self):
    self.evaluateRegex('(1.234[m**2/s])*5.678', '(1.234)*5.678')


class TestSymbolicRegex(unittest.TestCase):

  def setUp(self):
    self.tol = 1e-6

  def evaluateRegex(self, symbolicInput, expectedValue):
    result = re.sub(regexConfig.symbolic, symbolicMathRegexHandler, symbolicInput)
    self.assertEqual(result, expectedValue)

  def test_symbolic_regex_a(self):
    self.evaluateRegex('{1.234}', '1.234')

  def test_symbolic_regex_b(self):
    self.evaluateRegex('{1.234*2.0}', '2.468')

  def test_symbolic_regex_c(self):
    self.evaluateRegex('{10}', '10')

  def test_symbolic_regex_d(self):
    self.evaluateRegex('{10*2}', '20')

  def test_symbolic_regex_e(self):
    self.evaluateRegex('{1.0/2.0}', '0.5')

  def test_symbolic_regex_f(self):
    self.evaluateRegex('{2.0**2}', '4.0')

  def test_symbolic_regex_g(self):
    self.evaluateRegex('{1.0 + 2.0**2}', '5.0')

  def test_symbolic_regex_h(self):
    self.evaluateRegex('{(1.0 + 2.0)**2}', '9.0')

  def test_symbolic_regex_i(self):
    self.evaluateRegex('{((1.0 + 2.0)**2)**(0.5)}', '3.0')

  def test_symbolic_regex_j(self):
    self.evaluateRegex('{(1.2e3)*2}', '2400.0')

  def test_symbolic_regex_k(self):
    self.evaluateRegex('{1.2e3*2}', '2400.0')

  @unittest.expectedFailure
  def test_symbolic_regex_l(self):
    self.evaluateRegex('{2.0^2}', '4.0')

  @unittest.expectedFailure
  def test_symbolic_regex_m(self):
    self.evaluateRegex('{sqrt(4.0)}', '2.0')


class TestXMLProcessor(unittest.TestCase):
  def setUp(self):
    self.modPath = os.path.dirname(os.path.abspath(modPath))

  def diff_xml(self, source, target):
    self.assertTrue(filecmp.cmp(source, target))
    os.remove(source)

  def test_basic_xml(self):
    tmp = preprocessGEOSXML(self.modPath + '/tests/source_xml/basic.xml', verbose=0)
    self.diff_xml(tmp, self.modPath + '/tests/target_xml/raw_basic.xml')

  def test_includes_xml(self):
    tmp = preprocessGEOSXML(self.modPath + '/tests/source_xml/includes.xml', verbose=0)
    self.diff_xml(tmp, self.modPath + '/tests/target_xml/raw_includes.xml')

  def test_symbolic_xml(self):
    tmp = preprocessGEOSXML(self.modPath + '/tests/source_xml/symbolic.xml', verbose=0)
    self.diff_xml(tmp, self.modPath + '/tests/target_xml/raw_symbolic.xml')

  def test_formatting_xml(self):
    tmp = preprocessGEOSXML(self.modPath + '/tests/source_xml/symbolic.xml', verbose=0)
    format_xml_file(tmp)
    self.diff_xml(tmp, self.modPath + '/tests/target_xml/target_symbolic.xml')


def runUnitTests():
  print('\nRunning unit manager tests:')
  suite = unittest.TestLoader().loadTestsFromTestCase(TestUnitManager)
  unittest.TextTestRunner(verbosity=2).run(suite)

  print('\nRunning parameter regex handler tests:')
  suite = unittest.TestLoader().loadTestsFromTestCase(TestParameterRegex)
  unittest.TextTestRunner(verbosity=2).run(suite)

  print('\nRunning unit regex handler tests:')
  suite = unittest.TestLoader().loadTestsFromTestCase(TestUnitsRegex)
  unittest.TextTestRunner(verbosity=2).run(suite)

  print('\nRunning symbolic regex handler tests:')
  suite = unittest.TestLoader().loadTestsFromTestCase(TestSymbolicRegex)
  unittest.TextTestRunner(verbosity=2).run(suite)

  print('\nRunning xml processor tests:')
  suite = unittest.TestLoader().loadTestsFromTestCase(TestXMLProcessor)
  unittest.TextTestRunner(verbosity=2).run(suite)
