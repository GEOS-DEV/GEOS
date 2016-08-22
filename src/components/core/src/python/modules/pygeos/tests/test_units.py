
import pygeos
unitManager = pygeos.UnitManager()


# Setup tests
tol = 1e-6
tests = {'scaleTestA': {'unitStruct': ['1.234e1', 'meter'], 'expectedValue': 12.34},
         'scaleTestB': {'unitStruct': ['1.234E1', 'meter'], 'expectedValue': 12.34},
         'scaleTestC': {'unitStruct': ['1.234e+1', 'meter'], 'expectedValue': 12.34},
         'scaleTestD': {'unitStruct': ['1.234e-1', 'meter'], 'expectedValue': 0.1234},
         'scaleTestE': {'unitStruct': ['1.234', 'meter'], 'expectedValue': 1.234},
         'scaleTestF': {'unitStruct': ['2', 'meter'], 'expectedValue': 2.0},
         'prefixTestA': {'unitStruct': ['1.0', 'mumeter'], 'expectedValue': 1.0e-6},
         'prefixTestB': {'unitStruct': ['1.0', 'micrometer'], 'expectedValue': 1.0e-6},
         'prefixTestC': {'unitStruct': ['1.0', 'kilometer'], 'expectedValue': 1.0e3},
         'prefixTestD': {'unitStruct': ['1.0', 'meter'], 'expectedValue': 1.0},
         'prefixTestE': {'unitStruct': ['1.0', 'ms'], 'expectedValue': 1.0e-3},
         'prefixTestF': {'unitStruct': ['1.0', 'millisecond'], 'expectedValue': 1.0e-3},
         'prefixTestG': {'unitStruct': ['1.0', 'Ms'], 'expectedValue': 1.0e9},
         'compoundTestA': {'unitStruct': ['1.0', 'm/s'], 'expectedValue': 1.0},
         'compoundTestB': {'unitStruct': ['1.0', 'micrometer/s'], 'expectedValue': 1.0e-6},
         'compoundTestC': {'unitStruct': ['1.0', 'micrometer/ms'], 'expectedValue': 1.0e-3},
         'compoundTestD': {'unitStruct': ['1.0', 'micrometer/microsecond'], 'expectedValue': 1.0},
         'compoundTestE': {'unitStruct': ['1.0', 'm**2'], 'expectedValue': 1.0},
         'compoundTestF': {'unitStruct': ['1.0', 'km**2'], 'expectedValue': 1.0e6},
         'compoundTestG': {'unitStruct': ['1.0', 'kilometer**2'], 'expectedValue': 1.0e6},
         'compoundTestH': {'unitStruct': ['1.0', '(km*mm)'], 'expectedValue': 1.0},
         'compoundTestI': {'unitStruct': ['1.0', '(km*mm)**2'], 'expectedValue': 1.0},
         'etcTestA': {'unitStruct': ['1.0', 'bbl/day'], 'expectedValue': 0.000001840130728333},
         'etcTestB': {'unitStruct': ['1.0', 'cP'], 'expectedValue': 0.001}}


# Run tests
for k in sorted(tests.keys()):
  convertedValue = unitManager(tests[k]['unitStruct'])
  tmp = float(convertedValue)
  err = (tmp - tests[k]['expectedValue'])/(tmp + tests[k]['expectedValue'])
  if (err >= tol):
    print unitManager(tests[k]['unitStruct'])
    print('Warning: %s did not pass the unit conversion test!' % (k))
    print('         scale: %s, units: %s' % (tests[k]['unitStruct'][0], tests[k]['unitStruct'][1]))
    print('         value: %s, expected: %1.6e\n' % (convertedValue, tests[k]['expectedValue']))

print 'Done!'