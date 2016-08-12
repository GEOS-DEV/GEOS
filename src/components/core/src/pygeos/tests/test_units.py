
import sys
sys.path.append('../../')
import pygeos

unitManager = pygeos.UnitManager()
print unitManager(['1.234e1', 'meter**2/s'])