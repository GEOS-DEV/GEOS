"""Compare final material-average stresses to small-strain static references.

This reports differences, not a dynamics pass/fail. Run with the same
AXIS_PRESSURE_CASE used by runProblem. Oscillations and particle resolution
affect the comparison; the cylinder reference assumes axial plane strain.
"""
import csv
import math
import os

with open('boxAverageHistory.csv') as stream:
    reader=csv.reader(stream)
    next(reader)
    rows=[[float(value) for value in row] for row in reader if row]
if not rows:
    raise ValueError('boxAverageHistory.csv has no data rows')
row=rows[-1]
case=os.environ.get('AXIS_PRESSURE_CASE','borehole')
if case=='borehole':
    a,b,pi,po,nu=.01,.02,1e5,2e4,.25
    A=(pi*a*a-po*b*b)/(b*b-a*a)
    B=a*a*b*b*(pi-po)/(b*b-a*a)
    mean_term=2*B*math.log(b/a)/(b*b-a*a)
    reference=(A-mean_term,2*nu*A,A+mean_term)
elif case=='sphere':
    reference=(-2e4,)*3
else:
    raise ValueError('AXIS_PRESSURE_CASE must be borehole or sphere')
print(f'Final time: {row[0]:.8g} s; material volume: {row[17]:.8g} m^3')
for name,measured,expected in zip(('rr','yy','hoop'),row[1:4],reference):
    print(f'Mean sigma_{name}: {measured:.8g} Pa; static reference {expected:.8g} Pa; difference {measured-expected:.8g} Pa')
