"""Report final grid stretch, ring reaction and controller-compatible stress.

boxAverageHistory stresses are material-volume averages. The controller uses
full-domain volume including void: multiply by V_material / V_domain to compare.
This example starts with R=H=1, so V_domain=pi*Frr**2*Fyy.
"""
import csv
import math
from pathlib import Path

def last(path):
    with Path(path).open() as stream:
        reader=csv.reader(stream)
        next(reader)
        rows=[[float(value) for value in row] for row in reader if row]
    if not rows:
        raise ValueError(str(path)+' has no data rows')
    return rows[-1]

reaction=last('reactionHistory.csv')
box=last('boxAverageHistory.csv')
volume=math.pi*box[20]**2*box[21]
stress=box[1]*box[17]/volume
print(f'Final time: {reaction[0]:.8g}; grid Frr={reaction[1]:.8g}, Fyy={reaction[2]:.8g}')
print(f'Controller-compatible mean sigma_rr: {stress:.8g} Pa (compression negative)')
print(f'Material means sigma_rr={box[1]:.8g}, sigma_yy={box[2]:.8g}, sigma_hoop={box[3]:.8g} Pa')
print(f'Outer cylindrical wall reaction: {reaction[8]:.8g} N (full-ring radial scalar sum)')
