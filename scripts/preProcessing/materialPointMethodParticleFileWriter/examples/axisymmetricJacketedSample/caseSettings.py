"""Small, editable SI-unit test definition; no GEOS/MPI imports are needed.

These are illustrative elastic properties, not a material calibration. The
sample is an effective porous continuum, without explicit pores or pore fluid.
"""
import json
import math
import os
from pathlib import Path


def settings():
    variant = os.environ.get('JACKET_CASE', 'jacketed')
    if variant not in ('jacketed', 'homogeneous'):
        raise ValueError('JACKET_CASE must be jacketed or homogeneous')
    refine = int(os.environ.get('JACKET_REFINE', '1'))
    if refine < 1:
        raise ValueError('JACKET_REFINE must be a positive integer')
    pressure = float(os.environ.get('JACKET_PRESSURE', '10000'))
    ramp_time = float(os.environ.get('JACKET_RAMP_TIME', '0.001'))
    end_time = float(os.environ.get('JACKET_END_TIME', '0.002'))
    if not all(math.isfinite(v) for v in (pressure, ramp_time, end_time)):
        raise ValueError('Pressure and times must be finite')
    if pressure <= 0 or not 0 < ramp_time < end_time:
        raise ValueError('Need pressure > 0 and 0 < ramp time < end time')
    # density [kg/m3], Young modulus [Pa], Poisson ratio [-].
    materials = [dict(name='sample', density=800., E=50.e6, nu=0.20),
                 dict(name='steel', density=7850., E=200.e9, nu=0.30),
                 dict(name='rubber', density=1100., E=5.e6, nu=0.45)]
    if variant == 'homogeneous':
        for material in materials[1:]:
            material.update({key: materials[0][key] for key in ('density', 'E', 'nu')})
    return dict(variant=variant, refine=refine, pressure=pressure,
                ramp_time=ramp_time, end_time=end_time,
                radius=0.010, sample_half_height=0.010,
                plate_thickness=0.002, jacket_thickness=0.001,
                cell_width=0.001/refine, materials=materials,
                tracer_labels=['sample center', 'sample near plate', 'upper plate',
                               'side jacket', 'upper jacket'])


def pressure_at(time, case):
    """Cosine pressure ramp followed by a hold, matching interpType=1."""
    phase = min(1., max(0., float(time)/case['ramp_time']))
    return 0.5*case['pressure']*(1.-math.cos(math.pi*phase))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', type=Path, required=True)
    args = parser.parse_args()
    args.write.write_text(json.dumps(settings(), indent=2)+'\n')
