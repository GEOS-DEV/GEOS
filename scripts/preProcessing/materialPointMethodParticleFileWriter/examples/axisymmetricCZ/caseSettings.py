"""SI units; AXIS_CZ_CASE=elastic|debond and AXIS_CZ_UPDATE=FLIP|FMPM."""
import argparse
import json
import os
from pathlib import Path


def settings():
    mode = os.environ.get('AXIS_CZ_CASE', 'elastic')
    update = os.environ.get('AXIS_CZ_UPDATE', 'FLIP')
    displacement = os.environ.get('AXIS_CZ_DISPLACEMENT', 'Nodal')
    if mode not in ('elastic', 'debond') or update not in ('FLIP', 'FMPM'):
        raise ValueError('Use AXIS_CZ_CASE=elastic|debond and AXIS_CZ_UPDATE=FLIP|FMPM')
    if displacement not in ('Nodal', 'TypeA', 'TypeB'):
        raise ValueError('Use AXIS_CZ_DISPLACEMENT=Nodal|TypeA|TypeB')
    return dict(mode=mode, update=update, displacement=displacement,
                radius=.010, length=.008, speed=.1, density=1000., E=1.e7, nu=0.,
                stiffness=1.e8, failure_opening=5.e-6 if mode=='debond' else 1.,
                end_time=.002, tracer_labels=['lower interface', 'upper interface',
                                             'lower outer end', 'upper outer end'])


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', type=Path, required=True)
    parser.parse_args().write.write_text(json.dumps(settings(), indent=2)+'\n')
