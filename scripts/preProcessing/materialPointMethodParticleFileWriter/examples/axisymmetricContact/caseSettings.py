"""Small axial-impact case; all lengths are metres and stresses are pascals."""
import argparse
import json
import os
from pathlib import Path


def settings():
    update = os.environ.get('AXIS_CONTACT_UPDATE', 'FLIP')
    if update not in ('FLIP', 'FMPM'):
        raise ValueError('AXIS_CONTACT_UPDATE must be FLIP or FMPM')
    return dict(update=update, radius=0.010, length=0.008, gap=0.002,
                speed=1.0, density=1000.0, E=1.e6, nu=0.0,
                end_time=0.003, contact_time=0.001,
                tracer_labels=['lower interface', 'upper interface',
                               'lower outer end', 'upper outer end'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', type=Path, required=True)
    parser.parse_args().write.write_text(json.dumps(settings(), indent=2)+'\n')
