"""Axisymmetric contact/DFG isolation cases; SI units."""
import argparse
import json
import os


def settings():
    mode = os.environ.get('AXIS_DFG_CASE', 'damage')
    update = os.environ.get('AXIS_DFG_UPDATE', 'FLIP')
    if mode not in ('intact', 'damage', 'contact', 'combined'):
        raise ValueError('AXIS_DFG_CASE must be intact, damage, contact or combined')
    if update not in ('FLIP', 'FMPM'):
        raise ValueError('AXIS_DFG_UPDATE must be FLIP or FMPM')
    return dict(mode=mode, update=update, radius=.006, density=1000., E=1.e7,
                speed=.1, end_time=.004 if mode in ('contact', 'combined') else .001,
                dfg=int(mode != 'contact'), groups=3 if mode in ('contact', 'combined') else 1,
                tracer_labels=['below centre', 'above centre', 'lower end', 'upper end'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', required=True)
    with open(parser.parse_args().write, 'w') as stream:
        json.dump(settings(), stream, indent=2)
        stream.write('\n')
