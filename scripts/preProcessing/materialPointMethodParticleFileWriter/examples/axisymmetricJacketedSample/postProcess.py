#!/usr/bin/env python3
"""Plot GEOS histories. Homogeneous mode also checks a small-strain reference.

The bonded heterogeneous assembly has no simple uniform-stress exact solution.
Its plots are diagnostics, not a pass/fail material-model validation.
"""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from caseSettings import pressure_at


def read_csv(path):
    with path.open() as stream:
        reader = csv.reader(stream)
        names = [s.strip() for s in next(reader)]
        values = np.array([[float(v) for v in row] for row in reader if row], dtype=float)
    if values.ndim != 2 or len(values) < 5 or not np.all(np.isfinite(values)):
        raise ValueError('Missing, short or nonfinite history: '+str(path))
    return {name: values[:, i] for i, name in enumerate(names)}


def process(folder):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    case = json.loads((folder/'case.json').read_text())
    box = read_csv(folder/'boxAverageHistory.csv')
    traces = [read_csv(folder/f'jacketTracer_{i:06d}.csv') for i in range(5)]
    for time in [box['Time']]+[row['t'] for row in traces]:
        if time[-1] < 0.99*case['end_time'] or np.any(np.diff(time) <= 0):
            raise ValueError('History does not cover the run or contains repeated/out-of-order times')
    pressure = np.array([pressure_at(t, case) for t in box['Time']])
    sample = traces[0]
    volume_strain = case['materials'][0]['density']/sample['density']-1.
    fig, axs = plt.subplots(2, 2, figsize=(11, 7), layout='constrained')
    ax = axs[0, 0]
    ax.plot(box['Time']*1.e3, pressure*1.e-3, 'k--', label='bath pressure')
    for name, label in [('Sxx','radial'), ('Syy','axial'), ('Szz','hoop')]:
        ax.plot(box['Time']*1.e3, -box[name]*1.e-3, label=label)
    ax.set(ylabel='compression positive [kPa]', title='Sample-interior average stress')
    ax.legend(fontsize=8)
    ax = axs[0, 1]
    for row, label in zip(traces, case['tracer_labels']):
        ax.plot(row['t']*1.e3, -row['stressYY']*1.e-3, label=label)
    ax.set(ylabel='axial compression [kPa]', title='Lagrangian tracer stresses')
    ax.legend(fontsize=8)
    ax = axs[1, 0]
    ax.plot(sample['t']*1.e3, volume_strain*100., label='sample tracer')
    if case['variant'] == 'homogeneous':
        material = case['materials'][0]
        bulk = material['E']/(3.*(1.-2.*material['nu']))
        ax.plot(sample['t']*1.e3,
                [-100.*pressure_at(t, case)/bulk for t in sample['t']], 'k--', label='-p/K (quasistatic)')
    ax.set(ylabel='volume change [%]', title='Volume strain from conserved mass / density')
    ax.legend(fontsize=8)
    ax = axs[1, 1]
    for i in (0, 3):
        row = traces[i]
        ax.plot(row['t']*1.e3, (row['x']-row['x'][0])*1.e6, label=case['tracer_labels'][i])
    ax.set(ylabel='radial displacement from first record [um]', title='Radial contraction')
    ax.legend(fontsize=8)
    for ax in axs.flat:
        ax.set_xlabel('time [ms]'); ax.grid(alpha=0.2)
    fig.suptitle('Axisymmetric jacketed specimen: '+case['variant'])
    fig.savefig(folder/'histories.png', dpi=160); plt.close(fig)
    report = dict(variant=case['variant'], final_time=float(box['Time'][-1]),
                  final_sample_stress_Pa={name: float(box[name][-1]) for name in ('Sxx','Syy','Szz')},
                  final_sample_volume_strain=float(volume_strain[-1]))
    if case['variant'] == 'homogeneous':
        # Compare the final 20% of the hold. These smoke-test tolerances do not
        # replace a resolution/ramp-time convergence study.
        begin = case['end_time']-0.2*(case['end_time']-case['ramp_time'])
        tail = box['Time'] >= begin; tracer_tail = sample['t'] >= begin
        if tail.sum() < 5 or tracer_tail.sum() < 5:
            raise ValueError('Too few final-hold history samples')
        p = case['pressure']
        errors = {name: float(abs(np.mean(box[name][tail])+p)/p) for name in ('Sxx','Syy','Szz')}
        oscillation = max(float(np.ptp(box[name][tail])/p) for name in errors)
        volumetric_error = float(abs(np.mean(volume_strain[tracer_tail])+p/bulk)/(p/bulk))
        passed = max(errors.values()) < 0.05 and oscillation < 0.10 and volumetric_error < 0.20
        report.update(reference='sigma_rr=sigma_yy=sigma_tt=-p; dV/V=-p/K',
                      relative_stress_errors=errors, pressure_normalized_stress_range=oscillation,
                      relative_volume_strain_error=volumetric_error, homogeneous_check_passed=bool(passed))
    (folder/'historyChecks.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))
    if report.get('homogeneous_check_passed') is False:
        raise SystemExit('Homogeneous check failed: inspect histories.png; test longer settling and finer resolution.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', nargs='?', type=Path, default=Path('.'))
    process(parser.parse_args().directory)
