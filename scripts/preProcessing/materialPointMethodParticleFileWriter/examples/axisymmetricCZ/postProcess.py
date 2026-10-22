#!/usr/bin/env python3
"""Plot contact histories; diagnostics, not an analytical acceptance test."""
import csv
import json
from pathlib import Path
import numpy as np


def process(folder):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    case = json.loads((folder/'case.json').read_text())
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), layout='constrained')
    summaries = []
    for i, label in enumerate(case['tracer_labels']):
        with (folder/f'czTracer_{i:06d}.csv').open() as stream:
            reader = csv.reader(stream)
            names = [name.strip() for name in next(reader)]
            values = np.array([[float(x) for x in row] for row in reader if row])
        if values.ndim != 2 or len(values) < 5 or not np.all(np.isfinite(values)):
            raise ValueError('Missing, short or nonfinite tracer history')
        row = {name: values[:, j] for j, name in enumerate(names)}
        if row['t'][-1] < .99*case['end_time'] or np.any(np.diff(row['t']) <= 0):
            raise ValueError('Tracer history does not cover the run or has nonmonotonic times')
        if np.max(np.abs(row['velocityZ'])) > 1.e-12*case['speed']:
            raise ValueError('Unexpected circumferential velocity')
        for ax, data in zip(axes.flat, [row['velocityY'], row['stressYY']/1.e3,
                                        row['y']*1.e3, row['velocityX']]):
            ax.plot(row['t']*1.e3, data, label=label)
        summaries.append(dict(label=label, final_axial_velocity=float(row['velocityY'][-1]),
                              peak_axial_tension_Pa=float(row['stressYY'].max())))
    for ax, title, unit in zip(axes.flat,
                              ['Axial velocity', 'Axial tension', 'Tracer position', 'Radial velocity'],
                              ['m/s', 'kPa', 'mm', 'm/s']):
        ax.set(xlabel='time [ms]', ylabel=unit, title=title)
        ax.grid(alpha=.2)
    axes[0, 0].legend(fontsize=7)
    fig.suptitle('Axisymmetric CZ: '+case['mode']+' / '+case['update'])
    fig.savefig(folder/'czHistories.png', dpi=160)
    plt.close(fig)
    (folder/'czDiagnostics.json').write_text(json.dumps(summaries, indent=2)+'\n')
    print('Wrote czHistories.png and czDiagnostics.json. Inspect CZ force, separation and damage in the grid/particle outputs; no dynamics accuracy threshold is asserted.')


if __name__ == '__main__':
    process(Path('.'))
