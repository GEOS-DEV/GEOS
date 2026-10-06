import argparse
import csv

import matplotlib.pyplot as plt
import numpy as np


def read_history(path):
    """Read the blocks of ConvergenceHistory.txt: {name: (applications, error estimates)}."""
    histories, name, data = {}, None, []
    for line in open(path):
        line = line.strip()
        if line.startswith('## '):
            if name is not None:
                histories[name] = np.array(data)
            name, data = line[3:], []
        elif line and not line.startswith('#'):
            data.append([float(v) for v in line.split()])
    if name is not None:
        histories[name] = np.array(data)
    return histories


def main():

    parser = argparse.ArgumentParser(description="Script to generate the figure of the eigensolver comparison.")
    parser.add_argument('--outputDir', help='Path to the directory of the results', default='.')
    args = parser.parse_args()

    with open(args.outputDir + '/EigensolverComparison.csv') as f:
        rows = list(csv.DictReader(f))
    histories = read_history(args.outputDir + '/ConvergenceHistory.txt')

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))

    # Convergence of LOBPCG with and without deflation of the rigid-body modes
    styles = {'block_lobpcg': ('C0', '--', 'Block, no deflation'),
              'block_lobpcg_deflated': ('C0', '-', 'Block, deflation'),
              'beam10_lobpcg': ('C3', '--', 'Beam, no deflation'),
              'beam10_lobpcg_deflated': ('C3', '-', 'Beam, deflation')}
    for name, (color, style, label) in styles.items():
        h = histories[name]
        ax1.semilogy(h[:, 0], h[:, 1], style, color=color, label=label)
    ax1.axhline(1e-8, color='k', lw=0.8, ls=':')
    ax1.text(ax1.get_xlim()[1] * 0.98, 1.5e-8, 'tolerance', ha='right', fontsize=8)
    ax1.set_xlabel('Preconditioner applications')
    ax1.set_ylabel('Maximum error estimate')
    ax1.set_title('LOBPCG convergence', fontsize=10)
    ax1.grid(True, which='both', alpha=0.3)
    ax1.legend(fontsize=8)

    # Time of the four combinations for the two problems
    labels = ['Arnoldi', 'Arnoldi\n+ deflation', 'LOBPCG', 'LOBPCG\n+ deflation']
    width = 0.38
    positions = np.arange(4)
    for k, (problem, modes, color) in enumerate((('Block', '16', 'C0'), ('Beam', '10', 'C3'))):
        selected = [r for r in rows if r['Problem'] == problem and r['Modes'] == modes]
        times = [float(r['Time (s)']) for r in selected]
        converged = [r['Converged'] == 'yes' for r in selected]
        bars = ax2.bar(positions + (k - 0.5) * width, times, width, color=color, label=f'{problem}, {modes} modes')
        for bar, ok in zip(bars, converged):
            if not ok:
                bar.set_hatch('//')
                bar.set_alpha(0.5)
                ax2.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), 'not converged', rotation=90,
                         ha='center', va='bottom', fontsize=7)
    ax2.set_xticks(positions)
    ax2.set_xticklabels(labels, fontsize=8)
    ax2.set_ylabel('Eigensolve time on one core (s)')
    ax2.set_title('Cost of the four combinations', fontsize=10)
    ax2.grid(True, axis='y', alpha=0.3)
    ax2.legend(fontsize=8)

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
