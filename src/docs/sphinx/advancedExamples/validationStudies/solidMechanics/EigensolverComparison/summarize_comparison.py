"""Summarize the logs of run_comparison.sh into EigensolverComparison.csv and ConvergenceHistory.csv.

usage: python3 summarize_comparison.py LOGDIR [OUTPUTDIR]
"""
import os
import re
import sys

logdir = sys.argv[1]
outdir = sys.argv[2] if len(sys.argv) > 2 else '.'

RUNS = [
    ('block_arnoldi', 'Block', 16, 'Arnoldi', 'no'),
    ('block_arnoldi_deflated', 'Block', 16, 'Arnoldi', 'yes'),
    ('block_lobpcg', 'Block', 16, 'LOBPCG', 'no'),
    ('block_lobpcg_deflated', 'Block', 16, 'LOBPCG', 'yes'),
    ('beam10_arnoldi', 'Beam', 10, 'Arnoldi', 'no'),
    ('beam10_arnoldi_deflated', 'Beam', 10, 'Arnoldi', 'yes'),
    ('beam10_lobpcg', 'Beam', 10, 'LOBPCG', 'no'),
    ('beam10_lobpcg_deflated', 'Beam', 10, 'LOBPCG', 'yes'),
    ('beam20_arnoldi', 'Beam', 20, 'Arnoldi', 'no'),
    ('beam20_lobpcg_deflated', 'Beam', 20, 'LOBPCG', 'yes'),
]
summary = re.compile(r'(\d+) restarts/iterations, (\d+) operator applications, (\d+) linear solves '
                     r'\((\d+) linear iterations, ([\d.]+) s in linear solves\), ([\d.]+) s eigensolve')
row = re.compile(r'^\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)')
history = re.compile(r'LOBPCG: iteration\s+(\d+), converged\s+(\d+)/(\d+), max error estimate (\S+), operator applications (\d+)')

hist = {}
with open(os.path.join(outdir, 'EigensolverComparison.csv'), 'w') as out:
    out.write('Problem,Modes,Eigensolver,Rigid modes deflated,Iterations,Operator applications,Linear iterations,Time (s),Converged,Mode 7 (Hz)\n')
    for name, problem, modes, solver, deflated in RUNS:
        text = open(os.path.join(logdir, name + '.log')).read()
        s = summary.search(text)
        converged = 'converged to the tolerance' not in text
        f7 = None
        for line in text.splitlines():
            m = row.match(line)
            if m and int(m.group(1)) == 7 and ('e+' in m.group(2) or 'e-' in m.group(2)):
                f7 = float(m.group(2))
        linear = s.group(4) if int(s.group(3)) > 0 else '-'
        out.write(f'{problem},{modes},{solver},{deflated},{s.group(1)},{s.group(2)},{linear},{float(s.group(6)):.1f},'
                  f'{"yes" if converged else "no"},{f7:.4f}\n')
        if solver == 'LOBPCG':
            hist[name] = [(int(m.group(5)), float(m.group(4))) for m in history.finditer(text)]

with open(os.path.join(outdir, 'ConvergenceHistory.csv'), 'w') as out:
    # Maximum error estimate of the wanted modes of the LOBPCG runs, at each iteration
    out.write('Run,Operator applications,Maximum error estimate\n')
    for name, values in hist.items():
        for apps, err in values:
            out.write(f'{name},{apps},{err:.6e}\n')
