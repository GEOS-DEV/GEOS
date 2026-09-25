#!/usr/bin/env python3
"""Standalone checks of production mortar performance kernels (not a GEOS run).

Requires a C++14 compiler and NumPy. --lapack also requires LP64 liblapack/libblas.
The scheduler test explicitly emulates collectives with threads. --mpi builds
and runs the separate real-MPI test when mpicxx/mpiexec are available; it does not
substitute the emulator for MPI transport coverage.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import platform
import shutil
import struct
import subprocess
import sys
import time


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--cxx', default='g++')
    parser.add_argument('--lapack', action='store_true')
    parser.add_argument('--mpi', action='store_true')
    parser.add_argument('--sanitize', action='store_true')
    parser.add_argument('--benchmark', action='store_true')
    args = parser.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    here = Path(__file__).resolve().parent
    sm = here.parent.parent
    kernels = sm / 'kernels'
    env = dict(os.environ)
    for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
        env[name] = '1'
        os.environ[name] = '1'
    import numpy as np

    flags = ['-std=c++14', '-O3', '-Wall', '-Wextra', '-Werror', '-pthread']
    if args.sanitize:
        flags = ['-std=c++14', '-O1', '-g', '-fsanitize=address,undefined',
                 '-fno-omit-frame-pointer', '-fno-pie', '-no-pie', '-pthread']
    libs = ['-llapack', '-lblas'] if args.lapack else []
    if args.lapack:
        flags += ['-DGEOS_MPM_MORTAR_LAPACK=1']
    backend_names = ['reference', 'native'] + (['lapack'] if args.lapack else [])
    qr_source = kernels / 'MPMBoundaryMortarQR.cpp'
    results = {'scope': 'standalone kernels; no GEOS timestep or mechanics',
               'compiler': subprocess.check_output([args.cxx, '--version'], text=True).splitlines()[0],
               'platform': platform.platform(), 'lapack': args.lapack,
               'sanitizers': args.sanitize, 'tests': {}}
    command_log = []

    def run(command, **kwargs):
        command_log.append(' '.join(map(str, command)))
        return subprocess.run(list(map(str, command)), env=env, check=True, text=True,
                              capture_output=True, timeout=kwargs.pop('timeout', 120), **kwargs)

    performance = out / 'testPerformance'
    run([args.cxx, *flags, here / 'testPerformance.cpp', qr_source, *libs, '-o', performance])
    for backend in backend_names:
        log = run([performance, backend])
        results['tests'][f'search_scheduler_{backend}'] = log.stdout.strip()
        (out / f'scheduler_{backend}.log').write_text(log.stdout + log.stderr)

    # Re-run the existing production-header checks through each new backend.
    # The checked-in test remains unchanged. This generated wrapper changes only
    # which minimum-norm function project() calls, not the test inputs/assertions.
    source = (here.parent / 'testMPMBoundaryMortar.cpp').read_text()
    source = source.replace('#include "../kernels/MPMBoundaryMortar.hpp"',
                            f'#include "{kernels / "MPMBoundaryMortarQR.hpp"}"\n#include <cstdlib>')
    source = source.replace('bm::project(', 'projectBackend(')
    wrapper = '''
bm::Result projectBackend(std::vector<double> const & mass, std::vector<double> const & trial,
  std::vector<double> const & cutoff, std::vector<bm::Relation> const & boundary,
  std::vector<std::pair<bm::Index,bm::Index>> const & fallback, std::vector<bm::Row> const & rows,
  bm::Options const & options)
{
  char const * name=std::getenv("TEST_QR_BACKEND");
  return bm::project(mass,trial,cutoff,boundary,fallback,rows,options,
    bm::minimumNormSolver(bm::selectQRBackend(name?name:"reference")));
}
'''
    source = source.replace('int checks = 0;', wrapper + '\nint checks = 0;')
    generated = out / 'existing_checks_with_backend.cpp'
    generated.write_text(source)
    existing = out / 'existing_checks'
    run([args.cxx, *flags, generated, qr_source, *libs, '-o', existing])
    for backend in backend_names:
        env['TEST_QR_BACKEND'] = backend
        log = run([existing])
        results['tests'][f'existing_{backend}'] = log.stdout.strip()
        (out / f'existing_{backend}.log').write_text(log.stdout + log.stderr)

    rng = np.random.default_rng(72631)
    fixture = out / 'qr_input.bin'
    answer = out / 'qr_output.bin'

    def qr(B, b, backend, tolerance=1.e-12):
        m, n = B.shape
        with fixture.open('wb') as f:
            f.write(struct.pack('QQd', n, m, tolerance))
            f.write(np.asarray(B, dtype=np.float64).tobytes(order='C'))
            f.write(np.asarray(b, dtype=np.float64).tobytes())
        run([performance, '--qr', backend, fixture, answer])
        raw = answer.read_bytes()
        rank, seconds = struct.unpack('Qd', raw[:16])
        values = np.frombuffer(raw[16:], dtype=np.float64)
        return rank, seconds, values[:n].copy(), values[n:].copy()

    maxima = {name: {'relative_solution_error': 0., 'relative_equation_residual': 0.,
                     'impulse_identity_error': 0.} for name in backend_names}
    rank_cases = []
    for case in range(90):
        m, n = (int(rng.integers(4, 35)), int(rng.integers(5, 50)))
        rank = min(m, n) if case % 3 == 0 else max(1, min(m, n)-2)
        U = np.linalg.qr(rng.normal(size=(m, rank)))[0]
        V = np.linalg.qr(rng.normal(size=(n, rank)))[0]
        cond = (1., 1.e3, 1.e6)[case % 3]
        singular = np.geomspace(1., 1./cond, rank)
        B = (U * singular) @ V.T
        B /= np.linalg.norm(B, axis=1)[:, None]
        b = B @ rng.normal(size=n)
        u, s, vt = np.linalg.svd(B, full_matrices=False)
        expected = vt[:rank].T @ ((u[:, :rank].T @ b)/s[:rank])
        for backend in backend_names:
            found_rank, _, x, dual = qr(B, b, backend)
            if found_rank != rank:
                raise AssertionError((case, backend, 'rank', found_rank, rank))
            error = np.linalg.norm(x-expected)/max(1., np.linalg.norm(expected))
            residual = np.linalg.norm(B@x-b)/max(1., np.linalg.norm(b))
            # Dual recovery can amplify roundoff in nearly dependent rows. Track
            # the physical identity as well as the much better conditioned primal.
            impulse = np.linalg.norm(B.T@dual-x)/max(1., np.linalg.norm(x))
            if error > 2.e-7 or residual > 2.e-11 or impulse > 2.e-7:
                raise AssertionError((case, backend, error, residual, impulse))
            record = maxima[backend]
            record['relative_solution_error'] = max(record['relative_solution_error'], float(error))
            record['relative_equation_residual'] = max(record['relative_equation_residual'], float(residual))
            record['impulse_identity_error'] = max(record['impulse_identity_error'], float(impulse))
        rank_cases.append({'m': m, 'n': n, 'rank': rank, 'target_condition': cond})
    results['tests']['svd_cases'] = {'count': len(rank_cases), 'maxima': maxima, 'cases': rank_cases}

    if args.benchmark:
        records = []
        n, m = 826, 337
        B = rng.normal(size=(m, n))
        B /= np.linalg.norm(B, axis=1)[:, None]
        correlated = B.copy()
        correlated += 0.999*rng.normal(size=(1, n))
        correlated /= np.linalg.norm(correlated, axis=1)[:, None]
        for title, matrix in [('random', B), ('correlated', correlated)]:
            rhs = matrix@rng.normal(size=n)
            for backend in backend_names:
                times = []
                for repetition in range(5):
                    rank, seconds, x, _ = qr(matrix, rhs, backend)
                    if repetition:
                        times.append(seconds)
                records.append({'matrix': title, 'backend': backend, 'n': n, 'm': m,
                                'rank': rank, 'median_seconds': float(np.median(times)),
                                'min_seconds': min(times),
                                'relative_residual': float(np.linalg.norm(matrix@x-rhs)/np.linalg.norm(rhs))})
        results['benchmark'] = {'scope': 'synthetic same-size QR kernels only; NOT GEOS step timings',
                                'threads': 1, 'records': records}

    if args.mpi:
        compiler, launcher = shutil.which('mpicxx'), shutil.which('mpiexec')
        if compiler and launcher:
            binary = out / 'testMPI'
            run([compiler, *flags, here / 'testMPI.cpp', qr_source, *libs, '-o', binary])
            results['mpi'] = []
            for ranks in [1, 2, 3, 4, 8]:
                result = run([launcher, '-n', ranks, binary, 'lapack' if args.lapack else 'native'])
                results['mpi'].append({'ranks': ranks, 'output': result.stdout.strip()})
        else:
            results['mpi'] = 'NOT RUN: mpicxx/mpiexec unavailable; emulated tests do not cover MPI transport'
    else:
        results['mpi'] = 'NOT REQUESTED'
    (out / 'commands.log').write_text('\n'.join(command_log) + '\n')
    (out / 'results.json').write_text(json.dumps(results, indent=2) + '\n')
    print(json.dumps({k: v for k, v in results.items() if k != 'tests'}, indent=2))
    print('PASS: existing boundary tests, ordered-search/scheduler checks, and 90 SVD fixtures for each backend')


if __name__ == '__main__':
    main()
