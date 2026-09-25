#!/usr/bin/env python3
"""Compile/test the production CPU helper and compare it with an independent SVD.

This is an algebra test, not a GEOS, MPI transport, or particle mechanics run.
Requires a C++14 compiler, NumPy, and SciPy. No external packages are downloaded.
"""
import argparse
import json
from pathlib import Path
import subprocess
import numpy as np
import scipy.linalg


def solve(executable, mass, trial, boundary, ties, rows):
    n = len(mass)
    text = [f"{n} {len(rows)} {len(boundary)} {len(ties)}",
            " ".join(map(repr, mass.tolist())), " ".join(map(repr, trial.tolist()))]
    text += [f"{a} {b} {s:.17g} {d:.17g}" for a, b, s, d in boundary]
    text += [f"{a} {b}" for a, b in ties]
    for coefficients in rows:
        text.append(f"1 {len(coefficients)}")
        text += [f"{a} {b} {c:.17g}" for a, b, c in coefficients]
    proc = subprocess.run([str(executable), "--solve"], input="\n".join(text)+"\n",
                          text=True, capture_output=True, check=True)
    lines = proc.stdout.splitlines()
    return np.fromstring(lines[0], sep=" "), np.fromstring(lines[2], sep=" ")


def reference(mass, trial, boundary, ties, rows):
    n = len(mass)
    matrix, rhs = [], []
    for a, b, sign, offset in boundary:
        r = np.zeros(n); r[a] = 1
        if b >= 0: r[b] -= sign
        matrix.append(r); rhs.append(offset)
    for a, b in ties:
        r = np.zeros(n); r[a] = 1; r[b] -= 1
        matrix.append(r); rhs.append(0)
    for coefficients in rows:
        r = np.zeros(n)
        for a,b,c in coefficients: r[a] += c; r[b] -= c
        matrix.append(r); rhs.append(0)
    a = np.array(matrix)
    target = np.array(rhs)
    b = a / np.sqrt(mass)[None,:]
    scale = np.linalg.norm(b, axis=1)
    b = b / scale[:,None]
    g = (target-a@trial)/scale
    u,s,vt = scipy.linalg.svd(b, full_matrices=False)
    keep = s > s[0]*1e-12
    correction = vt[keep].T @ ((u[:,keep].T @ g)/s[keep])
    velocity = trial + correction / np.sqrt(mass)
    return velocity, np.max(np.abs(a@velocity-target))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--compiler', default='g++')
    args=parser.parse_args(); args.output.mkdir(parents=True, exist_ok=True)
    here=Path(__file__).resolve().parent
    executable=(args.output/'test_boundary_mortar').resolve()
    subprocess.run([args.compiler,'-std=c++14','-O2','-Wall','-Wextra','-Wpedantic','-Werror',
                    str(here/'testMPMBoundaryMortar.cpp'),'-o',str(executable)],check=True)
    basic=subprocess.check_output([str(executable)],text=True).strip()
    rng=np.random.default_rng(260923)
    max_difference=max_residual=0.0
    for k in range(120):
        nodes=6+k%10; n=nodes*2
        mass=10**rng.uniform(-2,2,n)
        trial=rng.normal(0,.1,n)
        wall=float(rng.normal(0,.2))
        # Field-wise mirror dependencies and a common prescribed boundary value.
        boundary=[(0,2,-1,2*wall),(1,3,-1,2*wall),
                  (4,-1,1,wall),(5,-1,1,wall)]
        ties=[(6,7)] if k%3==0 else []
        rows=[]
        for j in range(2+k%4):
            c=rng.uniform(.01,1,nodes);c/=c.sum()
            rows.append([(2*i,2*i+1,float(c[i])) for i in range(nodes)])
        if k%2==0: rows.append([(a,b,3.7*c) for a,b,c in rows[0]])
        actual, diagnostics=solve(executable,mass,trial,boundary,ties,rows)
        expected,residual=reference(mass,trial,boundary,ties,rows)
        delta=float(np.max(np.abs(actual-expected)))
        np.testing.assert_allclose(actual,expected,atol=3e-10,rtol=3e-10)
        max_difference=max(max_difference,delta);max_residual=max(max_residual,float(diagnostics[0]))
        # Changing the relation order cannot choose a different corner solution.
        reordered,_=solve(executable,mass,trial,list(reversed(boundary)),ties,rows)
        np.testing.assert_allclose(actual,reordered,atol=3e-10,rtol=3e-10)
    report=dict(basic=basic,randomized_cases=120,additional_reordered_cases=120,
                maximum_difference_from_svd=max_difference,
                maximum_final_average_mortar_jump=max_residual,
                scope='Production C++ helper only; not a full GEOS/MPI/physics test')
    (args.output/'reference_results.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__': main()
