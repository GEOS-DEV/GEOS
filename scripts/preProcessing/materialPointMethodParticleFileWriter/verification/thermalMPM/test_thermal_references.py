#!/usr/bin/env python3
"""Self-tests for analytical references and CSV post-processors (no GEOS required)."""
from __future__ import annotations
import csv
import json
import math
from pathlib import Path
import tempfile
import unittest

from thermal_cases import (CASES,L,Y,Z,T0,AMP,CP,RHO,Q_SPECIFIC,temperature,
                           contact_lumped_temperatures,three_field_lumped_temperatures,
                           contact_resistance,THERMAL_EXPANSION)
from thermal_post import run_case


def write_csv(path,rows):
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with path.open('w',newline='') as h:
        w=csv.DictWriter(h,fieldnames=keys);w.writeheader();w.writerows(rows)


def fixture(run_dir,case_id):
    case=CASES[case_id];kind=case['kind'];stop=case['stop']
    times=[stop*i/4 for i in range(5)]
    if kind=='contact':
        positions=[L/4,3*L/4]
    elif kind=='three_field':
        positions=[L/4,L/4,3*L/4]
    elif kind=='thermoelastic':
        positions=[L/5,4*L/5]
    elif kind=='plastic':
        positions=[L/3,2*L/3]
    else:
        positions=[L*(i+0.5)/32 for i in range(32)]
    for index,x in enumerate(positions):
        rows=[]
        for t in times:
            if kind=='contact':tmp=contact_lumped_temperatures(t,contact_resistance(case_id))[index]
            elif kind=='three_field':tmp=three_field_lumped_temperatures(t)[index]
            elif kind=='thermoelastic':tmp=T0+AMP
            elif kind=='plastic':tmp=T0+1000*0.01*(t/stop)/(RHO*CP)
            else:tmp=temperature(kind,x,t,case=case)
            xpos=x
            if kind=='thermoelastic':xpos=L/2+(x-L/2)*(1+THERMAL_EXPANSION*AMP*t/stop)
            epsp=0.01*t/stop if kind=='plastic' else 0.
            d=dict(t=t,x=xpos,y=0,z=0,temperature=tmp,internalEnergy=CP*(tmp-T0),
                   mass=RHO*L*Y*Z/len(positions),plasticStrainMagnitude=epsp)
            for suffix in ('XX','YY','ZZ','YZ','XZ','XY'):
                d['stress'+suffix]=1000. if (kind=='plastic' and suffix=='XX') else 0.
                d['plasticStrain'+suffix]=epsp if suffix=='XX' else 0.
            rows.append(d)
        write_csv(run_dir/('thermalVV_{:03d}.csv'.format(index)),rows)
    if kind=='contact':
        mean=T0+AMP/2
    elif kind=='three_field':
        mean=T0+AMP*0.375
    else:
        mean=T0
    rows=[]
    for t in times:
        effective_mean=mean
        if kind=='source':effective_mean=T0+Q_SPECIFIC*t/CP
        elif kind=='sink':effective_mean=T0-Q_SPECIFIC*t/CP
        # Mechanical internal energy may remain zero during pure conduction.
        rows.append(dict(Time=t,Density=RHO,Temperature=effective_mean,
                         **{'Internal Energy':0.},volume=L*Y*Z))
    write_csv(run_dir/'boxAverageHistory.csv',rows)


class TestReferenceModels(unittest.TestCase):
    def test_contact_total_energy(self):
        for resistance in (0.005,0.025,0.1,math.inf):
            for t in (0.,.01,.1,1.):
                a,b=contact_lumped_temperatures(t,resistance)
                self.assertAlmostEqual(a+b,2*T0+AMP,places=10)

    def test_three_field_total_energy(self):
        for t in (0.,.01,.1,1.):
            a,b,c=three_field_lumped_temperatures(t)
            self.assertAlmostEqual((a+b+2*c)/4,T0+0.375*AMP,places=10)

    def test_cosine_amplitude(self):
        a=CASES['cosine'];x=L/4
        self.assertAlmostEqual(temperature('cosine',x,1,case=a),T0,places=10)
        self.assertAlmostEqual(temperature('cosine',0,0,case=a),T0+AMP/2,places=10)

    def test_all_synthetic_cases_pass(self):
        with tempfile.TemporaryDirectory() as td:
            base=Path(td)
            for case in CASES:
                with self.subTest(case=case):
                    run=base/case;run.mkdir()
                    fixture(run,case)
                    result_dir=base/('results_'+case)
                    rc=run_case(case,['--run-dir',str(run),'--output-dir',str(result_dir),'--no-plot'])
                    data=json.loads((result_dir/(case+'_summary.json')).read_text())
                    self.assertEqual(rc,0,data)
                    self.assertEqual(data['status'],'PASS',data)

    def test_missing_output_fails(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)
            result=p/'out'
            rc=run_case('uniform',['--run-dir',td,'--output-dir',str(result),'--no-plot'])
            self.assertNotEqual(rc,0)
            data=json.loads((result/'uniform_summary.json').read_text())
            self.assertEqual(data['status'],'FAIL')


if __name__=='__main__':unittest.main()
