#!/usr/bin/env python3
"""GEOS MPM thermal V&V postprocessor. Python standard library required.

Reads GEOS tracer CSVs and boxAverageHistory.csv. The latter's `internalEnergy`
records mechanical stress-work rather than transported caloric thermal energy.
Conservation checks therefore use mass-weighted mean temperature at constant cp.
Produces detailed residuals,
quantitative acceptance metrics and optional plots. Never reports PASS with
missing output or missing physical observables.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path
import re
import sys

from thermal_cases import (CASES, L, Y, Z, T0, AMP, CP, RHO, ALPHA, K,
                           Q_SPECIFIC, THERMAL_EXPANSION, temperature,
                           contact_resistance, contact_lumped_temperatures,
                           three_field_lumped_temperatures)

TOL = {
    'uniform': 0.001,
    'cosine': 0.10,
    'gaussian': 0.13,
    'step': 0.17,
    'source_local': 0.035,
    'contact_perfect': 0.17,
    'source': 0.035,
    'capacity_linear':0.035,
    'conductivity_linear':0.10,
    'sink': 0.035,
    'dirichlet_steady': 0.10,
    'dirichlet_transient': 0.20,
    'contact': 0.15,
    'three_field': 0.20,
    'thermoelastic': 0.15,
    'plastic': 0.25,
}


def _number(v):
    try:
        f = float(v)
        return f if math.isfinite(f) else None
    except (ValueError, TypeError, OverflowError):
        return None


def _key(row, *names):
    have = {re.sub('[^a-z0-9]', '', k.lower()): k for k in row.keys()}
    for name in names:
        clean = re.sub('[^a-z0-9]', '', name.lower())
        if clean in have:
            return have[clean]
    return None


def _read_csv(p):
    with p.open(newline='') as h:
        dr = csv.DictReader(h)
        if not dr.fieldnames:
            return []
        return [{k.strip(): v for k, v in r.items() if k} for r in dr]


def _locate_run(case_id, source_dir, run_dir):
    if run_dir:
        p = Path(run_dir).resolve()
        if not p.is_dir():
            raise FileNotFoundError('Run directory does not exist: ' + str(p))
        return p
    # Avoid guessing one run among multiple results. Raise on ambiguity.
    prefix = 'thermalVV_' + case_id
    candidates = []
    for pattern in ('**/'+prefix+'/boxAverageHistory.csv', '**/'+prefix+'*/boxAverageHistory.csv',
                    '**/boxAverageHistory.csv'):
        candidates.extend(source_dir.glob(pattern))
        if candidates:
            break
    parents = sorted(set(f.parent for f in candidates))
    if len(parents) != 1:
        raise FileNotFoundError('Expected exactly one run directory for {} in {}; found {}. '
                                'Supply --run-dir explicitly.'.format(case_id, source_dir, len(parents)))
    return parents[0]


def read_tracers(run_dir):
    files = sorted(run_dir.glob('thermalVV_*.csv'))
    if not files:
        files = sorted(run_dir.glob('**/thermalVV_*.csv'))
    if not files:
        raise FileNotFoundError('No thermalVV_NNN.csv tracer files found in '+str(run_dir))
    traces = []
    for file in files:
        rows = _read_csv(file)
        if not rows:
            continue
        out = []
        for row in rows:
            tcol = _key(row,'t','time')
            xcol = _key(row,'x')
            vcol = _key(row,'temperature','particleTemperature')
            if not all((tcol,xcol,vcol)):
                raise ValueError('Missing time/x/temperature in tracer '+str(file))
            vals = dict(t=_number(row[tcol]),x=_number(row[xcol]),
                        y=_number(row.get(_key(row,'y'))),
                        z=_number(row.get(_key(row,'z'))),
                        temperature=_number(row[vcol]),
                        energy=_number(row.get(_key(row,'internalEnergy'))),
                        mass=_number(row.get(_key(row,'mass'))),
                        plastic_strain=_number(row.get(_key(row,'plasticStrainMagnitude'))),
                        sxx=_number(row.get(_key(row,'stressXX'))),
                        syy=_number(row.get(_key(row,'stressYY'))))
            for suffix in ('XX','YY','ZZ','YZ','XZ','XY'):
                vals['stress'+suffix] = _number(row.get(_key(row,'stress'+suffix)))
                vals['plasticStrain'+suffix] = _number(row.get(_key(row,'plasticStrain'+suffix)))
            if any(vals[k] is None for k in ('t','x','temperature')):
                raise ValueError('Invalid/NaN thermal tracer sample in '+str(file))
            out.append(vals)
        if out:
            if any(b['t'] <= a['t'] for a,b in zip(out, out[1:])):
                raise ValueError('Non-increasing output times in '+str(file))
            traces.append((file.stem, out))
    if not traces:
        raise ValueError('All tracer files were empty')
    return traces


def _box_history(run_dir):
    files = list(run_dir.glob('boxAverageHistory.csv')) or list(run_dir.glob('**/boxAverageHistory.csv'))
    if len(files) != 1:
        return []
    data = []
    for row in _read_csv(files[0]):
        t = _number(row.get(_key(row,'time')))
        u = _number(row.get(_key(row,'internalenergy')))
        rho = _number(row.get(_key(row,'density')))
        temp = _number(row.get(_key(row,'temperature')))
        vol = _number(row.get(_key(row,'volume')))
        if t is not None and temp is not None:
            # Only sensible for this suite's constant-cp, fixed-mass cases.
            # Do not mistake the existing GEOS internalEnergy for caloric energy.
            data.append(dict(t=t, mechanical_specific_energy=u, density=rho,
                             temperature=temp, volume=vol,
                             caloric_energy_J=CP*(temp-T0)*RHO*L*Y*Z))
    return data


def _score(samples, baseline, tolerance, *, relative=True):
    if not samples:
        return dict(status='FAIL', reason='no valid samples', samples=0)
    residual = [s['measured'] - s['expected'] for s in samples]
    rmse = math.sqrt(sum(r*r for r in residual)/len(residual))
    maxerr = max(abs(r) for r in residual)
    nrms = rmse / max(baseline, 1e-12) if relative else rmse
    return dict(status='PASS' if nrms <= tolerance else 'FAIL',
                samples=len(samples),rmse_K=rmse,max_abs_error_K=maxerr,
                normalized_rmse=nrms,tolerance=tolerance)


def _profile_ref(case_id, case, traces):
    kind = case['kind']
    metrics=[]
    tmin = 0.0
    # Avoid comparing discrete discontinuities as an infinite-domain exact
    # solution at t=0. Early-time solution affected by a finite boundary is
    # assessed only inside a central interval.
    for name, series in traces:
        for row in series:
            t,x = row['t'],row['x']
            if kind in ('gaussian','step','contact_perfect') and (x < .002 or x > L-.002):
                continue
            if kind in ('gaussian','step','contact_perfect') and t <= 0:
                continue
            if kind == 'source_local' and (abs(x-0.003)<0.0005 or abs(x-0.004)<0.0005):
                continue
            if kind == 'dirichlet_transient' and (x < L/20 or x > L-L/20):
                continue
            ref = temperature(kind,x,t,case=case)
            metrics.append(dict(trace=name,t=t,x=x,measured=row['temperature'],expected=ref,
                                residual=row['temperature']-ref))
    return metrics, _score(metrics, AMP if kind != 'uniform' else 1.0,
                           TOL[kind])


def _contact_ref(case_id,case,traces):
    kind=case['kind']
    metrics=[]
    if kind=='contact':
        if len(traces) < 2:
            raise ValueError('Contact analytic test needs two resolved tracers')
        tcurves = [dict((round(r['t'],11),r) for r in trace[1]) for trace in traces[:2]]
        times=sorted(set(tcurves[0])&set(tcurves[1]))
        rcontact=contact_resistance(case_id)
        for t in times:
            expected=contact_lumped_temperatures(t,rcontact)
            for i in range(2):
                measured=tcurves[i][t]['temperature']
                metrics.append(dict(trace=traces[i][0],t=t,x=tcurves[i][t]['x'],
                                    measured=measured,expected=expected[i],residual=measured-expected[i],
                                    field='AB'[i]))
        # The lumped reference assumes Bi<<1; report the approximation bound.
        bi = (L/2)/(float(case.get('conductivity',K))*rcontact) if math.isfinite(rcontact) else 0.0
    else:
        if len(traces) < 3:
            raise ValueError('Three-field test needs three field-center tracers')
        tcurves=[dict((round(r['t'],11),r) for r in trace[1]) for trace in traces[:3]]
        times=sorted(set.intersection(*(set(q) for q in tcurves)))
        for t in times:
            expected=three_field_lumped_temperatures(t)
            for i in range(3):
                measured=tcurves[i][t]['temperature']
                metrics.append(dict(trace=traces[i][0],t=t,x=tcurves[i][t]['x'],
                                    measured=measured,expected=expected[i],residual=measured-expected[i],
                                    field='ABC'[i]))
        bi = (L/2)/(K*contact_resistance('contact_global'))
    s=_score(metrics, AMP,TOL[kind])
    s.update(contact_resistance=contact_resistance(case_id), Bi=bi,
             note='Lumped analytical ODE limit; finite bulk gradients cause O(Bi) error. '
                  'Verify Bi small and separately run mesh/contact area refinement.')
    return metrics,s


def _energy_closed(case_id,case,traces,box):
    if len(box)<2:
        raise ValueError('Global boxAverageHistory.csv missing or insufficient for conservation test')
    e0=box[0]['caloric_energy_J']
    samples=[dict(t=r['t'],measured=r['caloric_energy_J'],expected=e0,
                  residual=r['caloric_energy_J']-e0) for r in box]
    scale = RHO * CP * L * Y * Z * AMP
    drift=max(abs(s['residual']) for s in samples)/max(scale,1e-12)
    return samples, dict(status='PASS' if drift<=1e-5 else 'FAIL',
                         samples=len(samples),relative_energy_drift=drift,tolerance=1e-5)


def _thermoelastic(case_id,case,traces):
    if len(traces)<2:
        raise ValueError('Free expansion needs two tracer histories to measure length change')
    lo=min(traces,key=lambda t: t[1][0]['x'])[1]
    hi=max(traces,key=lambda t: t[1][0]['x'])[1]
    if len(lo)<2 or len(hi)<2:
        raise ValueError('Need at least two times for thermal expansion check')
    measured=((hi[-1]['x']-lo[-1]['x'])/(hi[0]['x']-lo[0]['x'])-1.)
    expected=THERMAL_EXPANSION*AMP
    max_s=max(abs(r['sxx'] or 0.) for _,rr in traces for r in rr[-1:])
    metrics=[dict(t=hi[-1]['t'],measured=measured,expected=expected,residual=measured-expected)]
    rel=abs(measured-expected)/max(abs(expected),1e-12)
    return metrics,dict(status='PASS' if rel<=TOL['thermoelastic'] and max_s<=10.0 else 'FAIL',
                        strain_relative_error=rel,stress_magnitude_Pa=max_s,
                        tolerance=TOL['thermoelastic'])


def _plastic(case_id,case,traces):
    """Trapezoidal integration of measured sigma : d(epsilon_pl) for each tracer.

    For small strain and consistent stress/plastic strain conjugacy:
        Delta T = beta/(rho cp) integral sigma_ij d(epsp_ij)
    Tensor shear terms have factor two. This is a measurable energy oracle,
    not a rough end-point approximation. Plastic work must be nonzero.
    """
    beta=float(os.environ.get('MPM_THERMAL_PLASTIC_HEAT_FRACTION','1.0'))
    if not 0<=beta<=1:raise ValueError('plastic_heat_fraction must be in [0,1]')
    metrics=[]
    tensors=(('XX',1),('YY',1),('ZZ',1),('YZ',2),('XZ',2),('XY',2))
    for name,rr in traces:
        if len(rr)<2:continue
        work=0.0
        start=rr[0]['temperature']
        yielded=False
        for pre,post in zip(rr,rr[1:]):
            dwork=0.0
            for comp,factor in tensors:
                sigma_pre=pre['stress'+comp];sigma_post=post['stress'+comp]
                eps_pre=pre['plasticStrain'+comp];eps_post=post['plasticStrain'+comp]
                if any(v is None for v in (sigma_pre,sigma_post,eps_pre,eps_post)):
                    raise ValueError('All six stress and plastic strain components are needed: '+name)
                dwork+=factor*0.5*(sigma_pre+sigma_post)*(eps_post-eps_pre)
            work+=dwork
            yielded=yielded or abs(work)>1e-15
            predicted=beta*work/(RHO*CP)
            observed=post['temperature']-start
            metrics.append(dict(trace=name,t=post['t'],measured=observed,expected=predicted,
                                residual=observed-predicted,plastic_work_density_J_m3=work))
        if not yielded:
            raise ValueError('No plastic work detected; increase imposed strain or decrease yield stress')
    if not metrics:raise ValueError('No valid integrated plastic work data')
    rms=math.sqrt(sum(m['residual']**2 for m in metrics)/len(metrics))
    scale=max(max(abs(m['expected']) for m in metrics),1e-9)
    rel=rms/scale
    result=dict(status='PASS' if rel<=TOL['plastic'] else 'FAIL',
                samples=len(metrics), normalized_rmse=rel, tolerance=TOL['plastic'],
                plastic_heat_fraction=beta, note='Small-strain approximation; stress/plastic-strain tensors must be work-conjugate')
    return metrics,result


def _plot(output_dir,case_id,metrics):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        return None
    if not metrics:
        return None
    figure,ax=plt.subplots(figsize=(7,4))
    # Plot final-time measured and reference, or time histories for contact.
    max_t=max(m['t'] for m in metrics)
    last=[r for r in metrics if abs(r['t']-max_t)<1e-9*max(1.,max_t)]
    if len(last)>2:
        last.sort(key=lambda r:r.get('x',0.))
        ax.plot([r.get('x',0.) for r in last],[r['measured'] for r in last],'o',label='GEOS')
        ax.plot([r.get('x',0.) for r in last],[r['expected'] for r in last],'-',label='reference')
        ax.set_xlabel('x (m)')
    else:
        curves={}
        for m in metrics:
            curves.setdefault(m.get('trace',m.get('field','global')),[]).append(m)
        for name,curve in list(curves.items())[:5]:
            ax.plot([r['t'] for r in curve],[r['measured'] for r in curve],'-',label=name+' GEOS')
            ax.plot([r['t'] for r in curve],[r['expected'] for r in curve],'--',label=name+' reference')
        ax.set_xlabel('time (s)')
    ax.set_ylabel('temperature (K) / metric as labeled')
    ax.set_title(case_id)
    ax.grid(alpha=.3)
    ax.legend(fontsize=8)
    figfile=output_dir/(case_id+'_reference.png')
    figure.tight_layout();figure.savefig(figfile,dpi=160);plt.close(figure)
    return figfile


def _write_csv(filename,rows):
    if not rows:
        filename.write_text('')
        return
    keys=list(dict.fromkeys(k for r in rows for k in r.keys()))
    with filename.open('w',newline='') as h:
        writer=csv.DictWriter(h,fieldnames=keys)
        writer.writeheader();writer.writerows(rows)


def _evaluate(case_id,run_dir,outdir, *, generate_plots=True):
    case=CASES[case_id]; kind=case['kind']
    traces=read_tracers(run_dir)
    box=_box_history(run_dir)
    if kind in ('uniform','cosine','gaussian','step','source','source_local','sink','capacity_linear','conductivity_linear','contact_perfect','dirichlet_steady','dirichlet_transient'):
        metrics,summary=_profile_ref(case_id,case,traces)
    elif kind in ('contact','three_field'):
        metrics,summary=_contact_ref(case_id,case,traces)
    elif kind=='energy_closed':
        metrics,summary=_energy_closed(case_id,case,traces,box)
    elif kind=='thermoelastic':
        metrics,summary=_thermoelastic(case_id,case,traces)
    elif kind=='plastic':
        metrics,summary=_plastic(case_id,case,traces)
    else:
        raise RuntimeError('no post-processor for kind='+kind)
    # The box history measures mass-weighted temperature. With constant cp
    # and fixed mass, it is a rigorous caloric-energy conservation oracle.
    # The historical internalEnergy output is NOT transported thermal energy.
    if kind in ('source','sink'):
        if len(box) < 2:
            raise ValueError('boxAverageHistory required to verify source/sink energy')
        sign = 1 if kind == 'source' else -1
        residuals = [CP*(r['temperature']-T0) - sign * Q_SPECIFIC * r['t'] for r in box]
        energy_rmse = math.sqrt(sum(v*v for v in residuals)/len(residuals))
        summary['specific_energy_rmse_J_per_kg'] = energy_rmse
        summary['specific_energy_tolerance_J_per_kg'] = 0.025*Q_SPECIFIC*case['stop']
        if energy_rmse > summary['specific_energy_tolerance_J_per_kg']:
            summary['status'] = 'FAIL'
            summary['reason'] = 'source/sink energy gain inconsistent with integrated power'
    if kind in ('contact','three_field','contact_perfect'):
        if len(box) < 2:
            raise ValueError('boxAverageHistory required for contact energy conservation')
        ref = box[0]['caloric_energy_J']
        denom = max(RHO*CP*L*Y*Z*AMP,1e-12)
        drift = max(abs(b['caloric_energy_J']-ref) for b in box)/denom
        summary['relative_total_energy_drift'] = drift
        if drift > 1e-4:
            summary['status'] = 'FAIL'
            summary['reason'] = 'internal contact exchange not globally energy-conservative'
    # Optional solver-side caloric energy accounting, independent of the
    # temperature-profile oracle. Diagnostics are available with diagnostics=1
    # but runner environments do not always preserve solver stdout here.
    diag=[]
    for logfile in sorted(set(list(run_dir.glob('*.log')) + list(run_dir.glob('*.out')) +
                              list(run_dir.glob('**/run.log'))))[:20]:
        try:
            content=logfile.read_text(errors='replace')
        except OSError:
            continue
        for match in re.finditer(r'MPMThermal[^\n]*balanceResidual=([-+\d.eE]+)',content):
            value=_number(match.group(1))
            if value is None:
                summary['status']='FAIL'
                summary['reason']='Non-finite thermal balance diagnostic'
            else:
                diag.append(value)
    if diag:
        peak=max(abs(x) for x in diag)
        scale=max(1.,RHO*CP*L*Y*Z*AMP)
        summary['solver_caloric_balance_max_abs_J']=peak
        summary['solver_caloric_balance_samples']=len(diag)
        summary['solver_caloric_balance_threshold_J']=1.e-8*scale
        if peak > summary['solver_caloric_balance_threshold_J']:
            summary['status']='FAIL'
            summary['reason']='Thermal solver reports excessive caloric energy balance residual'
    else:
        summary['solver_caloric_balance']='UNAVAILABLE (no MPMThermal diagnostics in discoverable log)'
    # Machine-readable results and optional illustration.
    summary.update(case=case_id,kind=kind,description=case['description'],
                   run_dir=str(run_dir),tracer_count=len(traces),
                   last_output_time=max(r['t'] for _,rows in traces for r in rows))
    if abs(summary['last_output_time'] - case['stop']) > 0.02*case['stop']:
        summary['status']='FAIL'
        summary['reason']='run did not reach prescribed end time'
    _write_csv(outdir/(case_id+'_metrics.csv'),metrics)
    if generate_plots:
        p=_plot(outdir,case_id,metrics)
        if p:summary['plot']=str(p)
    return summary


def run_case(case_id,argv=None):
    parser=argparse.ArgumentParser(description='Analytical/energy V&V for '+case_id)
    parser.add_argument('--run-dir')
    parser.add_argument('--source-dir',default='.')
    parser.add_argument('--output-dir',default=None)
    parser.add_argument('--case-id',default=case_id)
    parser.add_argument('--no-plot',action='store_true')
    parser.add_argument('--raise-on-failure',action='store_true')
    args,_unknown=parser.parse_known_args(argv)
    outdir=Path(args.output_dir) if args.output_dir else Path(args.source_dir)/'analysis'
    outdir.mkdir(parents=True,exist_ok=True)
    try:
        run_dir=_locate_run(case_id,Path(args.source_dir),args.run_dir)
        summary=_evaluate(case_id,run_dir,outdir,generate_plots=not args.no_plot)
    except Exception as exc:
        summary=dict(case=case_id,status='FAIL',error=type(exc).__name__+': '+str(exc))
    summary_file=outdir/(case_id+'_summary.json')
    summary_file.write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n')
    tex=(r'\paragraph{'+case_id.replace('_',r'\_')+'} '+
         summary['status'].replace('_',r'\_')+'. '+
         (summary.get('reason') or summary.get('error') or '')+'\n')
    (outdir/(case_id+'_results.tex')).write_text(tex)
    print(json.dumps(summary,indent=2,sort_keys=True))
    if args.raise_on_failure and summary.get('status')!='PASS':return 1
    return 0 if summary.get('status')=='PASS' else 2


if __name__=='__main__':
    if len(sys.argv)<2 or sys.argv[1] not in CASES:
        raise SystemExit('usage: python thermal_post.py CASE_ID [--run-dir DIR] [--output-dir DIR]')
    raise SystemExit(run_case(sys.argv[1],sys.argv[2:]))
