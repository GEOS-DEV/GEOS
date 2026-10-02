#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-2.1-only
"""Real-solver checks for the fixed-composition EOS/capillary boundary-value solve.

The baseline is the preceding constant-density capillary-correctness checkpoint.
No mocked numerical states are used. Every solver invocation writes durable VTK
and logs. Table error, spatial TPFA error, and solver roundoff are separate checks.
"""
import argparse
import copy
import json
import math
import pathlib
import unittest
import xml.etree.ElementTree as ET
import testHydrostaticCapillary as base


def configure(tree, phases=('gas','oil'), entry=0., maximum=20000., increment=.1,
              tolerance=1e-6, solve=False, invariant=False, mass=True, chopping=False):
    constitutive=tree.find('Constitutive')
    fluid=tree.find('./Constitutive/InvariantImmiscibleFluid')
    constitutive.remove(fluid)
    n=len(phases)
    phases_xml='{'+','.join(phases)+'}'
    densities={'gas':500.,'oil':800.,'water':1000.}
    molecular={'gas':.044,'oil':.114,'water':.018}
    props=dict(name='fluid',phaseNames=phases_xml,componentNames=phases_xml,
               componentMolarWeight=base.vector(molecular[p] for p in phases))
    if invariant:
        props.update(densities=base.vector(densities[p] for p in phases),viscosities=base.vector([.001]*n))
        ET.SubElement(constitutive,'InvariantImmiscibleFluid',**props)
    else:
        hc=[p for p in phases if p!='water']
        # The inherited DeadOil reader expects two named tables if gas exists,
        # including gas/water. Its active hydrocarbon order still contains gas only.
        names=hc if len(hc)==2 or 'gas' not in phases else hc+hc
        props.update(surfaceDensities=base.vector(densities[p] for p in phases),
                     hydrocarbonFormationVolFactorTableNames='{'+','.join('B'+p for p in names)+'}',
                     hydrocarbonViscosityTableNames='{'+','.join('v'+p for p in names)+'}')
        if 'water' in phases:
            props.update(waterReferencePressure='10000000',waterFormationVolumeFactor='1',
                         waterCompressibility='1e-8',waterViscosity='.001')
        ET.SubElement(constitutive,'DeadOilFluid',**props)
    rp=tree.find('./Constitutive/BrooksCoreyRelativePermeability')
    rp.set('phaseNames',phases_xml)
    for key,value in [('phaseMinVolumeFraction',[0]*n),('phaseRelPermExponent',[2]*n),('phaseRelPermMaxValue',[1]*n)]: rp.set(key,base.vector(value))
    tree.find('./ElementRegions/CellElementRegion').set('materialList','{fluid,rock,relperm,capillary}')
    eq=tree.find('./FieldSpecifications/HydrostaticEquilibrium')
    offset=float(eq.get('datumElevation'))-10
    eq.set('componentNames',phases_xml)
    eq.set('componentFractionVsElevationTableNames','{'+','.join('z'+p for p in phases)+'}')
    eq.set('phaseContacts',base.vector(v+offset for v in ([3,7] if n==3 else [5])))
    eq.set('elevationIncrementInHydrostaticPressureTable',str(increment))
    eq.set('equilibrationTolerance',str(tolerance))
    eq.set('maxNumberOfEquilibrationIterations','100')
    functions=tree.find('Functions')
    for child in list(functions):
        if child.get('name')!='temperature': functions.remove(child)
    for phase in phases:
        ET.SubElement(functions,'TableFunction',name='z'+phase,coordinates=base.vector([offset,offset+10]),values=base.vector([1/n,1/n]))
    if not invariant:
        for phase in phases:
            if phase=='water': continue
            ET.SubElement(functions,'TableFunction',name='B'+phase,coordinates='{1000000,10000000,20000000}',
                          values='{1.27,1,.7}' if phase=='gas' else '{1.09,1,.9}')
            ET.SubElement(functions,'TableFunction',name='v'+phase,coordinates='{1000000,10000000,20000000}',values='{.001,.001,.001}')
    if n==3:
        options=dict(wettingIntermediateCapPressureTableName='pcwater',nonWettingIntermediateCapPressureTableName='pcgas')
        ET.SubElement(functions,'TableFunction',name='pcwater',coordinates='{0,1}',values=base.vector([maximum+entry,entry]))
        ET.SubElement(functions,'TableFunction',name='pcgas',coordinates='{0,1}',values=base.vector([entry,maximum+entry]))
    else:
        options=dict(wettingNonWettingCapPressureTableName='pc')
        values=[maximum+entry,entry] if 'water' in phases else [entry,maximum+entry]
        ET.SubElement(functions,'TableFunction',name='pc',coordinates='{0,1}',values=base.vector(values))
    ET.SubElement(constitutive,'TableCapillaryPressure',name='capillary',phaseNames=phases_xml,**options)
    tree.find('./Solvers/CompositionalMultiphaseFVM').set('useMass','1' if mass else '0')
    tree.find('./Solvers/CompositionalMultiphaseFVM').set('allowLocalCompDensityChopping','1' if chopping else '0')
    if solve:
        tree.find('./Constitutive/PressurePorosity').set('compressibility','1e-9')
        ET.SubElement(tree.find('Events'),'PeriodicEvent',name='solve',forceDt='1',target='/Solvers/flow')


def flux_metrics(values, gravity=9.81):
    zs=sorted({key[2] for key in values})
    dz=zs[1]-zs[0] if len(zs)>1 else 1.
    potential=flux=0.
    for (x,y,z),lower in values.items():
        upper=values.get((x,y,round(z+dz,6)))
        if upper is None: continue
        for ip in range(len(lower['phaseVolumeFraction'])):
            p0=lower['pressure'][0]-lower['capillary_phaseCapPressure'][ip]
            p1=upper['pressure'][0]-upper['capillary_phaseCapPressure'][ip]
            rho=.5*(lower['fluid_phaseMassDensity'][ip]+upper['fluid_phaseMassDensity'][ip])
            residual=p0-p1-rho*gravity*dz
            upstream=lower if residual>=0 else upper
            mobility=upstream['phaseMobility'][ip]
            if mobility>1e-12: potential=max(potential,abs(residual))
            flux=max(flux,abs(1e-13*mobility*residual/dz))
    return dict(max_mobile_potential_pa=potential,max_phase_mass_flux_per_area=flux)


def analytic_gas_oil_reference(z):
    """Independent integrated-EOS oracle while oil is present throughout.

    For B_o=1-c_o(P-P0), integrate B_o dP=-rho_o_surface g dz.
    Integrate dp_g/dP=rho_g(P)/rho_o(P), then use p_g(10)=P0
    and p_g(5)=p_o(5). These are closed-form primitives, not the
    production fixed-point/marching or capillary inversion algorithms.
    """
    c=1e-8
    def oil_primitive(delta): return delta-.5*c*delta*delta
    def inverse_oil(value): return 2*value/(1+math.sqrt(1-2*c*value))
    def gas_primitive(delta): return .625*(delta/3-(2/(9e-8))*math.log1p(-3e-8*delta))
    low,high=0.,100000.
    for _ in range(80):
        contact=.5*(low+high)
        datum=inverse_oil(oil_primitive(contact)-800*9.81*5)
        error=contact+gas_primitive(datum)-gas_primitive(contact)
        if error>0: high=contact
        else: low=contact
    contact=.5*(low+high)
    delta=inverse_oil(oil_primitive(contact)-800*9.81*(z-5))
    gas_delta=contact+gas_primitive(delta)-gas_primitive(contact)
    return 1e7+delta,max(0.,min(1.,(gas_delta-delta)/20000.))


class SelfConsistentCapillary(unittest.TestCase):
    case=base.HydrostaticCapillary.case
    compare_fields=base.HydrostaticCapillary.compare_fields

    def check_closure(self,values,phases=('gas','oil'),mass=True,invariant=False):
        molecular={'gas':.044,'oil':.114,'water':.018}
        surface={'gas':500.,'oil':800.,'water':1000.}
        for key,f in values.items():
            p=f['pressure'][0]
            self.assertGreater(p,0)
            densities=[]
            for ip,phase in enumerate(phases):
                expected=surface[phase] if invariant else (surface[phase]*math.exp(1e-8*(p-1e7)) if phase=='water' else surface[phase]/(1-(3e-8 if phase=='gas' else 1e-8)*(p-1e7)))
                self.assertAlmostEqual(f['fluid_phaseMassDensity'][ip],expected,delta=expected*2e-12)
                density=expected if mass else expected/molecular[phase]
                self.assertAlmostEqual(f['fluid_phaseDensity'][ip],density,delta=density*2e-12)
                densities.append(density)
            volumes=[z/rho for z,rho in zip(f['globalCompFraction'],densities)]
            total=sum(volumes)
            for actual,expected in zip(f['phaseVolumeFraction'],volumes):
                self.assertAlmostEqual(actual,expected/total,delta=2e-12)
                self.assertGreaterEqual(actual,-1e-12)
                self.assertLessEqual(actual,1+1e-12)
            self.assertAlmostEqual(sum(f['phaseVolumeFraction']),1,delta=2e-12)

    def test_counterexample_eos_and_actual_rest_state(self):
        _,old=self.case('old-compressible-counterexample',mode=None,nz=10,baseline=True,
                        mutate=lambda t:configure(t,tolerance=.001))
        old_metrics=flux_metrics(old)
        self.assertGreater(old_metrics['max_mobile_potential_pa'],.1)
        directory,current=self.case('consistent-counterexample',mode='gravityAligned',nz=10,
                                    mutate=lambda t:configure(t,solve=True))
        self.check_closure(current)
        current_metrics=flux_metrics(current)
        oracle_error=0.
        for (_,_,z),f in current.items():
            exact_pressure,exact_gas=analytic_gas_oil_reference(z)
            oracle_error=max(oracle_error,abs(f['pressure'][0]-exact_pressure))
            self.assertAlmostEqual(f['phaseVolumeFraction'][0],exact_gas,delta=2e-8)
        self.assertLess(oracle_error,.001)
        self.assertLess(current_metrics['max_mobile_potential_pa'],.001)
        self.assertLess(current_metrics['max_phase_mass_flux_per_area'],1e-10)
        final=base.records(directory,final=True)
        self.check_closure(final)
        for key in current:
            self.assertLess(abs(current[key]['pressure'][0]-final[key]['pressure'][0]),.01)
            for a,b in zip(current[key]['phaseVolumeFraction'],final[key]['phaseVolumeFraction']): self.assertLess(abs(a-b),1e-6)
        EVIDENCE['counterexample']={'before':old_metrics,'after':current_metrics,'after_step':flux_metrics(final),'analytic_pressure_error_pa':oracle_error}

    def test_phase_arrangements_permutations_and_molar_basis(self):
        for phases in [('gas','oil'),('oil','water'),('gas','water'),('gas','oil','water'),('water','oil','gas')]:
            with self.subTest(phases=phases):
                _,values=self.case('phases-'+'-'.join(phases),mode='gravityAligned',nz=10,
                                   mutate=lambda t:configure(t,phases=phases,maximum=10000 if len(phases)==3 else 20000,increment=.025))
                self.check_closure(values,phases)
                metrics=[flux_metrics(values)]
                if 'water' in phases:
                    # At appearance of the primary oil/gas reference phase,
                    # P is continuous but dP/dz can change. The actual TPFA
                    # interface flux is first-order on a coarse spatial mesh;
                    # establish convergence at a fixed finer integration table.
                    for nz in (20,40):
                        _,fine=self.case('interface-refinement-'+'-'.join(phases)+'-'+str(nz),mode='gravityAligned',nz=nz,
                                         mutate=lambda t:configure(t,phases=phases,maximum=10000 if len(phases)==3 else 20000,increment=.025))
                        self.check_closure(fine,phases)
                        metrics.append(flux_metrics(fine))
                    for coarse,fine in zip(metrics,metrics[1:]):
                        self.assertLess(fine['max_phase_mass_flux_per_area'],.7*coarse['max_phase_mass_flux_per_area'])
                else:
                    self.assertLess(metrics[0]['max_phase_mass_flux_per_area'],1e-10)
                EVIDENCE.setdefault('phase_arrangement_spatial_refinement',{})[','.join(phases)]=metrics
        _,molar=self.case('molar-basis',mode='gravityAligned',nz=10,mutate=lambda t:configure(t,mass=False))
        self.check_closure(molar,mass=False)

    def test_covariance_legacy_and_zero_gravity(self):
        _,reference=self.case('downward',mode='gravityAligned',nz=10,mutate=configure)
        _,legacy=self.case('legacy-corrected',mode=None,nz=10,mutate=configure)
        self.compare_fields(reference,legacy,1e-7)
        for name,rotation in [('upward',base.UPWARD),('oblique',base.oblique_rotation())]:
            _,values=self.case(name,mode='gravityAligned',rotation=rotation,translation=(13.,-8.,5.),nz=10,mutate=configure)
            self.compare_fields(reference,values,2e-6)
        _,zero=self.case('zero-gravity',mode='gravityAligned',gravity=(0.,0.,0.),nz=10,mutate=configure)
        self.check_closure(zero)
        for f in zero.values(): self.assertAlmostEqual(f['pressure'][0],1e7,delta=1e-7)
        self.assertLess(flux_metrics(zero,0.)['max_phase_mass_flux_per_area'],1e-14)

    def test_entry_presence_and_reference_switching(self):
        _,values=self.case('positive-entry-and-endpoint',mode='gravityAligned',nz=40,
                           mutate=lambda t:configure(t,entry=1000.,maximum=5000.,increment=.025))
        self.check_closure(values)
        saturations=[f['phaseVolumeFraction'][0] for f in values.values()]
        self.assertTrue(any(s==0 for s in saturations))
        self.assertTrue(any(s==1 for s in saturations))
        self.assertTrue(any(0<s<1 for s in saturations))
        self.assertLess(flux_metrics(values)['max_phase_mass_flux_per_area'],1e-9)
        # Datum lies above zero-Pc contact but below entry: anchor is an absent
        # continued gas phase, while present oil reconstructs primary pressure.
        def datum_in_entry(t):
            configure(t,entry=20000.)
            t.find('./FieldSpecifications/HydrostaticEquilibrium').set('datumElevation','5.1')
        _,continued=self.case('absent-datum-phase',mode='gravityAligned',nz=10,mutate=datum_in_entry)
        self.check_closure(continued)
        self.assertTrue(all(f['phaseVolumeFraction'][0]==0 for f in continued.values()))

    def test_table_and_spatial_convergence_are_separate(self):
        def stronger_compressibility(t,h):
            # Increase curvature so the convergence signal stays well above
            # floating-point accumulation noise; retain positive linear PVT.
            configure(t,increment=h,tolerance=1e-7)
            for name,values in [('Bgas','{1.3,1,.7}'),('Boil','{1.1,1,.9}')]:
                table=t.find("./Functions/TableFunction[@name='"+name+"']")
                table.set('coordinates','{9000000,10000000,11000000}')
                table.set('values',values)
        profiles=[]
        for h in [.5,.25,.125,.03125]:
            _,values=self.case('table-'+str(h),mode='gravityAligned',nz=8,
                               mutate=lambda t,h=h:stronger_compressibility(t,h))
            profiles.append(values)
        truth=profiles[-1]
        errors=[max(abs(f['pressure'][0]-truth[k]['pressure'][0]) for k,f in v.items()) for v in profiles[:-1]]
        self.assertGreater(errors[0],errors[1]*2.)
        self.assertGreater(errors[1],errors[2]*2.)
        spatial=[]
        for nz in [10,20,40]:
            _,values=self.case('spatial-'+str(nz),mode='gravityAligned',nz=nz,
                               mutate=lambda t:stronger_compressibility(t,.00625))
            spatial.append(flux_metrics(values))
        self.assertGreater(spatial[0]['max_phase_mass_flux_per_area'],spatial[1]['max_phase_mass_flux_per_area']*2.)
        self.assertGreater(spatial[1]['max_phase_mass_flux_per_area'],spatial[2]['max_phase_mass_flux_per_area']*2.)
        EVIDENCE['convergence']={'table_max_pressure_error_pa':errors,'spatial_tpfa':spatial}

    def test_invalid_pvt_endpoints_and_bounded_iterations(self):
        def mutate_table(t,values):
            configure(t)
            t.find("./Functions/TableFunction[@name='Bgas']").set('values',values)
        for name,mutation,diagnostic in [
            ('increasing-fvf',lambda t:mutate_table(t,'{.7,1,1.3}'),'nonincreasing formation-volume'),
            ('negative-fvf',lambda t:mutate_table(t,'{-1,-1,-1}'),'positive properties'),
            ('iterations',lambda t:(configure(t),t.find('./FieldSpecifications/HydrostaticEquilibrium').set('maxNumberOfEquilibrationIterations','1')),'did not converge'),
            ('incompatible-endpoint',lambda t:(configure(t),t.find("./Functions/TableFunction[@name='pc']").set('coordinates','{.1,1}')),'inconsistent mobile-phase'),
        ]:
            _,log=self.case(name,mode='gravityAligned',nz=10,mutate=mutation,success=False)
            self.assertIn(diagnostic,log)
        def surface_special(t):
            configure(t)
            table=t.find("./Functions/TableFunction[@name='Bgas']")
            table.set('coordinates','{0,1000000,10000000,20000000}')
            table.set('values','{.5,1.27,1,.7}')
        _,values=self.case('unused-surface-interval',mode='gravityAligned',nz=10,mutate=surface_special)
        self.check_closure(values)
        def skipped_bad_bulk(t):
            configure(t,increment=5.)
            table=t.find("./Functions/TableFunction[@name='Bgas']")
            table.set('coordinates','{0,1000000,2000000,3000000,10000000,20000000}')
            table.set('values','{1.3,1.27,1.28,1.28,1,.7}')
        _,log=self.case('invalid-unsampled-bulk-interval',mode='gravityAligned',nz=10,mutate=skipped_bad_bulk,success=False)
        self.assertIn('nonincreasing formation-volume',log)

    def test_legacy_corrected_controls_are_bounded(self):
        def tiny(t):
            configure(t,increment=1e-5)
        _,log=self.case('legacy-unresolvable-increment',mode=None,translation=(0.,0.,1e12),nz=10,mutate=tiny,success=False)
        self.assertIn('unresolvable',log)
        for key,value in [('equilibrationTolerance','0'),('maxNumberOfEquilibrationIterations','0'),('datumPressure','0')]:
            def invalid(t,key=key,value=value):
                configure(t);t.find('./FieldSpecifications/HydrostaticEquilibrium').set(key,value)
            _,log=self.case('legacy-invalid-'+key,mode=None,nz=10,mutate=invalid,success=False)
            self.assertIn('positive',log)

    def test_file_backed_pvt_and_piecewise_capillary(self):
        gas=base.ROOT/'gas.pvt'
        oil=base.ROOT/'oil.pvt'
        gas.write_text('1000000 1.27 .001\n10000000 1 .001\n20000000 .7 .001\n')
        oil.write_text('1000000 1.09 .001\n10000000 1 .001\n20000000 .9 .001\n')
        def files(t):
            configure(t)
            fluid=t.find('./Constitutive/DeadOilFluid')
            fluid.attrib.pop('hydrocarbonFormationVolFactorTableNames')
            fluid.attrib.pop('hydrocarbonViscosityTableNames')
            fluid.set('tableFiles','{'+str(gas)+','+str(oil)+'}')
            cap=t.find("./Functions/TableFunction[@name='pc']")
            cap.set('coordinates','{0,.1,.2,.5,1}')
            cap.set('values','{0,1000,2500,9000,20000}')
        _,values=self.case('file-backed-piecewise',mode='gravityAligned',nz=10,mutate=files)
        self.check_closure(values)
        self.assertLess(flux_metrics(values)['max_phase_mass_flux_per_area'],1e-10)
        gas.write_text('1000000 .7 .001\n10000000 1 .001\n20000000 1.3 .001\n')
        _,log=self.case('inadmissible-file-pvt',mode='gravityAligned',nz=10,mutate=files,success=False)
        self.assertIn('nonincreasing formation-volume',log)

    def test_tiny_saturation_is_not_zero_mobility(self):
        def mobile_endpoint(t):
            configure(t,invariant=True,entry=1000.)
            cap=t.find("./Functions/TableFunction[@name='pc']")
            cap.set('coordinates','{5e-9,1}')
            t.find('./Constitutive/BrooksCoreyRelativePermeability').set('phaseRelPermExponent','{.1,.1}')
        _,log=self.case('tiny-but-mobile-endpoint',mode='gravityAligned',nz=10,mutate=mobile_endpoint,success=False)
        self.assertIn('inconsistent mobile-phase',log)

    def test_realized_residual_roundtrip_and_chopping_refusal(self):
        def roundtrip(t):
            configure(t,invariant=True,entry=1000.,chopping=False)
            t.find("./Functions/TableFunction[@name='pc']").set('coordinates','{.01,1}')
            rp=t.find('./Constitutive/BrooksCoreyRelativePermeability')
            rp.set('phaseMinVolumeFraction','{.01,0}'); rp.set('phaseRelPermExponent','{.1,2}')
        _,log=self.case('residual-one-ulp-mobile',mode='gravityAligned',nz=10,mutate=roundtrip,success=False)
        self.assertIn('Realized hydrostatic EOS/capillary state',log)
        def chopping(t):
            configure(t,invariant=True,entry=1000.,chopping=True)
            t.find('./Constitutive/BrooksCoreyRelativePermeability').set('phaseRelPermExponent','{.1,2}')
        _,log=self.case('default-chopping-mobile-phase',mode='gravityAligned',nz=10,mutate=chopping,success=False)
        self.assertIn('allowLocalCompDensityChopping',log)

    def test_exact_entry_complementarity(self):
        def entry(t):
            configure(t,invariant=True,entry=1471.4999995,tolerance=1e-8)
            t.find('./FieldSpecifications/HydrostaticEquilibrium').set('datumPressure','100000')
        _,values=self.case('entry-without-fixed-snap-band',mode='gravityAligned',nz=10,mutate=entry)
        field=values[(1.25,1.25,5.5)]
        self.assertGreater(field['phaseVolumeFraction'][0],0.)
        self.assertAlmostEqual(field['phaseVolumeFraction'][0],2.5e-11,delta=3e-13)
        self.assertLess(flux_metrics(values)['max_mobile_potential_pa'],5e-8)

    def test_actual_three_phase_relative_permeability(self):
        def relperm(t,family,hidden=False):
            configure(t,invariant=True,phases=('gas','oil','water'),maximum=10000.)
            rp=t.find('./Constitutive/BrooksCoreyRelativePermeability')
            for key in ('phaseRelPermExponent','phaseRelPermMaxValue'): rp.attrib.pop(key)
            if family=='table':
                rp.tag='TableRelativePermeability'
                rp.attrib.pop('phaseMinVolumeFraction')
                rp.set('wettingIntermediateRelPermTableNames','{krw,krow}')
                rp.set('nonWettingIntermediateRelPermTableNames','{krg,krog}')
                rp.set('threePhaseInterpolator','BAKER')
                for name in ('krw','krow','krg','krog'):
                    ET.SubElement(t.find('Functions'),'TableFunction',name=name,coordinates='{.8,1}' if hidden and name=='krog' else '{0,1}',values='{0,1}')
                if hidden:
                    cap=t.find("./Functions/TableFunction[@name='pcwater']")
                    cap.set('coordinates','{.5,1}'); cap.set('values','{1000,0}')
                    t.find("./Functions/TableFunction[@name='pcgas']").set('values','{1000000,1001000}')
            else:
                rp.tag='BrooksCoreyStone2RelativePermeability'
                rp.set('phaseMinVolumeFraction','{0,.1,0}')
                for key,value in [('waterOilRelPermExponent','{2,.5}'),('gasOilRelPermExponent','{2,.5}'),('waterOilRelPermMaxValue','{1,1}'),('gasOilRelPermMaxValue','{1,1}')]: rp.set(key,value)
        for family in ('table','stone2'):
            _,values=self.case('actual-rp-'+family,mode='gravityAligned',nz=10,mutate=lambda t,family=family:relperm(t,family))
            self.check_closure(values,('gas','oil','water'),invariant=True)
            self.assertLess(flux_metrics(values)['max_phase_mass_flux_per_area'],1e-11)
        _,log=self.case('table-baker-hidden-mobile-oil',mode='gravityAligned',nz=10,mutate=lambda t:relperm(t,'table',True),success=False)
        self.assertIn('inconsistent mobile-phase',log)

    def test_uniform_constitutive_union_serial_and_mpi(self):
        def split(t,kind,selected=False):
            configure(t)
            mesh=t.find('Mesh'); mesh.clear()
            ET.SubElement(mesh,'InternalMesh',name='mesh',elementTypes='{C3D8}',xCoords='{0,5,10}',yCoords='{0,10}',zCoords='{0,10}',nx='{2,2}',ny='{4}',nz='{10}',cellBlockNames='{leftCells,rightCells}')
            regions=t.find('ElementRegions'); regions.clear()
            ET.SubElement(regions,'CellElementRegion',name='leftRegion',cellBlocks='{leftCells}',materialList='{fluid,rock,relperm,capillary}')
            right='{fluid,rock,relperm,capillaryB}' if kind=='capillary' else '{fluid,rock,relpermB,capillary}'
            ET.SubElement(regions,'CellElementRegion',name='rightRegion',cellBlocks='{rightCells}',materialList=right)
            t.find('./Solvers/CompositionalMultiphaseFVM').set('targetRegions','{leftRegion,rightRegion}')
            original=t.find('./Constitutive/TableCapillaryPressure' if kind=='capillary' else './Constitutive/BrooksCoreyRelativePermeability')
            other=copy.deepcopy(original); other.set('name','capillaryB' if kind=='capillary' else 'relpermB'); t.find('Constitutive').append(other)
            eq=t.find('./FieldSpecifications/HydrostaticEquilibrium'); eq.set('objectPath','ElementRegions/*')
            if selected:
                geometry=ET.SubElement(t,'Geometry'); ET.SubElement(geometry,'Box',name='left',xMin='{0,0,0}',xMax='{5,10,10}')
                eq.set('setNames','{left}')
                for name,field,component,value in [('pressure','pressure','-1','11000000'),('fraction0','globalCompFraction','0','.5'),('fraction1','globalCompFraction','1','.5')]:
                    ET.SubElement(t.find('FieldSpecifications'),'FieldSpecification',name=name,initialCondition='1',setNames='{all}',objectPath='ElementRegions/*',fieldName=field,component=component,scale=value)
        for kind in ('capillary','relperm'):
            for mode in (None,'gravityAligned'):
                for mpi in (False,True):
                    _,log=self.case('mixed-'+kind+'-'+str(mode)+'-'+str(mpi),mode=mode,nz=10,mpi=mpi,mutate=lambda t,kind=kind:split(t,kind),success=False)
                    self.assertIn('one spatially uniform fluid/capillary/relative-permeability',log)
        _,serial=self.case('uniform-selected-subregion',mode='gravityAligned',nz=10,mutate=lambda t:split(t,'capillary',True))
        _,parallel=self.case('uniform-selected-subregion-mpi',mode='gravityAligned',nz=10,mpi=True,mutate=lambda t:split(t,'capillary',True))
        self.compare_fields(serial,parallel,2e-6)

    def test_explicit_unsupported_model_refusals(self):
        def jfunction(t):
            configure(t)
            cap=t.find('./Constitutive/TableCapillaryPressure')
            cap.tag='JFunctionCapillaryPressure'
            cap.attrib.pop('wettingNonWettingCapPressureTableName')
            cap.set('wettingNonWettingJFunctionTableName','pc')
            cap.set('wettingNonWettingSurfaceTension','.02')
            cap.set('permeabilityDirection','X')
        _,log=self.case('rock-dependent-capillary-refused',mode='gravityAligned',nz=10,mutate=jfunction,success=False)
        self.assertIn('rock-dependent capillarity are not supported',log)
        def general_composition(t):
            configure(t,phases=('gas','water'))
            old=t.find('./Constitutive/DeadOilFluid')
            t.find('Constitutive').remove(old)
            folder=pathlib.Path(__file__).resolve().parents[1]/'inputFiles/compositionalMultiphaseFlow'
            ET.SubElement(t.find('Constitutive'),'CO2BrinePhillipsFluid',name='fluid',phaseNames='{gas,water}',componentNames='{co2,water}',
                          componentMolarWeight='{.044,.018}',phasePVTParaFiles='{'+str(folder/'pvtgas.txt')+','+str(folder/'pvtliquid.txt')+'}',
                          flashModelParaFile=str(folder/'co2flash.txt'))
            t.find('./FieldSpecifications/HydrostaticEquilibrium').set('componentNames','{co2,water}')
        _,log=self.case('general-composition-refused',mode='gravityAligned',nz=10,mutate=general_composition,success=False)
        self.assertIn('general compositional flash',log)

    def test_unselected_fields_match_initializer_free_control(self):
        def isolated(t,apply):
            configure(t,maximum=100000.,chopping=True)
            geometry=ET.SubElement(t,'Geometry')
            ET.SubElement(geometry,'Box',name='left',xMin='{0,0,0}',xMax='{2.5,10,10}')
            specs=t.find('FieldSpecifications')
            eq=specs.find('HydrostaticEquilibrium');eq.set('setNames','{left}');eq.set('phaseContacts','{-1}')
            if not apply: specs.remove(eq)
            for name,field,component,value in [('pressure','pressure','-1','11000000'),('fraction0','globalCompFraction','0','0'),('fraction1','globalCompFraction','1','1')]:
                ET.SubElement(specs,'FieldSpecification',name=name,initialCondition='1',setNames='{all}',objectPath='ElementRegions/reservoir',fieldName=field,component=component,scale=value)
        _,control=self.case('no-initializer-control',mode='gravityAligned',nz=10,mutate=lambda t:isolated(t,False))
        _,selected=self.case('selected-isolation-full-fields',mode='gravityAligned',nz=10,mutate=lambda t:isolated(t,True))
        for key,expected in control.items():
            if key[0]<=2.5: continue
            actual=selected[key]
            self.assertEqual(actual.keys(),expected.keys())
            for field in expected: self.assertEqual(actual[field],expected[field],(key,field,actual[field],expected[field]))

    def test_three_phase_oil_disappearance(self):
        def disappear(t,invariant):
            configure(t,invariant=invariant,phases=('gas','oil','water'),increment=.01)
            t.find("./Functions/TableFunction[@name='pcwater']").set('values','{20000,0}')
            t.find("./Functions/TableFunction[@name='pcgas']").set('values','{0,1000}')
        for invariant in (True,False):
            _,values=self.case('oil-disappearance-'+str(invariant),mode='gravityAligned',nz=40,mutate=lambda t,invariant=invariant:disappear(t,invariant))
            self.check_closure(values,('gas','oil','water'),invariant=invariant)
            reduced=[f for f in values.values() if f['phaseVolumeFraction'][1]==0 and min(f['phaseVolumeFraction'][0],f['phaseVolumeFraction'][2])>0]
            self.assertGreater(len(reduced),0)
            if invariant:
                for (_,_,z),f in values.items():
                    if f['phaseVolumeFraction'][1]==0 and min(f['phaseVolumeFraction'][0],f['phaseVolumeFraction'][2])>0:
                        self.assertAlmostEqual(f['phaseVolumeFraction'][0],(4905*z-26487)/21000,delta=1e-10)
            self.assertLess(flux_metrics(values)['max_phase_mass_flux_per_area'],1e-9)

    def test_native_sets_and_empty_mpi_rank(self):
        def native(t):
            configure(t)
            geometry=ET.SubElement(t,'Geometry')
            ET.SubElement(geometry,'Box',name='left',xMin='{0,0,0}',xMax='{2.5,10,10}')
            specs=t.find('FieldSpecifications')
            eq=specs.find('HydrostaticEquilibrium')
            eq.set('setNames','{left}')
            # Supply physically admissible explicit states only outside the
            # selected set, then check the selected native set is the sole target.
            ET.SubElement(specs,'FieldSpecification',name='pressure',initialCondition='1',setNames='{all}',objectPath='ElementRegions/reservoir',fieldName='pressure',scale='11000000')
            for ic in range(2):
                ET.SubElement(specs,'FieldSpecification',name='fraction'+str(ic),initialCondition='1',setNames='{all}',objectPath='ElementRegions/reservoir',fieldName='globalCompFraction',component=str(ic),scale='.5')
        _,serial=self.case('native-set',mode='gravityAligned',nz=10,mutate=native)
        _,parallel=self.case('native-set-empty-rank',mode='gravityAligned',nz=10,mutate=native,mpi=True)
        self.compare_fields(serial,parallel,2e-6)
        for (x,y,z),f in serial.items():
            if x>2.5: self.assertEqual(f['pressure'][0],11000000.)
            else: self.assertNotEqual(f['pressure'][0],11000000.)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--geos',required=True)
    parser.add_argument('--baseline',required=True)
    parser.add_argument('--mpiexec',required=True)
    parser.add_argument('--output',required=True)
    parser.add_argument('--test',action='append',help='Run named methods while developing; omit for complete acceptance')
    args=parser.parse_args()
    base.GEOS=str(pathlib.Path(args.geos).absolute())
    base.BASELINE=str(pathlib.Path(args.baseline).absolute())
    base.MPIEXEC=str(pathlib.Path(args.mpiexec).absolute())
    base.ROOT=pathlib.Path(args.output).absolute()
    base.ROOT.mkdir(parents=True,exist_ok=False)
    EVIDENCE={}
    suite=unittest.TestSuite(SelfConsistentCapillary(name) for name in args.test) if args.test else unittest.defaultTestLoader.loadTestsFromTestCase(SelfConsistentCapillary)
    result=unittest.TextTestRunner(verbosity=2).run(suite)
    (base.ROOT/'numerical-evidence.json').write_text(json.dumps(EVIDENCE,indent=2)+'\n')
    raise SystemExit(not result.wasSuccessful())
