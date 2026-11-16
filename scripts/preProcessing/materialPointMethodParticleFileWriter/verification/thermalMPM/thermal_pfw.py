"""GEOS MPM thermal verification inputs using the Thermal child XML interface.

This API is documented by geos_mpm_thermal_extension.patch.  The caller may
override the mesh size with MPM_THERMAL_NX; input creation requires the GEOS
PFW utilities on the Python import path (as for other verification inputs).
"""
from __future__ import annotations

import os
from thermal_cases import (CASES, L, Y, T0, AMP, CP, K, RHO,
                           Q_SPECIFIC, RC_GLOBAL, RC_MATERIAL, RC_GROUP,
                           THERMAL_EXPANSION, CAPACITY_SLOPE,
                           CONDUCTIVITY_SLOPE, RC_COMBINED, temperature)
import pfw_geometryObjects as geom
import pfw_tracerPoints as tracers


def _material_xml(name, conductivity, *, expand=False, plastic=False,
                  capacity_linear=False, conductivity_linear=False):
    tag = 'ThermalVonMisesJ' if plastic else (
        'LinearThermalElasticIsotropic' if (capacity_linear or conductivity_linear)
        else 'ThermalElasticIsotropic')
    attributes = {
        'name': name,
        'defaultDensity': RHO,
        'defaultBulkModulus': 666.6666666667,
        'defaultShearModulus': 400.0,
        'thermalSpecificHeat': CP,
        'thermalConductivity': conductivity,
        'thermalReferenceTemperature': T0,
    }
    if capacity_linear: attributes['thermalHeatCapacitySlope']=CAPACITY_SLOPE
    if conductivity_linear: attributes['thermalConductivitySlope']=CONDUCTIVITY_SLOPE
    if expand:
        attributes.update(solverThermalExpansion=1,
                          defaultDrainedLinearTEC=THERMAL_EXPANSION)
    if plastic:
        attributes['defaultYieldStrength'] = 1.0
    return '<{}\n  {}/>'.format(
        tag, '\n  '.join('{}="{}"'.format(k,v) for k,v in attributes.items()))


def _box(name, xa, xb, ya, yb, material, group, temp, *, surface=None, damage=0):
    surfaces = [False, False, False, False]
    if surface == 'left':
        surfaces[0] = True
    elif surface == 'right':
        surfaces[2] = True
    ob = geom.box(name, [float(xa),float(ya)], [float(xb),float(yb)],
                  dim=2, vel=[0.,0.,0.], mat=material, group=group,
                  flaggedSurfaces=surfaces, particleType=2)
    ob.temperature = temp
    ob.damage = damage
    return ob


def build(case_id):
    if case_id not in CASES:
        raise ValueError('unknown thermal verification case '+case_id)
    case = CASES[case_id]
    kind = case['kind']
    nx = int(os.environ.get('MPM_THERMAL_NX',case.get('nx',34)))
    if nx < 10:
        raise ValueError('MPM_THERMAL_NX must be >=10')
    dx = L/(nx-2)
    stop = float(case['stop'])
    dfg = case_id=='contact_dfg'
    contact = kind in ('contact','contact_perfect','three_field')
    # Physical units are SI. PFW does NOT automatically convert to mm/us/mg.
    pfw = {
        'caseName': 'thermalVV_'+case_id,
        'runDebug': False, 'mBatch': True, 'mSubmitJobs': False,
        'mCores': int(case.get('ranks',1)), 'mWallTime': '00:30:00',
        'autoRestart': False,
        'xmin':0.0, 'xmax':L, 'ymin':-Y/2, 'ymax':Y/2,
        'zmin':-dx/2, 'zmax':dx/2,
        'planeStrain':1,
        'periodic': [bool(case.get('periodic',False)),False,False],
        'xpar':int(case.get('xpar',1)), 'ypar':1, 'zpar':1,
        'nI':nx, 'nJ':6, 'nK':3, 'ppc':2,
        'endTime':stop, 'plotInterval':stop/5., 'restartInterval':stop*2.,
        'timeIntegrationOption':'ExplicitDynamic',
        'cflFactor':0.2, 'initialDt':1.e-9,
        'updateMethod':'PIC', 'cpdiDomainScaling':0,
        'subdivideParticles':0,
        'damageFieldPartitioning':int(dfg),
        'needsNeighborList':int(dfg),
        'enableContact':int(contact),
        'frictionCoefficient':0.0,
        'contactNormalType':'Difference',
        'contactGapCorrection':'Implicit',
        'contactGapActivationRelativeTolerance':1.e-9,
        'useSurfacePositionForContact':1,
        'computeNodalArea':int(contact),
        'boundaryConditionTypes': [0,0,0,0,1,1] if kind=='thermoelastic'
                                  else [1,1,1,1,1,1],
        'outputType':'silo', 'plotGridFields':1,
        'gridFieldNames':['gridMass','gridVelocity'],
        'plottableFields':'{ particleTemperature, particleThermalEnergy, particleMass, particleVelocity, particleStress }',
        'particleFileFields':['Velocity','MaterialType','ContactGroup','Temperature','SurfaceFlag',
                              'Damage','SurfaceNormal','SurfacePosition','RVector'],
        'writeStatistics':'all',
        'boxAverageHistory':1,
        'boxAverageMin':[0,-Y/2-1.e-9,-dx/2-1.e-9],
        'boxAverageMax':[L,Y/2+1.e-9,dx/2+1.e-9],
        'boxAverageWriteInterval':stop/30,
        'profileHistory':1,
        'profileDirection':'x','profileNumSlices':nx-2,
        'profileVariables':'{ temperature, internalEnergy, volumeFraction }',
        'profileWriteInterval':stop/30,
    }
    if dfg:
        pfw['neighborRadius']=2.5*dx
        pfw['maxSingleFieldStateFractionForSeparability']=-1
        pfw['thinFeatureDFGThreshold']=10
    if kind=='plastic':
        pfw.update(prescribedBoundaryFTable=1,
                   fTableInterpType='Linear',
                   fTable=[[0.,1.,1.,1.],[stop,1.02,1.,1.]],
                   boundaryConditionTypes=[2,2,2,2,1,1])
    else:
        pfw['prescribedBoundaryFTable']=0

    if kind=='three_field':
        boxes = [
            _box('bodyA',0,L/2,-Y/2,0,0,0,T0+AMP,surface='right'),
            _box('bodyB',0,L/2,0,Y/2,1,1,T0,surface='right'),
            _box('bodyC',L/2,L,-Y/2,Y/2,2,2,T0+AMP/4,surface='left')]
    elif contact:
        gap=dx/2 if case_id=='contact_gap' else 0.
        mb=1 if case_id in ('contact_material','contact_group','contact_combined','contact_same_group') else 0
        gb=0 if (dfg or case_id=='contact_same_group') else 1
        boxes=[
            _box('left',0,L/2-gap,-Y/2,Y/2,0,0,T0+AMP,
                 surface='right',damage=1. if dfg else 0.),
            _box('right',L/2+gap,L,-Y/2,Y/2,mb,gb,T0,
                 surface='left',damage=0.)]
    elif kind in ('source_local','dirichlet_steady','dirichlet_transient','conductivity_linear'):
        if kind=='source_local':
            # ParticleRegion2 (material 1) is the only heated material. Material
            # identities differ but share a contact group and mechanical field.
            boxes=[_box('before_heater',0,.003,-Y/2,Y/2,0,0,T0),
                   _box('heater',.003,.004,-Y/2,Y/2,1,0,T0),
                   _box('after_heater',.004,L,-Y/2,Y/2,0,0,T0)]
        else:
            # Three material regions allow unambiguous thermostat selection.
            # The analytical reference treats strips of one cell as reservoirs.
            init=lambda ob,pt: temperature(kind,float(pt[0]),0.,case=case)
            boxes=[_box('left_bath',0,dx,-Y/2,Y/2,0,0,init),
                   _box('interior',dx,L-dx,-Y/2,Y/2,1,0,init),
                   _box('right_bath',L-dx,L,-Y/2,Y/2,2,0,init)]
    else:
        if kind in ('gaussian','cosine','step','energy_closed'):
            init=lambda ob,pt: temperature(kind,float(pt[0]),0.,case=case)
        else:
            init=T0
        boxes=[_box('body',0,L,-Y/2,Y/2,0,0,init)]
    pfw['objects']=boxes
    nmaterials=max(ob.mat for ob in boxes)+1
    pfw['materials']=['thermalSolid'+str(i) for i in range(nmaterials)]
    pfw['materialPropertyString']='\n'.join(
        _material_xml(pfw['materials'][i],case.get('conductivity',K),
                      expand=(kind=='thermoelastic'),plastic=(kind=='plastic'),
                      capacity_linear=(kind=='capacity_linear'),
                      conductivity_linear=(kind=='conductivity_linear'))
        for i in range(nmaterials))

    thermal={'enabled':1,'defaultContactMode':'Resistance',
             'defaultContactResistance':RC_GLOBAL,
             'gapTolerance':1.e-6,
             'safetyFactor':0.45,'maxSubsteps':10000,
             'areaQuadrature':32,
             'diagnostics':1,
             'contactResistance':[], 'heatSources':[], 'isothermal':[]}
    if case_id=='contact_perfect':
        thermal['defaultContactResistance']=0.
        thermal['defaultContactMode']='Perfect'
    if case_id in ('contact_material','contact_group','contact_combined','contact_same_group'):
        thermal['contactResistance'].append(
            {'name':'materialPair','materials':['thermalSolid0','thermalSolid1'],
             'resistance':RC_MATERIAL})
    if case_id in ('contact_group','contact_combined','contact_same_material'):
        thermal['contactResistance'].append(
            {'name':'groupPair','contactGroups':[0,1],'resistance':RC_GROUP})
    if case_id=='contact_combined':
        thermal['contactResistance'].append(
            {'name':'combinedPair','contactGroups':[0,1],
             'materials':['thermalSolid0','thermalSolid1'],
             'resistance':RC_COMBINED})
    if case_id=='contact_insulating':
        thermal['contactResistance'].append(
            {'name':'insulatePair','contactGroups':[0,1],'mode':'Insulated'})
    if kind=='three_field':
        thermal['contactResistance'].append(
            {'name':'insulateAB','contactGroups':[0,1],'mode':'Insulated'})
    if kind in ('source','sink','source_local','capacity_linear'):
        source={'name':'selectedPower',
                'specificPower':Q_SPECIFIC*(1 if kind!='sink' else -1)}
        if kind=='source_local':
            source['targetRegions']=['ParticleRegion2']
        thermal['heatSources'].append(source)
    if kind in ('dirichlet_steady','dirichlet_transient','conductivity_linear'):
        thermal['isothermal']=[
            {'name':'hotBath','targetRegions':['ParticleRegion1'],
             'temperature':T0+AMP},
            {'name':'coldBath','targetRegions':['ParticleRegion3'],
             'temperature':T0}]
    if kind=='thermoelastic':
        # T(t)=300+100*t/stop, with no other thermal loading.
        thermal['isothermal'].append(
            {'name':'uniformRamp', 'targetRegions':['ParticleRegion1'],
             'temperature':T0, 'temperatureRate':AMP/stop})
    pfw['thermal']=thermal

    if kind=='three_field':
        points=[[L/4,-Y/4,0],[L/4,Y/4,0],[3*L/4,0,0]]
    elif contact:
        points=[[L/4,0,0],[3*L/4,0,0],
                [L/2-2*dx,0,0],[L/2+2*dx,0,0]]
    else:
        points=[[L*(i+.5)/32,0,0] for i in range(32)]
    tracers.set_tracers(pfw,points,
        ['particleID','mass','temperature','internalEnergy','plasticStrainMagnitude',
         'stressXX','stressYY','stressZZ','stressYZ','stressXZ','stressXY',
         'plasticStrainXX','plasticStrainYY','plasticStrainZZ',
         'plasticStrainYZ','plasticStrainXZ','plasticStrainXY'],
        write_interval=stop/30,output_prefix='thermalVV')
    return pfw
