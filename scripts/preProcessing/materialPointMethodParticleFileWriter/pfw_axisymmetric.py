"""Full-ring measures and input checks for no-swirl native axisymmetry."""
import math
import xml.etree.ElementTree as ET
import numpy as np


def ring_volume(point, widths):
    radius=float(point[0])
    if not math.isfinite(radius) or radius<=0:
        raise ValueError("Axisymmetric particles require finite positive radius")
    return 2*math.pi*radius*float(widths[0])*float(widths[1])


def validate_options(pfw):
    if not 0 <= float(pfw.get('xmin',0)) < float(pfw['xmax']):
        raise ValueError("Axisymmetric physical bounds require 0 <= xmin < xmax")
    if pfw.get('periodic',[False,False])[0]:
        raise ValueError("Axisymmetric radial periodicity is invalid")
    if pfw.get('planeStrainThickness') is not None:
        raise ValueError("Axisymmetric uses full rings; omit planeStrainThickness")
    unsupported=('damageFieldPartitioning','enableSurfaceTension','computeSPHJacobian','directionalOverlapCorrection',
                 'computeParticleSurfaceNormalsAndPositions','computeNodalArea','profileHistory',
                 'generalizedVortexMMS','enableBoreholePressure','enableConfiningPressure',
                 'useCrackTipDetection','prescribedBcTable','prescribedFTable','prescribedBoundaryFTable')
    for key in unsupported:
        if pfw.get(key,0):
            raise ValueError(key+' is not yet supported in Axisymmetric')
    if pfw.get('overlapCorrection','off') not in (0,'off','Off',None):
        raise ValueError('Axisymmetric requires overlapCorrection=off')
    if any(pfw.get('stressControl') or [0,0,0]):
        raise ValueError('Axisymmetric stress-controlled domain motion is not implemented')
    # Inspect event XML as well as solver attributes: the event activates the
    # pressure internally, so enableConfiningPressure is not a PFW control.
    events = ET.fromstring('<root>' + str(pfw.get('mpmEventsString', '') or '') + '</root>')
    for event in events.iter('ConfiningPressure'):
        if event.get('uniformPressure') != '1':
            raise ValueError('Axisymmetric ConfiningPressure requires uniformPressure="1"')
        if any(key in event.attrib for key in ('confiningPressureBoxMin', 'confiningPressureBoxMax',
                                               'confiningPressureSphereRadius')):
            raise ValueError('Uniform axisymmetric pressure must omit box and sphere exclusions')
    for event in events.iter('ReferenceCohesiveZones'):
        if pfw.get('useEvents') not in (1, '1'):
            raise ValueError('Axisymmetric ReferenceCohesiveZones requires useEvents=1; otherwise the interface never activates')
        if int(event.get('computeNormalsAndPositions', '0')):
            raise ValueError('Axisymmetric CZ requires explicit normals and positions')
    regions = pfw.get('cohesiveZoneRegions') or ''
    if isinstance(regions, (list, tuple)):
        regions = '\n'.join(regions)
    regions = ET.fromstring('<root>' + regions + '</root>')
    materials = ET.fromstring('<root>' + (pfw.get('materialPropertyString') or '') + '</root>')
    material_types = {item.get('name'): item.tag for item in materials.iter() if item.get('name')}
    for region in regions.iter('CohesiveZoneRegion'):
        if int(region.get('computeParticleSurfaceNormalsAndPositions', '0')):
            raise ValueError('Axisymmetric CZ requires explicit normals and positions')
        model = material_types.get(region.get('constitutiveModel'))
        if model is not None and model not in ('UncoupledCohesiveZone', 'CoupledCohesiveZone', 'BicrystalCohesiveZone'):
            raise ValueError('Axisymmetric CZ supports traction-separation laws; finite-thickness polymer CZ needs hoop-stretch kinematics')
        missing = set(('ContactGroup', 'SurfaceFlag', 'SurfaceNormal', 'SurfacePosition', 'CZTag')) - set(pfw.get('particleFileFields', []))
        if missing:
            raise ValueError('Axisymmetric CZ requires particleFileFields: ' + ', '.join(sorted(missing)))
    bc=pfw.get('boundaryConditionTypes')
    if bc is not None and any(value not in (0,1) for value in bc):
        raise ValueError('Axisymmetric supports Outflow and Symmetry boundaries')
    bc=list(bc or [0,0,0,0])
    if float(pfw.get('xmin',0))==0:
        bc[0]=1
    pfw['boundaryConditionTypes']=bc


def bind_objects(objects):
    """Do not pass a fictitious constant thickness to stochastic volume laws.

    A custom wrapper can implement setAxisymmetric() and integrate 2*pi*r*dA
    over each statistical cell. Explicit user StrengthScale data are untouched.
    Existing built-in area-only Weibull states lack radial first moments.
    """
    seen=set()
    def visit(obj):
        if obj is None or id(obj) in seen: return
        seen.add(id(obj))
        if isinstance(obj,(list,tuple)):
            for child in obj: visit(child)
            return
        fields=vars(obj) if hasattr(obj,'__dict__') else {}
        hook=getattr(type(obj),'setAxisymmetric',None)
        if hook is not None: hook(obj)
        elif '_nativeWeibull' in fields:
            raise ValueError('This Weibull object stores areas only; Axisymmetric needs ring-integrated statistical volumes or a setAxisymmetric hook')
        for name in ('subObject','subObjects','objects','object','objA','objB','subObjA','subObjB'):
            if name in fields: visit(fields[name])
    visit(objects)


def validate_particle(point, data, geometry):
    if data['particleType'] not in (0,2,3):
        raise ValueError('Axisymmetric supports Q1 point, CPDI and CPTI particles')
    # Ordinary body contact uses the existing ContactGroup fields. Geometry
    # stays meridional; mass and momentum already represent complete rings.
    group = data.get('group', 0)
    if not (0 <= group < 100 and int(group) == group):
        raise ValueError('Axisymmetric ContactGroup must be an integer from 0 to 99')
    tag = data.get('czTag', 0)
    if not (math.isfinite(tag) and tag >= 0 and int(tag) == tag):
        raise ValueError('Axisymmetric CZTag must be a nonnegative integer')
    if data.get('surfaceFlag',0) in (3,4):
        if group not in (0, 1):
            raise ValueError('Axisymmetric CZ requires contact groups 0 and 1')
        normal = np.asarray(data.get('surfaceNormal', [0,0,0]), dtype=float)
        position = np.asarray(data.get('surfacePosition', [0,0,0]), dtype=float)
        if (not np.all(np.isfinite(normal)) or not np.all(np.isfinite(position)) or
                np.linalg.norm(normal[:2]) <= 0 or abs(normal[2]) > 1.e-12 or abs(position[2]) > 1.e-12):
            raise ValueError('Axisymmetric CZ requires a nonzero meridional normal and meridional surface position')
    if np.any(np.asarray(data.get('surfaceTraction',[0,0,0]))):
        raise ValueError('Axisymmetric PFW surface-traction integration is not implemented; use body forces or initial velocity/stress')
    if abs(float(data.get('velocity',[0,0,0])[2]))>0:
        raise ValueError('Axisymmetric has no circumferential velocity (no swirl)')
    shape,g=geometry
    if shape==3:
        inward=0.5*(abs(g[0,0])+abs(g[1,0]))
    else:
        offsets=np.vstack((-(g[0]+g[1])/3,(2*g[0]-g[1])/3,(2*g[1]-g[0])/3))
        inward=-offsets[:,0].min()
    if point[0]-inward < -1e-12*max(1.0,abs(point[0])):
        raise ValueError('Initial axisymmetric material geometry crosses r=0')
