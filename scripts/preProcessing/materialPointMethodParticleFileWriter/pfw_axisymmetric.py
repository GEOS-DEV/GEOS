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
    unsupported=('computeSPHJacobian','directionalOverlapCorrection','profileHistory',
                 'generalizedVortexMMS','useCrackTipDetection')
    for key in unsupported:
        if pfw.get(key,0):
            raise ValueError(key+' is not yet supported in Axisymmetric')
    for key in ('enableBoreholePressure','enableConfiningPressure'):
        if pfw.get(key,0):
            raise ValueError(key+' is internal solver state; activate pressure through an MPM event with useEvents=1')
    if pfw.get('overlapCorrection','off') not in (0,'off','Off',None):
        raise ValueError('Axisymmetric requires overlapCorrection=off')
    # Inspect event XML as well as solver attributes: the event activates the
    # pressure internally, so enableConfiningPressure is not a PFW control.
    events = ET.fromstring('<root>' + str(pfw.get('mpmEventsString', '') or '') + '</root>')
    validate_pressure_events(pfw, events)
    for event in events.iter('ReferenceCohesiveZones'):
        if pfw.get('useEvents') not in (1, '1'):
            raise ValueError('Axisymmetric ReferenceCohesiveZones requires useEvents=1; otherwise the interface never activates')
    regions = pfw.get('cohesiveZoneRegions') or ''
    if isinstance(regions, (list, tuple)):
        regions = '\n'.join(regions)
    regions = ET.fromstring('<root>' + regions + '</root>')
    materials = ET.fromstring('<root>' + (pfw.get('materialPropertyString') or '') + '</root>')
    material_types = {item.get('name'): item.tag for item in materials.iter() if item.get('name')}
    for region in regions.iter('CohesiveZoneRegion'):
        model = material_types.get(region.get('constitutiveModel'))
        if model is not None and model not in ('UncoupledCohesiveZone', 'CoupledCohesiveZone', 'BicrystalCohesiveZone'):
            raise ValueError('Axisymmetric CZ supports traction-separation laws; finite-thickness polymer CZ needs hoop-stretch kinematics')
        missing = set(('ContactGroup', 'SurfaceFlag', 'SurfaceNormal', 'SurfacePosition', 'CZTag')) - set(pfw.get('particleFileFields', []))
        if missing:
            raise ValueError('Axisymmetric CZ requires particleFileFields: ' + ', '.join(sorted(missing)))
    if pfw.get('enableSurfaceTension'):
        pairs = _control_array(pfw, 'surfaceTensionPairs', [])
        if pairs.ndim != 2 or pairs.shape[1] != 3 or len(pairs) == 0:
            raise ValueError('surfaceTensionPairs requires rows [phaseA, phaseB, gamma]; phaseB=-1 is void')
        for a,b,gamma in pairs:
            if not (0 <= a < 100 and -1 <= b < 100 and a == int(a) and b == int(b) and a != b and gamma >= 0):
                raise ValueError('surfaceTensionPairs requires distinct integer groups, optional void -1, and nonnegative gamma')
        for key,default in (('surfaceTensionGradientTolerance',0),('surfaceTensionForceSign',1)):
            value = pfw.get(key)
            value = default if value is None else float(value)
            if not math.isfinite(value) or (key.endswith('Tolerance') and value < 0):
                raise ValueError(key + ' must be finite (tolerance also nonnegative)')
    validate_driving_options(pfw)


def validate_pressure_events(pfw, events):
    """Pressure selectors act on current grid coordinates in the (r,y) plane.

    Corners retain the existing three-value XML format; z is ignored by the
    axisymmetric selector. Confinement takes precedence in overlapping regions.
    """
    for event in list(events.iter('ConfiningPressure')) + list(events.iter('BoreholePressure')):
        if pfw.get('useEvents') not in (1, '1'):
            raise ValueError('Axisymmetric pressure events require useEvents=1')
        def number(key, default=None):
            raw = event.get(key, default)
            if raw is None:
                raise ValueError(event.tag + ' requires ' + key)
            value = float(raw)
            if not math.isfinite(value):
                raise ValueError(event.tag + ' requires finite ' + key)
            return value
        number('startPressure')
        number('endPressure')
        if number('interpType', 1) not in (0,1,2):
            raise ValueError('Pressure interpType must be 0, 1 or 2')
        if event.tag == 'BoreholePressure':
            if number('boreholeRadius') <= 0:
                raise ValueError('BoreholePressure requires positive boreholeRadius')
            continue
        uniform = number('uniformPressure', 0)
        sphere = number('confiningPressureSphereRadius', 0)
        if uniform not in (0,1) or sphere < 0:
            raise ValueError('uniformPressure must be 0/1 and sphere radius nonnegative')
        corners = []
        for key in ('confiningPressureBoxMin', 'confiningPressureBoxMax'):
            raw = event.get(key)
            if raw is None:
                corners.append(None)
                continue
            values = [float(value) for value in raw.strip().strip('{}').replace(',', ' ').split()]
            if len(values) != 3 or not all(math.isfinite(value) for value in values):
                raise ValueError(key + ' requires three finite values; Axisymmetric uses r,y and ignores z')
            corners.append(values)
        if uniform and (sphere > 0 or any(value is not None for value in corners)):
            raise ValueError('uniformPressure=1 must omit box and sphere exclusions')
        if not uniform and sphere <= 0 and any(value is None for value in corners):
            raise ValueError('ConfiningPressure requires uniform mode, a positive sphere radius, or both box corners')
        if all(value is not None for value in corners) and any(lo > hi for lo,hi in zip(*corners)):
            raise ValueError('Confining pressure box corners must be ordered (min <= max)')


def _control_array(pfw, key, default):
    """Accept Python/NumPy arrays and the existing GEOS-brace table strings."""
    import ast
    value = pfw.get(key)
    if value is None:
        value = default
    if isinstance(value, str):
        value = ast.literal_eval(value.replace('{', '[').replace('}', ']'))
    result = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(result)):
        raise ValueError(key + ' must contain finite numbers')
    return result


def validate_driving_options(pfw):
    """Wall driving is diagonal in (r,y), with a fixed revolution axis.

    The third grid stretch/rate remains 1/0; it is not particle hoop F/L.
    Stress control targets full-ring volume averages of sigma_rr and sigma_yy.
    """
    if pfw.get('prescribedFTable'):
        raise ValueError('Axisymmetric uses prescribedBoundaryFTable for wall motion; prescribedFTable periodic affine driving is not supported')
    control = _control_array(pfw, 'stressControl', [0,0,0])
    if control.shape != (3,) or np.any(~np.isin(control, [0,1,2])) or control[2] != 0:
        raise ValueError('Axisymmetric stressControl requires [radial, axial, 0] with entries 0, 1 or 2')
    if pfw.get('stressControlAdaptiveSeek') and np.any(control == 2):
        raise ValueError('stressControlAdaptiveSeek supports only control entries 0 or 1')
    stretch = _control_array(pfw, 'domainF', [1,1,1])
    rate = _control_array(pfw, 'domainL', [0,0,0])
    if stretch.shape != (3,) or rate.shape != (3,) or np.any(stretch <= 0) or stretch[2] != 1 or rate[2] != 0:
        raise ValueError('Axisymmetric grid controls require positive domainF, domainF[2]=1, domainL[2]=0')
    driven = (control[:2] > 0) | (rate[:2] != 0)
    if pfw.get('prescribedBoundaryFTable'):
        table = _control_array(pfw, 'fTable', [])
        if table.ndim != 2 or table.shape[0] == 0 or table.shape[1] != 4:
            raise ValueError('fTable requires rows [time, Frr, Fyy, 1]')
        if table[0,0] > 0 or np.any(np.diff(table[:,0]) <= 0) or np.any(table[:,3] != 1):
            raise ValueError('fTable must start at zero or earlier, increase in time and keep its last column 1')
        for d in range(2):
            if control[d] != 1 and (table[0,d+1] != 1 or np.any(table[:,d+1] <= 0)):
                raise ValueError('Prescribed F-table stretches must start at 1 and remain positive')
        driven |= np.any(table[:,1:3] != 1, axis=0)
    if np.any(control):
        table = _control_array(pfw, 'stressTable', [])
        if table.ndim != 2 or table.shape[0] == 0 or table.shape[1] != 4:
            raise ValueError('stressTable requires rows [time, sigma_rr, sigma_yy, unused]')
        if table[0,0] != 0 or np.any(np.diff(table[:,0]) <= 0):
            raise ValueError('stressTable must start at zero with strictly increasing times')
    if driven[1] and pfw.get('periodic', [False,False])[1]:
        raise ValueError('Axisymmetric axial wall driving requires nonperiodic axial boundaries')

    lower = [float(pfw.get('xmin',0)), float(pfw.get('ymin',0))]
    upper = [float(pfw['xmax']), float(pfw.get('ymax',1))]
    has_axis = lower[0] == 0
    types = _control_array(pfw, 'boundaryConditionTypes', [0,0,0,0])
    if types.shape not in ((4,), (6,)) or np.any(~np.isin(types,[0,1,2,3])):
        raise ValueError('Axisymmetric boundary types require 4 (or padded 6) entries from 0,1,2,3')
    if has_axis:
        types[0] = 1
    if len(types) == 6:
        types[4:] = 0
    pfw['boundaryConditionTypes'] = types.astype(int).tolist()
    schedules = [(types[:4], 0., float("inf"))]
    if pfw.get('prescribedBcTable'):
        table = _control_array(pfw, 'bcTable', [])
        if table.ndim != 2 or table.shape[0] == 0 or table.shape[1] not in (5,7):
            raise ValueError('bcTable requires time and four active face types')
        if table[0,0] != 0 or np.any(np.diff(table[:,0]) <= 0) or np.any(~np.isin(table[:,1:5],[0,1,2,3])):
            raise ValueError('bcTable must start at zero, increase in time and use face types 0,1,2,3')
        if has_axis and np.any(table[:,1] != 1):
            raise ValueError('The r=0 face must remain Symmetry (1) in every BC-table row')
        schedules = [(row[1:5], row[0], table[i+1,0] if i+1<len(table) else float('inf'))
                     for i,row in enumerate(table)]
    for types, begin, end in schedules:
        row_driven = (control[:2] > 0) | ((rate[:2] != 0) if not pfw.get('prescribedBoundaryFTable') else False)
        if pfw.get('prescribedBoundaryFTable'):
            f_table = _control_array(pfw, 'fTable', [])
            for i in range(1,len(f_table)):
                if f_table[i,0] > begin and f_table[i-1,0] < end:
                    row_driven |= f_table[i,1:3] != f_table[i-1,1:3]
        if pfw.get('updateMethod') == 'FMPM' and np.any(types == 3):
            raise ValueError('FMPM contact walls require an incremental law; use Moving walls or FLIP/PIC/XPIC')
        for d in range(2):
            if not row_driven[d]:
                continue
            if np.all(types[2*d:2*d+2] < 2):
                raise ValueError('Each driven direction requires a Moving or Contact face')
            if (types[2*d] == 1 and lower[d] != 0) or (types[2*d+1] == 1 and upper[d] != 0):
                raise ValueError('A driven direction cannot have a fixed Symmetry face away from the origin')
    enabled = _control_array(pfw, 'enablePrescribedBoundaryTransverseVelocities', [0]*6)
    velocities = _control_array(pfw, 'prescribedBoundaryTransverseVelocities', [[0,0]]*6)
    if enabled.shape not in ((4,), (6,)) or velocities.shape != (6,2):
        raise ValueError('Transverse wall controls require four or six flags and six velocity pairs')
    for face in range(4):
        if not enabled[face]:
            continue
        if velocities[face, 1 if face < 2 else 0] != 0:
            raise ValueError('Axisymmetric boundary velocities cannot prescribe swirl')
        if has_axis and face >= 2 and velocities[face,1] != 0:
            raise ValueError('A constant nonzero radial cap velocity is not regular at r=0; prescribe radial domain stretch instead')




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


def reconstructs_surfaces(pfw):
    """Global solver/event request consumed before cohesive initialization."""
    events = ET.fromstring('<root>' + str(pfw.get('mpmEventsString', '') or '') + '</root>')
    regions = pfw.get('cohesiveZoneRegions') or ''
    if isinstance(regions,(list,tuple)): regions = '\n'.join(regions)
    regions = ET.fromstring('<root>' + regions + '</root>')
    return any(int(region.get('computeParticleSurfaceNormalsAndPositions','0')) == 1
               for region in regions.iter('CohesiveZoneRegion')) or pfw.get('computeParticleSurfaceNormalsAndPositions') in (1,'1') or any(
        int(event.get('computeNormalsAndPositions','0')) == 1
        for event in events.iter('ReferenceCohesiveZones'))


def validate_particle(point, data, geometry, reconstruct_surfaces=False):
    if data['particleType'] not in (0,1,2,3):
        raise ValueError('Axisymmetric supports Q1 point, cubic B-spline, CPDI and CPTI particles')
    # Ordinary body contact uses the existing ContactGroup fields. Geometry
    # stays meridional; mass and momentum already represent complete rings.
    group = data.get('group', 0)
    if not (0 <= group < 100 and int(group) == group):
        raise ValueError('Axisymmetric ContactGroup must be an integer from 0 to 99')
    tag = data.get('czTag', 0)
    if not (math.isfinite(tag) and tag >= 0 and int(tag) == tag):
        raise ValueError('Axisymmetric CZTag must be a nonnegative integer')
    if data['particleType'] == 1 and (tag > 0 or data.get('surfaceFlag',0) in (3,4)):
        raise ValueError('Axisymmetric B-spline CZ particles need a B-spline cohesive area/reference operator; use Q1, CPDI or CPTI for CZ-tagged particles')
    if data.get('surfaceFlag',0) in (3,4):
        if group not in (0, 1):
            raise ValueError('Axisymmetric CZ requires contact groups 0 and 1')
        normal = np.asarray(data.get('surfaceNormal', [0,0,0]), dtype=float)
        position = np.asarray(data.get('surfacePosition', [0,0,0]), dtype=float)
        if (not np.all(np.isfinite(normal)) or not np.all(np.isfinite(position)) or
                (not reconstruct_surfaces and np.linalg.norm(normal[:2]) <= 0) or abs(normal[2]) > 1.e-12 or abs(position[2]) > 1.e-12):
            raise ValueError('Axisymmetric CZ requires finite meridional surface data and a nonzero normal unless reconstruction is requested')
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
