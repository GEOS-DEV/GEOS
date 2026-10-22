#!/usr/bin/env python3
"""Check generated geometry, ring measures, body labels and CZ options."""
import json
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np


def check(folder):
    case = json.loads((folder/'case.json').read_text())
    path = folder/'mpmParticleFile_axisymmetricCZ'
    with path.open() as stream:
        names = stream.readline().split()
    values = np.loadtxt(path, skiprows=1, ndmin=2)
    assert np.all(np.isfinite(values))
    def col(name):
        return values[:, names.index(name)]
    r, y, volume = col('PositionX'), col('PositionY'), col('ParticleVolume')
    groups, shapes = col('ContactGroup'), col('GeometryShape')
    g = values[:, [names.index('Geometry'+ij) for ij in ('XX', 'XY', 'YX', 'YY')]].reshape(-1, 2, 2)
    assert set(groups) == {0, 1} and set(shapes) == {3, 4}
    assert np.all(r > 0) and np.all(col('PositionZ') == 0)
    assert np.allclose(volume, 2*np.pi*r*np.linalg.det(g)*np.where(shapes == 4, .5, 1.), rtol=1.e-12, atol=0)
    assert np.all(col('VelocityX') == 0) and np.all(col('VelocityZ') == 0)
    assert np.all(col('SurfaceNormalZ') == 0)
    assert set(col('SurfaceFlag')).issubset({0, 3})
    expected = np.pi*case['radius']**2*case['length']
    surface_counts = {}
    for group in (0, 1):
        rows = groups == group
        assert set(shapes[rows]) == {3, 4}
        assert np.isclose(volume[rows].sum(), expected, rtol=1.e-12, atol=0)
        assert np.all(col('VelocityY')[rows] == (-1 if group == 0 else 1)*case['speed'])
        if group == 0:
            assert np.all(y[rows] < 0)
        else:
            assert np.all(y[rows] > 0)
        surface = rows & (col('SurfaceFlag') == 3)
        assert np.any(surface)
        surface_counts[str(group)] = int(np.count_nonzero(surface))
        assert np.all(col('SurfaceNormalX')[surface] == 0)
        assert np.all(col('SurfaceNormalY')[surface] == (1 if group == 0 else -1))
        assert np.allclose(y[surface]+col('SurfacePositionY')[surface],
                           (-1 if group == 0 else 1)*0, rtol=1.e-12, atol=1.e-15)
    assert np.all(col('CZTag') == 1)
    assert abs(np.sum(case['density']*volume*col('VelocityY'))) < 1.e-14
    xml = ET.parse(folder/'mpm_axisymmetricCZ.xml')
    solver = xml.find('.//SolidMechanics_MPM')
    for name, value in dict(formulation='Axisymmetric', useEvents='1', enableContact='1', damageFieldPartitioning='0',
                            contactNormalType='MassWeighted', contactGapCorrection='Implicit',
                            useSurfacePositionForContact='1', updateMethod=case['update']).items():
        assert solver.get(name) == value, (name, solver.get(name))
    assert solver.get('overlapCorrection') == 'Off', 'GEOS enums are case sensitive'
    assert float(solver.get('frictionCoefficient')) == 0
    assert xml.find('.//InternalMesh').get('axisymmetric') == '1'
    region=xml.find('.//CohesiveZoneRegion')
    assert region.get('constitutiveModel')=='interfaceLaw' and region.get('tag')=='1'
    event=xml.find('.//ReferenceCohesiveZones')
    assert event.get('czSurfaceDisplacementUpdate')==case['displacement']
    assert event.get('computeNormalsAndPositions')=='0'
    law=xml.find('.//UncoupledCohesiveZone')
    assert float(law.get('normalForceConstant'))==case['stiffness']
    assert float(law.get('maxNormalDisplacement'))==case['failure_opening']
    report = dict(particles=len(values), groups=2, volume_per_body=float(expected),
                  mass_per_body=float(case['density']*expected), update=case['update'],
                  useEvents=int(solver.get('useEvents')),
                  cz_surface_particles_per_group=surface_counts)
    (folder/'initialChecks.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    check(Path('.'))
