"""Validate emitted contact groups, damage, geometry and solver controls."""
import json
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np


def check(folder):
    case = json.loads((folder/'case.json').read_text())
    path = folder/'mpmParticleFile_axisymmetricDFG'
    names = path.open().readline().split()
    data = np.loadtxt(path, skiprows=1, ndmin=2)
    assert np.all(np.isfinite(data))
    def col(name):
        return data[:, names.index(name)]
    r, y, volume, groups = (col(k) for k in ('PositionX', 'PositionY', 'ParticleVolume', 'ContactGroup'))
    assert set(groups) == set(range(case['groups']))
    assert np.all(r > 0) and np.all(col('PositionZ') == 0)
    assert np.all(col('VelocityZ') == 0) and np.all(col('SurfaceFlag') == 0)
    assert np.all(col('SurfaceNormalX') == 0) and np.all(col('SurfaceNormalY') == 0)
    shape = col('GeometryShape')
    g = data[:, [names.index('Geometry'+ij) for ij in ('XX', 'XY', 'YX', 'YY')]].reshape(-1, 2, 2)
    assert np.allclose(volume, 2*np.pi*r*np.linalg.det(g)*np.where(shape == 4, .5, 1.), rtol=1e-12, atol=0)
    for group in range(case['groups']):
        assert set(shape[groups == group]) == {3, 4}
        length = .012 if case['groups'] == 1 else (.004 if group == 1 else .005)
        assert np.isclose(volume[groups == group].sum(), np.pi*case['radius']**2*length, rtol=1e-12)
    expected = (np.abs(y) < .00051) if case['mode'] in ('damage', 'combined') else np.zeros(len(y), dtype=bool)
    assert np.array_equal(col('Damage'), expected.astype(float))
    if case['mode'] in ('damage', 'combined'):
        assert np.any(expected & (y < 0)) and np.any(expected & (y > 0))
    assert abs(np.sum(volume*col('VelocityY'))) < 1e-14
    xml = ET.parse(folder/'mpm_axisymmetricDFG.xml')
    solver = xml.find('.//SolidMechanics_MPM')
    assert solver.get('formulation') == 'Axisymmetric'
    assert solver.get('damageFieldPartitioning') == str(case['dfg'])
    assert solver.get('enableContact') == '1' and solver.get('overlapCorrection') == 'Off'
    assert xml.find('.//CohesiveZoneRegion') is None
    report = dict(mode=case['mode'], particles=len(data), contact_groups=case['groups'],
                  velocity_fields=case['groups']*(2 if case['dfg'] else 1),
                  damaged_particles=int(expected.sum()), volume=float(volume.sum()))
    (folder/'initialChecks.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    check(Path('.'))
