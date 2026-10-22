#!/usr/bin/env python3
"""Check the actual generated particle/XML files and plot their r-y geometry."""
import argparse
import json
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np


def check(folder):
    case = json.loads((folder/'case.json').read_text())
    particle_file = folder/'mpmParticleFile_axisymmetricJacketedSample'
    with particle_file.open() as stream:
        names = stream.readline().split()
    data = np.loadtxt(particle_file, skiprows=1, ndmin=2)
    if not np.all(np.isfinite(data)):
        raise ValueError('Nonfinite particle data')
    col = lambda name: data[:, names.index(name)]
    centers = np.column_stack((col('PositionX'), col('PositionY')))
    shape = col('GeometryShape').astype(int)
    material = col('MaterialType').astype(int)
    volumes = col('ParticleVolume')
    edges = data[:, [names.index('Geometry'+ij) for ij in ('XX','XY','YX','YY')]].reshape(-1,2,2)
    assert set(shape) == {3, 4} and set(material) == {0, 1, 2}
    assert np.all(col('ContactGroup') == 0) and np.all(centers[:, 0] > 0)
    assert np.all(col('PositionZ') == 0)
    area = np.linalg.det(edges)*np.where(shape == 4, 0.5, 1.)
    np.testing.assert_allclose(volumes, 2*np.pi*centers[:, 0]*area, rtol=1.e-12, atol=0)
    assert np.all(area > 0)
    assert np.all(shape[material != 2] == 3) and np.all(shape[material == 2] == 4)
    radius = case['radius']; half = case['sample_half_height']
    plate = case['plate_thickness']; jacket = case['jacket_thickness']
    outer_radius = radius+jacket; outer_half = half+plate+jacket
    expected = np.array([np.pi*radius**2*2*half, np.pi*radius**2*2*plate,
                         np.pi*outer_radius**2*2*outer_half - np.pi*radius**2*2*(half+plate)])
    totals = np.array([volumes[material == i].sum() for i in range(3)])
    np.testing.assert_allclose(totals, expected, rtol=1.e-12, atol=0)
    # In addition to totals, check every material region and every domain vertex.
    tol = 1.e-12
    r, y = centers.T
    assert np.all((r[material == 0] < radius) & (abs(y[material == 0]) < half))
    assert np.all((r[material == 1] < radius) & (abs(y[material == 1]) > half) &
                  (abs(y[material == 1]) < half+plate))
    assert np.all((r[material == 2] > radius) | (abs(y[material == 2]) > half+plate))
    polygons = []
    for center, g, code in zip(centers, edges, shape):
        if code == 3:
            offsets = np.array([-g[0]-g[1], g[0]-g[1], g[0]+g[1], -g[0]+g[1]])/2
        else:
            offsets = np.array([-(g[0]+g[1]), 2*g[0]-g[1], 2*g[1]-g[0]])/3
        vertices = center+offsets
        assert np.all(vertices[:, 0] >= -tol) and np.all(vertices[:, 0] <= outer_radius+tol)
        assert np.all(abs(vertices[:, 1]) <= outer_half+tol)
        polygons.append(vertices)
    xml = ET.parse(folder/'mpm_axisymmetricJacketedSample.xml').getroot()
    solver = next(xml.iter('SolidMechanics_MPM'))
    mesh = next(xml.iter('InternalMesh'))
    event = next(xml.iter('ConfiningPressure'))
    assert solver.get('formulation') == 'Axisymmetric' and solver.get('useEvents') == '1'
    assert solver.get('timeIntegrationOption') == 'ExplicitDynamic' and solver.get('updateMethod') == 'PIC'
    assert mesh.get('axisymmetric') == '1' and mesh.get('spatialDimension') == '2'
    assert event.get('uniformPressure') == '1' and event.get('interpType') == '1'
    assert float(event.get('endPressure')) == case['pressure']
    assert float(event.get('endTime')) == case['ramp_time']
    material_xml = {node.get('name'): node for node in xml.iter('ElasticIsotropic')}
    for m in case['materials']:
        node = material_xml[m['name']]
        assert float(node.get('defaultDensity')) == m['density']
        assert float(node.get('defaultYoungModulus')) == m['E']
        assert float(node.get('defaultPoissonRatio')) == m['nu']
    rows = [dict(name=m['name'], particles=int(sum(material == i)), volume_m3=float(totals[i]),
                 mass_kg=float(totals[i]*m['density'])) for i, m in enumerate(case['materials'])]
    report = dict(variant=case['variant'], particles=len(data), components=rows,
                  total_volume_m3=float(totals.sum()), total_mass_kg=sum(row['mass_kg'] for row in rows))
    (folder/'initialChecks.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))
    draw(folder, polygons, material, case)
    return report


def draw(folder, polygons, material, case):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    from matplotlib.patches import Patch
    colors = ['#196482', '#888888', '#e8782d']
    fig, ax = plt.subplots(figsize=(6.3, 7.5), layout='constrained')
    collection = PolyCollection([1000*p for p in polygons], facecolors=[colors[i] for i in material],
                                edgecolors='#ffffff', linewidths=0.15)
    ax.add_collection(collection)
    ax.axvline(0., color='#cc3030', linestyle='--', linewidth=1.7)
    outer_r = 1000*(case['radius']+case['jacket_thickness'])
    outer_y = 1000*(case['sample_half_height']+case['plate_thickness']+case['jacket_thickness'])
    for y in np.linspace(-0.75*outer_y, 0.75*outer_y, 5):
        ax.annotate('', xy=(outer_r+0.2, y), xytext=(outer_r+2., y),
                    arrowprops=dict(arrowstyle='->', color='#23528c', lw=1.5))
    for r in (2., 5.5, 9.):
        for sign in (-1, 1):
            ax.annotate('', xy=(r, sign*(outer_y+0.2)), xytext=(r, sign*(outer_y+1.7)),
                        arrowprops=dict(arrowstyle='->', color='#23528c', lw=1.5))
    ax.set(xlim=(-0.8, outer_r+3.), ylim=(-outer_y-2., outer_y+2.), aspect='equal',
           xlabel='r [mm]', ylabel='axial position y [mm]',
           title='Axisymmetric jacketed sample\nActual initial CPDI/CPTI domains')
    labels = ['Sample (effective porous elastic)', 'Steel plates', 'Rubber jacket']
    if case['variant'] == 'homogeneous':
        labels = [label+' [same material]' for label in ('Core', 'Plates', 'Jacket')]
    ax.legend(handles=[Patch(color=c, label=label) for c, label in zip(colors, labels)],
              loc='upper center', bbox_to_anchor=(0.5, -0.10), frameon=False, fontsize=9)
    fig.savefig(folder/'initialGeometry.png', dpi=170)
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', nargs='?', type=Path, default=Path('.'))
    check(parser.parse_args().directory)
