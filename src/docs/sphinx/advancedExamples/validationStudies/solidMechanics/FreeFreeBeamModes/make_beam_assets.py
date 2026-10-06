"""Generate the static assets of the free-free beam example from a GEOS run.

This script is not run by the documentation build. Run it again if the results change.

1. Run GEOS on modalFreeFreeBeam_arnoldi_smoke.xml with modalNumModes="20" and a VTK output,
   for example by adding to modalFreeFreeBeam_base.xml:
       <Outputs><VTK name="vtkOutput" plotLevel="1"/> ... </Outputs>
       <PeriodicEvent name="vtkSnapshot" timeFrequency="1.0" target="/Outputs/vtkOutput"/>
2. python3 make_beam_assets.py --vtu <run>/vtkOutput/000001/beam/Level0/Region/rank_0.vtu \
                               --log <run>/run.log --geosDir <GEOS repository>

It writes FreeFreeBeamModeTable.csv, FreeFreeBeamModeShapes.txt, FreeFreeBeamMesh.png and FreeFreeBeamModes.png.
"""
import argparse
import re
import xml.etree.ElementTree as ElementTree

import numpy as np

from FreeFreeBeamModes_vs_EulerBernoulli import classify, euler_bernoulli_roots, timoshenko_estimate


def read_log(path):
    rows = []
    for line in open(path):
        m = re.match(r'^\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)', line)
        if m and ('e+' in m.group(2) or 'e-' in m.group(2)):
            rows.append((int(m.group(1)), float(m.group(2))))
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--vtu', required=True)
    parser.add_argument('--log', required=True)
    parser.add_argument('--geosDir', default='../../../../../../..')
    parser.add_argument('--outputDir', default='.')
    args = parser.parse_args()

    base = ElementTree.parse(args.geosDir + "/inputFiles/solidMechanics/modalFreeFreeBeam_base.xml")
    mesh_xml = base.find('Mesh/InternalMesh')
    steel = base.find('Constitutive/ElasticIsotropic')
    length = float(mesh_xml.get("xCoords").strip("{} ").split(",")[1])
    side = float(mesh_xml.get("yCoords").strip("{} ").split(",")[1])
    young = float(steel.get("defaultYoungModulus"))
    density = float(steel.get("defaultDensity"))
    poisson = float(steel.get("defaultPoissonRatio"))
    shear = young / (2.0 * (1.0 + poisson))

    # Theory: Euler-Bernoulli bending, Saint-Venant torsion of a square section (J = 0.1406 a^4, Ip = a^4 / 6), axial rod
    roots = euler_bernoulli_roots(8)
    bending = roots**2 / (2.0 * np.pi * length**2) * np.sqrt(young * (side**4 / 12.0) / (density * side**2))
    torsion = np.sqrt(shear * 0.1406 / (density * (1.0 / 6.0))) / (2.0 * length)
    axial = np.sqrt(young / density) / (2.0 * length)

    timoshenko = timoshenko_estimate(bending, roots, length, side / np.sqrt(12.0), young, shear, poisson)
    log = read_log(args.log)
    frequencies = [f for _, f in log]
    labels, theory = classify(frequencies[6:], bending, torsion, axial)

    with open(args.outputDir + '/FreeFreeBeamModeTable.csv', 'w') as out:
        out.write('Mode,Type,GEOS (Hz),Euler-Bernoulli / torsion / rod (Hz),Difference (%),Timoshenko estimate (Hz),Difference to Timoshenko (%)\n')
        for k in range(6):
            out.write(f'{k + 1},rigid body,{frequencies[k]:.2e},0,,,\n')
        for k, (label, f) in enumerate(zip(labels, frequencies[6:])):
            line = f'{k + 7},{label},{f:.4f},{theory[k]:.4f},{100.0 * (f / theory[k] - 1.0):+.1f}'
            if label.startswith('bending'):
                t = timoshenko[int(label.split()[1]) - 1]
                line += f',{t:.4f},{100.0 * (f / t - 1.0):+.1f}'
            else:
                line += ',,'
            out.write(line + '\n')

    # Mode shapes of the first four bending modes along the axis
    import pyvista as pv
    grid = pv.read(args.vtu)
    points = grid.points
    on_axis = np.where((np.abs(points[:, 1] - 0.5 * side) < 1e-6) & (np.abs(points[:, 2] - 0.5 * side) < 1e-6))[0]
    on_axis = on_axis[np.argsort(points[on_axis, 0])]
    modes = (7, 9, 11, 13)
    with open(args.outputDir + '/FreeFreeBeamModeShapes.txt', 'w') as out:
        out.write('# Mass-normalized transverse displacement of the nodes on the axis of the beam\n')
        out.write('# column 1 = x (m)\n')
        out.write('# columns 2 to 9 = y and z displacement of the modes 7, 9, 11 and 13 (first four bending modes)\n')
        for i in on_axis:
            row = [points[i, 0]] + [grid.point_data[f'modeShape{m}'][i, c] for m in modes for c in (1, 2)]
            out.write(' '.join(f'{v:.9e}' for v in row) + '\n')

    # Images
    pv.OFF_SCREEN = True
    outline = grid.outline()

    pl = pv.Plotter(shape=(2, 1), off_screen=True, window_size=(1500, 760))
    pl.subplot(0, 0)
    pl.add_mesh(grid, color='lightsteelblue', show_edges=True, line_width=0.6)
    pl.view_xy()
    pl.camera.parallel_projection = True
    pl.camera.focal_point = (0.5 * length, 0.5 * side, 0.0)
    pl.camera.position = (0.5 * length, 0.5 * side, 30.0)
    pl.camera.parallel_scale = 1.4
    pl.add_text('Side view of the mesh: 100 x 4 x 4 hexahedra, 10 m x 0.2 m x 0.2 m', font_size=15)
    pl.subplot(1, 0)
    pl.add_mesh(grid, color='lightsteelblue', show_edges=True, line_width=1.0)
    pl.camera.parallel_projection = False
    pl.camera.focal_point = (0.9, 0.1, 0.1)
    pl.camera.position = (-0.9, -0.8, 0.7)
    pl.camera.view_angle = 40.0
    pl.add_text('End of the beam: four hexahedra across the thickness', font_size=15)
    pl.screenshot(args.outputDir + '/FreeFreeBeamMesh.png')

    nrows = 7
    peak_display_strain = 0.02
    pl = pv.Plotter(shape=(nrows, 2), off_screen=True, window_size=(1500, 1250))
    for k in range(14):
        m = k + 7
        u = grid.point_data[f'modeShape{m}'].copy()
        # The polarization of a bending mode in the plane of the section is arbitrary: rotate it onto the y axis
        transverse = u[:, 1:3]
        _, _, vt = np.linalg.svd(transverse, full_matrices=False)
        axis = vt[0]
        rotated = np.stack([u[:, 0], transverse @ axis, transverse @ np.array([-axis[1], axis[0]])], axis=1)
        # The modes have no physical amplitude (they are mass-normalized). For the image, the displacement is
        # amplified so that the peak strain of every mode is the same (1 %), instead of the peak displacement.
        derivative = grid.compute_derivative(scalars=f'modeShape{m}')['gradient'].reshape(-1, 3, 3)
        strain = 0.5 * (derivative + np.transpose(derivative, (0, 2, 1)))
        peak_strain = np.percentile(np.max(np.abs(np.linalg.eigvalsh(strain)), axis=1), 99.0)
        scale = peak_display_strain / peak_strain
        shape = grid.copy()
        shape.points = grid.points + scale * rotated
        magnitude_raw = np.linalg.norm(rotated, axis=1)
        magnitude = magnitude_raw
        shape['displacement'] = magnitude / np.max(magnitude)
        pl.subplot(k % nrows, k // nrows)
        pl.add_mesh(shape, scalars='displacement', cmap='viridis', clim=(0.0, 1.0), show_scalar_bar=False)
        pl.add_mesh(outline, color='gray', line_width=1)
        pl.view_xy()
        pl.camera.parallel_projection = True
        pl.camera.focal_point = (0.5 * length, 0.0, 0.1)
        pl.camera.position = (0.5 * length, 0.0, 30.0)
        pl.camera.parallel_scale = 1.25
        pl.add_text(f'Mode {m}: {labels[k]}, {frequencies[6 + k]:.1f} Hz, peak {scale * np.max(magnitude_raw):.2f} m',
                    font_size=15, position='upper_left')
    gallery = pl.screenshot(return_img=True)

    # One colorbar for all the panels: the displacement of each mode is normalized by its maximum
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    height, width = gallery.shape[:2]
    dpi = 100
    fig = plt.figure(figsize=((width + 200) / dpi, (height + 30) / dpi), dpi=dpi)
    ax = fig.add_axes([0.0, 30.0 / (height + 30), width / (width + 200), height / (height + 30)])
    ax.imshow(gallery)
    ax.axis('off')
    cax = fig.add_axes([(width + 30) / (width + 200), 0.12, 24.0 / (width + 200), 0.76])
    fig.colorbar(matplotlib.cm.ScalarMappable(norm=matplotlib.colors.Normalize(0.0, 1.0), cmap='viridis'), cax=cax)
    cax.set_ylabel('Displacement magnitude / its maximum', fontsize=15)
    cax.tick_params(labelsize=13)
    fig.text(0.5, 0.006, 'The amplitude of a mode is arbitrary. Each mode is scaled so that its peak strain is 2 % (steel yields at about 0.1 %)', ha='center', fontsize=13)
    fig.savefig(args.outputDir + '/FreeFreeBeamModes.png', dpi=dpi)


if __name__ == "__main__":
    main()
