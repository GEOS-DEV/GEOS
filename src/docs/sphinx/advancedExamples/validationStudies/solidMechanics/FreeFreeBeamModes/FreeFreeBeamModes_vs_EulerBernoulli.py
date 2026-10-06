import argparse
import xml.etree.ElementTree as ElementTree

import matplotlib.pyplot as plt
import numpy as np


def euler_bernoulli_shape(x, beta_length, length):
    """Free-free Euler-Bernoulli mode shape, normalized to a mean square of one over the beam."""
    sigma = (np.cosh(beta_length) - np.cos(beta_length)) / (np.sinh(beta_length) - np.sin(beta_length))
    b = beta_length / length
    return np.cosh(b * x) + np.cos(b * x) - sigma * (np.sinh(b * x) + np.sin(b * x))


def euler_bernoulli_roots(count):
    """Roots of cos(x) cosh(x) = 1, found by bisection near (n + 1/2) pi."""
    roots = []
    for n in range(1, count + 1):
        lo, hi = (n + 0.5) * np.pi - 0.2, (n + 0.5) * np.pi + 0.2
        if n == 1:
            lo, hi = 4.5, 5.0
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            if (np.cos(lo) * np.cosh(lo) - 1.0) * (np.cos(mid) * np.cosh(mid) - 1.0) <= 0.0:
                hi = mid
            else:
                lo = mid
        roots.append(0.5 * (lo + hi))
    return np.array(roots)


def timoshenko_estimate(frequency, beta_length, length, gyration, young, shear, poisson):
    """Rayleigh-Timoshenko estimate of a bending frequency: shear deformation and rotary inertia lower the
    Euler-Bernoulli frequency by 1 / sqrt(1 + (beta r)^2 (1 + E / (kappa G))). It is first-order accurate for a
    free-free beam. kappa is the shear coefficient of a rectangular section (Cowper)."""
    kappa = 10.0 * (1.0 + poisson) / (12.0 + 11.0 * poisson)
    return frequency / np.sqrt(1.0 + (beta_length / length * gyration)**2 * (1.0 + young / (kappa * shear)))


def classify(frequencies, bending, torsion, axial):
    """Label the elastic modes and give their beam theory frequency.

    Pairs of equal frequencies are bending modes (the section is square). The single modes are the torsion mode and
    the axial mode, in the order of their frequencies.
    """
    labels, theory, singles = [], [], []
    i, n = 0, 0
    while i < len(frequencies):
        if i + 1 < len(frequencies) and abs(frequencies[i + 1] - frequencies[i]) < 1e-5 * frequencies[i]:
            for polarization in 'ab':
                labels.append(f'bending {n + 1} ({polarization})')
                theory.append(bending[n])
            n += 1
            i += 2
        else:
            singles.append(len(labels))
            labels.append(None)
            theory.append(None)
            i += 1
    for k, j in enumerate(singles):
        labels[j] = ['torsion 1', 'axial 1'][k]
        theory[j] = [torsion, axial][k]
    return labels, theory


def main():

    parser = argparse.ArgumentParser(description="Script to generate the figure of the free-free beam example.")
    parser.add_argument('--geosDir', help='Path to the GEOS repository', default='../../../../../../..')
    parser.add_argument('--outputDir', help='Path to the directory of the results', default='.')
    args = parser.parse_args()

    # Beam parameters, read from the GEOS input files
    base = ElementTree.parse(args.geosDir + "/inputFiles/solidMechanics/modalFreeFreeBeam_base.xml")
    mesh = base.find('Mesh/InternalMesh')
    steel = base.find('Constitutive/ElasticIsotropic')
    length = float(mesh.get("xCoords").strip("{} ").split(",")[1])
    side_y = float(mesh.get("yCoords").strip("{} ").split(",")[1])
    side_z = float(mesh.get("zCoords").strip("{} ").split(",")[1])
    young = float(steel.get("defaultYoungModulus"))
    density = float(steel.get("defaultDensity"))
    poisson = float(steel.get("defaultPoissonRatio"))
    area = side_y * side_z
    inertia = side_y * side_z**3 / 12.0
    shear = young / (2.0 * (1.0 + poisson))

    # GEOS results
    frequency = np.loadtxt(args.outputDir + "/FreeFreeBeamFrequencies.txt", usecols=(1,))
    elastic = frequency[6:]
    shapes = np.loadtxt(args.outputDir + "/FreeFreeBeamModeShapes.txt")
    x = shapes[:, 0]

    # Beam theory: Euler-Bernoulli bending (each frequency is double), Saint-Venant torsion of a square section
    # (J = 0.1406 a^4, polar moment a^4 / 6) and the axial frequency of a free-free rod
    roots = euler_bernoulli_roots(8)
    bending = roots**2 / (2.0 * np.pi * length**2) * np.sqrt(young * inertia / (density * area))
    torsion = np.sqrt(shear * 0.1406 / (density / 6.0 * 1.0)) / (2.0 * length)
    axial = np.sqrt(young / density) / (2.0 * length)
    labels, theory = classify(elastic, bending, torsion, axial)
    timoshenko = timoshenko_estimate(bending, roots, length, np.sqrt(inertia / area), young, shear, poisson)

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 4.2))

    # Frequencies: one position for each bending order, then torsion and axial
    names = ['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'T1', 'A1']
    theory_values = list(bending[:6]) + [torsion, axial]
    geos_values = []
    for n in range(6):
        geos_values.append(np.mean([f for f, label in zip(elastic, labels) if label.startswith(f'bending {n + 1} ')]))
    geos_values.append(elastic[labels.index('torsion 1')])
    geos_values.append(elastic[labels.index('axial 1')])
    position = np.arange(len(names))
    ax1.semilogy(position, theory_values, 'ko', mfc='none', ms=9, label='Euler-Bernoulli, torsion, rod')
    ax1.semilogy(position[:6], timoshenko[:6], 'b^', mfc='none', ms=7, label='Timoshenko estimate')
    ax1.semilogy(position, geos_values, 'rx', ms=9, label='GEOS')
    for p, g, t in zip(position, geos_values, theory_values):
        ax1.annotate('%+.1f %%' % (100.0 * (g / t - 1.0)), (p, g), textcoords='offset points', xytext=(0, 9),
                     ha='center', fontsize=7)
    ax1.margins(y=0.2)
    ax1.set_xticks(position)
    ax1.set_xticklabels(names)
    ax1.set_xlabel('Bending order (B), torsion (T), axial (A)')
    ax1.set_ylabel('Frequency (Hz)')
    ax1.set_title('Frequencies of the elastic modes', fontsize=10)
    ax1.grid(True, which='both', alpha=0.3)
    ax1.legend(loc='lower right')

    # Mode shapes: the polarization in the (y, z) plane is arbitrary, so project on the beam theory shape
    for axis, orders in ((ax2, (0, 1)), (ax3, (2, 3))):
        for order in orders:
            uy, uz = shapes[:, 1 + 2 * order], shapes[:, 2 + 2 * order]
            shape = euler_bernoulli_shape(x, roots[order], length) / np.sqrt(density * area * length)
            direction = np.array([np.sum(uy * shape), np.sum(uz * shape)])
            direction /= np.linalg.norm(direction)
            transverse = direction[0] * uy + direction[1] * uz
            color = 'C%d' % order
            axis.plot(x, shape, '-', color=color, label='Beam theory, bending %d' % (order + 1))
            axis.plot(x[::4], transverse[::4], 'o', color=color, mfc='none', label='GEOS, bending %d' % (order + 1))
        axis.set_xlabel('x (m)')
        axis.set_ylabel('Mass-normalized transverse displacement')
        axis.grid(True, alpha=0.3)
        axis.legend(fontsize=7)
    ax2.set_title('Mode shapes along the axis', fontsize=10)

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
