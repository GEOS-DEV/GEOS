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
    area = side_y * side_z
    inertia = side_y * side_z**3 / 12.0

    # GEOS results
    mode, frequency = np.loadtxt(args.outputDir + "/FreeFreeBeamFrequencies.txt", usecols=(0, 1), unpack=True)
    elastic = frequency[6:]
    x, uy1, uz1, uy2, uz2 = np.loadtxt(args.outputDir + "/FreeFreeBeamModeShapes.txt", unpack=True)

    # Euler-Bernoulli bending frequencies (each one is double: two polarizations) and the first axial frequency
    roots = euler_bernoulli_roots(6)
    bending = roots**2 / (2.0 * np.pi * length**2) * np.sqrt(young * inertia / (density * area))
    axial = np.sqrt(young / density) / (2.0 * length)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))

    # Frequencies: match each beam-theory frequency with the closest unused GEOS frequency
    used = np.zeros(elastic.size, dtype=bool)
    matched = []
    for f in list(bending) + [axial]:
        candidates = np.where(~used)[0]
        count = 2 if f in bending else 1
        for j in candidates[np.argsort(np.abs(elastic[candidates] - f))][:count]:
            used[j] = True
            matched.append((j, f))
    ax1.plot(np.arange(1, bending.size + 1), bending, 'k-o', mfc='none', label='Euler-Bernoulli bending')
    pairs = [(np.mean([elastic[j] for j, g in matched if g == f]), n + 1) for n, f in enumerate(bending)]
    ax1.plot([n for _, n in pairs], [f for f, _ in pairs], 'rx', ms=9, label='GEOS (mean of the pair)')
    ax1.set_xlabel('Bending mode order n')
    ax1.set_ylabel('Frequency (Hz)')
    ax1.set_yscale('log')
    ax1.grid(True, which='both', alpha=0.3)
    ax1.legend(loc='upper left')
    others = elastic[~used]
    ax1.set_title('Bending frequencies (torsion %.0f Hz and axial %.1f Hz are not shown)' % (others[0], axial), fontsize=9)

    # Mode shapes: the polarization in the (y, z) plane is arbitrary, so project on the beam-theory shape
    for (uy, uz), root, color, label in (((uy1, uz1), roots[0], 'C0', 'first bending'),
                                         ((uy2, uz2), roots[1], 'C1', 'second bending')):
        shape = euler_bernoulli_shape(x, root, length) / np.sqrt(density * area * length)
        direction = np.array([np.sum(uy * shape), np.sum(uz * shape)])
        direction /= np.linalg.norm(direction)
        transverse = direction[0] * uy + direction[1] * uz
        ax2.plot(x, shape, '-', color=color, label='Euler-Bernoulli, ' + label)
        ax2.plot(x[::4], transverse[::4], 'o', color=color, mfc='none', label='GEOS, ' + label)
    ax2.set_xlabel('x (m)')
    ax2.set_ylabel('Mass-normalized transverse displacement')
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=8)

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
