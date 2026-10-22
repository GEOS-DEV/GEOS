"""Plot a testRelperm_3_phase_*.txt table.

These tables hold (index, Sg, So, Sw, krg, kro, krw) rows, sampled along a
handful of fixed-Sw "legs" with Sg/So varying along each leg. This script
produces view of that data:

A filled contour of the sampled (Sw, Sg) points on a ternary diagram,
colored by one of krg/kro/krw (default kro), interpolated across the
sampled legs via Delaunay triangulation

Usage:
    python plot_relperm_table.py [path/to/testRelperm_3_phase_*.txt] [krg|kro|krw]

With no arguments, defaults to the example StoneII-model table and kro.
"""

import sys
import os
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

KR_LABELS = {'krg': '$k_{rg}$', 'kro': '$k_{ro}$', 'krw': '$k_{rw}$'}

def add_minor_lines(ax, z=0.0):
    """Add the dotted 10% gridlines to a ternary-diagram axes.

    `ax` may be a regular 2D Axes or an Axes3D - in the 3D case the lines
    are drawn in the z=`z` plane.
    """
    is_3d = hasattr(ax, "zaxis")

    for i in np.arange(0, 1.0 + 1e-9, 0.1):
        for Sw1, Sg1, Sw2, Sg2 in (
            (i, 1 - i, i, 0),
            (1 - i, i, 0, i),
            (i, 0, 0, i),
        ):
            x1, y1 = ternary_coord(Sw1, Sg1)
            x2, y2 = ternary_coord(Sw2, Sg2)
            if is_3d:
                ax.plot([x1, x2], [y1, y2], zs=z, zdir='z', linestyle=':', color='k', linewidth=1)
            else:
                ax.plot([x1, x2], [y1, y2], ':k', linewidth=1)


def ternaryfigure(num=None):
    """Create and set up a blank ternary-diagram figure.

    Opens/clears figure `num`, draws the triangle outline, the minor
    gridlines and the Sw / So / Sg vertex labels. 

    Returns (fig, ax).
    """
    fig = plt.figure()
    fig.clf()

    ax = fig.add_subplot()
    ax.set_aspect('equal')

    ax.set_axis_off()

    triangle_x = [0, 1 / np.sqrt(3), 2 / np.sqrt(3), 0]
    triangle_y = [0, 1, 0, 0]

    ax.plot(triangle_x, triangle_y, '-k', linewidth=2)
    ax.text(1 / np.sqrt(3) - 0.03, 1 + 0.03, '$S_w$')
    ax.text(2 / np.sqrt(3) + 0.01, -0.03, '$S_o$')
    ax.text(-0.05, -0.03, '$S_g$')

    add_minor_lines(ax)


    return fig, ax
def ternary_coord(Sw, Sg):
    """Convert (Sw, Sg) compositions to ternary-diagram cartesian coordinates.

    Sw, Sg may be scalars or arrays of the same shape (e.g. a meshgrid).
    The third component So = 1 - Sw - Sg is implicit. Returns (x, y) with
    the same shape as the inputs - y is just Sw, x places So at the right
    vertex (2/sqrt(3), 0), Sw at the top vertex (1/sqrt(3), 1) and Sg at the
    origin.
    """
    Sw = np.asarray(Sw, dtype=float)
    Sg = np.asarray(Sg, dtype=float)

    x = (1 - Sg) / np.cos(np.pi / 6) - Sw / np.tan(np.pi / 3)
    y = Sw

    return x, y


def load_relperm_table(path):
    """Load a testRelperm_3_phase_*.txt table.

    These files hold whitespace-separated columns (index, Sg, So, Sw, krg,
    kro, krw), with '#'-prefixed header comments. Returns a dict of 1D
    arrays keyed by column name.
    """
    data = np.loadtxt(path)
    columns = ['index', 'Sg', 'So', 'Sw', 'krg', 'kro', 'krw']
    return {name: data[:, i] for i, name in enumerate(columns)}

def plot_ternary_surface(table, value='kro'):
    """Filled contour of `value` (one of 'krg', 'kro', 'krw') over the
    sampled (Sw, Sg) points on a ternary diagram, interpolated across the
    sampled legs via Delaunay triangulation."""
    
    plt.rcParams.update({'font.size': 14})

    x, y = ternary_coord(table['Sw'], table['Sg'])

    fig, ax = ternaryfigure()
    cs = ax.tricontourf(x, y, table[value], levels=20, cmap='turbo')
    fig.colorbar(cs, ax=ax, label=KR_LABELS[value])

    return fig, ax


if __name__ == "__main__":
    data_directory = os.path.normpath(os.path.abspath("."))
    default_path = (Path(data_directory)
                    / "testRelperm_3_phase_example_3_phase_stoneii.txt")
    table_path = sys.argv[1] if len(sys.argv) > 1 else default_path
    value = sys.argv[2] if len(sys.argv) > 2 else 'kro'

    table = load_relperm_table(table_path)
    plot_ternary_surface(table, value=value)

    plt.show()
