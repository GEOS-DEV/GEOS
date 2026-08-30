"""Figure of the 1D tutorial: effect of cooling on failure.

Executed by Sphinx (the `.. plot::` directive of Example.rst) on every documentation build.
It reads a single versioned CSV file sitting next to it:

    cooling1D.csv

That CSV is produced offline by postprocess1DCooling.py from the TimeHistory (.hdf5) outputs
of both runs. The .hdf5 files are not versioned; this script therefore depends only on numpy
and matplotlib.

To regenerate the data after a change in the cases:
    geosx -i ThermoElastic_1DCooling_fim_smoke.xml        # in elastic/
    geosx -i ThermoDruckerPrager_1DCooling_fim_smoke.xml  # in druckerPrager/
    python3 postprocess1DCooling.py -e elastic -d druckerPrager -o <this_directory>
"""
import os

import numpy as np
import matplotlib.pyplot as plt

try:
    HERE = os.path.dirname(os.path.abspath(__file__))
except NameError:
    # the Sphinx `.. plot::` directive may execute this file without defining __file__
    HERE = os.getcwd()

fsize = 20
msize = 7
lw = 3
mew = 1.8
malpha = 0.6
lalpha = 0.8
N1 = 2          # one GEOS marker out of two, to keep the plot readable


def main():
    path = os.path.join(HERE, "cooling1D.csv")
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"cooling1D.csv not found in {HERE}. Regenerate it with postprocess1DCooling.py "
            f"(see the header of this file).")
    d = np.genfromtxt(path, delimiter=",", names=True)

    dT = d["dT"]
    sy = float(d["sigma_yield"][0])
    b = float(d["dp_b"][0])
    c = float(d["dp_c"][0])
    # cooling at which the yield surface is reached
    dT_yield = float(np.interp(sy, d["sigma_elastic_analytical"], dT))

    cmap = plt.get_cmap("tab10")
    fig, ax = plt.subplots(1, 3, figsize=(25, 7))

    # ---------------------------------------------------------------- panel 1
    # stress induced by the confined cooling, as a function of -dT
    a = ax[0]
    x = -dT
    a.plot(x, d["sigma_elastic_analytical"] / 1e3, '-', color=cmap(0), lw=lw, alpha=lalpha,
           label=r'Analytical elastic: $\sigma_{yy} = -E\,\alpha\,\Delta T$')
    a.plot(x[::N1], d["sigma_elastic_geos"][::N1] / 1e3, 'o', color=cmap(0), markersize=msize,
           fillstyle='none', mew=mew, alpha=malpha, label='GEOS: thermo-elastic')
    a.plot(x[::N1], d["sigma_druckerPrager_geos"][::N1] / 1e3, 's', color=cmap(1),
           markersize=msize, fillstyle='none', mew=mew, alpha=malpha,
           label='GEOS: thermo-plastic (Drucker-Prager)')
    a.axhline(sy / 1e3, color=cmap(1), ls='--', lw=lw * 0.7, alpha=0.9,
              label=r'Failure cap: $\sigma_{f} = c\,/\,(1 + b/3)$')
    a.axvline(-dT_yield, color='0.4', ls=':', lw=lw * 0.7)
    a.annotate(f'failure at $\\Delta T$ = {dT_yield:.1f} K',
               xy=(-dT_yield, sy / 1e3), xytext=(-dT_yield + 4, sy / 1e3 * 0.45),
               fontsize=fsize * 0.7, color='0.25',
               arrowprops=dict(arrowstyle='->', color='0.45', lw=1.4))
    a.set_xlabel(r'Cooling  $-\Delta T$  [K]', size=fsize, weight='bold')
    a.set_ylabel(r'Total stress  $\sigma_{yy}$  [kPa]', size=fsize, weight='bold')
    a.set_title('Thermally induced stress under confined cooling', size=fsize * 0.85)
    a.grid(True, alpha=0.35)
    a.legend(loc='upper left', fontsize=fsize * 0.62)
    a.tick_params(labelsize=fsize * 0.8)

    # ---------------------------------------------------------------- panel 2
    # loading path in the (P, Q) plane, together with the Drucker-Prager envelope
    a = ax[1]
    Pe, Qe = d["sigma_elastic_geos"] / 3.0, d["sigma_elastic_geos"]
    Pd, Qd = d["sigma_druckerPrager_geos"] / 3.0, d["sigma_druckerPrager_geos"]
    Pmax = max(Pe.max(), c / b) * 1.05
    Penv = np.linspace(0.0, c / b, 200)
    a.plot(Penv / 1e3, (c - b * Penv) / 1e3, 'k--', lw=lw * 0.8, alpha=lalpha,
           label=r'Drucker-Prager envelope: $Q = c - b\,P$')
    a.plot(Pe / 1e3, Qe / 1e3, '-', color=cmap(0), lw=lw, alpha=lalpha,
           label='thermo-elastic path')
    a.plot(Pd[::N1] / 1e3, Qd[::N1] / 1e3, 's', color=cmap(1), markersize=msize,
           fillstyle='none', mew=mew, alpha=malpha, label='thermo-plastic path')
    a.plot([sy / 3.0 / 1e3], [sy / 1e3], '*', color='k', markersize=msize * 2.4,
           label='failure point')
    a.set_xlim(0.0, Pmax / 1e3)
    a.set_xlabel(r'Mean stress  $P = \mathrm{tr}(\sigma)/3$  [kPa]', size=fsize, weight='bold')
    a.set_ylabel(r'Von Mises stress  $Q$  [kPa]', size=fsize, weight='bold')
    a.set_title('Stress path in the invariant plane', size=fsize * 0.85)
    a.grid(True, alpha=0.35)
    a.legend(loc='upper right', fontsize=fsize * 0.62)
    a.tick_params(labelsize=fsize * 0.8)

    # ---------------------------------------------------------------- panel 3
    # lateral displacement of the free face: yielding adds to it
    a = ax[2]
    a.plot(x, d["ux_elastic_analytical"] * 1e6, '-', color=cmap(0), lw=lw, alpha=lalpha,
           label='Analytical: elastic')
    a.plot(x[::N1], d["ux_elastic_geos"][::N1] * 1e6, 'o', color=cmap(0), markersize=msize,
           fillstyle='none', mew=mew, alpha=malpha, label='GEOS: thermo-elastic')
    a.plot(x, d["ux_druckerPrager_analytical"] * 1e6, '-', color=cmap(1), lw=lw, alpha=lalpha,
           label='Analytical: Drucker-Prager')
    a.plot(x[::N1], d["ux_druckerPrager_geos"][::N1] * 1e6, 's', color=cmap(1), markersize=msize,
           fillstyle='none', mew=mew, alpha=malpha, label='GEOS: thermo-plastic')
    a.axvline(-dT_yield, color='0.4', ls=':', lw=lw * 0.7)
    a.set_xlabel(r'Cooling  $-\Delta T$  [K]', size=fsize, weight='bold')
    a.set_ylabel(r'Lateral displacement  $u_x$  [$\mu$m]', size=fsize, weight='bold')
    a.set_title('Lateral displacement of the free face', size=fsize * 0.85)
    a.grid(True, alpha=0.35)
    a.legend(loc='lower left', fontsize=fsize * 0.62)
    a.tick_params(labelsize=fsize * 0.8)

    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
