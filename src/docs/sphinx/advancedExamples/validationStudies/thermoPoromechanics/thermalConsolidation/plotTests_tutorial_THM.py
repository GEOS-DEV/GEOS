"""Figure de verification du cas de consolidation thermo-poro-elastique.

Execute par Sphinx (directive `.. plot::` de Example.rst) a chaque build de la documentation.
Il ne lit que deux fichiers CSV versionnes a cote de lui :

    thermalConsolidation_geos.csv         resultats GEOS
    thermalConsolidation_analytical.csv   solution analytique

Ces CSV sont produits hors ligne par postprocessTests_tutorial_THM.py, a partir des sorties
TimeHistory (.hdf5) d'un run du cas benchmark. Les .hdf5 ne sont volontairement pas versionnes
(plusieurs Mo) ; ce script n'a donc besoin ni de h5py ni de mpmath, seulement de numpy et
matplotlib, ce qui garde le build de la doc leger et reproductible.

Pour regenerer les donnees apres une modification du cas :
    geosx -i ThermoPoroElastic_consolidation_benchmark_fim.xml
    python3 postprocessTests_tutorial_THM.py -o <ce_dossier>
"""
import os

import numpy as np
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))

# sous-echantillonnage des marqueurs GEOS (1 = tous les instants collectes)
N1 = 4
fsize = 32
msize = 8
lw = 4
mew = 2
malpha = 0.6
lalpha = 0.8


def _load(name):
    path = os.path.join(HERE, name)
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"{name} introuvable dans {HERE}. Regenerez-le avec postprocessTests_tutorial_THM.py "
            f"(voir l'en-tete de ce fichier).")
    return np.genfromtxt(path, delimiter=",", names=True)


def main():
    geos = _load("thermalConsolidation_geos.csv")
    ana = _load("thermalConsolidation_analytical.csv")

    tg = geos["time"][::N1]
    ta = ana["time"]

    fig, ax = plt.subplots(2, 2, figsize=(32, 18))
    cmap = plt.get_cmap("tab10")

    # (axe, colonnes analytiques, colonnes GEOS, labels, titre y)
    panels = [
        (ax[0, 0], ["p_zh0p0", "p_zh0p6", "p_zh0p8"], ["p_0m", "p_4p2m", "p_5p6m"],
         ["z/h = 0.0", "z/h = 0.6", "z/h = 0.8"], ["z = 0.0 m", "z = 4.2 m", "z = 5.6 m"],
         "Pore Pressure [Pa]", "lower left"),
        (ax[0, 1], ["T_zh0p0", "T_zh0p6", "T_zh0p8"], ["T_0m", "T_4p2m", "T_5p6m"],
         ["z/h = 0.0", "z/h = 0.6", "z/h = 0.8"], ["z = 0.0 m", "z = 4.2 m", "z = 5.6 m"],
         "Temperature [K]", "upper left"),
        (ax[1, 0], ["disp_1p4m", "disp_4p2m", "disp_7m"], ["disp_1p4m", "disp_4p2m", "disp_7m"],
         ["z = 1.4 m", "z = 4.2 m", "z = 7.0 m"], ["z = 1.4 m", "z = 4.2 m", "z = 7.0 m"],
         "Displacement [m]", "upper left"),
        (ax[1, 1], ["sxx_0m", "sxx_4p2m", "sxx_5p6m"], ["sxx_0m", "sxx_4p2m", "sxx_5p6m"],
         ["$\\sigma_{xx}$: z = 0.0 m", "$\\sigma_{xx}$: z = 4.2 m", "$\\sigma_{xx}$: z = 5.6 m"],
         ["$\\sigma_{xx}$: z = 0.0 m", "$\\sigma_{xx}$: z = 4.2 m", "$\\sigma_{xx}$: z = 5.6 m"],
         "Total stress [Pa]", "upper left"),
    ]

    for a, acols, gcols, alabels, glabels, ylabel, legloc in panels:
        for k in range(3):
            a.semilogx(ta, ana[acols[k]], color=cmap(k), lw=lw, alpha=lalpha,
                       label="Analytical: " + alabels[k])
            a.plot(tg, geos[gcols[k]][::N1], 'o', color=cmap(k), markersize=msize, alpha=malpha,
                   mec=cmap(k), fillstyle='none', mew=mew, label="GEOS: " + glabels[k])
        a.set_xlabel('Time [s]', size=fsize, weight="bold")
        a.set_ylabel(ylabel, size=fsize, weight="bold")
        a.grid(True)
        a.xaxis.set_tick_params(labelsize=fsize)
        a.yaxis.set_tick_params(labelsize=fsize)

    # contrainte verticale totale : constante, egale a la charge appliquee
    ax[1, 1].semilogx(ta, ana["syy"], 'k--', lw=lw * 0.6, alpha=0.7,
                      label='$\\sigma_{yy}$ total = -F (all z)')
    ax[1, 1].plot(tg, geos["syy_0m"][::N1], 'k+', markersize=msize, mew=mew, alpha=malpha,
                  label='GEOS $\\sigma_{yy}$ total')

    for a, *_rest in panels:
        a.legend(loc=_rest[-1], fontsize=fsize * 0.5)

    plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=0.4, hspace=0.4)
    plt.show()


if __name__ == "__main__":
    main()
