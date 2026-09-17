"""Post-traitement LOCAL du cas de consolidation thermo-poro-elastique.

A executer a la main, dans le repertoire ou GEOS a tourne (celui qui contient les .hdf5) :

    python3 postprocessTests_tutorial_THM.py [-o <dossier_de_sortie>]

Il lit les sorties TimeHistory de GEOS et calcule la solution analytique, puis ecrit deux
fichiers CSV legers :

    thermalConsolidation_geos.csv         resultats GEOS (une ligne par instant collecte)
    thermalConsolidation_analytical.csv   solution analytique (100 instants log-espaces)

Ce sont ces deux CSV qui sont versionnes dans le depot. Les .hdf5 (plusieurs Mo) ne le sont
pas : la figure de la documentation est produite par plotTests_tutorial_THM.py, qui ne lit
que ces CSV et n'a besoin ni de h5py ni de mpmath.

Entrees attendues (sorties du cas benchmark) :
    pressureHistory.hdf5  temperatureHistory.hdf5  displacementHistory.hdf5  stressHistory.hdf5
La sortie stress suppose que le XML active le Task 'stressCollection' et l'Output
'stressHistoryOutput' (voir ThermoPoroElastic_consolidation_base.xml).
"""
import argparse
import os

import numpy as np
import h5py
from mpmath import nsum, exp, sin, cos, inf
import math


# --------------------------------------------------------------------------------------
# Solution analytique (identique a la version historique du script de tracé)
# --------------------------------------------------------------------------------------
def computePP(x, t, h, phi12, phi11, r11, r22):
    p = nsum(lambda m: 1 / ((2 * m + 1) * math.pi / 2 / h) *
             (phi12 * exp(-r11 * ((2 * m + 1) * math.pi / 2 / h)**2 * t) +
              phi11 * exp(-r22 * ((2 * m + 1) * math.pi / 2 / h)**2 * t)) *
             sin(((2 * m + 1) * math.pi / 2 / h) * (h - x)), [0, inf])
    return 2 / h * p


def computeTemp(x, t, h, phi22, phi21, r11, r22):
    p = nsum(lambda m: 1 / ((2 * m + 1) * math.pi / 2 / h) *
             (phi22 * exp(-r11 * ((2 * m + 1) * math.pi / 2 / h)**2 * t) +
              phi21 * exp(-r22 * ((2 * m + 1) * math.pi / 2 / h)**2 * t)) *
             sin(((2 * m + 1) * math.pi / 2 / h) * (h - x)), [0, inf])
    return 2 / h * p


def computeDisp(x, t, h, t1, t2, r11, r22, A, M):
    p = nsum(lambda m: 1 / ((2 * m + 1) * math.pi / 2 / h)**2 *
             (t1 * exp(-r11 * ((2 * m + 1) * math.pi / 2 / h)**2 * t) +
              t2 * exp(-r22 * ((2 * m + 1) * math.pi / 2 / h)**2 * t)) *
             cos(((2 * m + 1) * math.pi / 2 / h) * (h - x)), [0, inf])
    return -A * x + 2 / h / M * p


def analyticalParameters():
    """Parametres du probleme et coefficients de la solution analytique."""
    E = 6.0e3
    nu = 0.4
    Ks = 1.0e27
    a_s = 3.0e-7
    phi = 0.2

    lamda = E * nu / (1.0 + nu) / (1.0 - 2.0 * nu)
    G = E / 2.0 / (1.0 + nu)
    K = E / 3.0 / (1.0 - 2.0 * nu)
    alpha = 1.0 - K / Ks
    M = lamda + 2.0 * G

    a_m = 9.0e-7
    m = 1.672e5
    beta = a_m * (3.0 * lamda + 2.0 * G) / 3.0
    perm = 4.0e-9
    mu = 1.0e-3
    km = perm / mu
    h = 7.0
    theta_0 = 0.0
    theta_a = 50.0
    F = 1.0
    Kt = 836

    nuu = 0.499
    M_biot = 2.0 * G * (nuu - nu) / (alpha**2 * (1.0 - 2.0 * nuu) * (1.0 - 2.0 * nu))
    a_p = 1.0 / M_biot
    p0 = alpha * M_biot / (lamda + 2.0 * G + alpha**2 * M_biot) * F
    A = (F - alpha * 0 - beta * theta_a) / M

    g11 = 1.0 / km * (a_p + alpha**2 / M)
    g12 = 1.0 / km * (alpha * beta / M - a_m)
    g21 = alpha * beta * theta_0 / (Kt * M)
    g22 = (m + beta**2 * theta_0 / M) / Kt
    gs = g11 * g22 - g12 * g21
    r11, r12 = g11 / gs, g12 / gs
    r21, r22 = g21 / gs, g22 / gs
    r1s = r12 * (theta_0 - theta_a) - r22 * p0
    r2s = r21 * p0 - r11 * (theta_0 - theta_a)
    rs = r11 * r22 - r12 * r21

    phi12 = -(r12 * r2s) / rs
    phi11 = -(r11 * r1s) / rs
    phi22 = -(r22 * r2s) / rs
    phi21 = -(r21 * r1s) / rs
    t1 = alpha * phi12 + beta * phi22
    t2 = alpha * phi11 + beta * phi21

    return dict(h=h, M=M, G=G, lamda=lamda, alpha=alpha, beta=beta, F=F, theta_a=theta_a, A=A,
                r11=r11, r22=r22, phi11=phi11, phi12=phi12, phi21=phi21, phi22=phi22,
                t1=t1, t2=t2)


# --------------------------------------------------------------------------------------
# Lecture des sorties GEOS
# --------------------------------------------------------------------------------------
def _valid_length(t):
    """Nombre d'echantillons reellement ecrits : le tampon HDF5 peut etre pre-alloue et
    complete par des zeros, auquel cas le temps repasse a 0."""
    for j in range(1, len(t)):
        if t[j] < 1e-12:
            return j
    return len(t)


def _pick(coords, targets):
    """Premier indice dont la cote y atteint chaque cible (meme regle que le script d'origine)."""
    out = []
    for c in targets:
        idx = next((j for j in range(len(coords)) if coords[j] >= c), len(coords) - 1)
        out.append(idx)
    return out


def readGeos(biot=1.0):
    for f in ("pressureHistory.hdf5", "temperatureHistory.hdf5",
              "displacementHistory.hdf5", "stressHistory.hdf5"):
        if not os.path.exists(f):
            raise FileNotFoundError(
                f"{f} introuvable. Lancez le cas benchmark (ThermoPoroElastic_consolidation_"
                f"benchmark_fim.xml) et executez ce script dans le repertoire de sortie.")

    hp = h5py.File("pressureHistory.hdf5", 'r')
    t = np.array(hp.get('pressure Time')).ravel()
    n = _valid_length(t)
    t = t[:n]
    pc = np.array(hp.get('pressure elementCenter'))[0, :, 1]
    P = np.array(hp.get('pressure'))[:n]
    ip = _pick(pc, (0.0, 4.2, 5.6))

    ht = h5py.File("temperatureHistory.hdf5", 'r')
    tc = np.array(ht.get('temperature elementCenter'))[0, :, 1]
    T = np.array(ht.get('temperature'))[:n]
    it = _pick(tc, (0.0, 4.2, 5.6))

    hd = h5py.File("displacementHistory.hdf5", 'r')
    dc = np.array(hd.get('totalDisplacement ReferencePosition'))[0, :, 1]
    U = np.array(hd.get('totalDisplacement'))[:n]
    iu = _pick(dc, (1.4, 4.2, 7.0))

    hs = h5py.File("stressHistory.hdf5", 'r')
    sc = np.array(hs.get('averageStress elementCenter'))[0, :, 1]
    S = np.array(hs.get('averageStress'))[:n]
    isg = [int(np.argmin(np.abs(sc - d))) for d in (0.0, 4.2, 5.6)]

    # contrainte TOTALE = averageStress - biot * pression, prises dans la MEME cellule
    # (maillages identiques -> l'indice de cellule stress est valable pour la pression)
    sxx = [0.5 * (S[:, i, 0] + S[:, i, 2]) - biot * P[:, i] for i in isg]
    syy = S[:, isg[0], 1] - biot * P[:, isg[0]]

    # positions reellement echantillonnees (utiles pour tracer des labels honnetes)
    pos = dict(p=[float(pc[i]) for i in ip], T=[float(tc[i]) for i in it],
               u=[float(dc[i]) for i in iu], s=[float(sc[i]) for i in isg])

    cols = [t]
    cols += [P[:, i] for i in ip]
    cols += [T[:, i] for i in it]
    cols += [-U[:, i, 1] for i in iu]          # signe du trace : tassement compte positif
    cols += list(sxx) + [syy]
    return np.column_stack(cols), pos


def computeAnalytical(par):
    h, M = par["h"], par["M"]
    lamda, G, alpha, beta, F = par["lamda"], par["G"], par["alpha"], par["beta"], par["F"]
    time = np.logspace(-2, 5, 100, endpoint=True)

    zp = [0.0, 0.6, 0.8]        # z/h pour pression, temperature et sigma_xx
    zu = [0.2, 0.6, 1.0]        # z/h pour le deplacement (z = 1.4, 4.2, 7.0 m)

    P = np.empty((len(time), 3)); T = np.empty((len(time), 3))
    U = np.empty((len(time), 3)); SXX = np.empty((len(time), 3))
    for k, zc in enumerate(zp):
        for i, ti in enumerate(time):
            p_an = float(computePP(zc * h, ti, h, par["phi12"], par["phi11"], par["r11"], par["r22"]))
            th_an = float(par["theta_a"] + computeTemp(zc * h, ti, h, par["phi22"], par["phi21"],
                                                       par["r11"], par["r22"]))
            P[i, k] = p_an
            T[i, k] = th_an + 273.15
            SXX[i, k] = -(lamda / M) * F - (2.0 * G / M) * (beta * th_an + alpha * p_an)
    for k, zc in enumerate(zu):
        for i, ti in enumerate(time):
            U[i, k] = -float(computeDisp(zc * h, ti, h, par["t1"], par["t2"],
                                         par["r11"], par["r22"], par["A"], M))
    SYY = np.full((len(time), 1), -F)
    return np.column_stack([time, P, T, U, SXX, SYY])


GEOS_HEADER = ("time,p_0m,p_4p2m,p_5p6m,T_0m,T_4p2m,T_5p6m,"
               "disp_1p4m,disp_4p2m,disp_7m,sxx_0m,sxx_4p2m,sxx_5p6m,syy_0m")
ANA_HEADER = ("time,p_zh0p0,p_zh0p6,p_zh0p8,T_zh0p0,T_zh0p6,T_zh0p8,"
              "disp_1p4m,disp_4p2m,disp_7m,sxx_0m,sxx_4p2m,sxx_5p6m,syy")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("-o", "--outdir", default=".", help="dossier ou ecrire les CSV")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    geos, pos = readGeos()
    fg = os.path.join(args.outdir, "thermalConsolidation_geos.csv")
    np.savetxt(fg, geos, delimiter=",", header=GEOS_HEADER, comments="", fmt="%.8e")
    print(f"{fg}  ({geos.shape[0]} instants)")
    print(f"  cotes echantillonnees  pression/temperature {pos['p']} m, "
          f"deplacement {pos['u']} m, contrainte {pos['s']} m")

    par = analyticalParameters()
    ana = computeAnalytical(par)
    fa = os.path.join(args.outdir, "thermalConsolidation_analytical.csv")
    np.savetxt(fa, ana, delimiter=",", header=ANA_HEADER, comments="", fmt="%.8e")
    print(f"{fa}  ({ana.shape[0]} instants)")


if __name__ == "__main__":
    main()
