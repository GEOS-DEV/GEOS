import h5py
import os
import argparse
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np

# some helpers
descr = 'Post processing for Class09 Cases \n' + \
        '\t use -f to specify a dest folder\n' + \
        '\t use -t to specify a drain or hyst\n'

parser = argparse.ArgumentParser(description=descr)
parser.add_argument('-f', nargs=1, help='specify dest folder')
parser.add_argument('-t', nargs=1, help='specify drain or hyst for ref loading')
args = parser.parse_args()


def form_block_value(suffix):
    # read hdf5 timeHistory
    fr = h5py.File(os.getcwd() + '/' + args.f[0] + '/rhoHistory' + suffix + '.hdf5', 'r')
    fs = h5py.File(os.getcwd() + '/' + args.f[0] + '/satHistory' + suffix + '.hdf5', 'r')
    fx = h5py.File(os.getcwd() + '/' + args.f[0] + '/xcompHistory' + suffix + '.hdf5', 'r')

    fporo = h5py.File(os.getcwd() + '/' + args.f[0] + '/poroHistory' + suffix + '.hdf5', 'r')
    fvol = h5py.File(os.getcwd() + '/' + args.f[0] + '/volHistory' + suffix + '.hdf5', 'r')

    # get the data
    dens = fr.get('fluid_phaseDensity')
    sat = fs.get('phaseVolumeFraction')
    xcomp = fx.get('fluid_phaseCompFraction')
    poro = fporo.get('rockPorosity_porosity')
    vol = fvol.get('elementVolume')

    # get some constant
    ntime, ncells, nphase = dens.shape
    _, _, ncomp = xcomp.shape
    times = np.arange(0, ntime, 1)
    nt = len(times)
    ti = 0
    mc = np.zeros((ncomp, nt, ncells))
    mp = np.zeros((nphase, nt, ncells))

    for time in times:
        for comp in range(ncomp):
            p = (comp + 1) % 2
            pv = poro[time, :, 0] * vol[0]
            # mass per comp
            mc[comp, ti, :] = pv * sat[time, :, p] * dens[time, :, p] * xcomp[time, :, comp]
            # mass per phase
            mp[p, ti, :] += pv * sat[time, :, p] * dens[time, :, p] * xcomp[time, :, comp]
        ti = ti + 1

    fr.close()
    fs.close()
    fx.close()
    fporo.close()
    fvol.close()

    return mc, mp


if (args.t[0] == 'hyst'):
    ecl = pd.read_csv(os.getcwd() + '/postt/ECLSLB_hyst.csv')
elif (args.t[0] == 'drain'):
    ecl = pd.read_csv(os.getcwd() + '/postt/ECLSLB_drain.csv')
else:
    raise NotImplemented

ref_co2g_ti = ecl['ECLSLB_co2g'].to_numpy()
ref_co2g = ecl['Unnamed: 1'].to_numpy()
ref_co2l_ti = ecl['ECLSLB_co2l'].to_numpy()
ref_co2l = ecl['Unnamed: 3'].to_numpy()
##

# get interior values
mi, mpi = form_block_value('')
# get interior values
mc, mpc = form_block_value('_bc')

##
id_co2_g = 3
id_co2_l = 2

# resize as time
nt = mi.shape[1]
sco2_l = np.zeros(nt)
sco2_g = np.zeros(nt)
stot = np.zeros(nt)

sl = np.zeros(nt)
sv = np.zeros(nt)

icount = 0
# some conversion
year_in_sec = 86400 * 365
dt = 1.0e8 / year_in_sec

for ti in range(mc.shape[1]):
    sco2_l[ti] = np.sum(mi[id_co2_l, ti, :])
    sco2_l[ti] += np.sum(mc[id_co2_l, ti, :])
    stot[ti] = sco2_l[ti]

    sco2_g[ti] = np.sum(mi[id_co2_g, ti, :])
    sco2_g[ti] += np.sum(mc[id_co2_g, ti, :])
    stot[ti] += sco2_g[ti]

plt.figure()
plt.plot(np.arange(0, nt) * dt , sco2_l, '-+')
plt.plot(np.arange(0, nt) * dt , sco2_g, '--+')
plt.plot(np.arange(0, nt) * dt , stot, '-r+')
#
plt.plot(ref_co2g_ti, ref_co2g * 1e9, 'k')
plt.plot(ref_co2l_ti, ref_co2l * 1e9, '--k')
plt.plot(ref_co2l_ti, (ref_co2g + ref_co2l) * 1e9, ':k')

plt.legend(['m_co2_l', 'm_co2_g', 'm_co2_tot', 'ref(ECLSLB) m_co2_l', 'ref(ECLSLB) m_co2_g', 'ref(ECLSLB) m_co2_tot'])
plt.xlabel('time[years]')
plt.ylabel('Co_2 mass [kg]')
plt.xlim([0, 50])
plt.ylim([0, 14e9])
plt.savefig('figure_compare_' + args.f[0] + '_' + args.t[0] + '.png')
