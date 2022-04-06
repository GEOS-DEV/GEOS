import h5py
import os
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np

print('In file ' + os.getcwd() )
pv_tbl = pd.read_csv('./postt/class09_pv.csv')
pv = (pv_tbl['vol']*pv_tbl['rockPorosity_referencePorosity']).to_numpy()

ecl = pd.read_csv('./post/ECLSLB_hyst.csv') 
ref_co2g_ti = ecl['ECLSLB_co2g'].to_numpy()
ref_co2g  = ecl['Unnamed: 1'].to_numpy()
ref_co2l_ti = ecl['ECLSLB_co2l'].to_numpy()
ref_co2l = ecl['Unnamed: 3'].to_numpy()

 
fr = h5py.File(os.getcwd() + '/rhoHistory.hdf5','r');
fs = h5py.File(os.getcwd() + '/satHistory.hdf5','r');
fx = h5py.File(os.getcwd() + '/xcompHistory.hdf5','r');

kr = list( fr.keys() )
ks = list( fs.keys() )
kx = list( fx.keys() )

dens = fr.get('fluid_phaseDensity')
sat = fs.get('phaseVolumeFraction')
xcomp = fx.get('fluid_phaseCompFraction')

ntime, ncells, nphase = dens.shape
_, _, ncomp = xcomp.shape

times = np.arange(1,ntime,1)
nt = len(times)
ti = 0
mc = np.zeros((ncomp,nt,ncells))

for time in times: 
	for comp in range(ncomp):
		mc[comp,ti,:] = pv*sat[time,:,int((comp-1)/2)]*dens[time,:,int((comp-1)/2)]*xcomp[time,:,comp]
	ti = ti + 1
	print(ti)


id_co2_g = 0
id_co2_l = 2 

sco2_l = np.zeros(nt)
sco2_g = np.zeros(nt)
icount = 0
year_in_sec = 86400*365
dt = 1.0e8/year_in_sec

for ti in range(mc.shape[1]):
	sco2_l[ti] = np.sum(mc[id_co2_l,ti,:])
	sco2_g[ti] = np.sum(mc[id_co2_g,ti,:])

plt.plot(np.arange(0,nt)*dt+dt, sco2_l)
plt.plot(np.arange(0,nt)*dt+dt, sco2_g)
#plt.plot(np.arange(0,nt)*dt, sco2_l + sco2_g)

plt.plot(ref_co2g_ti, ref_co2g*1e9,'k')
plt.plot(ref_co2l_ti, ref_co2l*1e9,'--k')

plt.legend(['m_co2_l','m_co2_g','ref(ECLSLB)','ref'])
plt.xlabel('time[years]')
plt.ylabel('Co_2 mass [kg]')
plt.xlim([0,50])
plt.ylim([0,14e9])
plt.savefig('figure_compare.png')

#plt.show()

