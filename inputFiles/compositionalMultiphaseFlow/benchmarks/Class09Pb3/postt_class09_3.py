import h5py
import os
from matplotlib import pyplot as plt
import numpy as np

print('In file ' + os.getcwd() )

fr = h5py.File(os.getcwd() + '/rhoHistory.hdf5','r');
fs = h5py.File(os.getcwd() + '/satHistory.hdf5','r');
fx = h5py.File(os.getcwd() + '/xcompHistory.hdf5','r');

kr = list( fr.keys() )
ks = list( fs.keys() )
kx = list( fx.keys() )

dens = fr[kr[-1]]
sat = fs[ks[0]] #weirdly enough
xcomp = fx[kx[-1]]

ntime, ncells, nphase = dens.shape
_, _, ncomp = xcomp.shape

times = np.arange(1,ntime,1)
nt = len(times)
ti = 0
mc = np.zeros((ncomp,nt,ncells))

for time in times: 
	for comp in range(ncomp):
		mc[comp,ti,:] = sat[time,:,int((comp-1)/2)]*dens[time,:,int((comp-1)/2)]*xcomp[time,:,comp]
	ti = ti + 1
	print(ti)

id_gas = 0
id_liq = 1

id_co2_g = 0
id_co2_l = 2 

sco2_l = np.zeros(nt)
sco2_g = np.zeros(nt)
icount = 0

for ti in range(mc.shape[1]):
	sco2_l[ti] = np.sum(mc[id_co2_l,ti,:])
	sco2_g[ti] = np.sum(mc[id_co2_g,ti,:])

plt.plot(np.arange(0,nt), sco2_l)
plt.plot(np.arange(0,nt), sco2_g)
#plt.plot(np.arange(0,nt), sco2_l + sco2_g)
plt.legend(['m_co2_l','m_co2_g','sum'])

