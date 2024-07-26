import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import h5py
import xml.etree.ElementTree as ElementTree
import math
from math import sin,cos,tan,exp,atan,asin
import csv


def main():
	poisson_ratio = 0.25 #calculated from shear and bulk modulus
	grainModulus = 1e27
	young_modulus = 100e6
	bulk_modulus = young_modulus / (3 * (1 - 2 * poisson_ratio))
	BiotCoefficient = 1- (bulk_modulus/grainModulus)


	file = open("simulation_result_0.csv")
	csvreader = csv.reader(file)
	header = next(csvreader)
	#print(header)
	rows = []
	for row in csvreader:
		rows.append(row)
	file.close() 

	rows = np.array(rows)
	zloc_0 = np.empty(len(rows[:,23]))
	pressure_0 = np.empty(len(rows[:,23]))
	tsxx_0 = np.empty(len(rows[:,23]))
	tsyy_0 = np.empty(len(rows[:,23]))
	tszz_0 = np.empty(len(rows[:,23]))
	for i in range(0,len(rows[:,23])):
		zloc_0[i]=-(float(rows[i,12]))
		pressure_0[i]=float(rows[i,15])		
		tsxx_0[i]=-(float(rows[i,22])-BiotCoefficient*pressure_0[i])/1.0e6
		tsyy_0[i]=-(float(rows[i,23])-BiotCoefficient*pressure_0[i])/1.0e6
		tszz_0[i]=-(float(rows[i,24])-BiotCoefficient*pressure_0[i])/1.0e6


	fluid_density = 1000
	solid_density = 2500
	porosity = 0.375
	bulk_density = (1-porosity)*solid_density + porosity*fluid_density
	gravity = 9.81
	poisson_ratio = 0.25 #calculated from shear and bulk modulus
	z_analytical= np.linspace(0, 1000, 100)
	#pp_analytical= 10.2*z_analytical/1000
	pp_analytical= fluid_density*gravity*z_analytical/1.0e6
	#szz_analtyical= 2405.5*9.8*(z_analytical-zloc_0[0])/1000000+18
	szz_analtyical= bulk_density*gravity*(z_analytical-zloc_0[0])/1.0e6

	#sxx_analtyical=0.15/(1-0.15)*(szz_analtyical-BiotCoefficient*pp_analytical)+BiotCoefficient*pp_analytical
	sxx_analtyical=poisson_ratio/(1-poisson_ratio)*(szz_analtyical-BiotCoefficient*pp_analytical)+BiotCoefficient*pp_analytical

	fsize = 20
	msize = 12
	lw = 6
	mew = 2
	malpha = 0.6
	lalpha = 0.8
	N1=1

	fig = plt.figure(figsize=(10,8))
	cmap = plt.get_cmap("tab10")

        
	plt.plot(tsxx_0[::N1], zloc_0[::N1], 'o', color=cmap(0), markersize=msize, alpha=malpha, mec=cmap(0), fillstyle='none', mew=mew, label= 'Sxx_Total_GEOS')
	plt.plot(sxx_analtyical, z_analytical, lw=lw, alpha=0.8, color='orange', linestyle= ':', label='Sxx_Total_Analytical')
	plt.plot(tsyy_0[::N1], zloc_0[::N1], 's', color=cmap(1), markersize=msize, alpha=malpha, mec=cmap(1), fillstyle='none', mew=mew, label= 'Syy_Total_GEOS')
	plt.plot(tszz_0[::N1], zloc_0[::N1], 'd', color=cmap(2), markersize=msize, alpha=malpha, mec=cmap(2), fillstyle='none', mew=mew, label= 'Szz_Total_GEOS')
	plt.plot(szz_analtyical, z_analytical, lw=lw, alpha=0.8, color='g', linestyle= ':', label='Szz_Total_Analytical')
	plt.plot(pressure_0[::N1]/1.0e6, zloc_0[::N1], 'x', color=cmap(3), markersize=msize, alpha=malpha, mec=cmap(3), fillstyle='none', mew=mew, label= 'Pore Pressure_GEOS')
	plt.plot(pp_analytical, z_analytical, lw=lw, alpha=0.8, color='r', linestyle= ':', label='Pore Pressure_Analytical')
	#plt.set_xlim(-1000, 1000)
	#plt.set_ylim(8000, 1000)
	plt.xlabel('Total Stresses [MPa]', size=fsize, weight="bold")
	plt.ylabel('Depth [m]', size=fsize, weight="bold")
	plt.legend(loc='upper right',fontsize=fsize*0.5)
	plt.grid(True)
	ax = plt.gca()
	ax.xaxis.set_tick_params(labelsize=fsize)
	ax.yaxis.set_tick_params(labelsize=fsize)
	ax.invert_yaxis()
	plt.savefig('Profile.png')
	plt.show() 


if __name__ == "__main__":
        main()



