import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import h5py
import xml.etree.ElementTree as ElementTree
from mpmath import *
import math
from math import sin,cos,tan,exp,atan,asin
import csv

def getHydromechanicalParametersFromXML(xmlFilePath):
	tree = ElementTree.parse(xmlFilePath)

	param1 = tree.find('Constitutive/ElasticIsotropic')
	param2 = tree.find('Constitutive/BiotPorosity')
	param3 = tree.find('Constitutive/CompressibleSinglePhaseFluid')

	hydromechanicalParameters = dict.fromkeys([
            "bulkModulus", "shearModulus", "youngModulus", "poissonRatio", "rockDensity", "poissonRatio", "biotCoefficient", "porosity", "fluidDensity", "traction"])

	hydromechanicalParameters["rockDensity"] = float(param1.get("defaultDensity"))
	hydromechanicalParameters["poissonRatio"] = float(param1.get("defaultPoissonRatio"))
	hydromechanicalParameters["youngModulus"] = float(param1.get("defaultYoungModulus"))

	E = hydromechanicalParameters["youngModulus"] 
	nu = hydromechanicalParameters["poissonRatio"]
	K = E / (3 * (1 - 2 * nu))
	G = E / (2 * (1 + nu))

	hydromechanicalParameters["poissonRatio"] = nu
	hydromechanicalParameters["bulkModulus"] = K
	hydromechanicalParameters["shearModulus"] = G

	Ks = float(param2.get("defaultGrainBulkModulus"))
	hydromechanicalParameters["biotCoefficient"] = 1.0 - K / Ks

	hydromechanicalParameters["porosity"] = float(param2.get("defaultReferencePorosity"))

	hydromechanicalParameters["fluidDensity"] = float(param3.get("defaultDensity"))

	param4 = tree.findall('FieldSpecifications/Traction')
	found_stress = False
	for elem in param4:
		if elem.get("name") == "tractionTop" and elem.get("tractionType") == "normal":
			traction = float(elem.get("scale")) * (-1)
			found_stress = True
		if found_stress: break 

	return hydromechanicalParameters

def inputStressGradientsMPa(stressXX=None, stressYY=None, stressZZ=None,porePressure=None):
    stress_gradients = {
        'stressXX': stressXX,
        'stressYY': stressYY,
        'stressZZ': stressZZ,
		'porePressure': porePressure
    }
    return stress_gradients 

def main():
	# File path
	xmlFile1Path = "../../../../../inputFiles/initialization/userdefinedStress_initialization_base.xml"
	xmlFile2Path = "../../../../../inputFiles/initialization/userdefinedStress_initialization_benchmark.xml"

	hydromechanicalParameters = getHydromechanicalParametersFromXML(xmlFile1Path)
	stress_gradients = inputStressGradientsMPa(0.17,0.27,0.24,0.1)
	sxx_grad = stress_gradients["stressXX"]
	syy_grad = stress_gradients["stressYY"]
	szz_grad = stress_gradients["stressZZ"]
	pp_grad = stress_gradients["porePressure"]
	BiotCoefficient = hydromechanicalParameters["biotCoefficient"]
	nu = hydromechanicalParameters["poissonRatio"]
	rhoF = hydromechanicalParameters["fluidDensity"]
	rhoR = hydromechanicalParameters["rockDensity"]
	phi = hydromechanicalParameters["porosity"]
	rhoB = (1-phi)*rhoR + phi*rhoF
	
	traction = hydromechanicalParameters["traction"]
	gravity = 9.8 
	
	file = open("simulation_result_0.csv")
	csvreader = csv.reader(file)
	header = next(csvreader)

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

	
	z_analytical= np.linspace(0, 1000, 100)

	pp_analytical= pp_grad*100000*z_analytical/1.0e6
	szz_analtyical= szz_grad*100000*(z_analytical-zloc_0[0])/1.0e6
	sxx_analtyical = sxx_grad*100000*(z_analytical-zloc_0[0])/1.0e6
	syy_analtyical = syy_grad*100000*(z_analytical-zloc_0[0])/1.0e6

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
	plt.plot(sxx_analtyical, z_analytical, lw=lw, alpha=0.8, color='b', linestyle= ':', label='Sxx_Total_Analytical')
	plt.plot(tsyy_0[::N1], zloc_0[::N1], 's', color=cmap(1), markersize=msize, alpha=malpha, mec=cmap(1), fillstyle='none', mew=mew, label= 'Syy_Total_GEOS')
	plt.plot(syy_analtyical, z_analytical, lw=lw, alpha=0.8, color='orange', linestyle= ':', label='Syy_Total_Analytical')
	plt.plot(tszz_0[::N1], zloc_0[::N1], 'd', color=cmap(2), markersize=msize, alpha=malpha, mec=cmap(2), fillstyle='none', mew=mew, label= 'Szz_Total_GEOS')
	plt.plot(szz_analtyical, z_analytical, lw=lw, alpha=0.8, color='y', linestyle= ':', label='Szz_Total_Analytical')
	plt.plot(pressure_0[::N1]/1.0e6, zloc_0[::N1], 'x', color=cmap(3), markersize=msize, alpha=malpha, mec=cmap(3), fillstyle='none', mew=mew, label= 'Pore Pressure_GEOS')
	plt.plot(pp_analytical, z_analytical, lw=lw, alpha=0.8, color='r', linestyle= ':', label='Pore Pressure_Analytical')
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
