import os
import sys
sys.path.append("/data/PLI/sytuan/GEOSX/Libs/matplotlib")

import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ElementTree

def getMechanicalParametersFromXML(xmlFilePath, xmlFilePath_case):
    tree = ElementTree.parse(xmlFilePath)
    tree_case = ElementTree.parse(xmlFilePath_case)

    model = tree_case.find('Tasks/TriaxialDriver')

    if model.get("material") == "ElasticIsotropic":
        param = tree.find('Constitutive/ElasticIsotropic')

        mechanicalParameters = dict.fromkeys(["bulkModulus", "shearModulus"])
        mechanicalParameters["bulkModulus"] = float(param.get("defaultBulkModulus"))
        mechanicalParameters["shearModulus"] = float(param.get("defaultShearModulus"))

    elif model.get("material") == "DruckerPrager":
        param = tree.find('Constitutive/DruckerPrager')

        mechanicalParameters = dict.fromkeys(["bulkModulus", "shearModulus", "cohesion", "frictionangle"])
        mechanicalParameters["bulkModulus"] = float(param.get("defaultBulkModulus"))
        mechanicalParameters["shearModulus"] = float(param.get("defaultShearModulus"))
        mechanicalParameters["cohesion"] = float(param.get("defaultCohesion"))
        mechanicalParameters["frictionAngle"] = float(param.get("defaultFrictionAngle"))
        mechanicalParameters["hardeningRate"] = float(param.get("defaultHardeningRate"))
        mechanicalParameters["initialStress"] = float(model.get("initialStress"))

    elif model.get("material") == "ExtendedDruckerPrager":
        param = tree.find('Constitutive/ExtendedDruckerPrager')

        mechanicalParameters = dict.fromkeys(
            ["bulkModulus", "shearModulus", "cohesion", "initialFrictionAngle", "residualFrictionAngle"])
        mechanicalParameters["bulkModulus"] = float(param.get("defaultBulkModulus"))
        mechanicalParameters["shearModulus"] = float(param.get("defaultShearModulus"))
        mechanicalParameters["cohesion"] = float(param.get("defaultCohesion"))
        mechanicalParameters["initialFrictionAngle"] = float(param.get("defaultInitialFrictionAngle"))
        mechanicalParameters["residualFrictionAngle"] = float(param.get("defaultResidualFrictionAngle"))

    elif model.get("material") == "DelftEgg":
        param = tree.find('Constitutive/DelftEgg')

        mechanicalParameters = dict.fromkeys(
            ["bulkModulus", "shearModulus", "preConsolidationPressure", "shapeParameter", "cslSlope"])
        mechanicalParameters["bulkModulus"] = float(param.get("defaultBulkModulus"))
        mechanicalParameters["shearModulus"] = float(param.get("defaultShearModulus"))
        mechanicalParameters["preConsolidationPressure"] = float(param.get("defaultPreConsolidationPressure"))
        mechanicalParameters["shapeParameter"] = float(param.get("defaultShapeParameter"))
        mechanicalParameters["cslSlope"] = float(param.get("defaultCslSlope"))

    elif model.get("material") == "ModifiedCamClay":
        param = tree.find('Constitutive/ModifiedCamClay')

        mechanicalParameters = dict.fromkeys(
            ["shearModulus", "preConsolidationPressure", "cslSlope", "recompressionIndex"])
        mechanicalParameters["shearModulus"] = float(param.get("defaultShearModulus"))
        mechanicalParameters["preConsolidationPressure"] = float(param.get("defaultPreConsolidationPressure"))
        mechanicalParameters["cslSlope"] = float(param.get("defaultCslSlope"))
        mechanicalParameters["recompressionIndex"] = float(param.get("defaultRecompressionIndex"))

    elif model.get("material") == "ViscoDruckerPrager":
        param = tree.find('Constitutive/ViscoDruckerPrager')

        mechanicalParameters = dict.fromkeys(["bulkModulus", "shearModulus", "cohesion", "frictionangle"])
        mechanicalParameters["bulkModulus"] = float(param.get("defaultBulkModulus"))
        mechanicalParameters["shearModulus"] = float(param.get("defaultShearModulus"))
        mechanicalParameters["cohesion"] = float(param.get("defaultCohesion"))
        mechanicalParameters["frictionAngle"] = float(param.get("defaultFrictionAngle"))
        mechanicalParameters["hardeningRate"] = float(param.get("defaultHardeningRate"))

    elif model.get("material") == "ViscoExtendedDruckerPrager":
        param = tree.find('Constitutive/ViscoExtendedDruckerPrager')

        mechanicalParameters = dict.fromkeys(
            ["bulkModulus", "shearModulus", "cohesion", "initialFrictionAngle", "residualFrictionAngle"])
        mechanicalParameters["bulkModulus"] = float(param.get("defaultBulkModulus"))
        mechanicalParameters["shearModulus"] = float(param.get("defaultShearModulus"))
        mechanicalParameters["cohesion"] = float(param.get("defaultCohesion"))
        mechanicalParameters["initialFrictionAngle"] = float(param.get("defaultInitialFrictionAngle"))
        mechanicalParameters["residualFrictionAngle"] = float(param.get("defaultResidualFrictionAngle"))

    elif model.get("material") == "ViscoDelftEgg":
        param = tree.find('Constitutive/ViscoDelftEgg')

        mechanicalParameters = dict.fromkeys(
            ["bulkModulus", "shearModulus", "preConsolidationPressure", "shapeParameter", "cslSlope"])
        mechanicalParameters["bulkModulus"] = float(param.get("defaultBulkModulus"))
        mechanicalParameters["shearModulus"] = float(param.get("defaultShearModulus"))
        mechanicalParameters["preConsolidationPressure"] = float(param.get("defaultPreConsolidationPressure"))
        mechanicalParameters["shapeParameter"] = float(param.get("defaultShapeParameter"))
        mechanicalParameters["cslSlope"] = float(param.get("defaultCslSlope"))

    elif model.get("material") == "ViscoModifiedCamClay":
        param = tree.find('Constitutive/ViscoModifiedCamClay')

        mechanicalParameters = dict.fromkeys(
            ["shearModulus", "preConsolidationPressure", "cslSlope", "recompressionIndex"])
        mechanicalParameters["shearModulus"] = float(param.get("defaultShearModulus"))
        mechanicalParameters["preConsolidationPressure"] = float(param.get("defaultPreConsolidationPressure"))
        mechanicalParameters["cslSlope"] = float(param.get("defaultCslSlope"))
        mechanicalParameters["recompressionIndex"] = float(param.get("defaultRecompressionIndex"))

    return mechanicalParameters

def main():
	# File paths
	path = "DruckerPragerResults.txt"
	xmlFilePath = "../../../../../../../inputFiles/triaxialDriver/triaxialDriver_base.xml"
	xmlFilePath_case = "../../../../../../../inputFiles/triaxialDriver/triaxialDriver_DruckerPrager.xml"
	imposedStrainFilePath = "../../../../../../../inputFiles/triaxialDriver/tables/axialStrain.geos"
	outputPath = "DruckerPrager.png"
		 
	# Load GEOSX results
	time, ax_strain, ra_strain1, ra_strain2, ax_stress, ra_stress1, ra_stress2, newton_iter, residual_norm = np.loadtxt(
		path, skiprows=5, unpack=True)

	# Extract mechanical parameters from XML file
	mechanicalParameters = getMechanicalParametersFromXML(xmlFilePath, xmlFilePath_case)
	bulkModulus = mechanicalParameters["bulkModulus"]
	shearModulus = mechanicalParameters["shearModulus"]
	cohesion = mechanicalParameters["cohesion"]
	frictionAngle = mechanicalParameters["frictionAngle"]
	hardeningRate = mechanicalParameters["hardeningRate"]
	initialStress = mechanicalParameters["initialStress"]

	# Compute Lamé modulus and Young modulus
	lameModulus = bulkModulus - 2.0/3.0*shearModulus
	youngModulus = 1.0 / ( 1.0/9.0/bulkModulus + 1.0/3.0/shearModulus )

	# Extract loading from input tables
	imp_strain = np.loadtxt(
		imposedStrainFilePath, skiprows=0, unpack=True)

	list_ax_strain_anal = []
	numStepPerLoadingPeriod = 1000

	for i in range(0,len(imp_strain)-1):
		dStrainPerStep = (imp_strain[i+1]-imp_strain[i])/numStepPerLoadingPeriod	
		loadingPeriod = np.arange(imp_strain[i],imp_strain[i+1]+dStrainPerStep,dStrainPerStep)
		list_ax_strain_anal = np.concatenate((list_ax_strain_anal, loadingPeriod), axis=0)

	list_ra_stress_anal = -10e6*np.ones(len(list_ax_strain_anal))
	list_ra_strain_anal = np.zeros(len(list_ax_strain_anal))
	list_ax_stress_anal = np.zeros(len(list_ax_strain_anal))

	# Friction and cohesion parameters
	frictionAngleRad = frictionAngle*3.1416/180.0
	cosFrictionAngle = np.cos(frictionAngleRad)
	sinFrictionAngle = np.sin(frictionAngleRad) 
	a = 6.0*cohesion/1.0e6*cosFrictionAngle/(3.0-sinFrictionAngle)
	b = 6.0*sinFrictionAngle/(3.0-sinFrictionAngle)

	# Elasto-plastic modulus
	plasticYoungModulus = 1.0/(1.0/youngModulus + (b-3.0)*(b-3.0)/9.0/hardeningRate)
	plasticModulusForRaStrain = 1.0/(1.0/2.0/shearModulus - (b-3.0)/2.0/hardeningRate)

	plasticYoungModulus_unload = 1.0/(1.0/youngModulus + (b+3.0)*(b+3.0)/9.0/hardeningRate)
	plasticModulusForRaStrain_unload = 1.0/(1.0/2.0/shearModulus + (b+3.0)/2.0/hardeningRate)
	
	list_ax_stress_anal[0] = initialStress
	list_ra_strain_anal[0] = 0

	for idx in range(1,len(list_ax_strain_anal)):
		delta_ax_strain_anal = list_ax_strain_anal[idx]-list_ax_strain_anal[idx-1]
		delta_ra_stress_anal = 0

		# Elastic trial
		delta_ra_strain_anal = (delta_ra_stress_anal-lameModulus*delta_ax_strain_anal)/(2.0*lameModulus+2.0*shearModulus)
		delta_ax_stress_anal = (lameModulus+2.0*shearModulus)*delta_ax_strain_anal + lameModulus/(lameModulus+shearModulus)*(delta_ra_stress_anal-lameModulus*delta_ax_strain_anal)

		ra_stress_anal = list_ra_stress_anal[idx]
		ax_stress_anal = list_ax_stress_anal[idx-1] + delta_ax_stress_anal
		ra_strain_anal = list_ra_strain_anal[idx-1] + delta_ra_strain_anal

		p_anal = (ax_stress_anal + 2.0 * ra_stress_anal) / 3.0 / 1.0e6
		q_anal = -(ax_stress_anal - ra_stress_anal) / 1.0e6

		# Plastic correction
		if(q_anal>=0): #loading
			F_anal = q_anal + b*p_anal - a

			if(F_anal>=0):
				delta_ax_strain_anal = list_ax_strain_anal[idx] - list_ax_strain_anal[idx-1]
				delta_ax_stress_anal = delta_ax_strain_anal*plasticYoungModulus 
			
				delta_ra_strain_anal = delta_ax_strain_anal	 - 	delta_ax_stress_anal / plasticModulusForRaStrain
		
				ax_stress_anal = list_ax_stress_anal[idx-1] + delta_ax_stress_anal
				ra_strain_anal = list_ra_strain_anal[idx-1] + delta_ra_strain_anal
		
				delta_a = (b-3.0)/3.0*delta_ax_stress_anal/1e6
		
				a += delta_a
		
		else: #unloading
			F_anal = -q_anal + b*p_anal - a

			if(F_anal>=0):
				delta_ax_strain_anal = list_ax_strain_anal[idx] - list_ax_strain_anal[idx-1]
				delta_ax_stress_anal = delta_ax_strain_anal*plasticYoungModulus_unload 
			
				delta_ra_strain_anal = delta_ax_strain_anal	 - 	delta_ax_stress_anal / plasticModulusForRaStrain_unload
		
				ax_stress_anal = list_ax_stress_anal[idx-1] + delta_ax_stress_anal
				ra_strain_anal = list_ra_strain_anal[idx-1] + delta_ra_strain_anal
		
				delta_a = (b+3.0)/3.0*delta_ax_stress_anal/1e6
		
				a += delta_a	

		list_ax_stress_anal[idx] = ax_stress_anal
		list_ra_strain_anal[idx] = ra_strain_anal
			
	list_p_anal = (list_ax_stress_anal + 2.0 * list_ra_stress_anal) / 3.0 / 1.0e6
	list_q_anal = -(list_ax_stress_anal - list_ra_stress_anal) / 1.0e6

	list_strain_vol_anal = list_ax_strain_anal + 2.0 * list_ra_strain_anal

	p_num = (ax_stress + 2.0 * ra_stress1) / 3.0 / 1.0e6
	q_num = -(ax_stress - ra_stress1) / 1.0e6

	strain_vol = ax_strain + 2.0 * ra_strain1

	#Visualization
	N1 = 1
	fsize = 30
	msize = 12
	lw = 6
	malpha = 0.5
	fig, ax = plt.subplots(1, 3, figsize=(37, 10))
	cmap = plt.get_cmap("tab10")

	ax[0].plot(-ax_strain * 100,
		       q_num,
		       'o',
		       color=cmap(0),
		       mec='b',
		       markersize=msize,
		       alpha=malpha,
		       label='Triaxial Driver')
	ax[0].plot(-ra_strain1 * 100, q_num, 'o', color=cmap(0), mec='b', markersize=msize, alpha=malpha)
	ax[0].plot(-list_ax_strain_anal* 100,
		       list_q_anal,
		       '-',
		       color='r',
		       mec='r',
		       markersize=msize,
		       alpha=malpha,
		       label='Semi-Analytic', linewidth=6)
	ax[0].plot(-list_ra_strain_anal * 100, 
			   list_q_anal, 
			   '-', 
			   color='r', 
			   mec='r', 
			   markersize=msize, 
			   alpha=malpha,
			   linewidth=6)
	ax[0].set_xlabel(r'Strain (%)', size=fsize, weight="bold")
	ax[0].set_ylabel(r'Deviatoric Stress (MPa)', size=fsize, weight="bold")
	#ax[0].legend(loc='lower right', fontsize=fsize)
	ax[0].xaxis.set_tick_params(labelsize=fsize)
	ax[0].yaxis.set_tick_params(labelsize=fsize)
	
	ax[1].plot(-ax_strain * 100,
		       -strain_vol * 100,
		       'o',
		       color=cmap(0),
		       mec='b',
		       markersize=msize,
		       alpha=malpha,
		       label='Triaxial Driver')
	ax[1].plot(-list_ax_strain_anal* 100,
		       -list_strain_vol_anal* 100,
		       '-',
		       color='r',
		       mec='r',
		       markersize=msize,
		       alpha=malpha,
		       label='Semi-Analytic', linewidth=6)
	ax[1].set_xlabel(r'Axial Strain (%)', size=fsize, weight="bold")
	ax[1].set_ylabel(r'Volumetric Strain (%)', size=fsize, weight="bold")
	#ax[1].legend(loc='lower right', fontsize=fsize)
	ax[1].xaxis.set_tick_params(labelsize=fsize)
	ax[1].yaxis.set_tick_params(labelsize=fsize)
	
	# Plan p-q
	ax[2].plot(-p_num,
		       q_num,
		       'o',
		       color=cmap(0),
		       mec='b',
		       markersize=msize,
		       alpha=malpha,
		       label='Triaxial Driver')
	ax[2].plot(-list_p_anal,
		       list_q_anal,
		       '-',
		       color='r',
		       mec='r',
		       markersize=msize,
		       alpha=malpha,
		       label='Semi-Analytic', linewidth=6)
	ax[2].set_xlabel(r'Mean stress (MPa)', size=fsize, weight="bold")
	ax[2].set_ylabel(r'Deviatoric Stress (MPa)', size=fsize, weight="bold")
	ax[2].legend(loc='lower right', fontsize=fsize)
	ax[2].xaxis.set_tick_params(labelsize=fsize)
	ax[2].yaxis.set_tick_params(labelsize=fsize)

	plt.subplots_adjust(left=0.2, bottom=0.1, right=0.9, top=0.9, wspace=0.4, hspace=0.4)

	plt.savefig(outputPath)
	os.system("xdg-open " + outputPath)

if __name__ == "__main__":
    main()
