import os
import sys
sys.path.append("/data/PLI/sytuan/GEOSX/Libs/matplotlib")
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ElementTree

def main():
	# File paths
	path = "ModifiedCamClayResults.txt"
	xmlFilePath = "../../../../../../../inputFiles/triaxialDriver/triaxialDriver_base.xml"
	xmlFilePath_case = "../../../../../../../inputFiles/triaxialDriver/triaxialDriver_ModifiedCamClay.xml"
	imposedStrainFilePath = "../../../../../../../inputFiles/triaxialDriver/tables/axialStrain.geos"
	outputPath = "ModifiedCamClay.png"
		 
	# Load GEOSX results
	time, ax_strain, ra_strain1, ra_strain2, ax_stress, ra_stress1, ra_stress2, newton_iter, residual_norm = np.loadtxt(
		path, skiprows=5, unpack=True)

	# Extract mechanical parameters from XML file
	tree = ElementTree.parse(xmlFilePath)
	tree_case = ElementTree.parse(xmlFilePath_case)
	model = tree_case.find('Tasks/TriaxialDriver')
	param = tree.find('Constitutive/ModifiedCamClay')

	refPressure = float(param.get("defaultRefPressure"))
	refStrainVol = float(param.get("defaultRefStrainVol"))
	shearModulus = float(param.get("defaultShearModulus"))
	preConsolidationPressure = float(param.get("defaultPreConsolidationPressure"))
	cslSlope = float(param.get("defaultCslSlope"))
	virginCompressionIndex = float(param.get("defaultVirginCompressionIndex"))
	recompressionIndex = float(param.get("defaultRecompressionIndex"))
	initialStress = float(model.get("initialStress"))

	# Extract loading from input tables
	imp_strain = np.loadtxt(
		imposedStrainFilePath, skiprows=0, unpack=True)

	list_strain_anal = []
	numStepPerLoadingPeriod = 10000

	for i in range(0,len(imp_strain)-1):
		dStrainPerStep = (imp_strain[i+1]-imp_strain[i])/numStepPerLoadingPeriod	
		loadingPeriod = np.arange(imp_strain[i],imp_strain[i+1]+dStrainPerStep,dStrainPerStep)
		list_strain_anal = np.concatenate((list_strain_anal, loadingPeriod), axis=0)

	
	list_stress_anal = np.zeros(len(list_strain_anal))
	
	p_anal = refPressure
	list_strain_anal += refStrainVol/3.0
	list_stress_anal += initialStress

	
	for idx in range(1,len(list_strain_anal)):
		delta_strain_anal = list_strain_anal[idx]-list_strain_anal[idx-1]
		
		# Compute elastic moduli
		bulkModulus = - p_anal/recompressionIndex
		
		# Elastic trial
		delta_stress_anal = 3.0*delta_strain_anal*bulkModulus
		
		stress_anal = list_stress_anal[idx-1] + delta_stress_anal

		p_anal = stress_anal
		q_anal = 0. # Isotropic loading
		
		# Plastic correction
		F_anal = q_anal*q_anal + cslSlope*cslSlope*p_anal*(p_anal-preConsolidationPressure)
		
		if(F_anal>=0):
			# Derivatives
			dF_dp = (2.0*p_anal-preConsolidationPressure)
			dF_dpc = -p_anal

			dG_dp = dF_dp # associated plastic rule was considered
			dpc_dlambda = -preConsolidationPressure/(virginCompressionIndex-recompressionIndex)*dG_dp
			hardeningRate = -dF_dpc*dpc_dlambda

			# Elasto-plastic modulus
			plasticBulkModulus = 1.0/(1.0/bulkModulus + dG_dp*dF_dp/hardeningRate)
			
			delta_strain_anal = list_strain_anal[idx] - list_strain_anal[idx-1]
			delta_stress_anal = 3.0*delta_strain_anal*plasticBulkModulus
			print(delta_stress_anal)
			break
			stress_anal = list_stress_anal[idx-1] + delta_stress_anal
			
			delta_lambda = dF_dp*delta_stress_anal/hardeningRate
			delta_pc = dpc_dlambda * delta_lambda
			preConsolidationPressure += delta_pc
			
			p_anal = stress_anal
			q_anal = 0. # Isotropic loading
		
				
		list_stress_anal[idx] = stress_anal
			
	#Visualization
	N1 = 1
	fsize = 30
	msize = 12
	lw = 6
	malpha = 0.5
	ax = plt.figure(figsize=(10, 10))
	cmap = plt.get_cmap("tab10")

	plt.plot(-ax_strain * 100,
		       -ax_stress*1e-3,
		       'o',
		       color=cmap(0),
		       mec='b',
		       markersize=msize,
		       alpha=malpha,
		       label='Triaxial Driver')
	

	plt.plot(-list_strain_anal* 100,
		       -list_stress_anal*1e-3,
		       '-',
		       color='r',
		       mec='r',
		       markersize=msize,
		       alpha=malpha,
		       label='Semi-Analytic', linewidth=6)
	
	#plt.yscale('log')
	plt.xlabel(r'Volumetric Strain (%)', size=fsize, weight="bold")
	plt.ylabel(r'Mean Stress (kPa)', size=fsize, weight="bold")	
	plt.xticks(fontsize=fsize)
	plt.yticks(fontsize=fsize, rotation=90)

	plt.subplots_adjust(left=0.2, bottom=0.1, right=0.9, top=0.9, wspace=0.4, hspace=0.4)

	plt.savefig(outputPath)
	os.system("xdg-open " + outputPath)

if __name__ == "__main__":
    main()
