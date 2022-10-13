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

	list_ax_strain_anal = []
	numStepPerLoadingPeriod = 10000

	for i in range(0,len(imp_strain)-1):
		dStrainPerStep = (imp_strain[i+1]-imp_strain[i])/numStepPerLoadingPeriod	
		loadingPeriod = np.arange(imp_strain[i],imp_strain[i+1]+dStrainPerStep,dStrainPerStep)
		list_ax_strain_anal = np.concatenate((list_ax_strain_anal, loadingPeriod), axis=0)
	
	list_ax_stress_anal = np.zeros(len(list_ax_strain_anal))
	list_ra_stress_anal = np.zeros(len(list_ax_strain_anal))
	list_ra_strain_anal = np.zeros(len(list_ax_strain_anal))
	
	p_anal = refPressure
	list_ax_strain_anal += refStrainVol # Oedometric compaction: zero lateral strain
	list_ax_stress_anal += initialStress # Assuming isotropic initial stress condition
	list_ra_stress_anal += initialStress

	for idx in range(1,len(list_ax_strain_anal)):
		delta_ax_strain_anal = list_ax_strain_anal[idx]-list_ax_strain_anal[idx-1]
		delta_ra_strain_anal = 0

		# Compute elastic moduli
		bulkModulus = - p_anal/recompressionIndex
		lameModulus = bulkModulus - 2.0/3.0*shearModulus
		
		# Elastic trial
		delta_ax_stress_anal = (lameModulus+2.0*shearModulus)*delta_ax_strain_anal + 2.0*lameModulus*delta_ra_strain_anal
		delta_ra_stress_anal = lameModulus*delta_ax_strain_anal + (2.0*lameModulus+2.0*shearModulus)*delta_ra_strain_anal

		ax_stress_anal = list_ax_stress_anal[idx-1] + delta_ax_stress_anal
		ra_stress_anal = list_ra_stress_anal[idx-1] + delta_ra_stress_anal

		p_anal = (ax_stress_anal + 2.0*ra_stress_anal)/3.0
		q_anal = ax_stress_anal - ra_stress_anal

		# Plastic correction
		F_anal = q_anal*q_anal + cslSlope*cslSlope*p_anal*(p_anal-preConsolidationPressure)		
		
		if(F_anal>=0):
			# Derivatives
			dF_dp = cslSlope*cslSlope*(2.0*p_anal-preConsolidationPressure)
			dF_dq = 2.0*q_anal
			dF_dpc = -cslSlope*cslSlope*p_anal

			dG_dp = dF_dp # associated plastic rule was considered
			dG_dq = dF_dq

			dpc_dlambda = -preConsolidationPressure/(virginCompressionIndex-recompressionIndex)*dG_dp
			hardeningRate = -dF_dpc*dpc_dlambda

			# Elasto-plastic coefficients
			coeff_1 = 1.0/bulkModulus + dG_dp*dF_dp/hardeningRate
			coeff_2 = dG_dp*dF_dq/hardeningRate
			coeff_3 = 3.0/2.0*dG_dq*dF_dp/hardeningRate
			coeff_4 = 1.0/2.0/shearModulus + 3.0/2.0*dG_dq*dF_dq/hardeningRate		
			denom = coeff_1*coeff_4 - coeff_2*coeff_3

			# Compute stress variations
			delta_ax_strain_anal = list_ax_strain_anal[idx] - list_ax_strain_anal[idx-1]
		
			delta_strain_vol = delta_ax_strain_anal
			delta_strain_shear = delta_ax_strain_anal
			delta_p_anal = (coeff_4*delta_strain_vol - coeff_2*delta_strain_shear)/denom
			delta_q_anal = (coeff_1*delta_strain_shear - coeff_3*delta_strain_vol)/denom
			
			delta_ax_stress_anal = (3.0*delta_p_anal+2*delta_q_anal)/3.0
			delta_ra_stress_anal = delta_ax_stress_anal-delta_q_anal
		
			# Update stress
			ax_stress_anal = list_ax_stress_anal[idx-1] + delta_ax_stress_anal
			ra_stress_anal = list_ra_stress_anal[idx-1] + delta_ra_stress_anal
			
			delta_lambda = (dF_dq*delta_q_anal + dF_dp*delta_p_anal)/hardeningRate
			delta_pc = dpc_dlambda * delta_lambda
			preConsolidationPressure += delta_pc
			
			p_anal = (ax_stress_anal + 2.0*ra_stress_anal)/3.0
			q_anal = ax_stress_anal - ra_stress_anal
						
		list_ax_stress_anal[idx] = ax_stress_anal
		list_ra_stress_anal[idx] = ra_stress_anal
		
	list_p_anal = (list_ax_stress_anal + 2.0 * list_ra_stress_anal) / 3.0
	list_q_anal = -(list_ax_stress_anal - list_ra_stress_anal)

	list_strain_vol_anal = list_ax_strain_anal + 2.0 * list_ra_strain_anal

	p_num = (ax_stress + 2.0 * ra_stress1) / 3.0
	q_num = -(ax_stress - ra_stress1)
	
	#Visualization
	N1 = 1
	fsize = 30
	msize = 12
	lw = 6
	malpha = 0.5
	fig, ax = plt.subplots(1, 3, figsize=(37, 10))
	cmap = plt.get_cmap("tab10")

	ax[0].plot(-ax_strain * 100,
		       -ax_stress*1e-3,
		       'o',
		       color=cmap(0),
		       mec='b',
		       markersize=msize,
		       alpha=malpha,
		       label='Triaxial Driver')
	ax[0].plot(-list_ax_strain_anal* 100,
		       -list_ax_stress_anal*1e-3,
		       '-',
		       color='r',
		       mec='r',
		       markersize=msize,
		       alpha=malpha,
		       label='Semi-Analytic', linewidth=6)
	
	ax[0].set_xlabel(r'Axial Strain (%)', size=fsize, weight="bold")
	ax[0].set_ylabel(r'Axial Stress (kPa)', size=fsize, weight="bold")
	#ax[0].legend(loc='lower right', fontsize=fsize)
	ax[0].xaxis.set_tick_params(labelsize=fsize)
	ax[0].yaxis.set_tick_params(labelsize=fsize)

	ax[1].plot(-ax_strain * 100,
		       -ra_stress1 * 1e-3,
		       'o',
		       color=cmap(0),
		       mec='b',
		       markersize=msize,
		       alpha=malpha,
		       label='Triaxial Driver')
	ax[1].plot(-list_ax_strain_anal* 100,
		       -list_ra_stress_anal* 1e-3,
		       '-',
		       color='r',
		       mec='r',
		       markersize=msize,
		       alpha=malpha,
		       label='Semi-Analytic', linewidth=6)
	ax[1].set_xlabel(r'Axial Strain (%)', size=fsize, weight="bold")
	ax[1].set_ylabel(r'Radial stress (kPa)', size=fsize, weight="bold")
	#ax[1].legend(loc='lower right', fontsize=fsize)
	ax[1].xaxis.set_tick_params(labelsize=fsize)
	ax[1].yaxis.set_tick_params(labelsize=fsize)
	
	# Plan p-q
	ax[2].plot(-p_num*1e-3,
		       q_num*1e-3,
		       'o',
		       color=cmap(0),
		       mec='b',
		       markersize=msize,
		       alpha=malpha,
		       label='Triaxial Driver')
	ax[2].plot(-list_p_anal*1e-3,
		       list_q_anal*1e-3,
		       '-',
		       color='r',
		       mec='r',
		       markersize=msize,
		       alpha=malpha,
		       label='Semi-Analytic', linewidth=6)
	ax[2].set_xlabel(r'Mean stress (kPa)', size=fsize, weight="bold")
	ax[2].set_ylabel(r'Deviatoric Stress (kPa)', size=fsize, weight="bold")
	ax[2].legend(loc='lower right', fontsize=fsize)
	ax[2].xaxis.set_tick_params(labelsize=fsize)
	ax[2].yaxis.set_tick_params(labelsize=fsize)
	
	plt.subplots_adjust(left=0.2, bottom=0.1, right=0.9, top=0.9, wspace=0.4, hspace=0.4)

	plt.savefig(outputPath)
	os.system("xdg-open " + outputPath)

if __name__ == "__main__":
    main()
