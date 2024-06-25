import os
import sys
sys.path.append("/data/PLI/sytuan/Libs")
sys.path.append("/data/PLI/sytuan/Libs/matplotlib")

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ElementTree
import numpy as np
import h5py
from numpy import genfromtxt

# This function allows rotating the stress matrix (with Voigt notation) around an angle on the plane x-y
def stressRotation(stress, phi_x):
    rotx = np.array([[np.cos(phi_x), np.sin(phi_x), 0.], [-np.sin(phi_x), np.cos(phi_x), 0.], [0., 0., 1.]])
    return np.dot(np.dot(np.transpose(rotx), stress), rotx)

def getLabData():
	# Read the data from the file
	filename = 'labData.txt'

	# Initialize lists to hold the data
	T_C = []
	Pconf_Bar = []

	# Read the data from the file
	with open(filename, 'r') as dataFile:
		# Skip the header line
		next(dataFile)
		for line in dataFile:
			# Split the line into columns and convert to float
			try:
				t_c, pconf_bar = map(float, line.split())
				# Append the data to the lists
				T_C.append(t_c)
				Pconf_Bar.append(pconf_bar)
			except:
				break
	return T_C, Pconf_Bar

def getHydromechanicalParametersFromXML(xmlFilePath):
    tree = ElementTree.parse(xmlFilePath)
    param1 = tree.find('Constitutive/ElasticIsotropic')
   
    hydromechanicalParameters = dict.fromkeys([
        "bulkModulus", "drainedLinearTEC"
    ])

    hydromechanicalParameters["bulkModulus"] = float(param1.get("defaultBulkModulus"))
    hydromechanicalParameters["drainedLinearTEC"] = float(param1.get("defaultDrainedLinearTEC"))

    return hydromechanicalParameters

def main():
	# File path
	xmlFilePath = "ThermoElasticOedometric_base.xml"

	# Extract info from XML
	hydromechanicalParameters = getHydromechanicalParametersFromXML(xmlFilePath)

	bulkModulus = hydromechanicalParameters["bulkModulus"]
	defaultThermalExpansionCoefficient = hydromechanicalParameters["drainedLinearTEC"]
	   

	# Extract stress
	hf_stress = h5py.File('stressHistory_rock.hdf5', 'r')
	time = np.array( hf_stress.get('rockSolid_stress Time') )
	center = np.array( hf_stress.get('rockSolid_stress elementCenter') )
	stress = np.array( hf_stress.get('rockSolid_stress') )

	hf_temperature = h5py.File('temperatureHistory_rock.hdf5', 'r')
	temperature = np.array( hf_temperature.get('temperature') )

	nElements = center.shape[1]
	nTimes = time.shape[0]   

	element_idx = 0

	# Compute total stress: 
	# With the actual version of GEOS, the output stress need to be combined with the temperature contribution to obtain the total tress as follows:
	dThermalExpansionCoefficient_dT = 0.0082e-5
	referenceTemperature = 0.0
	thermalExpansionCoefficient = defaultThermalExpansionCoefficient +  dThermalExpansionCoefficient_dT*(temperature[:,element_idx] - referenceTemperature)

	stress_xx_total = stress[:,element_idx,0] - 3.0 * bulkModulus * thermalExpansionCoefficient * temperature[:,element_idx]
	stress_yy_total = stress[:,element_idx,1] - 3.0 * bulkModulus * thermalExpansionCoefficient * temperature[:,element_idx]
	stress_zz_total = stress[:,element_idx,2] - 3.0 * bulkModulus * thermalExpansionCoefficient * temperature[:,element_idx]
	stress_yz_total = stress[:,element_idx,3]
	stress_xz_total = stress[:,element_idx,4]
	stress_xy_total = stress[:,element_idx,5]

	# Lab data
	T_lab, pconf_lab = getLabData()
			
	# Plot 
	plt.plot( temperature[:,element_idx],
			  -stress_xx_total * 1e-5, # convert to bar    
			  'k+',
			  label=f'GEOS with linear TEC')
	plt.plot( T_lab,
			  pconf_lab,  
			  'ro',
			  label='Lab data')
	
	

	plt.grid()
	plt.ylabel('Pconf [Bar]')
	plt.xlabel('Temperature [C]')
	plt.xlim(20.0, 90.0)
	plt.ylim(0.0, 200.0)   

	plt.legend(loc='lower right')


	plt.savefig('stress.png')

	os.system('xdg-open stress.png')  
	
if __name__ == "__main__":
	main()
