import os
import sys
import os
import argparse

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.linalg
from scipy import special

from scipy.sparse import diags
from scipy.sparse.linalg import spsolve

from xml.etree import ElementTree

# Analytical results for linear thermal behavior
def steadyState(Tin, Tout, Rin, Rout, radialCoordinate):
	return Tin + (Tout - Tin) * (np.log(radialCoordinate) - np.log(Rin)) / (np.log(Rout) - np.log(Rin))

def diffusionFunction(radialCoordinate, Rin, diffusionCoefficient, diffusionTime):
	return special.erfc(  (radialCoordinate - Rin) / 2.0 / np.sqrt( diffusionCoefficient * diffusionTime ) )

def computeTransientTemperature(Tin, Tout, Rin, radialCoordinate, thermalDiffusionCoefficient, diffusionTime):
	# Ref. Wang and Papamichos (1994), https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/94WR01774
	return Tout + (Tin-Tout) * np.sqrt(Rin/radialCoordinate) * diffusionFunction(radialCoordinate, Rin, thermalDiffusionCoefficient, diffusionTime)
	 
def computeThermalDiffusionCoefficient(thermalConductivity, volumetricHeatCapacity):
	return thermalConductivity / volumetricHeatCapacity
	
# Finite difference results for non-linear thermal behavior
def temperatureDependentThermalConductivity(lambda0, lambda_gradient, T, Treference):
    return lambda0 + lambda_gradient*(T-Treference)

def temperatureDependentVolumetricHeat(c0, c_gradient, T, Treference):
    return c0 + c_gradient*(T-Treference)

def coefficientMatrix(thermalConductivity, volumetricHeatCapacity, r, dt, N):
    # Coefficients for the matrix
    A = np.zeros((N+1, N+1))
    
    for i in range(1,N):
        dr = r[i] - r[i-1] 		
        r_i = r[i]

        diffusivity_i = thermalConductivity[i]/volumetricHeatCapacity[i]
        A[i, i-1] = - diffusivity_i * (dt/(dr**2) - dt/(2 * r_i * dr)) \
                    + (thermalConductivity[i+1] - thermalConductivity[i-1])/volumetricHeatCapacity[i]*dt/4/(dr**2)
        A[i, i] = 1.0 + 2.0 * diffusivity_i * dt / (dr**2)
        A[i, i+1] = - diffusivity_i * (dt/(dr**2) + dt/(2 * r_i * dr)) \
                    - (thermalConductivity[i+1] - thermalConductivity[i-1])/volumetricHeatCapacity[i]*dt/4/(dr**2)

    # Boundary conditions
    # No-flux at r=0 approximated by setting the flux between the first two cells to zero
    A[0, 0] = 1.0 
    A[N, N] = 1.0
    return A

def solve_radial_diffusion(r, tmax, dt, Tin, Tout, lambda0, lambda_gradient, c0, c_gradient, Treference):
    N = len(r)-1
    # Time setup
    n_steps = int(tmax / dt)
    
    # Time-stepping
    T = np.zeros(N+1) + Tout  # initial condition u(r, 0)
    T[0] = Tin
    for step in range(n_steps):
        thermalConductivity = temperatureDependentThermalConductivity(lambda0, lambda_gradient, T, Treference)
        
        volumetricHeatCapacity = temperatureDependentVolumetricHeat(c0, c_gradient, T, Treference)

        A = coefficientMatrix(thermalConductivity, volumetricHeatCapacity, r, dt, N)
        T = spsolve(A, T)
    
    return T


def extractDataFromXMLList(paramList):
	# Extract data from a list in XML such as "{ 1, 2, 3}"
	return paramList.replace('{', '').replace('}', '').strip().split(',')

def getWellboreGeometryFromXML(xmlFilePath):
	tree = ElementTree.parse(xmlFilePath)

	meshParam = tree.find('Mesh/InternalWellbore')
	radii = extractDataFromXMLList( meshParam.get("radius") )

	Rin = float(radii[0])
	Rout = float(radii[-1])

	return [Rin, Rout]

def getLoadingFromXML(xmlFilePath):
	tree = ElementTree.parse(xmlFilePath)
	fsParams = tree.findall('FieldSpecifications/FieldSpecification')

	for fsParam in fsParams:
		if ( (fsParam.get('fieldName') == "pressure") & (fsParam.get('initialCondition') != "1") ): 
			if fsParam.get('setNames') == "{ rneg }":
				Pin = float(fsParam.get('scale'))
			if fsParam.get('setNames') == "{ rpos }":
				Pout = float(fsParam.get('scale'))

	for fsParam in fsParams:
		if ( (fsParam.get('fieldName') == "temperature") & (fsParam.get('initialCondition') != "1") ): 
			if fsParam.get('setNames') == "{ rneg }":
				Tin = float(fsParam.get('scale'))
			if fsParam.get('setNames') == "{ rpos }":
				Tout = float(fsParam.get('scale'))


	tree_SinglePhaseThermalConductivities = tree.findall('Constitutive/SinglePhaseThermalConductivity')

	for tree_SinglePhaseThermalConductivity in tree_SinglePhaseThermalConductivities:
		if tree_SinglePhaseThermalConductivity.get('name') == "thermalCond_nonLinear":
			defaultThermalConductivity = float( extractDataFromXMLList( tree_SinglePhaseThermalConductivity.get('defaultThermalConductivityComponents') )[0] )

			thermalConductivityGradient = float( extractDataFromXMLList( tree_SinglePhaseThermalConductivity.get('thermalConductivityGradientComponents') )[0] )
			referenceTemperature = float( tree_SinglePhaseThermalConductivity.get('referenceTemperature') )


	tree_SolidInternalEnergies = tree.findall('Constitutive/SolidInternalEnergy')

	for tree_SolidInternalEnergy in tree_SolidInternalEnergies:
		if tree_SolidInternalEnergy.get('name') == "rockInternalEnergy_linear":
			volumetricHeatCapacity = float( tree_SolidInternalEnergy.get('referenceVolumetricHeatCapacity') )
			dVolumetricHeatCapacity_dTemperature = 0.0
	
	
	permeability = float( extractDataFromXMLList( tree.find('Constitutive/ConstantPermeability').get('permeabilityComponents') )[0] )

	porosity = float( tree.find('Constitutive/PressurePorosity').get('defaultReferencePorosity') )
	
	fluidViscosity = float( tree.find('Constitutive/ThermalCompressibleSinglePhaseFluid').get('defaultViscosity') )

	fluidCompressibility = float( tree.find('Constitutive/ThermalCompressibleSinglePhaseFluid').get('compressibility') )

	fluidThermalExpansionCoefficient = float( tree.find('Constitutive/ThermalCompressibleSinglePhaseFluid').get('thermalExpansionCoeff') )

	return [Pin, Pout, Tin, Tout, defaultThermalConductivity, thermalConductivityGradient, referenceTemperature, volumetricHeatCapacity, dVolumetricHeatCapacity_dTemperature, permeability, porosity, fluidViscosity, fluidCompressibility, fluidThermalExpansionCoefficient]


def main():

   # Initialize the argument parser
    parser = argparse.ArgumentParser(description="Script to generate figure from tutorial.")

    # Add arguments to accept individual file paths
    parser.add_argument('--geosDir', help='Path to the GEOS repository ', default='../../../../../../..')

    # Parse the command-line arguments
    args = parser.parse_args()

    geosDir = args.geosDir

	xmlFilePath = geosDir + "/inputFiles/singlePhaseFlow/"
	
	Rin, Rout = getWellboreGeometryFromXML(xmlFilePath+"thermalCompressible_temperatureDependentSinglePhaseThermalConductivity_benchmark.xml")

	Pin, Pout, Tin, Tout, defaultThermalConductivity, thermalConductivityGradient, referenceTemperature, volumetricHeatCapacity, dVolumetricHeatCapacity_dTemperature, permeability, porosity, fluidViscosity, fluidCompressibility, fluidThermalExpansionCoefficient = getLoadingFromXML(xmlFilePath+"thermalCompressible_2d_base.xml")

	plt.figure(figsize=(10,7))
	font = {'size'   : 16}
	plt.rc('font', **font)

	for chart_idx, idx in enumerate([1, 2, 5, 10]):
		# Numerical results
		data = pd.read_csv(f'data_{idx}.csv')
		data.dropna(inplace=True)
		data.drop_duplicates(inplace=True)
		data.reset_index(drop=True, inplace=True)
		
		radialCoordinate = (data['elementCenter:0']**2.0 + data['elementCenter:1']**2.0)**0.5
		temperature = data['temperature']
		#pressure = data['pressure']
		diffusionTime = data['Time'][0]

		# Analytical results for linear thermal behavior, for comparison
		radialCoordinate_anal = np.arange(Rin, Rout, (Rout-Rin)/100)   
       
		thermalDiffusionCoefficient = computeThermalDiffusionCoefficient(defaultThermalConductivity, volumetricHeatCapacity)

		T_transient_linear = computeTransientTemperature(Tin, Tout, Rin, radialCoordinate_anal, thermalDiffusionCoefficient, diffusionTime)
		
		# Analytical results of the steady state regime for comparison
		T_steadyState = steadyState(Tin, Tout, Rin, Rout, radialCoordinate_anal)

		# Finite different results for non-linear thermal behavior
		
		T_transient_nonLinear = solve_radial_diffusion(radialCoordinate_anal, diffusionTime, diffusionTime/1000, Tin, Tout, defaultThermalConductivity, thermalConductivityGradient, volumetricHeatCapacity, dVolumetricHeatCapacity_dTemperature, referenceTemperature)

		# Visualization
		# Temperature
		plt.subplot(2,2,chart_idx+1)
		plt.plot( radialCoordinate, temperature, 'k+' , label='GEOS' )
		plt.plot( radialCoordinate_anal, T_transient_nonLinear, 'g.' , label='FDM Non-Linear' )
		plt.plot( radialCoordinate_anal, T_transient_linear, 'r-' , label='Analytic Linear' )
		plt.plot( radialCoordinate_anal, T_steadyState, 'b-' , label='Steady State' )

		if chart_idx==1:
			plt.legend()

		if chart_idx in [2,3]:
			plt.xlabel('Radial distance from well center')

		if chart_idx in [0,2]:
			plt.ylabel('Temperature (°C)')

		plt.ylim(Tin-10,Tout+10)
		plt.xlim(0.1,0.3)
		plt.title('t = '+str(diffusionTime)+'(s)')
		plt.tight_layout()

	plt.show()

if __name__ == "__main__":
	main()

