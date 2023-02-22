import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ElementTree
import numpy as np
import h5py
from numpy import genfromtxt

def main():

	# File paths
	hdf5FilePath = "displacementHistory.hdf5"

	# Plot GEOSX results
	hf = h5py.File(hdf5FilePath, 'r')
	time = np.array( hf.get('totalDisplacement Time') )
	center = np.array( hf.get('totalDisplacement ReferencePosition') )
	displacement = np.array( hf.get('totalDisplacement') )

	nNodes = center.shape[1]
	xCoord = center[0, 0:nNodes, 0]
	yCoord = center[0, 0:nNodes, 1]
	
	ux_10000s = displacement[9, 0:nNodes, 0]
	uy_10000s = displacement[9, 0:nNodes, 1]
	
	ux_100000s = displacement[99, 0:nNodes, 0]
	uy_100000s = displacement[99, 0:nNodes, 1]

	rCoord, ur_10000s, ur_100000s = [], [], []

	# Extract radial data at theta = 0
	for idx in range(nNodes):
		if (yCoord[idx] < 1e-6):
			rCoord.append(xCoord[idx])
			ur_10000s.append(ux_10000s[idx]*1e6) # converted to um
			ur_100000s.append(ux_100000s[idx]*1e6) # converted to um
	
	# Reference results
	displacement_radial_analytic_1e4s = genfromtxt('displacement_radial_analytic_10000s.txt')
	displacement_radial_analytic_1e5s = genfromtxt('displacement_radial_analytic_100000s.txt')

	# Plot radial displacement at 1e4 (s)
	plt.plot( rCoord,
			  ur_10000s,        
			  'r+',
			  label='GEOSX: t = 1e4 (s)')

	plt.plot( displacement_radial_analytic_1e4s[:,0],
		      displacement_radial_analytic_1e4s[:,1]*1e6, # converted to um       
		      'r-',
		      label='Analytic: t = 1e4 (s)')

	# Plot radial displacement at 1e5 (s)
	plt.plot( rCoord,
			  ur_100000s,        
			  'b+',
			  label='GEOSX: t = 1e5 (s)')
	
	plt.plot( displacement_radial_analytic_1e5s[:,0],
		      displacement_radial_analytic_1e5s[:,1]*1e6, # converted to um       
		      'b-',
		      label='Analytic: t = 1e5 (s)')

	plt.grid()
	plt.ylabel(r'Displacement [$\mu$m]')
	plt.xlabel('Radial coordinate [m]')
	plt.xlim(0.15,0.4)

	plt.legend(loc='upper left')
	#plt.show()

	plt.savefig('displacement.png')

if __name__ == "__main__":
	main()
