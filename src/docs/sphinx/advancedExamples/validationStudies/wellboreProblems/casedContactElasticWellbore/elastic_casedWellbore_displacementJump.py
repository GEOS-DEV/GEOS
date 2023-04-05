import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ElementTree
import numpy as np
import h5py
import analyticalResults

def main():

	# File paths
	hdf5FilePath = "displacementJumpHistory.hdf5"

	# Plot GEOSX results
	hf = h5py.File(hdf5FilePath, 'r')
	time = np.array( hf.get('displacementJump Time') )
	displacementJump = np.array( hf.get('displacementJump') )

	nTime = time.shape[0]

	displacementJump_normal = displacementJump[0:nTime, 0, 0]
	
	
	# Reference results
	listTime = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
	listP0 = [0, -2e6, -4e6, -6e6, -8e6, -10e6, -6e6, -2e6, 2e6, 6e6, 10e6]

	list_displacementJump_anal = []
	for P0 in listP0:
		displacementJump_anal = analyticalResults.analyticDisplacementJump(P0) 
		list_displacementJump_anal.append(displacementJump_anal * 1e6) # converted to um

	# Plot displacement jump versus time
	plt.plot( time,
			  displacementJump_normal * 1e6,        # converted to um
			  'k+',
			  label='GEOSX')

	plt.plot( listTime,
			  list_displacementJump_anal,
			  'r-',
			  label='Analytic')
	
	plt.grid()
	plt.ylabel(r'Normal displacement jump [$\mu$m]')
	plt.xlabel('Time [s]')
	plt.xlim(0, 1)

	plt.legend(loc='upper left')
	plt.show()

if __name__ == "__main__":
	main()
