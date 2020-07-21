import sys
import matplotlib
import numpy as np

import matplotlib.pyplot as plt

if len(sys.argv) < 3:
    sys.exit('Usage: %s expFileName modelFileName' % sys.argv[0])


eTime, eArea = np.loadtxt(sys.argv[1], skiprows=0, unpack=True)
mTime, mArea = np.loadtxt(sys.argv[2], skiprows=1, unpack=True)

#total cell number of the domain (excluding the bounday cells)  

totalCellNum = 96 * 24

font = {'size'   : 12}

matplotlib.rc('font', **font)

plt.plot(mTime, mArea/totalCellNum, 'b', label='GEOSX model')
plt.plot(eTime, eArea, 'rX', label='Experiment')

plt.ylabel('Normalized suspended proppant area', multialignment='center')
plt.xlabel('Time (second)')

plt.xlim(0, 40)
plt.ylim(0.0, 0.2)

plt.legend(frameon=False)

plt.show()


