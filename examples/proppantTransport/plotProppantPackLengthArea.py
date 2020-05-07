import sys
import matplotlib
import numpy as np

import matplotlib.pyplot as plt

if len(sys.argv) < 3:
    sys.exit('Usage: %s expFileName modelFileName' % sys.argv[0])


eTime, eLength, eArea = np.loadtxt(sys.argv[1], skiprows=0, unpack=True)
mTime, mLength, mArea = np.loadtxt(sys.argv[2], skiprows=1, unpack=True)

# offSize is the size of the bounday cell which should be excluded from the model domain

offSize = 0.0127

# slot length 4 feet = 1.2192 m
totalLength = 1.2192

#total cell number of the domain (excluding the bounday cells)  

totalCellNum = 96 * 24

font = {'size'   : 12}

matplotlib.rc('font', **font)

fig = plt.figure(figsize=(8, 4))

plt.subplot(1, 2, 1)

plt.plot(mTime, (mLength-offSize)/totalLength, 'b', label='GEOSX model')
plt.plot(eTime, eLength, 'rX', label='Experiment')

plt.ylabel('Normalized proppant bed length', multialignment='center')
plt.xlabel('Time (second)')

plt.xlim(0, 40)
plt.ylim(0.5, 0.8)

plt.legend(frameon=False)

plt.subplot(1, 2, 2)

plt.plot(mTime, mArea/totalCellNum, 'b', label='GEOSX model')
plt.plot(eTime, eArea, 'rX', label='Experiment')

plt.ylabel('Normalized proppant bed area', multialignment='center')
plt.xlabel('Time (second)')

plt.xlim(0, 40)
plt.ylim(0.0, 0.5)

plt.legend(frameon=False)

plt.subplots_adjust(top=0.85, bottom=0.15, left=0.12, right=0.95, hspace=0.1, wspace=0.35)

plt.show()


