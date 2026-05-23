import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import os

# Set global font size to 12
plt.rcParams.update({'font.size': 12})

# Load data (ignore comment lines starting with '#')
data_dir = "../../../../../../../inputFiles/constitutiveDriver/testRelperm_data"
file_name = "test_relperm_sandstone_hysteresis.txt"
data_file = os.path.join(data_dir, "..", file_name)
if not os.path.exists(data_file):
    data_file = os.path.join(data_dir, file_name)
data = np.loadtxt(data_file)

# Columns:
# col 2 = gas saturation (Sg)
# col 4 = historical maximum gas saturation (Sghy)
# col 6 = gas relperm (krg)
Sg   = data[:, 1]
krg  = data[:, 5]
Sghy = data[:, 3]

# Create segments for the colored line
points = np.array([Sg, krg]).T.reshape(-1, 1, 2)
segments = np.concatenate([points[:-1], points[1:]], axis=1)

fig, ax = plt.subplots()

# Create the LineCollection
norm = plt.Normalize(Sghy.min(), Sghy.max())
lc = LineCollection(segments, cmap='jet', norm=norm)

# Set the colors based on the average Sghy of the two points making up each segment
lc.set_array(0.5 * (Sghy[:-1] + Sghy[1:]))
lc.set_linewidth(2) # 2pt line width

line = ax.add_collection(lc)

# Autoscale the axes since add_collection doesn't do it automatically
ax.autoscale()

ax.set_xlabel('Gas saturation (Sg)')
ax.set_ylabel('Gas relative permeability (krg)')
ax.set_title('Gas RelPerm vs Gas Saturation')
fig.colorbar(line, ax=ax, label='Sghy')
ax.grid(True, alpha=0.3)

plt.show()
