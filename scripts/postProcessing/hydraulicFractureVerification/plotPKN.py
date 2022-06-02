import sys
import matplotlib
import numpy as np

import math
import matplotlib.pyplot as plt


PI = 3.14159265359


def ReadPlotData(filename, headerlines, x, y):
    f = open(filename, "r")

    lines = f.readlines()

    for i in range(0, headerlines):
        lines.pop(0)

    #    print lines[0]
    for line in lines:
        vals = line.split(" ")
        x.append(float(vals[0]))
        y.append(float(vals[1]))


# if len(sys.argv) < 6:
#    sys.exit('Usage: %s prefix mu E\' q KI' % sys.argv[0])

prefix = sys.argv[1]
# mu = float(sys.argv[2])
# E  = float(sys.argv[3])
# q  = float(sys.argv[4])
# KI = float(sys.argv[5])
# h  = float(sys.argv[6])

mu = 0.001
E = 3.0e10
q = 4 * (0.001 * 10)
KI = 0.0
h = 20.0
nu = 0.25
G = E / (2 * (1 + nu))
Eprime = E / (1 - nu * nu)

t_sim, aper0 = np.loadtxt(prefix + "_aper0.curve", skiprows=1, unpack=True)
t_sim, Area = np.loadtxt(prefix + "_area.curve", skiprows=1, unpack=True)
t_sim, press0 = np.loadtxt(prefix + "_pressure0.curve", skiprows=1, unpack=True)


Length = 2 * Area / h


soln_time = np.zeros(1000)
soln_aper0 = np.zeros_like(soln_time)
soln_Length = np.zeros_like(soln_time)
soln_press0 = np.zeros_like(soln_time)

for i in range(0, len(soln_time)):
    soln_time[i] = t_sim[-1] / 1000.0 * i
    soln_aper0[i] = 2.18 * pow(mu * q * q / (Eprime * h) * soln_time[i], 0.2)
    soln_Length[i] = (
        0.39 * pow(Eprime * pow(q, 3) / (mu * pow(h, 4)), 0.2) * pow(soln_time[i], 0.8)
    )
    # soln_press0[i] = soln_aper0[i] * Eprime / ( 2 * h )
    # soln_press0[i] = 2.5 * pow( pow(G,4)*mu*q*q/(4*pow(1-nu,4)*pow(h,6))*soln_time[i],0.2 )
    soln_press0[i] = 2 / h * pow(mu * q * pow(Eprime, 3) * soln_Length[i] / PI, 0.25)


N1 = 4

fig1 = plt.figure()  # a new figure window
ax1 = fig1.add_subplot(3, 1, 1)  # specify (nrows, ncols, axnum)
plt.plot(
    t_sim[1::N1], aper0[1::N1] * 1000, "o", fillstyle="none", markersize=5, label="GEOS"
)
plt.plot(soln_time, soln_aper0 * 1000, "k", label="Nordgren (1972)")
# plt.xlabel('time (s)')
plt.ylabel("Wellbore\n Aperture (mm)", multialignment="center")
plt.legend(loc="lower right")

ax2 = fig1.add_subplot(3, 1, 2)  # specify (nrows, ncols, axnum)
plt.plot(t_sim[1::N1], Length[1::N1], "o", fillstyle="none", markersize=5, label="GEOS")
plt.plot(soln_time, soln_Length, "k", label="Nordgren (1972)")
# plt.xlabel('time (s)')
plt.ylabel("Fracture\n Half-Length (m)", multialignment="center")
# plt.xlim([ 0.0  , 40.0 ])
# plt.ylim([ 0.0  , 50.0 ])


ax2 = fig1.add_subplot(3, 1, 3)  # specify (nrows, ncols, axnum)
plt.plot(
    t_sim[1::N1],
    press0[1::N1] / 1.0e6,
    "o",
    fillstyle="none",
    markersize=5,
    label="GEOS",
)
plt.plot(soln_time, soln_press0 / 1.0e6, "k", label="Nordgren (1972)")
plt.xlabel("time (s)")
plt.ylabel("Wellbore\nPressure (MPa)", multialignment="center")
plt.ylim([0.0, 2.0])

# plt.savefig(prefix + '_vsTime.eps')
plt.show()
