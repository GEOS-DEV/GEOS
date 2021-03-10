import numpy as np
import matplotlib.pyplot as plt

a = 0.1
p0 = -10.0

K = 5.5556e9
G = 4.16667e9

listr = np.arange(a,10*a,0.1*a)
listSigrr = []
listSigtt = []
for r in listr:
  sigrr = p0*a*a/r/r
  sigtt = -sigrr
  listSigrr.append(sigrr)
  listSigtt.append(sigtt)


fig = plt.figure(figsize=[13,5])

plt.subplot(121)
plt.plot(listr,listSigrr, 'k', linewidth=2, label='Analytic')


x, y = [], []
for line in open('radialStress.curve', 'r'):
	if not (line.strip().startswith("#") or line.strip()==''):
		values = [float(s) for s in line.split()]
		x.append(values[0])
		y.append(values[1]*1e-6)
plt.plot(x,y,'ko', label='GEOSX result')
plt.xlabel('r (m)')
plt.ylabel('Radial stress (MPa)')
plt.xlim(a,3*a)
plt.subplot(122)
plt.plot(listr,listSigtt, 'k', linewidth=2, label='Analytic')


x, y = [], []
for line in open('tangentStress.curve', 'r'):
	if not (line.strip().startswith("#") or line.strip()==''):
		values = [float(s) for s in line.split()]
		x.append(values[0])
		y.append(values[1]*1e-6)
plt.plot(x,y,'ko', label='GEOSX result')
plt.xlabel('r (m)')
plt.ylabel('Tangent stress (MPa)')
plt.xlim(a,3*a)
plt.legend()
plt.show()
