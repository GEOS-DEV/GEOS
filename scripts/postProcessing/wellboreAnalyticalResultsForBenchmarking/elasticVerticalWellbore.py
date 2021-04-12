import numpy as np
import matplotlib.pyplot as plt

a = 0.1
p0 = -10.0

K = 5.5556e9
G = 4.16667e9

fig = plt.figure(figsize=[13,10])

# Radial stress
r, sig_geosx, sig_anal, dsig, dsigRelative = [], [], [], [], []
for line in open('radialStress.curve', 'r'):
	if not (line.strip().startswith("#") or line.strip()==''):
		values = [float(s) for s in line.split()]
		rval = values[0]
		sigVal = values[1]*1e-6 # convert to MPa
		sigValAnal = p0*a*a/rval/rval
		dsigVal = sigVal - sigValAnal

		r.append( rval )
		sig_geosx.append( sigVal )
		sig_anal.append( sigValAnal )
		dsig.append( dsigVal )
		dsigRelative.append( (dsigVal / sigValAnal) * 100 )

plt.subplot(221)
plt.plot(r, sig_geosx, 'ko', label='GEOSX result')
plt.plot(r, sig_anal,  'k', linewidth=2, label='Analytic')
plt.ylabel('Stress (MPa)')
plt.xlim(a,3*a)
plt.title("Radial stress")

ax = plt.subplot(223)
plt.plot(r, dsigRelative, 'k', linewidth=2, label='Relative error')
plt.xlabel('r (m)')
plt.ylabel('Relative error (%)')
plt.ylim(-10,10)
ax2 = ax.twinx()
plt.plot(r, dsig,  'k--', linewidth=2, label='Absolute error')
plt.ylim(-1,1)
plt.xlim(a,3*a)
plt.title('Radial stress error')

# Tangent stress
r, sig_geosx, sig_anal, dsig, dsigRelative = [], [], [], [], []
for line in open('tangentStress.curve', 'r'):
	if not (line.strip().startswith("#") or line.strip()==''):
		values = [float(s) for s in line.split()]
		rval = values[0]
		sigVal = values[1]*1e-6 # convert to MPa
		sigValAnal = -p0*a*a/rval/rval
		dsigVal = sigVal - sigValAnal

		r.append( rval )
		sig_geosx.append( sigVal )
		sig_anal.append( sigValAnal )
		dsig.append( dsigVal )
		dsigRelative.append( (dsigVal / sigValAnal) * 100 )

plt.subplot(222)
plt.plot(r, sig_geosx, 'ko', label='GEOSX result')
plt.plot(r, sig_anal,  'k', linewidth=2, label='Analytic')
plt.ylabel('Stress (MPa)')
plt.xlim(a,3*a)
plt.legend()
plt.title('Tangent stress')

ax = plt.subplot(224)
l1 = plt.plot(r, dsigRelative, 'k', linewidth=2, label='Relative error')
plt.xlabel('r (m)')
plt.ylim(-10,10)
ax.legend(loc=0);
ax2 = ax.twinx()
l2 = plt.plot(r, dsig,  'k--', linewidth=2, label='Absolute error')
plt.ylabel('Absolute error (MPa)')
plt.ylim(-1,1)
plt.xlim(a,3*a)
plt.title('Tangent stress error')
ls = l1+l2
ax.legend(ls, [ ls[0].get_label(), ls[1].get_label() ])

plt.show()
