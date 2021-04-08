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
plt.plot(r, dsig,  'k--', linewidth=2, label='Error')
plt.ylabel('Radial stress (MPa)')
plt.xlim(a,3*a)

plt.subplot(223)
plt.plot(r, dsigRelative, 'k', linewidth=2)
plt.xlabel('r (m)')
plt.ylabel('Relative error (%)')
plt.xlim(a,3*a)
plt.ylim(-10,10)


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
plt.plot(r, dsig,  'k--', linewidth=2, label='Error')
plt.ylabel('Tangent stress (MPa)')
plt.xlim(a,3*a)
plt.legend()

plt.subplot(224)
plt.plot(r, dsigRelative, 'k', linewidth=2)
plt.xlabel('r (m)')
plt.xlim(a,3*a)
plt.ylim(-10,10)

plt.show()
