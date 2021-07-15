from scipy.special import kn
from math import factorial, floor
import numpy as np
import matplotlib.pyplot as plt


def FFunction(s, R):
    P = kn(0,R * (s ** 0.5)) / (s * kn(0,s ** 0.5))
    U = -(R * kn(1,R * (s ** 0.5)) - kn(1,s ** 0.5)) / (s * R * (s ** 0.5) * kn(0,s ** 0.5))
    Srr = -(-R * kn(1,R * (s ** 0.5)) + kn(1,s ** 0.5)) / (R ** 2. * (s ** 1.5) * kn(0,s ** 0.5))
    Stt = -P - Srr
   
    return [P, U, Srr, Stt]


def Vfunction(i, N):
    
	sum1 = 0.
	kmin = int(floor((i + 1.) / 2.))

	kmax = min(i, N)
    
	for k in range(kmin, kmax+1):
		sum1 = sum1 + (1.*(k**N) * factorial(2 * k) / (factorial(N - k) * factorial(k) * factorial(k - 1) * factorial(i - k) * factorial(2 * k - i)))


	return ((-1.)**(N + i)) * sum1

def StehfestTransform(t, R):
	N = 5
	sum1 = 0.
	sum2 = 0.
	sum3 = 0.
	sum4 = 0.
	for j in range(1,2*N+1):
		Lresult = FFunction(j * np.log(2.) / t, R)
		sum1 = sum1 + Vfunction(j, N) * Lresult[0]
		sum2 = sum2 + Vfunction(j, N) * Lresult[1]
		sum3 = sum3 + Vfunction(j, N) * Lresult[2]
		sum4 = sum4 + Vfunction(j, N) * Lresult[3]

	return [sum1* np.log(2.) / t, sum2* np.log(2.) / t, sum3* np.log(2.) / t, sum4* np.log(2.) / t] 

def parameters(t):
	eta =bBiot * (1.-2.*nu)/2./(1.-nu)
	kappa = permeability/fluidViscosity
	
	Ku = K + MBiot*(bBiot**2.)
	nuu = (3.*Ku - 2.*G) / (6.*Ku + 2.*G)

	BSkempton = 3.*(nuu-nu)/( bBiot*(1.-2.*nu)*(1.+nuu) )
	coefc = 2.*kappa*(BSkempton**2.)*G*(1.-nu)*( (1.+nuu)**2.)/9./(1.-nuu)/(nuu-nu)

	
	tstar=t*coefc/(a**2.)
	return tstar

t = 78.
a = 0.1
p0 = 10e6

phi = 0.3
MBiot = 15.8e9
bBiot = 0.771
E = 20.6e9
nu = 0.189

K = E/3./(1.-2.*nu)
G = E/2./(1.+nu)

permeability = 1e-17
fluidViscosity = 0.001

T = parameters(t)

listR = np.arange(1.,10.,0.01)
listP = []
listU = []
listSrr = []
listStt = []
for R in listR:

	result = StehfestTransform(T, R)
	listP.append(result[0])
	listU.append(result[1])
	listSrr.append(result[2])
	listStt.append(result[3])

listr = [R*a for R in listR]
listp = [1e-6*P*p0 for P in listP]
listUr = [1e6*U*a*p0*bBiot*(1.-2.*nu)/(2.*G*(1.-nu)) for U in listU]
listSigrrTot = [1e-6*Srr*bBiot*(1.-2.*nu)/(1.-nu)*p0 for Srr in listSrr]
listSigttTot = [1e-6*Stt*bBiot*(1.-2.*nu)/(1.-nu)*p0 for Stt in listStt]
listSigrr = [bBiot*val1+val2 for val1,val2 in zip(listp, listSigrrTot)]
listSigtt = [bBiot*val1+val2 for val1,val2 in zip(listp, listSigttTot)]

fig = plt.figure(figsize=[13,10])

plt.subplot(221)
plt.plot(listr,listp, 'k', linewidth=2, label='Analytic')


x, y = [], []
for line in open('pressure.curve', 'r'):
	if not (line.strip().startswith("#") or line.strip()==''):
		values = [float(s) for s in line.split()]
		x.append(values[0])
		y.append(values[1]*1e-6)
plt.plot(x,y,'ko', label='GEOSX result')
plt.ylabel('Pore pressure (MPa)')
plt.xlim(a, 10*a)
plt.legend()


plt.subplot(222)
plt.plot(listr,listUr, 'k', linewidth=2, label='Analytic')

x, y = [], []
for line in open('radialDisplacement.curve', 'r'):
	if not (line.strip().startswith("#") or line.strip()==''):
		values = [float(s) for s in line.split()]
		x.append(values[0])
		y.append(values[1]*1e6)
plt.plot(x,y,'ko', label='GEOSX result')
plt.ylabel('Radial displacement (um)')
plt.xlim(a, 10*a)


plt.subplot(223)
plt.plot(listr,listSigrr, 'k', linewidth=2, label='Analytic')

x, y = [], []
for line in open('radialStress.curve', 'r'):
	if not (line.strip().startswith("#") or line.strip()==''):
		values = [float(s) for s in line.split()]
		x.append(values[0])
		y.append(values[1]*1e-6)
plt.plot(x,y,'ko', label='GEOSX result')
plt.xlabel('r (m)')
plt.ylabel('Effective radial stress (MPa)')
plt.xlim(a, 10*a)


plt.subplot(224)
plt.plot(listr,listSigtt, 'k', linewidth=2, label='Analytic')

x, y = [], []
for line in open('tangentStress.curve', 'r'):
	if not (line.strip().startswith("#") or line.strip()==''):
		values = [float(s) for s in line.split()]
		x.append(values[0])
		y.append(values[1]*1e-6)
plt.plot(x,y,'ko', label='GEOSX result')
plt.xlabel('r (m)')
plt.ylabel('Effective tangent stress (MPa)')
plt.xlim(a, 10*a)

plt.show()
