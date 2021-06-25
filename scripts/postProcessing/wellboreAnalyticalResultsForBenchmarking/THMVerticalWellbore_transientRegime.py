from scipy.special import kv, iv, kn, erfc
from math import factorial, floor
import numpy as np
import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings("ignore")



def FFunction(s, R):

    
    T = T1/s* kn(0, r*(s/kappaT)** 0.5) / kn(0, ri*(s/kappaT)** 0.5)
    p = alpha_e*T1/(1-c/kappaT)/s * ( kn(0, r*(s/kappaT)** 0.5) / kn(0, ri*(s/kappaT)** 0.5) - kn(0, r*(s/c)** 0.5) / kn(0, ri*(s/c)** 0.5) )

    int_rT = -( T1*r*kn(1, r*(s/kappaT)**0.5) ) / ( s*(s/kappaT)**0.5 * kn(0, ri*(s/kappaT)**0.5) )
    int_rp = -alpha_e/(1-c/kappaT) * T1*r/s * ( kn(1, r*(s/kappaT)**0.5) / ((s/kappaT)**0.5 * kn(0, ri*(s/kappaT)**0.5)) - kn(1, r*(s/c)**0.5) / ((s/c)**0.5 * kn(0, ri*(s/c)**0.5)) )
    A2 = T1*ri/(1-c/kappaT)/G/s * ( ( eta*alpha_e + eta_d*(1-c/kappaT) ) * kn(1, ri*(s/kappaT)**0.5) / ((s/kappaT)**0.5 * kn(0, ri*(s/kappaT)**0.5)) - eta*alpha_e * kn(1, ri*(s/c)**0.5) / ((s/c)**0.5 * kn(0, ri*(s/c)**0.5)) )
    
    ur = eta/G/r*int_rp + eta_d/G/r*int_rT + A2/r
    sig_rr = -2.*eta/r**2. * int_rp - 2.*eta_d/r**2. * int_rT - 2.*G/r**2.*A2
    sig_tt = 2.*eta/r**2. * int_rp + 2.*eta_d/r**2. * int_rT + 2.*G/r**2.*A2 - 2.*eta*p - 2.*eta_d*T
    
    return [T,p,ur,sig_rr,sig_tt]

def StehfestTransform(t,r):
    N = 2
    sum1 = 0.
    sum2 = 0.
    sum3 = 0.
    sum4 = 0.
    sum5 = 0.
    for j in range(1,2*N+1):
        Lresult = FFunction(j * np.log(2.) / t, r)

        sum1 = sum1 + Vfunction(j, N) * Lresult[0]
        sum2 = sum2 + Vfunction(j, N) * Lresult[1]
        sum3 = sum3 + Vfunction(j, N) * Lresult[2]
        sum4 = sum4 + Vfunction(j, N) * Lresult[3]
        sum5 = sum5 + Vfunction(j, N) * Lresult[4]
       
    return [sum1* np.log(2.) / t, sum2* np.log(2.) / t, sum3* np.log(2.) / t, sum4* np.log(2.) / t, sum5* np.log(2.) / t] 

def Vfunction(i, N):
    
    sum1 = 0.
    kmin = int(floor((i + 1.) / 2.))

    kmax = min(i, N)

    for k in range(kmin, kmax+1):
        sum1 = sum1 + (1.*(k**N) * factorial(2 * k) / (factorial(N - k) * factorial(k) * factorial(k - 1) * factorial(i - k) * factorial(2 * k - i)))
    return ((-1.)**(N + i)) * sum1


# Early time approximation
def I1(r,x,t):
    XX = (r-ri)/(2.*(x*t)**0.5)
    return erfc(XX)#( (ri/r)**0.5 - (r-ri)/4.*(c*t/(ri*r**3.))**0.5*XX + (9.*ri**2-2.*ri*r-7.*r**2.)/(32.*ri**(3./2.)*r**(5./2.))*c*t*(1./4.+2.*XX**2.) )*erfc(XX) + ( (r-ri)/4.*(c*t/(ri*r**3.*np.pi))**0.5 - (9.*ri**2-2.*ri*r-7.*r**2.)/(32.*ri**(3./2.)*r**(5./2.)*np.pi**0.5)*2.*c*t*XX )*np.exp(-XX**2.) 


# Input parameters
ri=0.1
T1 = 60. 
#p1 = T1/(6.)*1e6

#Westerly Granite, see Tab 1 in Wang and Papamichos (1994)

# Notations by Cheng (2016)
alpha_e = 0.39e6 #equivalent to parameter cp,PaC
c =0.22e-4 #m2/s
kappaT = 0.5*c/0.23 #equivalent to parameter c0, the ratio given in Tab1 of Wang and Papamichos (1994) is for c/c0, not c0/c; Also, a factor 0.5 was added to match the results given by the Figs of Wang and Papamichos.

# Equivalent notations by Wang and Papamichos (1994)
cp = alpha_e
c0 = kappaT
Tw = T1
#pw = p1

# Additional parameters to calculate stresses and displacement
#Westerly Granite, see Tab 3.2 in Cheng (2016)
#alpha_d = 6e5 #N/m2/K
Ku = 39.3e9
nuu = 0.331
nu = 0.25
G = 3.*Ku*(1.-2.*nuu)/2./(1.+nuu)

bBiot = 1.0 #eta / ( (1.-2.*nu)/2./(1.-nu) ) Note: actual GEOSX version does only accept b=1
eta = bBiot * (1.-2.*nu)/2./(1.-nu) # 0.15
#eta_d = alpha_d * (1.-2.*nu)/2./(1.-nu)

# Prepare input parameter for GEOSX
tstar = 0.1
fluidViscosity = 1e-3
porosity = 0.01

a = ri 
t = tstar * a * a / c 

#bBiot = eta / ( (1.-2.*nu)/2./(1.-nu) )
BSkempton = 3.*(nuu-nu)/( bBiot*(1.-2.*nu)*(1.+nuu) )
kappa = c / ( 2.*(BSkempton**2.)*G*(1.-nu)*( (1.+nuu)**2.)/9./(1.-nuu)/(nuu-nu) )

permeability = kappa * fluidViscosity

K = 2.*G*(1.+nu) / 3. / (1.-2.*nu)
MBiot = (Ku - K) / bBiot / bBiot
fluidCompressibility = 1./MBiot/porosity

#alpha_m = alpha_d * 3. * (1.-2.*nu) /2./G/(1.+nu)
#alpha_f = ( ( cp*kappa/alpha_m/c - 2.*(nuu-nu)/BSkempton/(1.+nuu)/(1.-nu) ) / porosity + 1. ) * alpha_m
alpha_m = cp*kappa/c / ( 2.*(nuu-nu)/BSkempton/(1.+nuu)/(1.-nu) )
alpha_d = alpha_m / ( 3. * (1.-2.*nu) /2./G/(1.+nu) )
alpha_f = alpha_m
eta_d = alpha_d * (1.-2.*nu)/2./(1.-nu)

print("Bulk modulus: " + str(K))
print("Shear modulus: " + str(G))
print("Biot's coefficient: " + str(bBiot))

print("Permeability: " + str(permeability))
print("Porosity: " + str(porosity))
print("Fluid compressibility: " + str(fluidCompressibility))

print("Thermal diffusion coefficient: " + str(kappaT))
print("Thermal stress coefficient: " + str(alpha_d))
print("Thermal porosity coefficient: " + str( (alpha_f-alpha_m)*porosity ))
print("Diffusion time: " + str(t))

fig = plt.figure(figsize=[10,10])
list_tstar = [0.01,0.1,1.]
for iTime in range(len(list_tstar)):
    tstar = list_tstar[iTime]
    t = tstar*ri**2./c
    
    listr = np.arange(ri,4.*ri,0.01*ri)
    listT = []
    listp = []
    listUr = []
    listSigRR = []
    listSigTT = []
    for r in listr:
        result = StehfestTransform(t, r)
        listT.append(result[0])
        listp.append(result[1])
        listUr.append(result[2])
        listSigRR.append(result[3])
        listSigTT.append(result[4])
        
    
    if (iTime  == 0):
        line_type = 'k'
    if (iTime  == 1):
        line_type = 'k--'
    if (iTime  == 2):
        line_type = 'k-.'
    if (iTime  == 3):
        line_type = 'k:'
    if (iTime  == 4):
        line_type = 'b'
        

    #plt.subplot(111)
    plt.plot(listr/ri,listp, line_type, linewidth=2, label=r'$tstar$='+str(tstar))

    plt.ylabel('Ppore (Pa)')
    plt.xlabel('Normalized Radial Distance')
    plt.xlim(1, 4)
    #plt.ylim(-0.8,1.2)
    plt.legend()
	

# Geosx results
r_geosx, p_geosx = [], []
for line in open('pressure.curve', 'r'):
	if not (line.strip().startswith("#") or line.strip()==''):
		values = [float(s) for s in line.split()]
		rVal = values[0]
		pVal = values[1]
		r_geosx.append( rVal/ri )
		p_geosx.append( pVal )

plt.plot(r_geosx, p_geosx, 'ko', label='GEOSX result')

plt.show()
