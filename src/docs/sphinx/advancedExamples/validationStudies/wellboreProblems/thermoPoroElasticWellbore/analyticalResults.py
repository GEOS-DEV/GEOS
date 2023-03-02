import sys
sys.path.append('../')
import numpy as np
import h5py
import wellboreAnalyticalSolutions

def analyticalResults(t):
	ri = 0.1 # inner radius of the wellbore
	Ti = 100. # temperature applied on the inner surface of the wellbore
	
	# Notations and data by Cheng (2016)
	# Data of Rock salt

	K = 20.7e9 # drained bulk modulus
	G = 12.4e9 # drained shear modulus
	beta_d = 3*4e-5 # drained volumetric thermal expansion coefficient of the porous medium
	Ks = 23.5e9	# bulk modulus of the solid skeleton
	cf = 5e-10 # fluid compressibility
	muf = 1e-3 # fluid viscosity
	beta_f = 1e-5 # volumetric thermal expansion coefficient of fluid   

	porosity = 0.0010
	permeability = 1e-21 # intrinsic permeability 
	kT = 6.6 # thermal conducitivity
	volumetricHeatCapacity = 1.89e6 # product of specific heat and density
	
	# Compute parameters for the analytical solutions
	nu = (3.0*K-2.0*G)/(6.0*K+2.0*G)
	E = 2. * G * (1. + nu)
	alpha = 1.0 - K / Ks
	M = 1.0/( porosity*cf + (alpha-porosity)/Ks )
	M11 = K + 4. / 3. * G
	kappa = permeability/muf
	kappaT = kT/volumetricHeatCapacity
	c =  kappa * (M * M11 / (M11 + alpha**2. * M))	
	alpha_d = K*beta_d

	Ku = K + M*alpha*alpha
	S = (3.0*Ku + 4.0*G) /M /(3.0*K+4.0*G)

	beta_s = beta_d
	beta_v = porosity*(beta_f - beta_s)
	beta_e = beta_d*alpha + beta_v
	
	beta_c = beta_e - 3.0*alpha *K /(3.0*K+4.0*G) *beta_d
	alpha_e = beta_c/S
	
	# Compute analytical results
	r = np.arange(ri, 10. * ri, 0.01 * ri)

	T_thermal, p_thermal, ur_thermal, sig_rr_thermal, sig_tt_thermal, sig_zz_thermal = wellboreAnalyticalSolutions.inTime_thermal(t, r, ri, Ti, G, nu, alpha, alpha_d, c, kappaT, alpha_e)


	return [ r, T_thermal, p_thermal, ur_thermal, sig_rr_thermal, sig_tt_thermal, sig_zz_thermal ]

