import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import h5py
from numpy import genfromtxt
import edpWellbore_analyticalResults as edpAnal
import elastoPlasticWellboreAnalyticalSolutions as epwAnal

def stressRotation(stress, phi_x):
    rotx = np.array([[np.cos(phi_x), np.sin(phi_x), 0.], [-np.sin(phi_x), np.cos(phi_x), 0.], [0., 0., 1.]])
    return np.dot(np.dot(np.transpose(rotx), stress), rotx)

def main():
	# Get stress, time and element center from GEOSX results
	stress_field_name = 'rock_stress'
	hf_stress = h5py.File('stressHistory_rock.hdf5', 'r')
	time = np.array( hf_stress.get(stress_field_name + ' Time') )
	center = np.array( hf_stress.get(stress_field_name + ' elementCenter') )
	stress = np.array( hf_stress.get(stress_field_name) )
	
	# Get the deformed wellbore radius
	hf_disp = h5py.File("displacementHistory.hdf5", 'r')
	displacement = np.array( hf_disp.get('totalDisplacement') )
	node_position = np.array( hf_disp.get('totalDisplacement ReferencePosition') )
	
	a0 = node_position[0, 0, 0]
	da = displacement[:, 0, 0]
	a = a0 + da
	
	# Compute total stress
	stress_xx = ( sum(stress[26,:,6*i+0] for i in range(8)) )/8.0
	stress_yy = ( sum(stress[26,:,6*i+1] for i in range(8)) )/8.0
	stress_zz = ( sum(stress[26,:,6*i+2] for i in range(8)) )/8.0
	stress_yz = ( sum(stress[26,:,6*i+3] for i in range(8)) )/8.0
	stress_xz = ( sum(stress[26,:,6*i+4] for i in range(8)) )/8.0
	stress_xy = ( sum(stress[26,:,6*i+5] for i in range(8)) )/8.0
	
	# Coordinate of elemnt center
	nElements = center.shape[1]
	xCoord = center[0, :, 0]
	yCoord = center[0, :, 1]

	rCoord = np.sqrt( xCoord*xCoord + yCoord*yCoord )

	# Compute stress components in cylindrical coordinate system
	stress_rr = np.zeros(stress_xx.shape)
	stress_tt = np.zeros(stress_xx.shape)

	for idx_elem in range(stress.shape[1]):
		stressMatrix_cartesian = np.array([[stress_xx[idx_elem],stress_xy[idx_elem],stress_xz[idx_elem]],\
                         				   [stress_xy[idx_elem],stress_yy[idx_elem],stress_yz[idx_elem]],\
                         				   [stress_xz[idx_elem],stress_yz[idx_elem],stress_zz[idx_elem]]])

		if(yCoord[idx_elem] != 0):
			phi_x = np.arctan( xCoord[idx_elem]/yCoord[idx_elem] )
		else:
			phi_x = 0

		stressMatrix_cylindirical = stressRotation(stressMatrix_cartesian, phi_x)
		stress_rr[idx_elem] = stressMatrix_cylindirical[1][1]
		stress_tt[idx_elem] = stressMatrix_cylindirical[0][0]

	# Stress invariants at wellbore surface
	stress_xx_wellboresurface = ( sum(stress[:,0,6*i+0] for i in range(8)) )/8.0
	stress_yy_wellboresurface = ( sum(stress[:,0,6*i+1] for i in range(8)) )/8.0
	stress_zz_wellboresurface = ( sum(stress[:,0,6*i+2] for i in range(8)) )/8.0
	stress_yz_wellboresurface = ( sum(stress[:,0,6*i+3] for i in range(8)) )/8.0
	stress_xz_wellboresurface = ( sum(stress[:,0,6*i+4] for i in range(8)) )/8.0
	stress_xy_wellboresurface = ( sum(stress[:,0,6*i+5] for i in range(8)) )/8.0
	p = (stress_xx_wellboresurface + stress_yy_wellboresurface + stress_zz_wellboresurface)/3.0
	q = np.sqrt(1.5) * np.sqrt( (stress_xx_wellboresurface-p)**2.0 + (stress_yy_wellboresurface-p)**2.0 + (stress_zz_wellboresurface-p)**2.0 + \
        stress_xy_wellboresurface**2.0 + stress_xz_wellboresurface**2.0 + stress_yz_wellboresurface**2.0 )

	# Analytical results 

	# At pw = 2 (MPa)
	a0_a_ratio = 1.075 # This ratio corresponds to pw = 2 (MPa)
	r_ep_analytic,srr_ep_analytic,stt_ep_analytic,szz_ep_analytic,r_elas_analytic,srr_elas_analytic,stt_elas_analytic,szz_elas_analytic,\
    pw_analytic,p_wellsurface_analytic,q_wellsurface_analytic = edpAnal.EDP(a0_a_ratio)
		
	# Compute pressure at wellbore surface
	list_a0_a_ratio = np.arange(1.005,1.1,0.001)
	list_pw_analytic = []
	list_p_wellsurface_analytic = []
	list_q_wellsurface_analytic = []
	
	for a0_a_ratio_val in list_a0_a_ratio:
		tmp = edpAnal.EDP(a0_a_ratio_val)
		pw_analytic = tmp[8]
		p_wellsurface_analytic = tmp[9]
		q_wellsurface_analytic = tmp[10]
		list_pw_analytic.append(pw_analytic)
		list_p_wellsurface_analytic.append(p_wellsurface_analytic)
		list_q_wellsurface_analytic.append(q_wellsurface_analytic)

	# The case of elastic material for comparision
	sv = 15e6
	sh = 11.25e6
	pw = 2e6
	r_elas_assumed,srr_elas_assumed,stt_elas_assumed,szz_elas_assumed = epwAnal.solution_elastic(sh,sv,1.0,pw)

	# Plots	
	plt.figure(figsize=(15,10))

	# Plot GEOS results versus analytical results for elasto-plastic material
	plt.subplot(2,2,1)
	plt.plot( rCoord/rCoord[0],-stress_rr,'r+',label=r'$\sigma_{rr}$: GEOS')
	plt.plot(r_ep_analytic,srr_ep_analytic, 'r', label=r'$\sigma_{rr}$: Analytic Plastic')
	plt.plot(r_elas_analytic,srr_elas_analytic, 'r')

	plt.plot( rCoord/rCoord[0],-stress_tt,'b+',label=r'$\sigma_{\theta\theta}$: GEOS')
	plt.plot(r_ep_analytic,stt_ep_analytic, 'b', label=r'$\sigma_{\theta\theta}$: Analytic Plastic')
	plt.plot(r_elas_analytic,stt_elas_analytic, 'b')

	plt.plot( rCoord/rCoord[0],-stress_zz,'g+',label=r'$\sigma_{zz}$: GEOS')
	plt.plot(r_ep_analytic,szz_ep_analytic, 'g', label=r'$\sigma_{zz}$: Analytic Plastic')
	plt.plot(r_elas_analytic,szz_elas_analytic, 'g')

	plt.plot([r_ep_analytic[0],r_ep_analytic[0]],[0.0,20e6],'k--', label='Elastic-Plastic boundary')

	plt.ylabel('Stresses [Pa]')
	plt.xlabel(r'Normalized radial coordinate, $r/a$')
	plt.xlim(1.0,10.0)
	plt.ylim(0.0,20e6)
	plt.yticks(np.linspace(0.0,20e6,11))
	plt.xscale('log', subs=range(2,10))
	plt.legend()

	# Plot GEOS results for elasto-plastic material versus analytical results for elastic material
	plt.subplot(2,2,2)
	plt.plot( rCoord/rCoord[0],-stress_rr,'r+',label=r'$\sigma_{rr}$: GEOS')
	plt.plot(r_elas_assumed,srr_elas_assumed, 'r', label=r'$\sigma_{rr}$: Analytic Elastic')

	plt.plot( rCoord/rCoord[0],-stress_tt,'b+',label=r'$\sigma_{\theta\theta}$: GEOS')
	plt.plot(r_elas_assumed,stt_elas_assumed, 'b', label=r'$\sigma_{\theta\theta}$: Analytic Elastic')

	plt.plot( rCoord/rCoord[0],-stress_zz,'g+',label=r'$\sigma_{zz}$: GEOS')
	plt.plot(r_elas_assumed,szz_elas_assumed, 'g', label=r'$\sigma_{zz}$: Analytic Elastic')
	
	plt.ylabel('Stresses [Pa]')
	plt.xlabel(r'Normalized radial coordinate, $r/a$')
	plt.xlim(1.0,10.0)
	plt.ylim(0.0,20e6)
	plt.yticks(np.linspace(0.0,20e6,11))
	plt.xscale('log', subs=range(2,10))
	plt.legend()


	# Plot p-q plan
	initialFrictionAngle = 15.27 # deg
	finalFrictionAngle = 23.05   # deg
	initialFrictionAngle *= np.pi/180.0 # converted to rad
	finalFrictionAngle *= np.pi/180.0 # converted to rad
	param_b_i = edpAnal.compute_param_b(initialFrictionAngle)
	param_b_f = edpAnal.compute_param_b(finalFrictionAngle)
	p_yieldSurface = np.arange(0,15e6,1)
	q_yieldSurface_i = p_yieldSurface*param_b_i
	q_yieldSurface_f = p_yieldSurface*param_b_f

	plt.subplot(2,2,3)
	plt.plot(-p,q,'ko',label='GEOS')
	plt.plot(list_p_wellsurface_analytic,list_q_wellsurface_analytic,'b',label='Analytic')
	plt.plot(p_yieldSurface,q_yieldSurface_i,'g--',label='Initial yield surface')
	plt.plot(p_yieldSurface,q_yieldSurface_f,'r--',label='Final yield surface')

	plt.xlabel('p (Pa)')
	plt.ylabel('q (Pa)')
	plt.xlim(0,15e6)
	plt.ylim(0,10e6)
	plt.legend()

	# Plot a0/a versus pw
	pw = stress_xx_wellboresurface
	plt.subplot(2,2,4)
	plt.plot(a0/a, -pw,'ko',label='GEOS')
	plt.plot(list_a0_a_ratio,list_pw_analytic,'b',label='Analytic')
	plt.xlabel('a0/a')
	plt.ylabel('pw (Pa)')
	plt.xlim(1.0,1.1)
	plt.ylim(0,12e6)
	plt.legend()

	plt.show()

if __name__ == "__main__":
	main()
	
