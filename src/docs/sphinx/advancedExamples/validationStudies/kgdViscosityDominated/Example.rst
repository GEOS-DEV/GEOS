.. _kgdViscosityDominated:


#######################################################
Viscosity dominated KGD hydraulic fracture
#######################################################

------------------------------------------------------------------
Description of the case
------------------------------------------------------------------

This example is similar to the example :ref:`kgdViscosityDominated` except that the viscosity dominated regime is now considered. More work is spent to move the fracturing fluid than to split the intact rock. To make the case identical to the viscosity dominated asymptotic solution, a negligible rock toughness :math:`K_{Ic}` will be defined such as fracture is always propagating following fluid front. Asymptotic solutions of the fracture length :math:`\ell`, the net pressure :math:`p_0` and the fracture aperture :math:`w_0` at the injection point for the KGD fracture are provided by `(Adachi and Detournay, 2002) <https://onlinelibrary.wiley.com/doi/abs/10.1002/nag.213?casa_token=Jp094LF1vgQAAAAA:cTFllaSy9ze1t6Cf2-PzUq2k51ZkM0dChlMRcsvkKzB1ILTxCboU4-0LrKv8Jao-Sx5t3O1hLiX38UEk>`__: 

.. math::
   \ell = 0.6152 (\frac{ E_p Q_0^3 }{ 12\mu })^{ 1/6 } t^{ 2/3 }

   w_0^2 = 2.1 (\frac{ 12\mu Q_0 }{ E_p })^{ 1/2 } \ell

   w_0 p_0 = 0.62 (12\mu Q_0 E_p)^{ 1/2 }

where the plane modulus :math:`E_p` is defined by

.. math:: E_p = \frac{ E }{ 1-\nu^2 }

and the term :math:`X` is given as:

.. math:: X = \frac{ 256 }{ 3 \pi^2 } \frac{ K_{Ic}^4 }{ \mu Q_0 {E_p}^3 }


**Input file**

The input xml files for this test case are located at:

.. code-block:: console

  inputFiles/multiphysics/kgdViscosityDominated_Base.xml

and

.. code-block:: console

  inputFiles/multiphysics/kgdViscosityDominated_Example.xml

The corresponding integrated test with coarser mesh and smaller injection duration is also prepared:

.. code-block:: console

  inputFiles/multiphysics/kgdViscosityDominated_Smoke.xml

Fluid viscosity and rock toughness are defined by the xml blocks below. We note that setting an absolute zero value for the toughness may lead to a floating point issue. Therefore, a small value of :math:`K_{Ic}=10^{-6}` MPa should be considered.

.. literalinclude:: ../../../../../../inputFiles/multiphysics/kgdViscosityDominated_Base.xml
  :language: xml
  :start-after: <!-- Sphinx_Solvers_SurfaceGenerator -->
  :end-before:  <!-- Sphinx_Solvers_SurfaceGenerator_End -->

.. literalinclude:: ../../../../../../inputFiles/multiphysics/kgdViscosityDominated_Base.xml
  :language: xml
  :start-after: <!-- Sphinx_Constitutive_CompressibleSinglePhaseFluid -->
  :end-before:  <!-- Sphinx_Constitutive_CompressibleSinglePhaseFluid_End -->

A good agreement between GEOSX results and analytical solutions is shown in the comparison below:

.. plot::

	import numpy as np
	import matplotlib.pyplot as plt
	import xml.etree.ElementTree as ElementTree

	class kgd:
		def __init__(self, E, nu, KIc, mu, Q0, t):
			Ep = E / ( 1.0 - nu**2.0 )

			self.t  = t
			self.Q0 = Q0
			self.mu = mu
			self.Ep = Ep
			self.X  = 256.0 / ( 3.0 * np.pi**2.0 ) * ( KIc**4.0 ) / ( mu * Q0 * Ep**3.0 ) 


		def analyticalSolution(self):
			t  = self.t
			mu = self.mu
			Q0 = self.Q0
			Ep = self.Ep
			X  = self.X
		
			halfLength = 0.6152 * ( ( Ep * Q0**3.0 ) / ( 12.0 * mu ) )**( 1.0/6.0 ) * t**( 2.0/3.0 )
			
			inletAperture = np.sqrt( 2.1 * ( 12.0 * mu * Q0 / Ep )**( 0.5 ) * halfLength )

			inletPressure = 0.62 * (12.0 * mu * Q0 * Ep)**( 0.5 ) / inletAperture
	
			return [ halfLength, inletAperture , inletPressure ]

	def getParametersFromXML( xmlFilePath ):
		tree = ElementTree.parse(xmlFilePath + "_Example.xml")

		maxTime = float(tree.find('Events').get('maxTime'))

		tree = ElementTree.parse(xmlFilePath + "_Base.xml")

		elasticParam = tree.find('Constitutive/ElasticIsotropic')

		youngModulus = float(elasticParam.get('defaultYoungModulus'))
		poissonRatio = float(elasticParam.get('defaultPoissonRatio'))

		viscosity = float( tree.find('Constitutive/CompressibleSinglePhaseFluid').get('defaultViscosity') )

		toughness = float( tree.find('Solvers/SurfaceGenerator').get('rockToughness') )

		fluidDensity = float( tree.find('Constitutive/CompressibleSinglePhaseFluid').get('defaultDensity') )

		injectionRate = -2.0 * float( tree.find('FieldSpecifications/SourceFlux').get('scale') ) / fluidDensity

		return [ maxTime, youngModulus, poissonRatio, toughness, viscosity, injectionRate ]


	def main():
		xmlFilePathPrefix = "../../../../../../inputFiles/multiphysics/kgdViscosityDominated"

		tMax, E, nu, KIc, mu, Q0 = getParametersFromXML( xmlFilePathPrefix )

		t    = np.arange(0.01*tMax,tMax,0.01*tMax)

		model = kgd( E, nu, KIc, mu, Q0, t )
		halfLength, inletAperture , inletPressure = model.analyticalSolution()


		# GEOSX results
		t_geosx, halfLength_geosx = [], []
		for line in open('surfaceArea.curve', 'r'):
			if not (line.strip().startswith("#") or line.strip()==''):
				values = [float(s) for s in line.split()]
				t_geosx.append( values[0] )
				halfLength_geosx.append( values[1] / 2.0 )

		inletAperture_geosx = []
		for line in open('inletAperture.curve', 'r'):
			if not (line.strip().startswith("#") or line.strip()==''):
				values = [float(s) for s in line.split()]
				inletAperture_geosx.append( values[1] * 1e3 )

		inletPressure_geosx = []
		for line in open('inletPressure.curve', 'r'):
			if not (line.strip().startswith("#") or line.strip()==''):
				values = [float(s) for s in line.split()]
				inletPressure_geosx.append( values[1] / 1e6 )


		fig = plt.figure(figsize=[15,10])

		plt.subplot(221)
		plt.plot(t_geosx, halfLength_geosx, 'ko', label='GEOSX result')
		plt.plot(t, halfLength,  'k', linewidth=2, label='Analytic')
		plt.ylabel('Fracture half-length (m)')
		plt.xlabel('Injection time (s)')

		plt.subplot(222)
		plt.plot(t_geosx, inletAperture_geosx, 'ko', label='GEOSX result')
		plt.plot(t, inletAperture * 1e3,  'k', linewidth=2, label='Analytic')
		plt.ylabel('Inlet aperture (mm)')
		plt.xlabel('Injection time (s)')

		plt.subplot(223)
		plt.plot(t_geosx, inletPressure_geosx, 'ko', label='GEOSX result')
		plt.plot(t, inletPressure / 1e6,  'k', linewidth=2, label='Analytic')
		plt.ylabel('Inlet fluid pressure (MPa)')
		plt.xlabel('Injection time (s)')

		plt.legend()
		plt.show()

	if __name__ == "__main__":
		main()

------------------------------------------------------------------
To go further
------------------------------------------------------------------

**Feedback on this example**

This concludes the viscosity dominated KGD example.
For any feedback on this example, please submit a `GitHub issue on the project's GitHub page <https://github.com/GEOSX/GEOSX/issues>`_.
