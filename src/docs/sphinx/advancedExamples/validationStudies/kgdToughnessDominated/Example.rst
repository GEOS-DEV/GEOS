.. _kgdToughnessDominated:


#######################################################
Toughness dominated KGD problem
#######################################################

------------------------------------------------------------------
Description of the case
------------------------------------------------------------------

We consider a plane-strain hydraulic fracture propagating in an infinite homogeneous medium due to the injection of a fluid into the fracture with a rate :math:`Q_0` in a period from 0 to :math:`t_{max}`. The injected fluid flows in the fracture with respect to the lubrication equation resulting from the mass conservation and the Poiseuille law. It relates the fracture aperture to the fluid pressure via the fluid viscosity :math:`\mu`. On the other side, the fluid pressure is linked to the fracture aperture through the mechanical deformation of the solid matrix that is characterized by rock elastic properties: the Young modulus :math:`E` and Poisson ratio :math:`\nu` or the anisotropic elastic stiffnesses. The fracture growth is controlled by comparing the stress intensity factor and rock toughness :math:`K_{Ic}`.

Analytical results exist for some cases. We consider for example the toughness dominated regime and isotropic homogeneous impermeable rock, the exact analytical result of the fracture length :math:`\ell`, the net pressure :math:`p_0` and the fracture aperture :math:`w_0` at the injection point can be obtained by following closed-form solutions `(Bunger et al., 2005) <https://link.springer.com/article/10.1007%2Fs10704-005-0154-0>`__: 

.. math::
   \ell = 0.9324 X^{ -1/6 } (\frac{ E_p Q_0^3 }{ 12\mu })^{ 1/6 } t^{ 2/3 }

   w_0^2 = 0.5 X^{ 1/2 } (\frac{ 12\mu Q_0 }{ E_p })^{ 1/2 } \ell

   w_0 p_0 = 0.125 X^{ 1/2 } (12\mu Q_0 E_p)^{ 1/2 }

where the plane modulus :math:`E_p` is defined by

.. math:: E_p = \frac{ E }{ 1-\nu^2 }

We noted also:

.. math:: X = \frac{ 256 }{ 3 \pi^2 } \frac{ K_{Ic}^4 }{ \mu Q_0 {E_p}^3 }


**Input file**

The xml input files for this test case are located at:

.. code-block:: console

  inputFiles/multiphysics/kgdToughnessDominated_Base.xml

and

.. code-block:: console

  inputFiles/multiphysics/kgdToughnessDominated_Example.xml

the corresponding integrated test with coarser mesh and smaller injection time is defined by the xml:

.. code-block:: console

  inputFiles/multiphysics/kgdToughnessDominated_Smoke.xml

-----------------------------------------------------------
The physic solvers
-----------------------------------------------------------

The solver ``SurfaceGenerator`` define rock toughness :math:`K_{Ic}` as: 

.. literalinclude:: ../../../../../../inputFiles/multiphysics/kgdToughnessDominated_Base.xml
  :language: xml
  :start-after: <!-- Sphinx_Solvers_SurfaceGenerator -->
  :end-before:  <!-- Sphinx_Solvers_SurfaceGenerator_End -->

Rock and fracture deformation are modeled by the solid mechanics solver ``SolidMechanicsLagrangianSSLE``. In this solver we define a ``targetRegions`` that includes both the continuum region and the fracture region. The name of the contact constitutive behavior is also declared in this solver by the ``contactRelationName`` beside the ``solidMaterialNames``.

.. literalinclude:: ../../../../../../inputFiles/multiphysics/kgdToughnessDominated_Base.xml
  :language: xml
  :start-after: <!-- Sphinx_Solvers_SolidMechanicsLagrangianSSLE -->
  :end-before:  <!-- Sphinx_Solvers_SolidMechanicsLagrangianSSLE_End -->

The single phase fluid flow inside the fracture is modeled by the finite volume method implemented in the solver ``SinglePhaseFVM`` as:   

.. literalinclude:: ../../../../../../inputFiles/multiphysics/kgdToughnessDominated_Base.xml
  :language: xml
  :start-after: <!-- Sphinx_Solvers_SinglePhaseFVM -->
  :end-before:  <!-- Sphinx_Solvers_SinglePhaseFVM_End -->

All these elementary solvers are wrapped up in the solver ``Hydrofracture`` to model the coupling between fluid flow inside the fracture, rock deformation and fracture opening/closer and propagation. A fully coupled scheme is defined by setting a flag ``FIM`` for ``couplingTypeOption``.  

.. literalinclude:: ../../../../../../inputFiles/multiphysics/kgdToughnessDominated_Base.xml
  :language: xml
  :start-after: <!-- Sphinx_Solvers_Hydrofracture -->
  :end-before:  <!-- Sphinx_Solvers_Hydrofracture_End -->

-----------------------------------------------------------
The constitutive laws
-----------------------------------------------------------

The constitutive law ``CompressibleSinglePhaseFluid`` defines default and reference fluid's viscosity, compressibility and density. For this toughness dominated example, a negligible fluid viscosity is defined:
 
.. literalinclude:: ../../../../../../inputFiles/multiphysics/kgdToughnessDominated_Base.xml
  :language: xml
  :start-after: <!-- Sphinx_Constitutive_CompressibleSinglePhaseFluid -->
  :end-before:  <!-- Sphinx_Constitutive_CompressibleSinglePhaseFluid_End -->

The isotropic elastic Young modulus and Poisson ratio are defined by the ``ElasticIsotropic`` block. We note that the density of rock defined in this block is useless because gravity effect is ignored in this example.

.. literalinclude:: ../../../../../../inputFiles/multiphysics/kgdToughnessDominated_Base.xml
  :language: xml
  :start-after: <!-- Sphinx_Constitutive_ElasticIsotropic -->
  :end-before:  <!-- Sphinx_Constitutive_ElasticIsotropic_End -->

-----------------------------------------------------------
Mesh
-----------------------------------------------------------

Internal mesh generator is used to generate the geometry of this example. The domain size is large enough comparing to the final size of the fracture. A sensitivity analysis has shown that the domain size in the direction perpendicular to the fracture plane, i.e. x-axis, must be at least ten times of the final fracture half-legnth to minimize the numerical error. However, smaller size along the fracture plane, i.e. y-axis, of only two times the fracture half-length is good enough for the numerical convergence. It is also important to note that at least to layer is required in z-axis to ensure a good match between the numerical results and analytical solutions. This is because of the consideration of the node based fracture propagation criterion. Also in x-axis, bias parameter ``xBias`` is added for optimizing the mesh by reducing the size of the element nearby the fracture plane.

.. literalinclude:: ../../../../../../inputFiles/multiphysics/kgdToughnessDominated_Base.xml
  :language: xml
  :start-after: <!-- Sphinx_Mesh_InternalMesh -->
  :end-before:  <!-- Sphinx_Mesh_InternalMesh_End -->

.. figure:: mesh.png
   :align: center
   :width: 1000
   :figclass: align-center

-----------------------------------------------------------
Defining the initial fracture
-----------------------------------------------------------

The initial fracture is defined by a nodeset occupying a small area where the KGD fracture starts to propagate as:

.. literalinclude:: ../../../../../../inputFiles/multiphysics/kgdToughnessDominated_Example.xml
  :language: xml
  :start-after: <!-- Sphinx_Geometry_InitFracture -->
  :end-before:  <!-- Sphinx_Geometry_InitFracture_End -->
 
This initial ``ruptureState`` condition is defined for this area by the following ``FieldSpecification`` block:

.. literalinclude:: ../../../../../../inputFiles/multiphysics/kgdToughnessDominated_Base.xml
  :language: xml
  :start-after: <!-- Sphinx_FieldSpecifications_InitFracture -->
  :end-before:  <!-- Sphinx_FieldSpecifications_InitFracture_End -->

-----------------------------------------------------------
Defining the fracture plane
-----------------------------------------------------------

The plane within it the KGD fracture propagate is known, so it is convenient to verify the fracture propagation criterion only on that given plane to reduce the computational cost. The fracture plane is defined by a separable nodeset by the following initial ``FieldSpecification`` condition: 

.. literalinclude:: ../../../../../../inputFiles/multiphysics/kgdToughnessDominated_Example.xml
  :language: xml
  :start-after: <!-- Sphinx_Geometry_FracturePlane -->
  :end-before:  <!-- Sphinx_Geometry_FracturePlane_End -->

.. literalinclude:: ../../../../../../inputFiles/multiphysics/kgdToughnessDominated_Base.xml
  :language: xml
  :start-after: <!-- Sphinx_FieldSpecifications_FracturePlane -->
  :end-before:  <!-- Sphinx_FieldSpecifications_FracturePlane_End -->

-----------------------------------------------------------
Defining the injection rate
-----------------------------------------------------------

Fluid is injected into a sub-area of the initial fracture. It should be remarked that only half of the injection rate is defined in this boundary condition because only half-wing of the KGD fracture is modeled regarding its symmetry. It is important to note also that the mass injection rate is actually defined, not the volume injection rate. More precisely, the value given for ``scale`` is :math:`Q_0 \rho_f/2` not :math:`Q_0 /2`.

.. literalinclude:: ../../../../../../inputFiles/multiphysics/kgdToughnessDominated_Base.xml
  :language: xml
  :start-after: <!-- Sphinx_FieldSpecifications_InjSource -->
  :end-before:  <!-- Sphinx_FieldSpecifications_InjSource_End -->
	

The parameters used in the simulation are summarized in the following table.

  +----------------+-----------------------+------------------+-------------------+
  | Symbol         | Parameter             | Units            | Value             |
  +================+=======================+==================+===================+
  | :math:`Q_0`    | Injection rate        | [m:sup:`3`/s]    | 10:sup:`-4`       |
  +----------------+-----------------------+------------------+-------------------+
  | :math:`E`      | Young's modulus       | [GPa]            | 30                |
  +----------------+-----------------------+------------------+-------------------+
  | :math:`\nu`    | Poisson's ratio       | [ - ]            | 0.25              |
  +----------------+-----------------------+------------------+-------------------+
  | :math:`\mu`    | Fluid viscosity       | [Pa.s]           | 10:sup:`-6`       |
  +----------------+-----------------------+------------------+-------------------+
  | :math:`K_{Ic}` | Rock toughness        | [MPa.m:sup:`1/2`]| 1                 |
  +----------------+-----------------------+------------------+-------------------+


---------------------------------
Inspecting results
---------------------------------

Fracture propagation during time is shown in the figure below. 

.. figure:: propagation.gif
   :align: center
   :width: 1000
   :figclass: align-center

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
		
			halfLength = 0.9324 * X**( -1.0/6.0 ) * ( ( Ep * Q0**3.0 ) / ( 12.0 * mu ) )**( 1.0/6.0 ) * t**( 2.0/3.0 )
			print( 0.9324 * X**( -1.0/6.0 ) * ( ( Ep * Q0**3.0 ) / ( 12.0 * mu ) )**( 1.0/6.0 ) )
			inletAperture = np.sqrt( 0.5 * X**( 0.5 ) * ( 12.0 * mu * Q0 / Ep )**( 0.5 ) * halfLength )

			inletPressure = 0.125 * X**( 0.5 ) * (12.0 * mu * Q0 * Ep)**( 0.5 ) / inletAperture
	
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
		xmlFilePathPrefix = "../../../../../../inputFiles/multiphysics/kgdToughnessDominated"

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

This concludes the toughness dominated KGD example.
For any feedback on this example, please submit a `GitHub issue on the project's GitHub page <https://github.com/GEOSX/GEOSX/issues>`_.
