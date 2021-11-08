.. _kgdToughnessDominated:


#######################################################
Toughness dominated KGD hydraulic fracture
#######################################################

------------------------------------------------------------------
Description of the case
------------------------------------------------------------------

In this example, we consider a plane-strain hydraulic fracture propagating in an infinite, homogeneous and elastic medium,  due to fluid injection at a rate :math:`Q_0` during a period from 0 to :math:`t_{max}`. Two dimensional KGD fracture is characterized as a vertical fracture with a rectangle-shaped cross section. For verification purpose, the presented numerical model is restricted to the assumptions used to analytically solve this problem `(Bunger et al., 2005) <https://link.springer.com/article/10.1007%2Fs10704-005-0154-0>`__. Vertical and impermeable fracture surface is assumed, which eliminate the effect of fracture plane inclination and fluid leakoff. The injected fluid flows within the fracture, which is assumed to be governed by the lubrication equation resulting from the mass conservation and the Poiseuille law. Fracture profile is related to fluid pressure distribution, which is mainly dictated by fluid viscosity :math:`\mu`. In addition, fluid pressure contributes to the fracture development through the mechanical deformation of the solid matrix, which is characterized by rock elastic properties, including the Young modulus :math:`E`, and the Poisson ratio :math:`\nu`.

For toughness-dominated fractures, more work is spent to split the intact rock than that applied to move the fracturing fluid. To make the case identical to the toughness dominated asymptotic solution, incompressible fluid with an ultra-low viscosity of 0.001 cp and medium rock toughness should be defined. Fracture is propagating with the creation of new surface if the stress intensity factor exceeds rock toughness :math:`K_{Ic}`.

In toughness-storage dominated regime, asymptotic solutions of the fracture length :math:`\ell`, the net pressure :math:`p_0` and the fracture aperture :math:`w_0` at the injection point for the KGD fracture are provided by `(Bunger et al., 2005) <https://link.springer.com/article/10.1007%2Fs10704-005-0154-0>`__:

.. math::
   \ell = 0.9324 X^{ -1/6 } (\frac{ E_p Q_0^3 }{ 12\mu })^{ 1/6 } t^{ 2/3 }

   w_0^2 = 0.5 X^{ 1/2 } (\frac{ 12\mu Q_0 }{ E_p })^{ 1/2 } \ell

   w_0 p_0 = 0.125 X^{ 1/2 } (12\mu Q_0 E_p)^{ 1/2 }

where the plane modulus :math:`E_p` is defined by

.. math:: E_p = \frac{ E }{ 1-\nu^2 }

and the term :math:`X` is given as:

.. math:: X = \frac{ 256 }{ 3 \pi^2 } \frac{ K_{Ic}^4 }{ \mu Q_0 {E_p}^3 }


**Input file**

The input xml files for this test case are located at:

.. code-block:: console

  inputFiles/hydraulicFracturing/kgdToughnessDominated_base.xml

and

.. code-block:: console

  inputFiles/hydraulicFracturing/kgdToughnessDominated_benchmark.xml

The corresponding integrated test with coarser mesh and smaller injection duration is also prepared:

.. code-block:: console

  inputFiles/hydraulicFracturing/kgdToughnessDominated_Smoke.xml

-----------------------------------------------------------
Mechanics solvers
-----------------------------------------------------------

The solver ``SurfaceGenerator`` defines rock toughness :math:`K_{Ic}` as:

.. literalinclude:: ../../../../../../inputFiles/hydraulicFracturing/kgdToughnessDominated_base.xml
  :language: xml
  :start-after: <!-- Sphinx_Solvers_SurfaceGenerator -->
  :end-before:  <!-- Sphinx_Solvers_SurfaceGenerator_End -->

Rock and fracture deformation are modeled by the solid mechanics solver ``SolidMechanicsLagrangianSSLE``. In this solver, we define ``targetRegions`` that includes both the continuum region and the fracture region. The name of the contact constitutive behavior is also specified in this solver by the ``contactRelationName``, besides the ``solidMaterialNames``.

.. literalinclude:: ../../../../../../inputFiles/hydraulicFracturing/kgdToughnessDominated_base.xml
  :language: xml
  :start-after: <!-- Sphinx_Solvers_SolidMechanicsLagrangianSSLE -->
  :end-before:  <!-- Sphinx_Solvers_SolidMechanicsLagrangianSSLE_End -->

The single phase fluid flow inside the fracture is solved by the finite volume method in the solver ``SinglePhaseFVM`` as:

.. literalinclude:: ../../../../../../inputFiles/hydraulicFracturing/kgdToughnessDominated_base.xml
  :language: xml
  :start-after: <!-- Sphinx_Solvers_SinglePhaseFVM -->
  :end-before:  <!-- Sphinx_Solvers_SinglePhaseFVM_End -->

All these elementary solvers are combined in the solver ``Hydrofracture`` to model the coupling between fluid flow within the fracture, rock deformation, fracture opening/closure and propagation. A fully coupled scheme is defined by setting a flag ``FIM`` for ``couplingTypeOption``.

.. literalinclude:: ../../../../../../inputFiles/hydraulicFracturing/kgdToughnessDominated_base.xml
  :language: xml
  :start-after: <!-- Sphinx_Solvers_Hydrofracture -->
  :end-before:  <!-- Sphinx_Solvers_Hydrofracture_End -->

-----------------------------------------------------------
The constitutive laws
-----------------------------------------------------------

The constitutive law ``CompressibleSinglePhaseFluid`` defines the default and reference fluid viscosity, compressibility and density. For this toughness dominated example, ultra low fluid viscosity is used:

.. literalinclude:: ../../../../../../inputFiles/hydraulicFracturing/kgdToughnessDominated_base.xml
  :language: xml
  :start-after: <!-- Sphinx_Constitutive_CompressibleSinglePhaseFluid -->
  :end-before:  <!-- Sphinx_Constitutive_CompressibleSinglePhaseFluid_End -->

The isotropic elastic Young modulus and Poisson ratio are defined in the ``ElasticIsotropic`` block. The density of rock defined in this block is useless, as gravity effect is ignored in this example.

.. literalinclude:: ../../../../../../inputFiles/hydraulicFracturing/kgdToughnessDominated_base.xml
  :language: xml
  :start-after: <!-- Sphinx_Constitutive_ElasticIsotropic -->
  :end-before:  <!-- Sphinx_Constitutive_ElasticIsotropic_End -->

-----------------------------------------------------------
Mesh
-----------------------------------------------------------

Internal mesh generator is used to generate the geometry of this example. The domain size is large enough comparing to the final size of the fracture. A sensitivity analysis has shown that the domain size in the direction perpendicular to the fracture plane, i.e. x-axis, must be at least ten times of the final fracture half-length to minimize the boundary effect. However, smaller size along the fracture plane, i.e. y-axis, of only two times the fracture half-length is good enough. It is also important to note that at least two layers are required in z-axis to ensure a good match between the numerical results and analytical solutions, due to the node based fracture propagation criterion. Also in x-axis, bias parameter ``xBias`` is added for optimizing the mesh by refining the elements near the fracture plane.

.. literalinclude:: ../../../../../../inputFiles/hydraulicFracturing/kgdToughnessDominated_base.xml
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

The initial fracture is defined by a nodeset occupying a small area where the KGD fracture starts to propagate:

.. literalinclude:: ../../../../../../inputFiles/hydraulicFracturing/kgdToughnessDominated_benchmark.xml
  :language: xml
  :start-after: <!-- Sphinx_Geometry_InitFracture -->
  :end-before:  <!-- Sphinx_Geometry_InitFracture_End -->

This initial ``ruptureState`` condition must be specified for this area in the following ``FieldSpecification`` block:

.. literalinclude:: ../../../../../../inputFiles/hydraulicFracturing/kgdToughnessDominated_base.xml
  :language: xml
  :start-after: <!-- Sphinx_FieldSpecifications_InitFracture -->
  :end-before:  <!-- Sphinx_FieldSpecifications_InitFracture_End -->

-----------------------------------------------------------
Defining the fracture plane
-----------------------------------------------------------

The plane within which the KGD fracture propagates is predefined to reduce the computational cost. The fracture plane is outlined by a separable nodeset by the following initial ``FieldSpecification`` condition:

.. literalinclude:: ../../../../../../inputFiles/hydraulicFracturing/kgdToughnessDominated_benchmark.xml
  :language: xml
  :start-after: <!-- Sphinx_Geometry_FracturePlane -->
  :end-before:  <!-- Sphinx_Geometry_FracturePlane_End -->

.. literalinclude:: ../../../../../../inputFiles/hydraulicFracturing/kgdToughnessDominated_base.xml
  :language: xml
  :start-after: <!-- Sphinx_FieldSpecifications_FracturePlane -->
  :end-before:  <!-- Sphinx_FieldSpecifications_FracturePlane_End -->

-----------------------------------------------------------
Defining the injection rate
-----------------------------------------------------------

Fluid is injected into a sub-area of the initial fracture. Only half of the injection rate is defined in this boundary condition because only half-wing of the KGD fracture is modeled regarding its symmetry. Hereby, the mass injection rate is actually defined, instead of the volume injection rate. More precisely, the value given for ``scale`` is :math:`Q_0 \rho_f/2` (not :math:`Q_0 /2`).

.. literalinclude:: ../../../../../../inputFiles/hydraulicFracturing/kgdToughnessDominated_base.xml
  :language: xml
  :start-after: <!-- Sphinx_SourceFlux_InjSource -->
  :end-before:  <!-- Sphinx_SourceFlux_InjSource_End -->

The parameters used in the simulation are summarized in the following table.

  +----------------+-----------------------+--------------------+-------------------+
  | Symbol         | Parameter             | Units              | Value             |
  +================+=======================+====================+===================+
  | :math:`Q_0`    | Injection rate        | [m\ :sup:`3`/s]    | 10\ :sup:`-4`     |
  +----------------+-----------------------+--------------------+-------------------+
  | :math:`E`      | Young's modulus       | [GPa]              | 30                |
  +----------------+-----------------------+--------------------+-------------------+
  | :math:`\nu`    | Poisson's ratio       | [ - ]              | 0.25              |
  +----------------+-----------------------+--------------------+-------------------+
  | :math:`\mu`    | Fluid viscosity       | [Pa.s]             | 10\ :sup:`-6`     |
  +----------------+-----------------------+--------------------+-------------------+
  | :math:`K_{Ic}` | Rock toughness        | [MPa.m\ :sup:`1/2`]| 1                 |
  +----------------+-----------------------+--------------------+-------------------+

---------------------------------
Inspecting results
---------------------------------

Fracture propagation during the fluid injection period is shown in the figure below.

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

			inletAperture = np.sqrt( 0.5 * X**( 0.5 ) * ( 12.0 * mu * Q0 / Ep )**( 0.5 ) * halfLength )

			inletPressure = 0.125 * X**( 0.5 ) * (12.0 * mu * Q0 * Ep)**( 0.5 ) / inletAperture

			return [ halfLength, inletAperture , inletPressure ]

	def getParametersFromXML( xmlFilePath ):
		tree = ElementTree.parse(xmlFilePath + "_benchmark.xml")

		maxTime = float(tree.find('Events').get('maxTime'))

		tree = ElementTree.parse(xmlFilePath + "_base.xml")

		elasticParam = tree.find('Constitutive/ElasticIsotropic')

		youngModulus = float(elasticParam.get('defaultYoungModulus'))
		poissonRatio = float(elasticParam.get('defaultPoissonRatio'))

		viscosity = float( tree.find('Constitutive/CompressibleSinglePhaseFluid').get('defaultViscosity') )

		toughness = float( tree.find('Solvers/SurfaceGenerator').get('rockToughness') )

		fluidDensity = float( tree.find('Constitutive/CompressibleSinglePhaseFluid').get('defaultDensity') )

		injectionRate = -2.0 * float( tree.find('FieldSpecifications/SourceFlux').get('scale') ) / fluidDensity

		return [ maxTime, youngModulus, poissonRatio, toughness, viscosity, injectionRate ]


	def main():
		xmlFilePathPrefix = "../../../../../../inputFiles/hydraulicFracturing/kgdToughnessDominated"

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
