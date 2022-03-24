.. _pennyFracToughnessDominated:


###########################################################
Toughness-Storage-Dominated Penny Shaped Hydraulic Fracture
###########################################################


**Context**

In this example, we simulate the growth of a radial hydraulic fracture in toughness-storage-dominated regime, a classic benchmark in hydraulic fracturing. The developed fracture is characterized as a planar fracture with an elliptical cross-section perpendicular to the fracture plane and a circular fracture tip. This problem is solved using the hydrofracture solver in GEOSX. The modeling predictions on the temporal evolutions of the fracture characteristics (length, aperture, and pressure) are verified against the analytical solutions `(Savitski and Detournay, 2002)  <https://www.sciencedirect.com/science/article/pii/S0020768302004924>`__. 


**Input file**

This example uses no external input files. Everything we need is contained within two GEOSX input files:

.. code-block:: console

  inputFiles/hydraulicFracturing/pennyShapedToughnessDominated_base.xml

.. code-block:: console

  inputFiles/hydraulicFracturing/pennyShapedToughnessDominated_benchmark.xml


A python script for post-processing the simulation results is also prepared:

.. code-block:: console

  src/docs/sphinx/advancedExamples/validationStudies/pennyFracToughnessDominated/pennyFracToughnessDominatedFigure.py


------------------------------------------------------------------
Description of the case
------------------------------------------------------------------

We model a radial fracture emerging from a point source and forming a perfect circular shape in an infinite, isotropic, and homogenous elastic domain. As with the KGD problem, we simplify the model to a radial fracture in a toughness-storage-dominated propagation regime. For toughness-dominated fractures, more work is spent on splitting the intact rock than on moving the fracturing fluid. Storage-dominated propagation occurs if most of the fracturing fluid is contained within the propagating fracture. In this analysis, incompressible fluid with ultra-low viscosity (:math:`0.001 cp`) and medium rock toughness (:math:`3.0 MPa { \sqrt{m} }`) are specified. In addition, an impermeable fracture surface is assumed to eliminate the effect of fluid leak-off. This way, the GEOSX simulations represent cases within the valid range of the toughness-storage-dominated assumptions. 

In this model, the injected fluid within the fracture follows the lubrication equation resulting from mass conservation and Poiseuille's law. The fracture propagates by creating new surfaces if the stress intensity factor exceeds the local rock toughness :math:`K_{IC}`. By symmetry, the simulation is reduced to a quarter-scale to save computational cost. For verification purposes, a plane strain deformation is considered in the numerical model. 

In this example, we set up and solve a hydraulic fracture model to obtain the temporal solutions of the fracture radius :math:`R`, the net pressure :math:`p_0` and the fracture aperture :math:`w_0` at the injection point for the penny-shaped fracture developed in this toughness-storage-dominated regime. The numerical predictions from GEOSX are then compared with the corresponding asymptotic solutions `(Savitski and Detournay, 2002)  <https://www.sciencedirect.com/science/article/pii/S0020768302004924>`__:

.. math:: R(t) = 0.8546 (\frac{ E_p^2 Q_0^2 t^2 }{ K_p^2 })^{ 1/5 }

.. math:: w_0(t) = 0.6537 (\frac{ K_p^4 Q_0 t }{ E_p^4 })^{ 1/5 } 

.. math:: p_0(t) = 0.3004 (\frac{ K_p^6 }{ E_p Q_0 t })^{ 1/5 } 

where the plane modulus :math:`E_p` is related to Young's modulus :math:`E` and Poisson's ratio :math:`\nu`:

.. math:: E_p = \frac{ E }{ 1-\nu^2 }

The term :math:`K_p` is proportional to the rock toughness :math:`K_{IC}`:

.. math:: K_p = \frac{ 8 }{ \sqrt{2 \pi} } K_{IC}
 


For this example, we focus on the ``Mesh``,
the ``Constitutive``, and the ``FieldSpecifications`` tags.

------------------------------------------------------------------
Mesh
------------------------------------------------------------------

The following figure shows the mesh used in this problem.


.. _problemSketchPennyShapedToughnessDominatedFig:
.. figure:: mesh.png
   :align: center
   :width: 500
   :figclass: align-center

   Generated mesh

We use the internal mesh generator to create a computational domain
(:math:`400\, m \, \times 400 \,  m \, \times 800 \, m`), as parametrized in the ``InternalMesh`` XML tag. 
The structured mesh contains 80 x 80 x 60 eight-node brick elements in the x, y, and z directions respectively. 
Such eight-node hexahedral elements are defined as ``C3D8`` elementTypes, and their collection forms a mesh
with one group of cell blocks named here ``cb1``. Local refinement is performed for the elements in the vicinity of the fracture plane. 

Note that the domain size in the direction perpendicular to the fracture plane, i.e. z-axis, must be at least ten times of the final fracture radius to minimize possible boundary effects. 


.. literalinclude:: ../../../../../../inputFiles/hydraulicFracturing/pennyShapedToughnessDominated_benchmark.xml
    :language: xml
    :start-after: <!-- SPHINX_MESH -->
    :end-before: <!-- SPHINX_MESH_END -->


The fracture plane is defined by a nodeset occupying a small region within the computation domain, where the fracture tends to open and propagate upon fluid injection:


.. literalinclude:: ../../../../../../inputFiles/hydraulicFracturing/pennyShapedToughnessDominated_benchmark.xml
  :language: xml
  :start-after: <!-- SPHINX_FRACPLANE -->
  :end-before: <!-- SPHINX_FRACPLANE_END -->


------------------------
Solid mechanics solver
------------------------

GEOSX is a multi-physics platform. Different combinations of
physics solvers available in the code can be applied
in different regions of the domain and be functional at different stages of the simulation.
The ``Solvers`` tag in the XML file is used to list and parameterize these solvers.

Three elementary solvers are combined in the solver ``Hydrofracture`` to model the coupling between fluid flow within the fracture, rock deformation, fracture deformation and propagation:


.. literalinclude:: ../../../../../../inputFiles/hydraulicFracturing/pennyShapedToughnessDominated_benchmark.xml
  :language: xml
  :start-after: <!-- SPHINX_HYDROFRACSOLVER -->
  :end-before: <!-- SPHINX_HYDROFRACSOLVER_END -->


- Rock and fracture deformation are modeled by the solid mechanics solver ``SolidMechanicsLagrangianSSLE``. In this solver, we define ``targetRegions`` that includes both the continuum region and the fracture region. The name of the contact constitutive behavior is specified in this solver by the ``contactRelationName``.


.. literalinclude:: ../../../../../../inputFiles/hydraulicFracturing/pennyShapedToughnessDominated_benchmark.xml
  :language: xml
  :start-after: <!-- SPHINX_MECHANICALSOLVER -->
  :end-before: <!-- SPHINX_MECHANICALSOLVER_END -->


- The single-phase fluid flow inside the fracture is solved by the finite volume method in the solver ``SinglePhaseFVM``.


.. literalinclude:: ../../../../../../inputFiles/hydraulicFracturing/pennyShapedToughnessDominated_benchmark.xml
  :language: xml
  :start-after: <!-- SPHINX_SINGLEPHASEFVM -->
  :end-before: <!-- SPHINX_SINGLEPHASEFVM_END -->


- The solver ``SurfaceGenerator`` defines the fracture region and rock toughness ``rockToughness="3.0e6"``. With ``nodeBasedSIF="1"``, a node-based Stress Intensity Factor (SIF) calculation is chosen for the fracture propagation criterion. 


.. literalinclude:: ../../../../../../inputFiles/hydraulicFracturing/pennyShapedToughnessDominated_benchmark.xml
  :language: xml
  :start-after: <!-- SPHINX_SURFACEGENERATOR -->
  :end-before: <!-- SPHINX_SURFACEGENERATOR_END -->


------------------------------
Constitutive laws
------------------------------

For this problem, a homogeneous and isotropic domain with one solid material is assumed. Its mechanical properties and associated fluid rheology are specified in the ``Constitutive`` section. 
``ElasticIsotropic`` model is used to describe the mechanical behavior of ``rock`` when subjected to fluid injection.
The single-phase fluid model ``CompressibleSinglePhaseFluid`` is selected to simulate the response of ``water`` upon fracture propagation.


.. literalinclude:: ../../../../../../inputFiles/hydraulicFracturing/pennyShapedToughnessDominated_base.xml
    :language: xml
    :start-after: <!-- SPHINX_MATERIAL -->
    :end-before: <!-- SPHINX_MATERIAL_END -->


All constitutive parameters such as density, viscosity, bulk modulus, and shear modulus are specified in the International System of Units.


-----------------------------------------------------------
Initial and boundary conditions
-----------------------------------------------------------

The next step is to specify fields, including:

  - The initial values: the ``waterDensity``, ``separableFace`` and the ``ruptureState`` of the propagating fracture have to be initialized,
  - The boundary conditions: fluid injection rates and the constraints of the outer boundaries have to be set.

In this example, a mass injection rate ``SourceFlux`` (``scale="-6.625"``) is applied at the surfaces of the initial fracture. Only one fourth of the total injection rate is defined in this boundary condition because only a quarter of the fracture is modeled (the problem is symmetric). The value given for ``scale`` is :math:`Q_0 \rho_f/4` (not :math:`Q_0 /4`). 
All the outer boundaries are subject to roller constraints. 
These boundary conditions are set through the ``FieldSpecifications`` section.


.. literalinclude:: ../../../../../../inputFiles/hydraulicFracturing/pennyShapedToughnessDominated_base.xml
    :language: xml
    :start-after: <!-- SPHINX_BC -->
    :end-before: <!-- SPHINX_BC_END -->


The parameters used in the simulation are summarized in the following table.

+------------------+-------------------------+--------------------+--------------------+
| Symbol           | Parameter               | Unit               | Value              |
+==================+=========================+====================+====================+
| :math:`K`        | Bulk Modulus            | [GPa]              | 20.0               |
+------------------+-------------------------+--------------------+--------------------+
| :math:`G`        | Shear Modulus           | [GPa]              | 12.0               |
+------------------+-------------------------+--------------------+--------------------+
| :math:`K_{IC}`   | Rock Toughness          | [MPa.m\ :sup:`1/2`]| 3.0                |
+------------------+-------------------------+--------------------+--------------------+
| :math:`\mu`      | Fluid Viscosity         | [Pa.s]             | 1.0x10\ :sup:`-6`  |
+------------------+-------------------------+--------------------+--------------------+
| :math:`Q_0`      | Injection Rate          | [m\ :sup:`3`/s]    | 0.0265             |
+------------------+-------------------------+--------------------+--------------------+
| :math:`t_{inj}`  | Injection Time          | [s]                | 400                |
+------------------+-------------------------+--------------------+--------------------+



---------------------------------
Inspecting results
---------------------------------

The following figure shows the distribution of :math:`\sigma_{zz}` at :math:`t=400 s` within the computational domain..

.. _problemVerificationPennyShapedToughnessDominatedFig1:
.. figure:: szz.png
   :align: center
   :width: 500
   :figclass: align-center

   Simulation result of :math:`\sigma_{zz}` at :math:`t=400 s`

 
The figure below compares the asymptotic solutions (curves) and the GEOSX simulation results (markers) for this analysis. The time history plots of fracture radius, fracture aperture and fluid pressure at the point source match the asymptotic solutions, confirming the accuracy of GEOSX simulations. 

.. plot:: docs/sphinx/advancedExamples/validationStudies/pennyFracToughnessDominated/pennyFracToughnessDominatedFigure.py


------------------------------------------------------------------
To go further
------------------------------------------------------------------


**Feedback on this example**

For any feedback on this example, please submit a `GitHub issue on the project's GitHub page <https://github.com/GEOSX/GEOSX/issues>`_.
