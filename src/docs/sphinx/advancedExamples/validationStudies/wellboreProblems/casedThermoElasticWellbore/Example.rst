.. _AdvancedExampleCasedThermoElasticWellbore:


####################################################
Cased ThemoElastic Wellbore Problem
####################################################

------------------------------------------------------------------
Problem description
------------------------------------------------------------------

This example uses the ``ThermalSinglePhasePoromechanics`` solver to handle a cased wellbore problem subjected to an uniform temperature change on the inner surface of the casing. The completed wellbore is composed of a steel casing, a cement sheath and rock formation. Isotropic linear ThermoElastic behavior is assumed for all the three materials. No separation or thermal barrier is allowed for the casing-cement and cement-rock contact interfaces. Plane strain condition is also asumed.

.. _problemSketchCasedThermalElasticWellboreFig:
.. figure:: sketch.png
   :align: center
   :width: 500
   :figclass: align-center

   Sketch of a cased thermoelastic wellbore 

Solution to this axisymetric problem can be obtained in the cylindrical coordinate system by an implicit 1D finite difference method `(Jane and Lee 1999) <https://www.sciencedirect.com/science/article/abs/pii/S0093641399000828>`__. Results of such analysis will be considered as reference resutls for validating GEOSX results.


**Input file**

This benchmark example uses no external input files and everything required is
contained within two GEOSX xml files that are located at:

.. code-block:: console

  inputFiles/wellbore/CasedThermoElasticWellbore_base.xml

and

.. code-block:: console

  inputFiles/wellbore/CasedThermoElasticWellbore_benchmark.xml

The corresponding integrated test is

.. code-block:: console

  inputFiles/wellbore/CasedThermoElasticWellbore_smoke.xml

-----------------------------------------------------------
Geometry and Mesh
-----------------------------------------------------------

The internal wellbore mesh generator ``InternalWellbore`` is employed to create the mesh of this wellbore problem. The radii of the casing cylinder, the cement sheath cylinder and the far-field boundary of the surrounding rock formation are defined by a vector ``radius``. In the tangent direction, ``theta`` angle is specified from 0 to 90 degree for a quarter of the geometry regarding the axisymmetry of the problem. The problem is plane strain and radial thermal diffusion on therefore only one single element in z-axis is needed. The trajectory of the well is defined by ``trajectory``, which is vertical in this case. The ``autoSpaceRadialElems`` parameters allow optimally increasing the element size from local zone around the wellbore to the far-field zone. In this example, the auto spacing option is only applied for the rock formation. The ``useCartesianOuterBoundary`` with a value 3 specified for the rock layer transforms the far-field boundary to a circular shape. The ``cellBlockNames`` and ``elementTypes`` define the regions and related element types associated to casing, cement sheath and rock. 
 
.. literalinclude:: ../../../../../../../inputFiles/wellbore/CasedThermoElasticWellbore_benchmark.xml
  :language: xml
  :start-after: <!-- SPHINX_WellboreMesh -->
  :end-before: <!-- SPHINX_WellboreMeshEnd -->

.. _meshCasedThermalElasticWellboreFig:
.. figure:: mesh.png
   :align: center
   :width: 500
   :figclass: align-center

   An optimized mesh for the considering cased wellbore.

-----------------------------------------------------------
Material properties
-----------------------------------------------------------

The bulk and shear drained elastic moduli of the materials as well as its drained linear thermal expansion coefficient relating temperature change to stress and displacement changes are defined within the ``constitution`` tag as follows:
 
.. literalinclude:: ../../../../../../../inputFiles/wellbore/CasedThermoElasticWellbore_base.xml
  :language: xml
  :start-after: <!-- SPHINX_ThermoElasticProperties -->
  :end-before: <!-- SPHINX_ThermoElasticPropertiesEnd -->

Here the solid density is also defined but it is not used because the gravitational effect is ignored in this example. To mimic a thermoelastic coupling without fluid flow, a negligible porosity and a zero Biot's coefficient are defined as:

.. literalinclude:: ../../../../../../../inputFiles/wellbore/CasedThermoElasticWellbore_base.xml
  :language: xml
  :start-after: <!-- SPHINX_PoroElasticProperties -->
  :end-before: <!-- SPHINX_PoroElasticPropertiesEnd -->

In this XML block, the Biot's coefficient is not directly defined but via the elastic bulk modulus :math:`K_{s}` of the solid skeleton as :math:`b_{Biot} = 1 - K/K_{s}`. In this example, we define a skeleton bulk modulus that is identical to the drained bulk modulus :math:`K` defined above to enforce the Biot's coefficient to zero.

The thermal conductivities and the volumetric heat capacities of casing, cement and rock are defined by following XML blocks:

.. literalinclude:: ../../../../../../../inputFiles/wellbore/CasedThermoElasticWellbore_base.xml
  :language: xml
  :start-after: <!-- SPHINX_ThermalProperties -->
  :end-before: <!-- SPHINX_ThermalPropertiesEnd -->

and

.. literalinclude:: ../../../../../../../inputFiles/wellbore/CasedThermoElasticWellbore_base.xml
  :language: xml
  :start-after: <!-- SPHINX_HeatCapacityProperties -->
  :end-before: <!-- SPHINX_HeatCapacityPropertiesEnd -->

A negligible permeability is defined for all the three layers to simulate a thermoelastic problem without fluid flow.

.. literalinclude:: ../../../../../../../inputFiles/wellbore/CasedThermoElasticWellbore_base.xml
  :language: xml
  :start-after: <!-- SPHINX_PermeabilityProperties -->
  :end-before: <!-- SPHINX_PermeabilityPropertiesEnd -->

Also, a neglibile volumetric heat capacity is defined for fluid to completely remove the thermal convection effect such that only thermal transfer via the diffusion phenomenon is considered.

.. literalinclude:: ../../../../../../../inputFiles/wellbore/CasedThermoElasticWellbore_base.xml
  :language: xml
  :start-after: <!-- SPHINX_FluidProperties -->
  :end-before: <!-- SPHINX_FluidPropertiesEnd -->

Other fluid properties such as viscosity, thermal expansion coefficient, etc. are not used because fluid flow is not simulated and pore pressure is zero everywhere.

-----------------------------------------------------------
Boundary conditions
-----------------------------------------------------------

The mechanical boundary conditions are applied to ensure the axisymmetric plane strain conditions such as:

.. literalinclude:: ../../../../../../../inputFiles/wellbore/CasedThermoElasticWellbore_base.xml
  :language: xml
  :start-after: <!-- SPHINX_PlaneStrainAxisymmetryBC -->
  :end-before: <!-- SPHINX_PlaneStrainAxisymmetryBCEnd -->

Besides, the far-field boundray is assumed to be fixed because the local changes on the wellbore must have neglibile effect on the far-field boundary. 

.. literalinclude:: ../../../../../../../inputFiles/wellbore/CasedThermoElasticWellbore_base.xml
  :language: xml
  :start-after: <!-- SPHINX_FixedFarFieldBC -->
  :end-before: <!-- SPHINX_FixedFarFieldBCEnd -->

The stress free condition on the inner surface of the casing is defined by:

.. literalinclude:: ../../../../../../../inputFiles/wellbore/CasedThermoElasticWellbore_base.xml
  :language: xml
  :start-after: <!-- SPHINX_ZeroInnerTractionBC -->
  :end-before: <!-- SPHINX_ZeroInnerTractionBCEnd -->

Initial resevoir temperature, that is also the far-field temperature and the temperature of a cold fluid applied on the inner surface of the casing are defined as

.. literalinclude:: ../../../../../../../inputFiles/wellbore/CasedThermoElasticWellbore_base.xml
  :language: xml
  :start-after: <!-- SPHINX_TemperatureBC -->
  :end-before: <!-- SPHINX_TemperatureBCEnd -->

It is important to remark that the initial stress of each layers must be set with accordance to the initial temperature: :math:`\sigma_{0} = 3K\alpha T_{0}` where :math:`\sigma_{0}` is the initial principal stress, :math:`T_{0}` is the initial temperature, :math:`K` is the drained bulk modulus and :math:`alpha` is the drained linear thermal expansion coefficient of the materials.

.. literalinclude:: ../../../../../../../inputFiles/wellbore/CasedThermoElasticWellbore_base.xml
  :language: xml
  :start-after: <!-- SPHINX_StressCasingInit -->
  :end-before: <!-- SPHINX_StressCasingInitEnd -->

Zero pore pressure is set everywhere to simulate a thermoelastic problem only:

.. literalinclude:: ../../../../../../../inputFiles/wellbore/CasedThermoElasticWellbore_base.xml
  :language: xml
  :start-after: <!-- SPHINX_ZeroPorePressureBC -->
  :end-before: <!-- SPHINX_ZeroPorePressureBCEnd -->

-----------------------------------------------------------
Collecting output data
-----------------------------------------------------------

It is convenient to collect data in hdf5 format that can be easily post-processed using Python. To collect the temperature field in the three layers for all the time steps, following XML blocks need to be defined:

.. literalinclude:: ../../../../../../../inputFiles/wellbore/CasedThermoElasticWellbore_base.xml
  :language: xml
  :start-after: <!-- SPHINX_TemperatureCollection -->
  :end-before: <!-- SPHINX_TemperatureCollectionEnd -->

.. literalinclude:: ../../../../../../../inputFiles/wellbore/CasedThermoElasticWellbore_base.xml
  :language: xml
  :start-after: <!-- SPHINX_TemperatureTimeHistory -->
  :end-before: <!-- SPHINX_TemperatureTimeHistoryEnd -->

Similarly, following blocks are needed to collect the solid stress:

.. literalinclude:: ../../../../../../../inputFiles/wellbore/CasedThermoElasticWellbore_base.xml
  :language: xml
  :start-after: <!-- SPHINX_SolidStressCollection -->
  :end-before: <!-- SPHINX_SolidStressCollectionEnd -->

.. literalinclude:: ../../../../../../../inputFiles/wellbore/CasedThermoElasticWellbore_base.xml
  :language: xml
  :start-after: <!-- SPHINX_SolidStressTimeHistory -->
  :end-before: <!-- SPHINX_SolidStressTimeHistoryEnd -->

The displacement field can be collected for the whole domain using ``nodeManager`` as follows

.. literalinclude:: ../../../../../../../inputFiles/wellbore/CasedThermoElasticWellbore_base.xml
  :language: xml
  :start-after: <!-- SPHINX_DisplacementCollection -->
  :end-before: <!-- SPHINX_DisplacementCollectionEnd -->

.. literalinclude:: ../../../../../../../inputFiles/wellbore/CasedThermoElasticWellbore_base.xml
  :language: xml
  :start-after: <!-- SPHINX_DisplacementTimeHistory -->
  :end-before: <!-- SPHINX_DisplacementTimeHistoryEnd -->

Also, periodic events are required to allow collecting these data. For example, periodic events for collecting the displacement field are defined as:

.. literalinclude:: ../../../../../../../inputFiles/wellbore/CasedThermoElasticWellbore_benchmark.xml
  :language: xml
  :start-after: <!-- SPHINX_DisplacementPeriodicEvent -->
  :end-before: <!-- SPHINX_DisplacementPeriodicEventEnd -->


---------------------------------
Results and benchmark
---------------------------------

A good agreement between the GEOSX results and analytical results for temperature distribution around the cased wellbore is shown in the figures below:

.. plot:: docs/sphinx/advancedExamples/validationStudies/wellboreProblems/CasedThermoElasticWellbore/thermoElastic_casedWellbore_temperature.py

.. _resultsCasedThermoElasticWellbore_temperatureFig:
.. figure:: temperature.png
   :align: center
   :width: 500
   :figclass: align-center

   Temperature distribution around the cased wellbore

and the validation for the radial displacement around the cased wellbore is shown below:

.. plot:: docs/sphinx/advancedExamples/validationStudies/wellboreProblems/CasedThermoElasticWellbore/thermoElastic_casedWellbore_displacement.py

.. _resultsCasedThermoElasticWellbore_displacementFig:
.. figure:: displacement.png
   :align: center
   :width: 500
   :figclass: align-center

   Radial displacement around the cased wellbore

The validations of the total radial and hoop stress components computed by GEOSX againt reference results are shown in the figure below:

.. plot:: docs/sphinx/advancedExamples/validationStudies/wellboreProblems/CasedThermoElasticWellbore/thermoElastic_casedWellbore_stress.py

.. _resultsCasedThermoElasticWellbore_stressFig:
.. figure:: stress.png
   :align: center
   :width: 1000
   :figclass: align-center

   Distribution of radial and hoop stress around the cased wellbore

------------------------------------------------------------------
To go further
------------------------------------------------------------------

**Feedback on this example**

This concludes the cased wellbore example.
For any feedback on this example, please submit a `GitHub issue on the project's GitHub page <https://github.com/GEOSX/GEOSX/issues>`_.
