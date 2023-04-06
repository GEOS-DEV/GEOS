.. _AdvancedExampleCasedElasticWellbore_ImperfectInterfaces:


###############################################
Cased Elastic Wellbore with Impefect Interfaces
###############################################

------------------------------------------------------------------
Problem description
------------------------------------------------------------------

This example uses the ``LagrangianContact`` solver to handle a cased wellbore problem with imperfect contact interfaces. The completed wellbore is composed of a steel casing, a cement sheath and rock formation. Isotropic linear elastic behavior is assumed for all the three materials. The contact surfaces between these materials are supposed to obey the Lagrangian contact law. 

.. _problemSketchCasedElasticWellboreInterfacesFig:
.. figure:: sketch.png
   :align: center
   :width: 600
   :figclass: align-center

   A cased wellbore with imperfect casing-cement and cement-rock interfaces 

With a radial compression loading, the imperfect contact intefaces behave just like the perfect one (see :ref:`AdvancedExampleCasedElasticWellbore`). With a radial tension loading on the inner face of the wellbore, the casing is debonded from the cement layer. Analytical results of the radial displacement :math:`u_{r}`, in the casing is expressed as `(Hervé and Zaoui, 1995) <https://link.springer.com/chapter/10.1007%2F978-94-015-8494-4_55>`__ :

.. math::
   u_{r} = Ar - \frac{B}{r}

where :math:`r` is the radial coordinate, :math:`A` and :math:`B` are constants that are obtained by solving the boundary conditions, as detailed in the post-processing script. The outer face of the casing as well as the inner face of the cement layer are free of stress because of the debonding at the casing-cement interface. Thefore, the displacement jump at the cement-rock interface is nil and the displacement jump accross the casing-cement interface is equal to :math:`u_r(r=r_{out_casing})` where :math:`r_{out_casing}` is the outer raidus of the casing. 


**Input file**

This benchmark example uses no external input files and everything required is
contained within two GEOSX xml files that are located at:

.. code-block:: console

  inputFiles/wellbore/CasedElasticWellbore_ImperfectInterfaces_base.xml

and

.. code-block:: console

  inputFiles/wellbore/CasedElasticWellbore_ImperfectInterfaces_benchmark.xml

The corresponding integrated test is

.. code-block:: console

  inputFiles/wellbore/CasedElasticWellbore_ImperfectInterfaces_smoke.xml

In this example, we should focus on following XML blocks


--------------------------------------------------------------------
Cylinder geometry
--------------------------------------------------------------------

The nodesets defining the casing-cement and cement-rock interfaces are curved. In this example, we use the ``Cylinder`` geometry to select these nodesets as shown below. This geometry is defined by the centers of its faces, ``firstFaceCenter`` and ``secondFaceCenter``, and its inner and outer radii, ``innerRadius`` and ``outerRadius``. We should note that the inner radius is optional as it is only needed for defining a hollow cylinder (i.e. an annulus). In this example, the inner radius is needed because we are selecting only the nodes on the casing-cement and the cement-rock interfaces.
 
.. literalinclude:: ../../../../../../../inputFiles/wellbore/CasedElasticWellbore_ImperfectInterfaces_benchmark.xml
  :language: xml
  :start-after: <!-- SPHINX_CementSheathInterfaces -->
  :end-before: <!-- SPHINX_CementSheathInterfacesEnd -->


--------------------------------------------------------------------   
Events
--------------------------------------------------------------------

In this example, we need to define a solo event for generating the imperfect contact surfaces as shown below:

.. literalinclude:: ../../../../../../../inputFiles/wellbore/CasedElasticWellbore_ImperfectInterfaces_benchmark.xml
  :language: xml
  :start-after: <!-- SPHINX_EventSurfaceGen -->
  :end-before: <!-- SPHINX_EventSurfaceGenEnd -->

where the surface generation solver is defined by:

.. literalinclude:: ../../../../../../../inputFiles/wellbore/CasedElasticWellbore_ImperfectInterfaces_base.xml
  :language: xml
  :start-after: <!-- SPHINX_SolversSurfaceGenerator -->
  :end-before: <!-- SPHINX_SolversSurfaceGeneratorEnd -->

Here the ``rockToughness`` is defined by default but it is omitted by this simulation.

To collect the displacement jump across the imperfect interfaces, we define also two periodic events as shown below:

.. literalinclude:: ../../../../../../../inputFiles/wellbore/CasedElasticWellbore_ImperfectInterfaces_benchmark.xml
  :language: xml
  :start-after: <!-- SPHINX_EventDisplacementJumpHistoryCollection -->
  :end-before: <!-- SPHINX_EventDisplacementJumpHistoryCollectionEnd -->

The corresponding ``Tasks`` and ``Outputs`` targets need to be defined together with these events.

--------------------------------------------------------------------	       
NumericalMethods
--------------------------------------------------------------------

The ``stabilizationName`` that is required in the ``LagrangianContact`` solver is defined by: 

.. literalinclude:: ../../../../../../../inputFiles/wellbore/CasedElasticWellbore_ImperfectInterfaces_base.xml
  :language: xml
  :start-after: <!-- SPHINX_NumericalMethods -->
  :end-before: <!-- SPHINX_NumericalMethodsEnd -->

--------------------------------------------------------------------	       
Contact region and material
--------------------------------------------------------------------

The imperfect contact surfaces between casing, cement and rock layers are defined as ``Fracture`` as shown below:

.. literalinclude:: ../../../../../../../inputFiles/wellbore/CasedElasticWellbore_ImperfectInterfaces_base.xml
  :language: xml
  :start-after: <!-- SPHINX_SurfaceElementRegion -->
  :end-before: <!-- SPHINX_SurfaceElementRegion -->

Here the ``faceBlock`` name, ``faceElementSubRegion``, is needed for defining Tasks for collecting displacement Jump across the contact surfaces. The ``defaultAperture`` defined in this block is the default hydraulic aperture, that should not be confused with the mechanical aperture. For this pure mechanical problem, this default hydraulic aperture parameter is omitted. The fracture material given in the ``materialList`` is defined as follows:

.. literalinclude:: ../../../../../../../inputFiles/wellbore/CasedElasticWellbore_ImperfectInterfaces_base.xml
  :language: xml
  :start-after: <!-- SPHINX_MaterialContact -->
  :end-before: <!-- SPHINX_MaterialContactEnd -->

For this pure mechanical problem without fluid flow and shearing stress acting on the contact surface, all the parameters defined in this block are omitted.

---------------------------------
Results and benchmark
---------------------------------

The GEOSX results of displacement jump accross the casing-cement and cement-rock interfaces are shown in the figure below: 

.. _CasedElasticWellboreInterfacesDisplacementJump3DFig:
.. figure:: displacementJump3D.png
   :align: center
   :width: 600
   :figclass: align-center

   Displacement jumps across the casing-cement and cement-rock interfaces

We can observe an expected zero displacement jump at the cement-rock interface under a tension stress on the inner surface of the casing. This is because the applied stress do not cause any strain on the cement and rock layers after the debonding occured at the casing-cement interface. The displacement jump at the casing-cement interface is homogenous and varies with time because the tension stress on the inner surface of the casing varies with time as defined in the XML file. A perfect validation of GEOSX results versus theoretical results is shown in the figure below:

.. plot:: docs/sphinx/advancedExamples/validationStudies/wellboreProblems/casedContactElasticWellbore/elastic_casedWellbore_displacementJump.py

.. _CasedElasticWellboreInterfacesDisplacementJumpVsTimeFig:
.. figure:: displacementJump.png
   :align: center
   :width: 600
   :figclass: align-center

   Displacement jump versus time at the casing-cement interface

------------------------------------------------------------------
To go further
------------------------------------------------------------------

**Feedback on this example**

This concludes the cased wellbore example.
For any feedback on this example, please submit a `GitHub issue on the project's GitHub page <https://github.com/GEOSX/GEOSX/issues>`_.
