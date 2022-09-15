.. _ViscoplasticModelsComparison:


#############################################################
A comparison of viscoplastic models in GEOSX
#############################################################


**Context**

This example aims to compare results obtained by the ``TriaxialDriver`` for different viscoplastic models implemented in GEOSX. The ``TriaxialDriver`` that mimics the boundary conditions of the laboratory triaxial compression tests allows us to verify the correctness as well as the numerical convergence of the implemented constitutive models without passing through the finite element or the finite volume step.

**Objectives**

At the end of this example, you will know:
 - how to define different viscoplastic models,
 - how to set up a triaxial test scenario with the ``TriaxialDriver``,
 - how to postprocess and compare results of different viscoplastic models.


**Input file**

The base XML file for this comparison is located at:

.. code-block:: console

  inputFiles/triaxialDriver/triaxialDriver_base.xml

This base file is accompanied by different xml files that are defined for different viscoplastic models. The name of each file is clearly chosen to highlight the corresponding viscoplastic model. For example, the Visco Extended Drucker-Prager model is defined by the file:

.. code-block:: console

  inputFiles/triaxialDriver/triaxialDriver_ViscoExtendedDruckerPrager.xml

This example also uses a set of table files located at:

.. code-block:: console

  inputFiles/triaxialDriver/tables/


Last, a Python script for post-processing the results is provided:

.. code-block:: console

  src/docs/sphinx/advancedExamples/validationStudies/viscoplasticity/modelComparisons/VPModelsComparison.py


------------------------------------------------------------------
Description of the case
------------------------------------------------------------------

The ``TriaxialDriver`` in GEOSX is employed to verify the correctness of the numerical results obtained by different viscoplastic models in GEOS including:

 - the Drucker-Prager model (see :ref:`DruckerPrager`),
 - the Visco Drucker-Prager model (see :ref:`ViscoDruckerPrager`),
 - the Extended Drucker-Prager model (see :ref:`DruckerPragerExtended`),
 - the Visco Extended Drucker-Prager model (see :ref:`ViscoDruckerPragerExtended`),

To replicate a typical triaxial test, we use a table function to specify loading conditions in axial and radial directions. The resulting strains and stresses in both directions are numerically calculated and saved into a simple ASCII output file. A python script is given for interpreting and plotting the comparison between the considered constitutive models.


For this example, we focus on the ``Constitutive``, and the ``Task`` tags.


------------------------------
Constitutive laws
------------------------------

The material definition of the considered viscoplastic models are given in the base file. Each model is defined by a separate xml block as:


.. literalinclude:: ../../../../../inputFiles/triaxialDriver/triaxialDriver_base.xml
    :language: xml
    :start-after: <!-- SPHINX_MATERIAL -->
    :end-before: <!-- SPHINX_MATERIAL_END -->


All constitutive parameters such as density, viscosity, bulk modulus, and shear modulus are specified in the International System of Units.



------------------------------------------------------------------
Task
------------------------------------------------------------------

One additional xml file per viscoplastic model is given to define the corresponding task of that model. For example the task of the Visco Extended Drucker-Prager model is defined as:

.. literalinclude:: ../../../../../inputFiles/triaxialDriver/triaxialDriver_ViscoExtendedDruckerPrager.xml
    :language: xml
    :start-after: <!-- SPHINX_TASK -->
    :end-before: <!-- SPHINX_TASK_END -->


In addition to triaxial tests, volumetric and oedometer tests can be simulated by changing the ``strainControl`` mode, and by defining loading controls in axial and radial direction accordingly. A volumetric test can be modelled by setting the axial and radial control functions to the same strain function, whereas an oedometer test runs by setting the radial strain to zero (see :ref:`TriaxialDriver`).


---------------------------------
Inspecting results
---------------------------------

The simulation results are saved in a text files, named after the corresponding viscoplastic model. For example results of the Visco Extended Drucker-Prager model is stored in the file ``ViscoExtendedDruckerPragerResults.txt``.
These output files has a brief header explaining the meaning of each column. Each row corresponds to one timestep of the driver, starting from initial conditions in the first row. They can be interpreted and plotted by the provided python script for visualizing the comparison between the models:


.. plot:: docs/sphinx/advancedExamples/validationStudies/viscoplasticity/modelComparisons/VPModelsComparison.py


------------------------------------------------------------------
To go further
------------------------------------------------------------------


**Feedback on this example**

For any feedback on this example, please submit a `GitHub issue on the project's GitHub page <https://github.com/GEOSX/GEOSX/issues>`_.
