.. _AdvancedExampleModifiedCamClay:


#######################################################################
Modified CamClay model: Triaxial Driver versus Semi-analytical solution 
#######################################################################

------------------------------------------------------------------
Problem description
------------------------------------------------------------------

This example uses the Triaxial Driver to simulate an elasto-plastic triaxial compression test of a Modified CamClay solid. Oedometric condition with zero lateral strain together with loading/unloading axial strain periods are imposed. Semi analytical results for the mean and shear stress variations :math:`\delta p` and :math:`\delta q` can be obtained from the imposed vertical strain variation by solving the following equation system:

.. math::
   \delta \varepsilon_V = \delta p(\frac{1}{K} + \frac{1}{h}\frac{\partial F}{\partial p}\frac{\partial G}{\partial p}) + \delta q \frac{1}{h}\frac{\partial F}{\partial q}\frac{\partial G}{\partial p}

.. math::
   \delta \varepsilon_V = \delta p \frac{3}{2h}\frac{\partial F}{\partial p}\frac{\partial G}{\partial q} + \delta q(\frac{1}{2\mu} + \frac{1}{h}\frac{\partial F}{\partial q}\frac{\partial G}{\partial q})

where :math:`K` and :math:`\mu` are elastic bulk and shear moduli, :math:`F` and :math:`G` the plastic yield surface and the plastic potential and :math:`h` the hardening rate that is defined by:

.. math::
   h = -\frac{\partial F}{\partial \varepsilon^{vp}_{vol}}\frac{\partial G}{\partial p}

in which :math:`\varepsilon^{vp}_{vol}` is the volumetric visco-plastic strain. These solutions are implemented in a Python script associated to this DOC for verifying GEOSX results.


**Input files**

This benchmark example uses no external input files and everything required is
contained within two GEOSX xml files that are located at:

.. code-block:: console

  inputFiles/triaxialDriver/triaxialDriver_base.xml

and

.. code-block:: console

  inputFiles/triaxialDriver/triaxialDriver_ModifiedCamClay.xml

This example also uses a set of table files located at:

.. code-block:: console

  inputFiles/triaxialDriver/tables/


A Python script for the semi-analytical solutions presented above as well as for post-processing the GEOSX results is provided at:

.. code-block:: console

  src/docs/sphinx/advancedExamples/validationStudies/viscoplasticity/ModifiedCamClay/TriaxialDriver_vs_SemiAnalytic_ModifiedCamClay.py

For this example, we focus on the ``Task`` and the ``Constitutive`` tags.

------------------------------------------------------------------
Task
------------------------------------------------------------------

The imposed axial strain loading/unloading periods, the constant lateral confining stress as well as the initial stress are defined in the ``Task`` block as 

.. literalinclude:: ../../../../../inputFiles/triaxialDriver/triaxialDriver_ModifiedCamClay.xml
    :language: xml
    :start-after: <!-- SPHINX_TASK -->
    :end-before: <!-- SPHINX_TASK_END -->

------------------------------
Constitutive laws
------------------------------

The elasto-plastic parameters are defined as


.. literalinclude:: ../../../../../inputFiles/triaxialDriver/triaxialDriver_base.xml
    :language: xml
    :start-after: <!-- SPHINX_MATERIAL_MODIFIED_CAMCLAY -->
    :end-before: <!-- SPHINX_MATERIAL_MODIFIED_CAMCLAY_END -->


All constitutive parameters such as density, viscosity, bulk modulus, and shear modulus are specified in the International System of Units.

--------------------------------------------------------------
A comparison between GEOSX results and semi-analytical results
--------------------------------------------------------------

The simulation results are saved in a text file, named ``ModifiedCamClayResults.txt``. A perfect comparison between the results given by the TriaxialDriver solver in GEOSX and the semi-analytical results presented above is show below 


.. plot:: docs/sphinx/advancedExamples/validationStudies/viscoplasticity/ModifiedCamClay/TriaxialDriver_vs_SemiAnalytic_ModifiedCamClay.py


------------------------------------------------------------------
To go further
------------------------------------------------------------------


**Feedback on this example**

For any feedback on this example, please submit a `GitHub issue on the project's GitHub page <https://github.com/GEOSX/GEOSX/issues>`_.
