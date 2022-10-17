.. _AdvancedExampleViscoModifiedCamClay:


#############################################################################
Visco Modified CamClay model: Triaxial Driver versus Semi-analytical solution 
#############################################################################

------------------------------------------------------------------
Problem description
------------------------------------------------------------------

This example uses the Triaxial Driver to simulate a visco-elasto-plastic oedometric compression test of a Visco Modified CamClay solid. Oedometric condition with zero lateral strain together with loading/unloading axial strain periods are imposed. Semi analytical results for the mean and shear stress variations :math:`\Delta p` and :math:`\Delta q` can be etablished, considering the Perzyna approach, for the considered oedometric boundary conditions as `(Runesson et al. 1999) <https://onlinelibrary.wiley.com/doi/abs/10.1002/(SICI)1099-1484(199901)4:1%3C75::AID-CFM60%3E3.0.CO;2-4>`__ :

.. math::
   \Delta p = K(\Delta\varepsilon_{V} - \Delta\lambda \frac{\partial G}{\partial p})

.. math::
   \Delta q = 2\mu(\Delta\varepsilon_{V} - \Delta\lambda \frac{3}{2}\frac{\partial G}{\partial q})

where :math:`K` and :math:`\mu` are elastic bulk and shear moduli, :math:`G` the plastic potential and :math:`\Delta\lambda` the visco-plastic multiplier that can be approximated by:

.. math::
   \Delta\lambda = \frac{\Delta t}{t_*} \frac{F}{3\mu\frac{\partial F}{\partial q}\frac{\partial G}{\partial q} + K\frac{\partial F}{\partial p}\frac{\partial G}{\partial p} + h}

in which :math:`\Delta t` is the time increment, :math:`t_*` the relaxation time, :math:`F` the stress function defining the visco-plastic yield surface and :math:`h` the hardening rate that is defined by:

.. math::
   h = -\frac{\partial F}{\partial \lambda}


These solutions are implemented in a Python script associated to this DOC for verifying GEOSX results.


**Input files**

This validation example uses two GEOSX xml files that are located at:

.. code-block:: console

  inputFiles/triaxialDriver/triaxialDriver_base.xml

and

.. code-block:: console

  inputFiles/triaxialDriver/triaxialDriver_ViscoModifiedCamClay.xml

It also uses a set of table files located at:

.. code-block:: console

  inputFiles/triaxialDriver/tables/


A Python script for the semi-analytical solutions presented above as well as for post-processing the GEOSX results is provided at:

.. code-block:: console

  src/docs/sphinx/advancedExamples/validationStudies/viscoplasticity/ViscoModifiedCamClay/TriaxialDriver_vs_SemiAnalytic_ViscoModifiedCamClay.py

For this example, we focus on the ``Task`` and the ``Constitutive`` tags.

------------------------------------------------------------------
Task
------------------------------------------------------------------

The imposed axial strain loading/unloading periods, the lateral zero strain as well as the initial stress are defined in the ``Task`` block as 

.. literalinclude:: ../../../../../../../inputFiles/triaxialDriver/triaxialDriver_ViscoModifiedCamClay.xml
    :language: xml
    :start-after: <!-- SPHINX_TASK -->
    :end-before: <!-- SPHINX_TASK_END -->

------------------------------
Constitutive laws
------------------------------

The elasto-visco-plastic parameters are defined as


.. literalinclude:: ../../../../../../../inputFiles/triaxialDriver/triaxialDriver_base.xml
    :language: xml
    :start-after: <!-- SPHINX_MATERIAL_VISCO_MODIFIED_CAMCLAY -->
    :end-before: <!-- SPHINX_MATERIAL_VISCO_MODIFIED_CAMCLAY_END -->


All constitutive parameters such as density, viscosity, and the bulk and shear moduli are specified in the International System of Units.

--------------------------------------------------------------
A comparison between GEOSX results and semi-analytical results
--------------------------------------------------------------

The simulation results are saved in a text file, named ``ViscoModifiedCamClayResults.txt``. A comparison between the results given by the TriaxialDriver solver in GEOSX and the semi-analytical results presented above is show below. This discrepancy between these results may due to the discrepancy between the Duvaut-Lions approach and the Perzyna approach for time dependant behavior when applying for the Modified CamClay model as discussed by `Runesson et al. (1999) <https://onlinelibrary.wiley.com/doi/abs/10.1002/(SICI)1099-1484(199901)4:1%3C75::AID-CFM60%3E3.0.CO;2-4>`__. 


.. plot:: docs/sphinx/advancedExamples/validationStudies/viscoplasticity/ViscoModifiedCamClay/TriaxialDriver_vs_SemiAnalytic_ViscoModifiedCamClay.py


------------------------------------------------------------------
To go further
------------------------------------------------------------------


**Feedback on this example**

For any feedback on this example, please submit a `GitHub issue on the project's GitHub page <https://github.com/GEOSX/GEOSX/issues>`_.
