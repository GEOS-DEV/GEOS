.. _AdvancedExampleViscoExtendedDruckerPrager:


####################################################################################
Visco Extended Drucker-Prager model: Triaxial Driver versus Semi-analytical solution 
####################################################################################

------------------------------------------------------------------
Problem description
------------------------------------------------------------------

This example uses the Triaxial Driver to simulate a triaxial compression test of a Visco Extended Drucker-Prager solid. Constant lateral confining stress together with loading/unloading axial strain periods are imposed. Imposed axial strain range are high enough for allowing plastic yield in both loading and unloading period. This complicated senario is optimal for verifying the numerical convergence of the Visco Extended Drucker-Prager constitutive model implimented in GEOSX. 

Semi analytical result for axial stress variation :math:`\Delta\sigma_{V}` and lateral strain variation :math:`\Delta\varepsilon_{V}` can be etablished for the considered triaxial boundary conditions following the theoretical basis of the Perzyna time dependent approach presented by `(Runesson et al. 1999) <https://onlinelibrary.wiley.com/doi/abs/10.1002/(SICI)1099-1484(199901)4:1%3C75::AID-CFM60%3E3.0.CO;2-4>`__ as:

.. math::
   \Delta\sigma_{V} = (\Delta\varepsilon_{V}-\Delta\lambda\frac{\theta b-3}{3}) E

.. math::
   \Delta\varepsilon_{H} = \Delta\varepsilon_{V} - \frac{\Delta\sigma_{V}}{2\mu} + \frac{3}{2}\Delta\lambda

where :math:`E` and :math:`\mu` are the elastic Young and shear moduli. The visco-plastic multiplier :math:`\Delta\lambda` can be approximated by:

.. math::
   \Delta\lambda = \frac{\Delta t}{t_*} \frac{F}{3\mu + K\theta b^2 + h}

in which :math:`\Delta t` is the time increment, :math:`t_*` the relaxation time, :math:`F` the stress function defining the visco-plastic yield surface, :math:`K` the elastic bulk modulus, :math:`b` the frictional parameter defining the visco-plastic yield surface, :math:`\theta` the dilation ratio defining the plastic potential and :math:`h` the hardening rate. The hardening rate :math:`h` is defined by

.. math::
   h = \frac{\partial F}{\partial b} \frac{\partial b}{\partial \lambda}

These solutions were etablished for a positive shear stress :math:`q = -(\sigma_{V} - \sigma_{H})` (negative sign convention for compression stress). For the case when the plastic yield occurs at a negative shear stress, we have

.. math::
   \Delta\sigma_{V} = (\Delta\varepsilon_{V}-\Delta\lambda\frac{\theta b+3}{3}) E

.. math::
   \Delta\varepsilon_{H} = \Delta\varepsilon_{V} - \frac{\Delta\sigma_{V}}{2\mu} - \frac{3}{2}\Delta\lambda

These solutions are implemented in a Python script associated to this DOC for verifying GEOSX results.


**Input files**

This validation example uses two GEOSX xml files that are located at:

.. code-block:: console

  inputFiles/triaxialDriver/triaxialDriver_base.xml

and

.. code-block:: console

  inputFiles/triaxialDriver/triaxialDriver_ViscoExtendedDruckerPrager.xml

It also uses a set of table files located at:

.. code-block:: console

  inputFiles/triaxialDriver/tables/


A Python script for the semi-analytical solutions presented above as well as for post-processing the GEOSX results is provided at:

.. code-block:: console

  src/docs/sphinx/advancedExamples/validationStudies/viscoplasticity/ViscoExtendedDruckerPrager/TriaxialDriver_vs_SemiAnalytic_ViscoExtendedDruckerPrager.py

For this example, we focus on the ``Task`` and the ``Constitutive`` tags.

------------------------------------------------------------------
Task
------------------------------------------------------------------

The imposed axial strain loading/unloading periods, the constant lateral confining stress as well as the initial stress are defined in the ``Task`` block as 

.. literalinclude:: ../../../../../../../inputFiles/triaxialDriver/triaxialDriver_ViscoExtendedDruckerPrager.xml
    :language: xml
    :start-after: <!-- SPHINX_TASK -->
    :end-before: <!-- SPHINX_TASK_END -->

------------------------------
Constitutive laws
------------------------------

The elasto-visco-plastic parameters are defined as


.. literalinclude:: ../../../../../../../inputFiles/triaxialDriver/triaxialDriver_base.xml
    :language: xml
    :start-after: <!-- SPHINX_MATERIAL_VISCO_EXTENDED_DRUCKER_PRAGER -->
    :end-before: <!-- SPHINX_MATERIAL_VISCO_EXTENDED_DRUCKER_PRAGER_END -->


All constitutive parameters such as density, viscosity, and bulk and shear moduli are specified in the International System of Units.

--------------------------------------------------------------
A comparison between GEOSX results and semi-analytical results
--------------------------------------------------------------

The simulation results are saved in a text file, named ``ViscoExtendedDruckerPragerResults.txt``. A comparison between the results given by the TriaxialDriver solver in GEOSX and the approximated semi-analytical results presented above is show below. Interestingly we observed that the Duvaut-Lions approache implemented in GEOSX can fit perfectly with the Perzyna approach that was considered for deriving the analytical results. This consistency between these time dependence approahes is because of the linear hardening law of the considered constitutive model as already discussed by `(Runesson et al. 1999) <https://onlinelibrary.wiley.com/doi/abs/10.1002/(SICI)1099-1484(199901)4:1%3C75::AID-CFM60%3E3.0.CO;2-4>`__ . 


.. plot:: docs/sphinx/advancedExamples/validationStudies/viscoplasticity/ViscoExtendedDruckerPrager/TriaxialDriver_vs_SemiAnalytic_ViscoExtendedDruckerPrager.py


------------------------------------------------------------------
To go further
------------------------------------------------------------------


**Feedback on this example**

For any feedback on this example, please submit a `GitHub issue on the project's GitHub page <https://github.com/GEOSX/GEOSX/issues>`_.
