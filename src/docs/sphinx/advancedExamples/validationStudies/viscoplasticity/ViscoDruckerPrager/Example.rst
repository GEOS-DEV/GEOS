.. _AdvancedExampleViscoDruckerPrager:


###########################################################################
Visco Drucker-Prager model: Triaxial Driver versus Semi-analytical solution 
###########################################################################

------------------------------------------------------------------
Problem description
------------------------------------------------------------------

This example uses the Triaxial Driver to simulate a triaxial compression test of a Visco Drucker-Prager solid. Constant lateral confining stress together with loading/unloading axial strain periods are imposed. Imposed axial strain range are high enough for allowing visco-plastic yield in both loading and unloading period. This complicated senario is optimal for verifying the numerical convergence of the Visco Drucker-Prager constitutive model implimented in GEOSX. 

Semi analytical result for axial stress variation :math:`\delta\sigma_{V}` and lateral strain variation :math:`\delta\varepsilon_{V}` can be etablished for the considered triaxial boundary conditions as `(Runesson et al. 1999) <https://onlinelibrary.wiley.com/doi/abs/10.1002/(SICI)1099-1484(199901)4:1%3C75::AID-CFM60%3E3.0.CO;2-4>`__ :

.. math::
   \delta\sigma_{V} = (\delta\varepsilon_{V}-\delta\lambda\frac{b^\prime-3}{3}) E

.. math::
   \delta\varepsilon_{H} = \delta\varepsilon_{V} - \frac{\delta\sigma_{V}}{2\mu} + \frac{3}{2}\delta\lambda

where :math:`E` and :math:`\mu` are the elastic Young and shear moduli. The visco-plastic multiplier :math:`\delta\lambda` can be approximated by:

.. math::
   \delta\lambda = \frac{\delta t}{t_*} \frac{F}{3\mu + Kbb^\prime + h}

in which :math:`\delta t` is the time increment, :math:`t_*` the relaxation time, :math:`F` the stress function defining the visco-plastic yield surface, :math:`K` the elastic bulk modulus, :math:`b` the frictional parameter defining the visco-plastic yield surface, :math:`b^\prime` the dilation parameter defining the plastic potential and :math:`h` the hardening rate. These solutions are applied only when plastic yield condition is staisfied. The cohesion parameter defining the plastic yield surface is updated with stress change as

.. math::
   \delta a = h \delta\lambda

These solutions were etablished for a positive shear stress `q = -(\sigma_{V} - \sigma_{H})` (negative sign convention for compression stress). For the case when the plastic yield occurs at a negative shear stress, we have

.. math::
   \delta\sigma_{V} = (\delta\varepsilon_{V}-\delta\lambda\frac{b^\prime+3}{3}) E

.. math::
   \delta\varepsilon_{H} = \delta\varepsilon_{V} - \frac{\delta\sigma_{V}}{2\mu} - \frac{3}{2}\delta\lambda

These solutions are implemented in a Python script associated to this DOC for verifying GEOSX results.


**Input files**

This benchmark example uses no external input files and everything required is
contained within two GEOSX xml files that are located at:

.. code-block:: console

  inputFiles/triaxialDriver/triaxialDriver_base.xml

and

.. code-block:: console

  inputFiles/triaxialDriver/triaxialDriver_ViscoDruckerPrager.xml

This example also uses a set of table files located at:

.. code-block:: console

  inputFiles/triaxialDriver/tables/


A Python script for the semi-analytical solutions presented above as well as for post-processing the GEOSX results is provided at:

.. code-block:: console

  src/docs/sphinx/advancedExamples/validationStudies/viscoplasticity/ViscoDruckerPrager/TriaxialDriver_vs_SemiAnalytic_ViscoDruckerPrager.py

For this example, we focus on the ``Task`` and the ``Constitutive`` tags.

------------------------------------------------------------------
Task
------------------------------------------------------------------

The imposed axial strain loading/unloading periods, the constant lateral confining stress as well as the initial stress are defined in the ``Task`` block as 

.. literalinclude:: ../../../../../inputFiles/triaxialDriver/triaxialDriver_ViscoDruckerPrager.xml
    :language: xml
    :start-after: <!-- SPHINX_TASK -->
    :end-before: <!-- SPHINX_TASK_END -->

------------------------------
Constitutive laws
------------------------------

The elasto-plastic parameters such as the elastic moduli, the friction angle and the cohesion defining the plastic yield surface as well as the hardening rate are defined as


.. literalinclude:: ../../../../../inputFiles/triaxialDriver/triaxialDriver_base.xml
    :language: xml
    :start-after: <!-- SPHINX_MATERIAL_VISCO_DRUCKER_PRAGER -->
    :end-before: <!-- SPHINX_MATERIAL_VISCO_DRUCKER_PRAGER_END -->


All constitutive parameters such as density, viscosity, bulk modulus, and shear modulus are specified in the International System of Units.

--------------------------------------------------------------
A comparison between GEOSX results and semi-analytical results
--------------------------------------------------------------

The simulation results are saved in a text file, named ``ViscoDruckerPragerResults.txt``. A very well comparison between the results given by the TriaxialDriver solver in GEOSX and the approximated semi-analytical results presented above is show below 


.. plot:: docs/sphinx/advancedExamples/validationStudies/viscoplasticity/ViscoDruckerPrager/TriaxialDriver_vs_SemiAnalytic_ViscoDruckerPrager.py


------------------------------------------------------------------
To go further
------------------------------------------------------------------


**Feedback on this example**

For any feedback on this example, please submit a `GitHub issue on the project's GitHub page <https://github.com/GEOSX/GEOSX/issues>`_.
