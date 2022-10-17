.. _AdvancedExampleDruckerPrager:


#####################################################################
Drucker-Prager model: Triaxial Driver versus Semi-analytical solution 
#####################################################################

------------------------------------------------------------------
Problem description
------------------------------------------------------------------

This example uses the Triaxial Driver to simulate an elasto-plastic triaxial compression test of a Drucker-Prager solid. Constant lateral confining stress together with loading/unloading axial strain periods are imposed. Imposed axial strain range are high enough for allowing plastic yield in both loading and unloading period. This complicated senario is optimal for verifying the numerical convergence of the Drucker-Prager constitutive model implimented in GEOSX. 

Semi analytical result for axial stress variation :math:`\Delta\sigma_{V}` and lateral strain variation :math:`\Delta\varepsilon_{V}` can be etablished for the considered triaxial boundary conditions as:

.. math::
   \Delta\sigma_{V} = \Delta\varepsilon_{V} E_{ep}

.. math::
   \Delta\varepsilon_{H} = \Delta\varepsilon_{V} - \frac{\Delta\sigma_{V}}{E^\prime_{ep}}

where :math:`E_{ep}` and :math:`E^\prime_{ep}` are elasto-plastic Young moduli that can be obtained from the elastic Young and shear moduli (:math:`E` and :math:`\mu`), the frictional parameter :math:`b`, the dilation parameter :math:`b^\prime` and the hardening rate :math:`h` of the Drucker-Prager model by:

.. math::
   \frac{1}{E_{ep}} = \frac{1}{E} + \frac{(b^\prime-3)(b-3)}{9h}

.. math::
   \frac{1}{E^\prime_{ep}} = \frac{1}{2\mu} - \frac{b-3}{2h}

These solutions are applied only when plastic yield condition is statisfied. The cohesion parameter defining the plastic yield surface is updated with stress change as

.. math::
   \Delta a = \frac{b-3}{3}\Delta\sigma_{V}


These solutions were etablished for a positive shear stress :math:`q = -(\sigma_{V} - \sigma_{H})` (negative sign convention for compression stress). For the case when the plastic yield occurs at a negative shear stress, we have

.. math::
   \frac{1}{E_{ep}} = \frac{1}{E} + \frac{(b^\prime+3)(b+3)}{9h}

.. math::
   \frac{1}{E^\prime_{ep}} = \frac{1}{2\mu} + \frac{b+3}{2h}

and

.. math::
   \Delta a = \frac{b+3}{3}\Delta\sigma_{V} 

These solutions are implemented in a Python script associated to this DOC for verifying GEOSX results.


**Input files**

This validation example uses two GEOSX xml files that are located at:

.. code-block:: console

  inputFiles/triaxialDriver/triaxialDriver_base.xml

and

.. code-block:: console

  inputFiles/triaxialDriver/triaxialDriver_DruckerPrager.xml

It also uses a set of table files located at:

.. code-block:: console

  inputFiles/triaxialDriver/tables/


A Python script for the semi-analytical solutions presented above as well as for post-processing the GEOSX results is provided at:

.. code-block:: console

  src/docs/sphinx/advancedExamples/validationStudies/viscoplasticity/DruckerPrager/TriaxialDriver_vs_SemiAnalytic_DruckerPrager.py

For this example, we focus on the ``Task`` and the ``Constitutive`` tags.

------------------------------------------------------------------
Task
------------------------------------------------------------------

The imposed axial strain loading/unloading periods, the constant lateral confining stress as well as the initial stress are defined in the ``Task`` block as 

.. literalinclude:: ../../../../../../../inputFiles/triaxialDriver/triaxialDriver_DruckerPrager.xml
    :language: xml
    :start-after: <!-- SPHINX_TASK -->
    :end-before: <!-- SPHINX_TASK_END -->

------------------------------
Constitutive laws
------------------------------

The elasto-plastic parameters are defined as


.. literalinclude:: ../../../../../../../inputFiles/triaxialDriver/triaxialDriver_base.xml
    :language: xml
    :start-after: <!-- SPHINX_MATERIAL_DRUCKER_PRAGER -->
    :end-before: <!-- SPHINX_MATERIAL_DRUCKER_PRAGER_END -->


All constitutive parameters such as density, viscosity, and bulk and shear moduli are specified in the International System of Units.

--------------------------------------------------------------
A comparison between GEOSX results and semi-analytical results
--------------------------------------------------------------

The simulation results are saved in a text file, named ``DruckerPragerResults.txt``. A perfect comparison between the results given by the TriaxialDriver solver in GEOSX and the semi-analytical results presented above is show below 


.. plot:: docs/sphinx/advancedExamples/validationStudies/viscoplasticity/DruckerPrager/TriaxialDriver_vs_SemiAnalytic_DruckerPrager.py


------------------------------------------------------------------
To go further
------------------------------------------------------------------


**Feedback on this example**

For any feedback on this example, please submit a `GitHub issue on the project's GitHub page <https://github.com/GEOSX/GEOSX/issues>`_.
