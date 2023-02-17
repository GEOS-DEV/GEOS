.. _AdvancedExampleCasedThermoElasticWellbore:


####################################################
Cased ThemoElastic Wellbore Problem
####################################################

------------------------------------------------------------------
Problem description
------------------------------------------------------------------

This example uses the ``ThermalSinglePhasePoromechanics`` solver to handle a cased wellbore problem subjected to a temperature change on the inner surface of the casing. The completed wellbore is composed of a steel casing, a cement sheath and rock formation. Isotropic linear ThermoElastic behavior is assumed for all the three materials. No separation or thermal barrier is allowed for the casing-cement and cement-rock contact interfaces.

Solution to this axisymetric problem can be obtained in the cylindrical coordinate system by an implicit 1D finite difference method (see e.g. `(Jane and Lee 1999) <https://www.sciencedirect.com/science/article/abs/pii/S0093641399000828>`__ ). Results of such analysis will be considered as reference resutls for validating GEOSX results.


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

---------------------------------
Results and benchmark
---------------------------------

A good agreement between the GEOSX results and analytical results is shown in the figure below:

.. plot:: docs/sphinx/advancedExamples/validationStudies/wellboreProblems/CasedThermoElasticWellbore/thermoElastic_casedWellbore_temperature.py

.. plot:: docs/sphinx/advancedExamples/validationStudies/wellboreProblems/CasedThermoElasticWellbore/thermoElastic_casedWellbore_displacement.py

.. plot:: docs/sphinx/advancedExamples/validationStudies/wellboreProblems/CasedThermoElasticWellbore/thermoElastic_casedWellbore_stress.py

------------------------------------------------------------------
To go further
------------------------------------------------------------------

**Feedback on this example**

This concludes the cased wellbore example.
For any feedback on this example, please submit a `GitHub issue on the project's GitHub page <https://github.com/GEOSX/GEOSX/issues>`_.
