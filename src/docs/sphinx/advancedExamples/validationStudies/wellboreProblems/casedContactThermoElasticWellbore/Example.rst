.. _AdvancedExampleCasedContactThermoElasticWellbore:


######################################################################
Cased ThermoElastic Wellbore Problem with Imperfect Contact Interfaces
######################################################################

**Input file**

This benchmark example uses no external input files and everything required is
contained within two GEOS xml files that are located at:

.. code-block:: console

  inputFiles/wellbore/CasedThermoElasticWellbore_ImperfectInterfaces_base.xml

and

.. code-block:: console

  inputFiles/wellbore/CasedThermoElasticWellbore_ImperfectInterfaces_benchmark.xml

The corresponding integrated test is

.. code-block:: console

  inputFiles/wellbore/CasedThermoElasticWellbore_ImperfectInterfaces_smoke.xml

---------------------------------
Results and benchmark
---------------------------------

The GEOS results of displacement jump across the casing-cement and cement-rock interfaces are shown in the figure below: 

.. _CasedThermoElasticWellboreInterfacesThermalDebondingFig:
.. figure:: thermalDebonding.png
   :align: center
   :width: 600
   :figclass: align-center

   Displacement jumps across the casing-cement and cement-rock interfaces

The GEOS results and analytical results for temperature distribution around the cased wellbore is shown in the figures below:

.. plot:: docs/sphinx/advancedExamples/validationStudies/wellboreProblems/casedContactThermoElasticWellbore/thermoElastic_casedContactWellbore_temperature.py

.. _problemCasedContactThermoElasticWellbore_Temperature_Fig:
.. figure:: temperature.png
   :align: center
   :width: 800
   :figclass: align-center

   Temperature field.

and the radial displacement around the wellbore is shown below:

.. plot:: docs/sphinx/advancedExamples/validationStudies/wellboreProblems/casedContactThermoElasticWellbore/thermoElastic_casedContactWellbore_displacement.py

.. _problemCasedContactThermoElasticWellbore_Displacement_Fig:
.. figure:: displacement.png
   :align: center
   :width: 800
   :figclass: align-center

   The displacement field.

The total radial and hoop stress (tangent stress) components computed by GEOS and reference results are shown in the figure below:

.. plot:: docs/sphinx/advancedExamples/validationStudies/wellboreProblems/casedContactThermoElasticWellbore/thermoElastic_casedContactWellbore_stress.py

.. _problemCasedContactThermoElasticWellbore_Stresses_Fig:
.. figure:: stress.png
   :align: center
   :width: 800
   :figclass: align-center

   The stress field.

------------------------------------------------------------------
To go further
------------------------------------------------------------------

**Feedback on this example**

This concludes the cased wellbore example.
For any feedback on this example, please submit a `GitHub issue on the project's GitHub page <https://github.com/GEOS-DEV/GEOS/issues>`_.
