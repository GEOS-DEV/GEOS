.. _AdvancedExampleWellboreNonLinearThermalDiffusionTemperatureDependentSinglePhaseThermalConductivity:


#####################################################################################################################
Non-Linear Thermal Diffusion Around a Wellbore: The Case with Temperature Dependent Single Phase Thermal Conductivity
#####################################################################################################################

------------------------------------------------------------------
Problem description
------------------------------------------------------------------

This example is an extension of the linear thermal diffusion problem presented in :ref:`AdvancedExamplePureThermalDiffusionWellbore` to model wellbore cooling upon CO2 injection. It uses the thermal single-phase flow solver to model a non-linear thermal diffusion problem around a wellbore where the single phase thermal conductivity of the porous rock depends linearly on the temperature.


**Input file**

This benchmark example uses no external input file and everything required is
contained within two GEOS xml files that are located at:

.. code-block:: console

  inputFiles/singlePhaseFlow/thermalCompressible_2d_base.xml

and

.. code-block:: console

  inputFiles/singlePhaseFlow/thermalCompressible_temperatureDependentSinglePhaseThermalConductivity_benchmark.xml


In this example, we focus on the ``Constitutive`` tag.

-----------------------------------------------------------
Constitutive
-----------------------------------------------------------

The reference value of the single phase thermal conductivity of the porous medium around the wellbore and its derivative with respect to temperature are defined in the ``SinglePhaseThermalConductivity`` XML block:  

.. literalinclude:: ../../../../../../../inputFiles/singlePhaseFlow/thermalCompressible_2d_base.xml
  :language: xml
  :start-after: <!-- SPHINX_SinglePhaseThermalConductivity_nonLinear -->
  :end-before: <!-- SPHINX_SinglePhaseThermalConductivity_nonLinearEnd -->


---------------------------------
Results and benchmark
---------------------------------

A good agreement between the results obtained using GEOS and the reference results that are obtained by the classical finite difference method is shown in the figure below:


.. plot:: docs/sphinx/advancedExamples/validationStudies/wellboreProblems/nonLinearThermalDiffusion_TemperatureDependentSinglePhaseThermalConductivity/temperatureDependentSinglePhaseThermalConductivity_plot.py


------------------------------------------------------------------
To go further
------------------------------------------------------------------

**Feedback on this example**

This concludes the example of a non-linear thermal diffusion problem around a wellbore due to temperature dependent single phase thermal conductivity of porous rock.
For any feedback on this example, please submit a `GitHub issue on the project's GitHub page <https://github.com/GEOS-DEV/GEOS/issues>`_.
