.. _AdvancedExampleNonLinearThermalOedometric:


#################################################################################
Thermo-elastic oedometric test with a temperature-dependent expansion coefficient
#################################################################################

------------------------------------------------------------------
Problem description
------------------------------------------------------------------

This example validates the temperature-dependent thermal expansion coefficient model
:ref:`SolidThermalExpansion <TemperatureDependentThermalExpansionCoefficient>` on a drained
thermo-elastic oedometric test. A cube is fixed in the lateral directions :math:`x` and :math:`y`,
fixed at its base and free on its top face. The pore pressure is held at zero and gravity is
neglected. A uniform temperature is imposed, increasing linearly from 0 to 90 over
:math:`[0, 10^5]` s (heating), or decreasing from 90 to 0 (cooling).

Since the lateral strains vanish and the top face is free, the lateral stress is

.. math::
   \sigma_{xx} = \sigma_{yy} = - \frac{E}{1-\nu} \int_{T_{init}}^{T} \alpha(s) \, ds

where :math:`E` and :math:`\nu` are the Young's modulus and Poisson's ratio. For the linear law
:math:`\alpha(T) = \alpha_0 + \frac{d\alpha}{dT} (T - T_0)`, the mean value of :math:`\alpha` over
:math:`[T_{init}, T]` is its value at the midpoint, so that

.. math::
   \sigma_{xx} = - \frac{E}{1-\nu} \, \alpha\!\left( \frac{T_{init} + T}{2} \right) (T - T_{init}).

This analytical solution is used to validate the GEOS results.

**Input files**

Each case uses a base file and a smoke or benchmark file that defines the mesh and the time
stepping:

.. code-block:: console

  inputFiles/thermoPoromechanics/ThermoElasticOedometric_linearTEC_base.xml
  inputFiles/thermoPoromechanics/ThermoElasticOedometric_linearTEC_benchmark.xml
  inputFiles/thermoPoromechanics/ThermoElasticOedometric_linearTEC_smoke.xml
  inputFiles/thermoPoromechanics/ThermoElasticOedometric_tecTable_base.xml
  inputFiles/thermoPoromechanics/ThermoElasticOedometric_tecTable_benchmark.xml
  inputFiles/thermoPoromechanics/ThermoElasticOedometric_tecTable_smoke.xml
  inputFiles/thermoPoromechanics/ThermoElasticOedometric_linearTEC_cooling_benchmark.xml
  inputFiles/thermoPoromechanics/ThermoElasticOedometric_linearTEC_cooling_smoke.xml

-----------------------------------------------------------
Constitutive
-----------------------------------------------------------

The porous solid refers to the thermal expansion model through ``solidThermalExpansionModelName``:

.. literalinclude:: ../../../../../../../inputFiles/thermoPoromechanics/ThermoElasticOedometric_linearTEC_base.xml
  :language: xml
  :start-after: <!-- SPHINX_PorousSolidProperties -->
  :end-before: <!-- SPHINX_PorousSolidPropertiesEnd -->

The linear law is defined as:

.. literalinclude:: ../../../../../../../inputFiles/thermoPoromechanics/ThermoElasticOedometric_linearTEC_base.xml
  :language: xml
  :start-after: <!-- SPHINX_SolidThermalExpansionLinear -->
  :end-before: <!-- SPHINX_SolidThermalExpansionLinearEnd -->

The tabulated variant refers to a ``TableFunction`` of temperature:

.. literalinclude:: ../../../../../../../inputFiles/thermoPoromechanics/ThermoElasticOedometric_tecTable_base.xml
  :language: xml
  :start-after: <!-- SPHINX_SolidThermalExpansionTabulated -->
  :end-before: <!-- SPHINX_SolidThermalExpansionTabulatedEnd -->

.. literalinclude:: ../../../../../../../inputFiles/thermoPoromechanics/ThermoElasticOedometric_tecTable_base.xml
  :language: xml
  :start-after: <!-- SPHINX_TECTable -->
  :end-before: <!-- SPHINX_TECTableEnd -->

The table is not linear: its slope changes at 20 and 50, and it is constant above its last point (70),
up to the final temperature of 90. A time step that crosses one of these points is integrated exactly,
segment by segment. The analytical solution integrates the same piecewise-linear coefficient.

--------------------------------------------------------------------
FieldSpecifications
--------------------------------------------------------------------

The uniform temperature is imposed with a time-dependent table:

.. literalinclude:: ../../../../../../../inputFiles/thermoPoromechanics/ThermoElasticOedometric_linearTEC_base.xml
  :language: xml
  :start-after: <!-- SPHINX_TemperatureBC -->
  :end-before: <!-- SPHINX_TemperatureBCEnd -->

.. literalinclude:: ../../../../../../../inputFiles/thermoPoromechanics/ThermoElasticOedometric_linearTEC_base.xml
  :language: xml
  :start-after: <!-- SPHINX_TemperatureTable -->
  :end-before: <!-- SPHINX_TemperatureTableEnd -->

---------------------------------
Results
---------------------------------

The figure below compares the GEOS lateral stress with the analytical solution for the three cases.
The two agree to round-off: the thermal strain increment of each time step is the exact integral of
:math:`\alpha` over the step, so the result does not depend on the time step. In the cooling case
the lateral stress becomes tensile.

.. plot:: docs/sphinx/advancedExamples/validationStudies/thermoPoromechanics/nonLinearThermalOedometric/plotThermoElasticOedometric.py

The data of the figure (``thermoElasticOedometric.csv``) is produced from the benchmark runs by
``postprocessThermoElasticOedometric.py``, located next to this page.

------------------------------------------------------------------
To go further
------------------------------------------------------------------

**Feedback on this example**

For any feedback on this example, please submit a
`GitHub issue on the project's GitHub page <https://github.com/GEOS-DEV/GEOS/issues>`_.
