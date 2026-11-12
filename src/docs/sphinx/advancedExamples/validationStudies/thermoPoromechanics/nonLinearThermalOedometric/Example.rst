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
Sealed heating: pore pressure and porosity
------------------------------------------------------------------

The drained test above has no pore pressure, so it does not exercise the porosity. In the sealed variant
the same oedometric sample cannot drain: no fluid leaves it, and the uniform temperature rises from 0 to
90 with the kinked table above. The porosity then uses the coefficient of the ``SolidThermalExpansion``
model, like the stress; ``defaultPorosityTEC`` is not given (it is rejected when the model is attached):

.. literalinclude:: ../../../../../../../inputFiles/thermoPoromechanics/ThermoPoroElasticSealedHeating_base.xml
  :language: xml
  :start-after: <!-- SPHINX_SealedPorousSolid -->
  :end-before: <!-- SPHINX_SealedPorousSolidEnd -->

With :math:`I(T) = \int_{T_{init}}^{T} \alpha(s)\,ds`, Biot coefficient :math:`b`, grain bulk modulus
:math:`K_s`, and a fluid of density :math:`\rho_0 \exp( c_f p - \beta_f (T - T_{init}) )`, the vertical
strain, the porosity and the fluid mass balance read

.. math::
   \varepsilon_{zz} = \frac{3 K I + b\,p}{K + \frac{4}{3} G}, \qquad
   \phi = \phi_0 + b\,\varepsilon_{zz} + \frac{b - \phi_0}{K_s}\,p - 3 (b - \phi_0)\,I, \qquad
   \phi\, e^{c_f p - \beta_f (T - T_{init})} = \phi_0 .

The last equation gives the pressure at each temperature (one scalar root, solved to round-off), and
the effective stresses follow. GEOS matches this solution to :math:`10^{-13}` with the fully implicit
coupling and to :math:`2 \times 10^{-10}` with the sequential coupling, whose iterations stop at
``newtonTol``. The time step does not matter: the ten steps of the smoke decks give the same result as
any finer stepping. The dashed curve is the exact pressure when the porosity uses a constant coefficient
:math:`\alpha(0)` instead: the final pressure would be 30 % lower.

.. plot:: docs/sphinx/advancedExamples/validationStudies/thermoPoromechanics/nonLinearThermalOedometric/plotSealedHeating.py

.. code-block:: console

  inputFiles/thermoPoromechanics/ThermoPoroElasticSealedHeating_base.xml
  inputFiles/thermoPoromechanics/ThermoPoroElasticSealedHeating_fim_smoke.xml
  inputFiles/thermoPoromechanics/ThermoPoroElasticSealedHeating_sequential_smoke.xml

The data (``sealedHeating.csv``) is produced by ``postprocessSealedHeating.py``.

------------------------------------------------------------------
Mesh convergence: cooling of a column
------------------------------------------------------------------

The previous tests are uniform. Here a column of length 1 along :math:`x`, confined in :math:`y` and
:math:`z`, fixed at :math:`x = 0` and free at :math:`x = 1`, starts at 90; its end :math:`x = 0` is held
at 0 and the other end is insulated. The pressure is held at zero. The fluid and the solid have the
same volumetric heat capacity :math:`C`, so that the energy balance does not depend on porosity:

.. literalinclude:: ../../../../../../../inputFiles/thermoPoromechanics/ThermoElasticColumnCooling_base.xml
  :language: xml
  :start-after: <!-- SPHINX_ColumnHeatCapacities -->
  :end-before: <!-- SPHINX_ColumnHeatCapacitiesEnd -->

The temperature is then the series solution of the heat equation with diffusivity :math:`D = k / C`,

.. math::
   T(x,t) = 90 \sum_{n \ge 0} \frac{4}{(2n+1)\pi} \sin(\lambda_n x)\, e^{-\lambda_n^2 D t},
   \qquad \lambda_n = \frac{(2n+1)\pi}{2},

and, since :math:`\sigma_{xx} = 0` everywhere, the lateral stress at each point is
:math:`\sigma_{yy} = -\frac{E}{1-\nu} \int_{90}^{T(x,t)} \alpha(s)\,ds`. The column crosses all the
kinks of the table. ``convergenceThermoElasticColumn.py`` runs GEOS on 10 to 160 cells, with a time
step proportional to the square of the cell size, and fits the order of the error on the temperature,
the lateral stress and the end displacement at three times. All three converge at order 2.0.

.. plot:: docs/sphinx/advancedExamples/validationStudies/thermoPoromechanics/nonLinearThermalOedometric/plotConvergenceStudies.py column

------------------------------------------------------------------
Mesh convergence: steady cooling of a wellbore wall
------------------------------------------------------------------

A quarter of an annulus :math:`0.1 < r < 1`, in plane strain and drained, is held at 0 on its inner wall
(free) and at 90 on its outer boundary (fixed). At steady state the temperature is
:math:`T(r) = 90 \ln(r/a) / \ln(b/a)`, and with the thermal strain :math:`e(r) = \int_{90}^{T(r)} \alpha`,
the displacement is

.. math::
   u(r) = \beta \frac{J(r)}{r} + C_1 r + \frac{C_2}{r}, \qquad J(r) = \int_a^r e(s)\, s\, ds, \qquad
   \beta = \frac{3K}{\lambda + 2G},

with :math:`C_1, C_2` given by :math:`\sigma_{rr}(a) = 0` and :math:`u(b) = 0`; the stresses follow.
``convergenceThermoElasticWellbore.py`` runs GEOS on 8 x 8 to 64 x 64 cells: the temperature, the three
stresses and the wall displacement converge at order 2.0.

.. plot:: docs/sphinx/advancedExamples/validationStudies/thermoPoromechanics/nonLinearThermalOedometric/plotConvergenceStudies.py wellbore

.. code-block:: console

  inputFiles/thermoPoromechanics/ThermoElasticColumnCooling_base.xml
  inputFiles/thermoPoromechanics/ThermoElasticColumnCooling_smoke.xml
  inputFiles/thermoPoromechanics/ThermoElasticWellboreSteadyCooling_base.xml
  inputFiles/thermoPoromechanics/ThermoElasticWellboreSteadyCooling_smoke.xml

The two scripts write the CSV files of the figures and exit with an error if an observed order falls
below 1.8.

------------------------------------------------------------------
To go further
------------------------------------------------------------------

**Feedback on this example**

For any feedback on this example, please submit a
`GitHub issue on the project's GitHub page <https://github.com/GEOS-DEV/GEOS/issues>`_.
