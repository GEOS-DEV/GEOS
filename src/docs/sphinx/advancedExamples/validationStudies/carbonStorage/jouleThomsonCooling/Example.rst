.. |oslash| unicode:: U+00F8

.. |degree| unicode:: U+00B0

.. _jouleThomsonCoolingExample:

##########################################################################
Cooling during CO2 Injection into a Gas Reservoir
##########################################################################


**Context**

When a real gas expands from a high to a low pressure at constant enthalpy, its temperature changes. This is the Joule-Thomson (JT) effect. For CO2 at reservoir conditions the gas cools when it expands. During Carbon Sequestration with Enhanced Gas Recovery, CO2 is injected into a depleted natural gas reservoir and the pressure drops steeply between the wellbore and the reservoir, so JT cooling may occur close to the injection well. If it were large, it could lead to CO2 or CH4 hydrate formation, to freezing of the residual water (with a loss of injectivity), or to thermal stresses able to fracture the formation.

`(Oldenburg, 2007) <https://doi.org/10.1016/j.enconman.2007.01.010>`__ uses the TOUGH2/EOS7C simulator to quantify the magnitude of this cooling. The model is first verified against the classical JT expansion experiment, and is then applied to an idealized one-dimensional radial gas reservoir with constant-rate CO2 injection. The study shows that the cooling grows when the permeability decreases (a larger pressure drop is required to inject the same mass of CO2), and that it is partly compensated by the heat released by the rock matrix and the residual water.

This example reproduces the lowest-permeability case of the paper, :math:`k = 5\times 10^{-15}` m\ :sup:`2` (5 mD), whose results are shown in the paper. It is the case with the strongest cooling among the constant-rate injections with an injection temperature equal to the reservoir temperature, and therefore the most demanding one for the thermal formulation. It can serve as a validation of the thermal compositional multiphase flow solver of GEOS, in particular of the coupling between

- a real-gas equation of state, which provides the enthalpy of the CO2-CH4-H2O mixture and thus the JT coefficient,
- the energy balance, including the advection of enthalpy, the heat stored in the rock and in the immobile water, and
  conduction,
- thermal wells with prescribed mass and volumetric rates.

------------------------------------------------------------------------
Brief case description
------------------------------------------------------------------------

**Joule-Thomson cooling**

In the classical JT experiment, a gas flows through a porous plug in a thermally insulated system. The expansion is isenthalpic and the temperature drop is proportional to the pressure drop through the JT coefficient :math:`\mu_{JT}`,

.. math::

    \Delta T \approx \mu_{JT} \, \Delta P, \qquad \mu_{JT} = \left(\frac{\partial T}{\partial P}\right)_H = \frac{1}{c_p}\left[ T \left(\frac{\partial v}{\partial T}\right)_P - v \right],

where :math:`c_p` and :math:`v` are the specific heat capacity and the specific volume of the gas. For CO2 at 45 \ |degree|\  C (318.15 K) and 50 bar, :math:`\mu_{JT}` is approximately 0.9 K/bar according to the NIST reference equation of state (`Lemmon et al., 1998 <https://webbook.nist.gov/chemistry/fluid/>`__), and it decreases when the pressure increases (approximately 0.4 K/bar at 100 bar). The product :math:`\mu_{JT}\,\Delta P` is thus the maximum cooling that can be obtained for a given pressure drop. In a reservoir, this maximum is never reached: the rock grains and the residual water, with a volumetric heat capacity much larger than that of the gas, release heat as the gas cools down. This is a transient effect: the matrix progressively loses its heat during a prolonged expansion. The cooling therefore depends not only on the pressure drop, but also on the porosity, the residual water saturation, the heat capacity of the grains and the thermal conductivity.

**Reservoir model of the paper**

The paper considers an idealized depleted gas reservoir, with the geometry and the properties given in the paper. It is represented by a one-dimensional radially symmetric system, 50 m thick and 1130 m in radius, as illustrated below.

.. code-block:: console

    injection well (r = 0)  --->  CO2-rich gas  | CH4-rich gas  --->  production (r ~ 1100 m)

The reservoir is initially filled with methane at 50 bar and 45 \ |degree|\ C. The water is immobile (the water saturation is equal to the residual saturation :math:`S_{lr} = 0.2`), and the only mobile phase is the gas. Pure CO2 is injected at a constant rate of 3 kg/s and at the reservoir temperature. Methane is produced at the outer boundary, at 0.56 kg/s, to represent the production of four surrounding wells located at a distance of 1100 m in a five-spot pattern. The system is closed: no fluid or heat crosses the outer boundary. Because the injection temperature is equal to the initial reservoir temperature, every temperature change in the reservoir is caused by the JT effect (attenuated by the heat provided by the rock and the water), which makes this setup well suited to isolate and validate this effect.

The parameters of the problem are summarized in the following table.

.. list-table:: Parameters of the problem
    :widths: 30 70
    :header-rows: 1
    :class: align-top

    * - Parameter
      - Value
    * - Geometry
      - 1D radial, :math:`R = 1130` m, thickness 50 m

        GEOS is modeled as a 1\ |degree|\  sector, :math:`r \in [0.1, 1130]` m.
    * - Permeability
      - :math:`5\times 10^{-15}` m\ :sup:`2` (5 mD)
    * - Porosity
      - 0.30
    * - Rock volumetric heat capacity
      - :math:`2.6\times 10^{6}` J/m\ :sup:`3`/K
    * - Thermal conductivity
      - 2.51 W/m/K
    * - Residual saturations
      - :math:`S_{lr} = 0.2` (water), :math:`S_{gr} = 0.01` (gas)
    * - Relative permeability
      - Corey curves
    * - Initial pressure
      - 50 bar (:math:`5\times 10^{6}` Pa)
    * - Initial temperature
      - 45 \ |degree|\ C (318.15 K)
    * - Initial fluids
      - CH4 gas, water at :math:`S_l = S_{lr} = 0.2`

        GEOS initializes this via mass fractions: :math:`z_{CO_2} = 0`, :math:`z_{CH_4} = 0.115401`, :math:`z_{H_2O} = 0.884599`.
    * - CO2 injection
      - 3 kg/s of CO2 at 45 \ |degree|\ C
    * - CH4 production
      - 0.56 kg/s at the outer boundary (1100 m)
    * - Simulated time
      - 15 years
    * - Equation of state
      - Peng-Robinson for the gas phase

        The paper uses EOS7C/GasEOS, while GEOS explicitly uses S\ |oslash|\ reide-Whitson for the aqueous phase.

..
  **What the paper finds for this case**

..
  For a constant injection rate, the pressure build-up near the well scales inversely with the permeability. The paper reports a cooling of about 2 \ |degree|\ C for 100 mD, 8 \ |degree|\ C for 10 mD, and about 14 \ |degree|\ C for 5 mD (Figure 8) after 15 years of injection. For 5 mD (Figure 7), the cooled region extends over more than 100 m after 15 years and it remains inside the region occupied by nearly pure CO2. The minimum temperature remains around 31 \ |degree|\ C, well above the temperatures at which CO2 hydrates are stable at these pressures (below approximately 10 \ |degree|\ C), so that the problem is not expected to be a concern for CSEGR in this setting.

------------------------------------------------------------------------
Input base
------------------------------------------------------------------------

This benchmark test is based on the XML file located below:

.. code-block:: console

  inputFiles/compositionalMultiphaseFlow/benchmarks/jouleThomsonCooling/jouleThomsonCooling.xml

The mesh is a thin wedge of cylinder (``InternalWellbore``), 1\ |degree|\  wide (``theta`` from -0.5\ |degree|\  to 0.5\ |degree|\ ). The 109 radial cells are defined explicitly through the ``radius`` and ``nr`` arrays: the cell size grows geometrically with the distance to the injection well (the ratio between successive radii decreases from 1.2 close to the well to about 1.08 at the outer boundary), which refines the mesh where the pressure and temperature gradients are the strongest. The 360\ |degree|\  of the paper are accounted for by dividing all the well rates by 360.

.. literalinclude:: ../../../../../../../inputFiles/compositionalMultiphaseFlow/benchmarks/jouleThomsonCooling/jouleThomsonCooling.xml
    :language: xml
    :start-after: <!-- SPHINX_MESH -->
    :end-before: <!-- SPHINX_MESH_END -->

The mesh is only one cell thick in the angular and vertical directions, so the flow is purely radial. The two wells are vertical and are located in the mid-plane of the wedge (:math:`\theta = 0^{\circ}`): the injector is at :math:`r = 0.101` m, in the first cell of the mesh, and the producer is at :math:`r = 1129.5` m, in the outermost cell, which ends at the outer boundary at 1130 m. As in the paper, the production is therefore taken from the outermost grid block. Since there is a single cell over the thickness, the position of the perforation along the well (``distanceFromHead``) has no influence on the solution. The well-to-reservoir transmissibility of each perforation is prescribed explicitly (``transmissibility``). The two wells are listed in the ``ElementRegions`` together with the reservoir region. All the regions share the same fluid, relative permeability and thermal conductivity models:

.. literalinclude:: ../../../../../../../inputFiles/compositionalMultiphaseFlow/benchmarks/jouleThomsonCooling/jouleThomsonCooling.xml
    :language: xml
    :start-after: <!-- SPHINX_REGIONS -->
    :end-before: <!-- SPHINX_REGIONS_END -->

------------------------------------------------------------------------
Constitutive models
------------------------------------------------------------------------

**Fluid model**

The fluid is a mixture of three components (CO2, CH4 and H2O) distributed in two phases (an aqueous phase named ``water`` and a ``gas`` phase) through a two-phase flash. It is defined by the ``CompositionalTwoPhaseThermalFluidPhillipsBrine`` model:

.. literalinclude:: ../../../../../../../inputFiles/compositionalMultiphaseFlow/benchmarks/jouleThomsonCooling/jouleThomsonCooling.xml
    :language: xml
    :start-after: <!-- SPHINX_FLUID -->
    :end-before: <!-- SPHINX_FLUID_END -->
    :emphasize-lines: 4,10-19

The main ingredients of this model are the following:

- ``equationsOfState`` lists one equation of state for each phase, in the order of ``phaseNames``. The gas phase uses the   Peng-Robinson (PR) cubic equation of state, as in the paper (EOS7C computes the properties of the CO2-CH4-H2O gas mixtures with PR). The aqueous phase uses the S\ |oslash|\ reide-Whitson variant of PR, which modifies the attractive term for water and uses salinity-dependent interaction parameters, in order to represent the mutual solubilities of the gases and of water. The aqueous phase stays at its residual saturation in this problem and acts mostly as a heat reservoir, so that the solubility effects are weak, but the model remains fully compositional.
- the component data (``componentMolarWeight`` in kg/mol, critical pressures and temperatures, acentric factors) are the standard values of CO2, CH4 and H2O, given in this order, which is the order of ``componentNames``. The same order is used for the ``injectionStream``, the ``globalCompFraction`` components and the binary and heat capacity coefficients. 
- the aqueous phase density and viscosity follow the brine correlations of Phillips et al., with ``salinity`` equal to 0 (fresh water).
- ``pressureCoordinates`` and ``temperatureCoordinates`` define the pressure-temperature grid on which the fluid properties are evaluated. They must bracket all the states visited during the simulation: here from 1 bar (the minimum bottom-hole pressure of the producer) to 300 bar (the maximum bottom-hole pressure of the injector), and from -10 \ |degree|\ C to 110 \ |degree|\ C. The states actually visited are between about 5 and 7 MPa, and between 31 and 45 \ |degree|\ C.

.. _jouleThomsonCoolingHeatCapacity:

**Heat capacity coefficients**

The enthalpy of a phase is the sum of an ideal-gas contribution and of a departure (residual) contribution given by the equation of state, which accounts for the non-ideal behavior of the real gas. The ideal-gas part is obtained by integrating the isobaric ideal-gas heat capacity of each component, which is provided by ``componentHeatCapacityCoefficients`` (one row per component, in the order of ``componentNames``):

.. literalinclude:: ../../../../../../../inputFiles/compositionalMultiphaseFlow/benchmarks/jouleThomsonCooling/jouleThomsonCooling.xml
    :language: xml
    :start-after: <!-- SPHINX_FLUID -->
    :end-before: <!-- SPHINX_FLUID_END -->
    :lines: 15-19
    :dedent: 6

Each row contains the five coefficients :math:`a_0, \ldots, a_4` of a polynomial in the temperature difference
with respect to 25 \ |degree|\ C:

.. math::

    c_p^{0}(T) = \sum_{i=0}^{4} a_i \, (T - 298.15)^i \qquad [\mathrm{J/mol/K}], \quad T \text{ in K}.

The first coefficient :math:`a_0` is therefore the ideal-gas molar heat capacity at 25 \ |degree|\ C (37.14 J/mol/K for CO2, 35.71 for CH4 and 33.59 for H2O). The following table gives the resulting values at the reservoir temperature, which are in agreement with the reference ideal-gas heat capacities of the NIST equations of state within 0.01 J/mol/K over the whole temperature range of the fluid tables (-10 to 110 \ |degree|\ C).

.. list-table:: Ideal-gas heat capacity :math:`c_p^0` at 318.15 K
    :widths: 20 20 20 20
    :header-rows: 1

    * - Component
      - CO2
      - CH4
      - H2O
    * - :math:`c_p^0` (J/mol/K)
      - 38.04
      - 36.51
      - 33.69

These coefficients are important for two reasons. First, the JT coefficient is inversely proportional to the heat capacity (see the formula above), and the heat capacity of the gas, which is the sum of the ideal-gas part and of the real-gas departure, controls how much the gas cools for a given enthalpy variation. Second, the heat capacity sets the thermal inertia of the fluids against that of the rock, which determines how fast the cooling is attenuated by the matrix.

.. _jouleThomsonBinaryCoeff:

**Binary interaction parameters**

The Peng-Robinson equation of state uses a mixing rule for the attractive term of the mixture in which the cross term of each pair of components :math:`(i,j)` is scaled by :math:`(1 - k_{ij})`. The binary interaction parameters :math:`k_{ij}` are given by ``componentBinaryCoeff`` as a symmetric matrix with a zero diagonal, in the order CO2, CH4, H2O (the same order as ``componentNames``), so that :math:`k_{CO_2-CH_4} = 0.1`, :math:`k_{CO_2-H_2O} = 0.1896` and :math:`k_{CH_4-H_2O} = 0.4850`. The CO2-CH4 coefficient is the one that matters most for this problem, since the gas phase is a CO2-CH4 mixture in the region where the CO2 front progresses and displaces the methane. The coefficients involving water are larger, and they only affect the small amount of water vapor in the gas phase at 45 \ |degree|\ C and the solubility of the gases in the aqueous phase, which stays at its residual saturation.

**Rock, relative permeability and thermal conductivity**

The rock is modeled as an incompressible porous medium (``compressibility`` = 0 in ``rockPorosity``) with a constant permeability and a constant porosity of 0.3. Its heat capacity is given as
a volumetric value, :math:`\rho_R C_R = 2600 \times 1000 = 2.6 \times 10^{6}` J/m\ :sup:`3`/K, which is the product
of the grain density and the grain heat capacity of the paper. The matrix contributes :math:`(1-\phi)\rho_R C_R` per unit of bulk volume to the
energy accumulation, which is the term that attenuates the JT cooling in the paper.

.. literalinclude:: ../../../../../../../inputFiles/compositionalMultiphaseFlow/benchmarks/jouleThomsonCooling/jouleThomsonCooling.xml
    :language: xml
    :start-after: <!-- SPHINX_ROCK -->
    :end-before: <!-- SPHINX_ROCK_END -->

The relative permeabilities are Brooks-Corey power laws with the residual saturations of the paper. With a water residual saturation of 0.2
(``phaseMinVolumeFraction``) and an initial water saturation equal to 0.2, the water relative permeability is zero and the water does not move:

.. literalinclude:: ../../../../../../../inputFiles/compositionalMultiphaseFlow/benchmarks/jouleThomsonCooling/jouleThomsonCooling.xml
    :language: xml
    :start-after: <!-- SPHINX_RELPERM -->
    :end-before: <!-- SPHINX_RELPERM_END -->

The thermal conductivity of the formation is constant and isotropic:

.. literalinclude:: ../../../../../../../inputFiles/compositionalMultiphaseFlow/benchmarks/jouleThomsonCooling/jouleThomsonCooling.xml
    :language: xml
    :start-after: <!-- SPHINX_THC -->
    :end-before: <!-- SPHINX_THC_END -->

.. note::
    The *GEOS* model is not strictly identical to the TOUGH2/EOS7C model of the paper. The main differences are:

    - the Corey relative permeabilities of the paper are approximated by Brooks-Corey laws with exponents 4 and 2 (the water curve
      is the same, the gas curve differs slightly);
    - the capillary pressure (van Genuchten) and the molecular diffusion of the paper (Table 2) are not included in the
      *GEOS* input. They are expected to have little influence, since the water remains at its residual saturation, but
      they are a difference with the published model;
    - the aqueous phase uses the S\ |oslash|\ reide-Whitson equation of state, whereas EOS7C uses PR for the whole system;
    - the producer is a well with a prescribed surface gas rate, located at 1129.5 m in the outermost grid block, instead of a source term in that block.

------------------------------------------------------------------------
Flow solver
------------------------------------------------------------------------

The flow is solved by ``CompositionalMultiphaseFVM``, with a TPFA discretization, on the reservoir region only:

.. literalinclude:: ../../../../../../../inputFiles/compositionalMultiphaseFlow/benchmarks/jouleThomsonCooling/jouleThomsonCooling.xml
    :language: xml
    :start-after: <!-- SPHINX_TPFA -->
    :end-before: <!-- SPHINX_TPFA_END -->

.. literalinclude:: ../../../../../../../inputFiles/compositionalMultiphaseFlow/benchmarks/jouleThomsonCooling/jouleThomsonCooling.xml
    :language: xml
    :start-after: <!-- SPHINX_FLOW_SOL -->
    :end-before: <!-- SPHINX_FLOW_SOL_END -->

The ``isThermal`` flag activates the energy balance equation, so that the temperature is an unknown of the problem. The
``useMass`` flag selects a mass-based formulation: the compositions and the well rates are expressed as mass
fractions and mass rates (kg/s). The time step is controlled by the maximum relative changes of pressure (20 %), temperature (2 %)
and phase volume fraction (20 %) per time step, and ``maxCompFractionChange`` limits the update of the compositions during
the Newton iterations.

The reservoir and the wells are coupled in a fully implicit way by ``CompositionalMultiphaseReservoir``. The problem is small
(one-dimensional), so a direct linear solver is used:

.. literalinclude:: ../../../../../../../inputFiles/compositionalMultiphaseFlow/benchmarks/jouleThomsonCooling/jouleThomsonCooling.xml
    :language: xml
    :start-after: <!-- SPHINX_RES_SOL -->
    :end-before: <!-- SPHINX_RES_SOL_END -->

------------------------------------------------------------------------
Initial and boundary conditions
------------------------------------------------------------------------

.. _jouleThomsonInitialState:

**Initial state**

The initial state is uniform over the whole reservoir and is specified with ``FieldSpecification`` tags:

.. literalinclude:: ../../../../../../../inputFiles/compositionalMultiphaseFlow/benchmarks/jouleThomsonCooling/jouleThomsonCooling.xml
    :language: xml
    :start-after: <!-- SPHINX_INIT -->
    :end-before: <!-- SPHINX_INIT_END -->

It is chosen to match the initial conditions of the paper, for the following reasons.

- **Pressure, 50 bar.** This is the initial pressure of the depleted reservoir. It is low compared to the pressure of
  the CO2 supply, which is what creates a large pressure drop around the injector. The values are uniform (no
  ``HydrostaticEquilibrium`` is required) since the model is one cell thick and the problem is
  controlled by the radial flow.
- **Temperature, 45 °C (318.15 K).** This is the reservoir temperature of the paper, and it is the same as the injection
  temperature. The initial thermal state is therefore in equilibrium with the injected fluid, and **any temperature
  deviation observed afterwards is the result of the JT expansion** (and of its attenuation by the heat of the rock),
  not of a thermal contrast between the injected fluid and the reservoir.
- **Composition.** The reservoir contains methane and water, with no CO2. In *GEOS*, the saturations are not primary inputs: the initial
  phase distribution is computed by a flash from the global composition, pressure and temperature. The paper
  prescribes an initial water saturation equal to the residual saturation, :math:`S_w = 0.2`, so that the water is immobile and the only mobile
  phase is gas. The global mass fractions are calibrated to obtain exactly this saturation. At 50 bar and 45 \ |degree|\ C, the density of CH4 gas
  is :math:`\rho_g = 32.36` kg/m\ :sup:`3` and the density of water is :math:`\rho_w = 992.3` kg/m\ :sup:`3`, so that

  .. math::

      z_{CH_4} = \frac{\rho_g (1 - S_w)}{\rho_g (1 - S_w) + \rho_w S_w} = \frac{32.36 \times 0.8}{32.36 \times 0.8 + 992.3 \times 0.2} \approx 0.1154,
      \qquad z_{H_2O} = 1 - z_{CH_4} \approx 0.8846,

  which are the values of ``initialCH4`` and ``initialH2O`` (``initialCO2`` is the CO2 mass fraction and is equal to zero). The mutual
  solubility between the phases is neglected in this estimate, and its effect on the saturation is very small.

.. note::
    These mass fractions are specific to the initial pressure and temperature. If these values are modified (for example
    to reproduce the cases of the paper with a reservoir at 75 \ |degree|\ C), the densities change and the compositions must be
    computed again to keep :math:`S_w = 0.2`.

**Boundary conditions**

No boundary condition is imposed: all the outer boundaries of the domain are closed for the mass and for the heat
(no-flow and no heat flux). This corresponds to the closed reservoir of the paper, in which the pressure
increases with the injection and the heat is exchanged only between the fluids and the rock. The only sources and sinks
are the two wells.

.. _jouleThomsonCoolingWellControls:

**Well controls**

The paper specifies constant-rate injection of CO2 and constant-rate production of CH4. These are the controls that govern the
pressure drop between the wells and the reservoir and therefore the JT cooling. The injector is defined as follows:

.. literalinclude:: ../../../../../../../inputFiles/compositionalMultiphaseFlow/benchmarks/jouleThomsonCooling/jouleThomsonCooling.xml
    :language: xml
    :start-after: <!-- SPHINX_INJ_CONTROL -->
    :end-before: <!-- SPHINX_INJ_CONTROL_END -->

and the producer as:

.. literalinclude:: ../../../../../../../inputFiles/compositionalMultiphaseFlow/benchmarks/jouleThomsonCooling/jouleThomsonCooling.xml
    :language: xml
    :start-after: <!-- SPHINX_PROD_CONTROL -->
    :end-before: <!-- SPHINX_PROD_CONTROL_END -->

The correspondence with the paper is the following:

.. list-table:: Well controls and their relation to the paper
    :widths: 20 35 45
    :header-rows: 1

    * - Item
      - *GEOS* input
      - Relation to the paper
    * - Injector control
      - ``control="massRate"`` with ``InjectionMassRateConstraint``, ``massRate`` = 8.3333e-3 kg/s
      - Constant-rate injection of 3 kg/s. The rate is divided by 360 because the model is a 1\ |degree|\  sector
        (:math:`3 \times 1/360`).
    * - Injected fluid
      - ``injectionStream`` = {0.9999, 1e-4, 0}, ``injectionTemperature`` = 318.15 K
      - Pure CO2 at :math:`T_{inj}` = 45 \ |degree|\ C, equal to the reservoir temperature. The stream is
        99.99 % CO2 (in mass), with a trace of CH4 and no H2O, which is equivalent to pure CO2 for this problem.
    * - Injector safeguard
      - ``MaximumBHPConstraint``, 300 bar at the mid-height of the reservoir (``referenceElevation`` = 25 m)
      - Not in the paper. This limit is far from the pressures reached (about 70 bar near the well), so
        it never becomes active, and the injection remains at constant rate during the whole simulation. It only
        prevents an unphysical pressure if the model is modified (for example with a lower permeability).
    * - Producer control
      - ``control="phaseVolRate"`` for the ``gas`` phase, :math:`2.28872\times 10^{-3}` m\ :sup:`3`/s at surface conditions
        (101325 Pa, 288.15 K)
      - Constant production of 0.56 kg/s of CH4. In the sector, this is :math:`0.56/360 = 1.556 \times 10^{-3}` kg/s,
        which is converted to a volumetric rate using the density of methane at surface conditions,
        :math:`0.68` kg/m\ :sup:`3`.
    * - Producer safeguard
      - ``MinimumBHPConstraint``, 1 bar (``referenceElevation`` = 25 m)
      - Not in the paper. It is inactive here since the reservoir pressure is above 50 bar.
    * - Well status
      - ``statusTableName="wellStatus"``
      - The ``wellStatus`` table function keeps both wells closed for :math:`t < 0` and opens them at :math:`t = 0`, which is the start of the injection.
    * - Crossflow
      - ``enableCrossflow="0"`` for both wells
      - The injector can only inject and the producer can only produce.

**Role of the** ``wellStatus`` **table.** Both wells refer to the same table function, ``wellStatus``, through the ``statusTableName`` attribute. The table gives the status of the well as a function of time: a value of 1 means that the well is open and a value of 0 that it is closed. Here the table is defined by three points (:math:`t = -10^8`, 0 and :math:`10^8` s, with the values 0, 1 and 1) and a ``lower`` (piecewise constant) interpolation, so that the wells are closed for :math:`t < 0` and open from :math:`t = 0` until the end of the simulation. The injection and the production therefore start simultaneously at the beginning of the solver applications, with the reservoir at its uniform initial state. The same mechanism can be used to schedule shut-ins or delayed starts of a well without modifying the well controls, for instance to run a variant of the case in which the producer is opened later.

The mass rate of the producer in the paper and the surface volumetric rate used here are equivalent as long as the
produced gas is methane. A simple volumetric balance indicates that, after 15 years, the CO2 front remains at about 500-600 m
from the injector, so the producer, at 1129.5 m, only produces methane during the simulation, and the rate is
constant in mass. The ``writeCSV`` flag exports the time history of the rates, of the bottom-hole pressure and of the
temperature of each well, which is useful to monitor the injection overpressure.

------------------------------------------------------------------------
Events and output
------------------------------------------------------------------------

The simulation covers 15 years (a year is 3.15576e7 s). The events are defined in the ``Events`` block, whose ``minTime`` and ``maxTime`` are -1 day and slightly more than 15 years, so that the last output event is reached. They are of five kinds:

- **Solver applications.** The coupled reservoir and well solver is applied by two periodic events. During the first 10 s the maximum time step is 1 s, to handle the sudden start of the injection and the initial pressure transient near the well. After that, the maximum time step is 30 days (2.592e6 s). In both cases the actual time step is further controlled by the targets defined in the flow solver (maximum relative changes of pressure, temperature and phase volume fraction).
- **VTK outputs.** Eight solo events write the full fields at the times of the curves shown in Figure 7 of the paper: the initial state, then 1, 10, 20, 30 and 90 days, 1 year and 15 years. They use ``targetExactTimestep`` so that the time step is cut to land exactly on the output time, and they are listed before the solver events so that they record the state at their target time.
- **Statistics.** The ``CompositionalMultiphaseStatistics`` task (pressure, temperature, phase and component masses, CFL numbers) is executed at the initial time and then once per year.
- **Time history.** Around each of the seven output times after the initial state, a periodic event triggers the ``PackCollection`` tasks that collect the cell centers, the pressure, the temperature and the phase compositions. The collection is hourly, over a 3-hour window starting one hour before the output time, and the collected data are written once by the ``timeHistory`` output at 15 years. These files provide the radial profiles used for the comparison with the paper.
- **Restart.** Restart files are written every 5 years.

.. note::
    The time tags in *GEOS* are in seconds. The ``maxTime`` of the events is slightly larger than 15 years, so that the
    last output event is reached.

------------------------------------------------------------------------
Post-processing and reporting
------------------------------------------------------------------------

The profiles of the CO2 mass fraction in the gas phase, of the pressure and of the temperature along the radial direction are post-processed with the Python script ``plot_joule_thomson_profile.py``, located next to the input file. It reads the time history written by the ``timeHistory`` output (``jouleThomsonCooling.hdf5``), selects the output times of the simulation (1, 10, 20, 30 and 90 days, 1 year and 15 years, with 1 year = 365.25 days), and plots the profiles of the first row of the wedge on a logarithmic radial axis from 0.1 to 1000 m. For the temperature, the script also overlays the points digitized from Figure 7 of the paper (``oldenburg.json``), shown as dots with the same color as the corresponding *GEOS* curve. There is no reference data for the pressure and the CO2 fraction. The paper has no curve for 1 day, which is an additional output of the *GEOS* run.

The main features to retrieve are:

- a cold region around the injector whose temperature decreases and whose extent increases with time. After 15 years the cooling is
  about 14 \ |degree|\ C with respect to the initial 45 \ |degree|\ C (minimum temperature around 31 \ |degree|\ C), and the cooled region extends over
  more than 100 m;
- a pressure increase near the well, of about 22 bar after 15 years, which decreases quickly with the distance to the well;
- a cooled and pressurized region that remains inside the region where the gas is almost pure CO2;
- a temperature that remains at 45 \ |degree|\ C farther from the well.

**CO2 front**

The CO2 mass fraction in the gas phase is the quantity that locates the injected CO2: it is equal to 1 behind the front and to 0 ahead of it, where the gas is methane. The front is already at about 9 m from the well after 1 day. It reaches about 25 m after 10 days, 40 m after 30 days, 70 m after 90 days, 150 m after 1 year and 550-600 m after 15 years, which is about half of the domain. The producer, at 1129.5 m, therefore only produces methane during the whole simulation, as assumed in the definition of the well controls. The front radius grows approximately as the square root of time, as expected for a constant-rate injection in a radial geometry (the injected volume is proportional to :math:`r^2`), and slightly more slowly at late times because the CO2 is compressed when the pressure increases. A volume balance gives the same order of magnitude: 3 kg/s during 15 years is :math:`1.4 \times 10^{9}` kg of CO2, which fills the gas-filled pore volume (:math:`\phi\,(1-S_{lr})\,\pi r^2 h`) up to 500-600 m for CO2 densities between 120 and 150 kg/m\ :sup:`3`, which are typical values at 50-70 bar and 31-45 \ |degree|\ C.

The front is sharp, with a transition from 1 to 0 over a factor of about 1.5 in radius. In the absence of molecular diffusion and dispersion in the *GEOS* model, this width is controlled by the numerical dispersion of the scheme and by the mesh size, and the paper, which includes molecular diffusion, can be expected to show a somewhat wider transition. The CO2 profile is important for the interpretation of the two other profiles, since the viscosity and the density of the gas change across the front.

.. plot:: ../inputFiles/compositionalMultiphaseFlow/benchmarks/jouleThomsonCooling/plot_joule_thomson_profile.py plot_co2_fraction
    :align: center
    :caption: Radial profiles of the CO2 mass fraction in the gas phase at the output times. The sharp front moves from about 9 m (1 day) to 550-600 m (15 years).

**Pressure**

The pressure build-up is largest at the injector and is concentrated in the first tens of meters at early times. At the first cell (:math:`r \approx 0.1` m), the pressure is about 64 bar after 1 day, which is an increase of about 14 bar with respect to the initial pressure of 50 bar, and it rises slowly to about 66 bar after 10 days, 68 bar after 90 days, 69 bar after 1 year and 72 bar after 15 years (an increase of about 22 bar). In the region filled by CO2, the pressure decreases almost linearly with :math:`\log r`, with a drop of about 4-5 bar per decade of radius at 15 years. This is the expected signature of a quasi-steady radial Darcy flow, :math:`\Delta p \propto \frac{q\,\mu}{2\pi k h}\ln r`, and a rough estimate with a reservoir flow rate of :math:`3/140 \approx 0.02` m\ :sup:`3`/s and a CO2 viscosity of about :math:`2\times 10^{-5}` Pa.s gives about 6 bar per decade, which is of the same order. The slope changes slightly across the CO2 front, where the gas mobility is different.

The pressure perturbation travels much faster than the CO2 front, since it propagates by diffusion of the pressure (the compressibility of the gas is the only storage mechanism) while the front is transported with the fluid. After 1 day, the pressure is already perturbed up to about 30-50 m while the front is at 9 m, and after 1 year the whole domain is affected. At the outer boundary, the pressure first decreases slightly below the initial value (about 49.7 bar after 1 year) because of the production, before the pressure wave of the injector arrives. It then increases to about 53 bar after 15 years: the system is closed and the injected mass (3 kg/s) is larger than the produced mass (0.56 kg/s), so that the average pressure of the reservoir increases. The pressure at the well is far below the 300 bar limit of the injector, which confirms that the injection stays at a constant rate.

.. plot:: ../inputFiles/compositionalMultiphaseFlow/benchmarks/jouleThomsonCooling/plot_joule_thomson_profile.py plot_pressure
    :align: center
    :caption: Radial profiles of the pressure at the output times. The pressure is highest at the injector (about 72 bar after 15 years) and tends to the initial value of 50 bar far from the well at early times.

**Temperature and comparison with the paper**

The temperature profiles show the three zones that are expected for this problem. Far from the injector, the temperature is unchanged (45 \ |degree|\ C). Behind the cold front, the gas is cooled by the JT expansion between the well and the rest of the reservoir, and the temperature goes through a minimum. Near the well, the temperature increases again toward the injection temperature (about 43 \ |degree|\ C in the first cell), because the injected gas has not yet expanded much when it enters the first cells. The minimum is deeper and further from the well as time goes on: it is about 38 \ |degree|\ C at 0.7 m after 1 day, 36 \ |degree|\ C at about 2 m after 10 days, 34-35 \ |degree|\ C at 5-8 m after 90 days, 33 \ |degree|\ C at about 12 m after 1 year, and 31 \ |degree|\ C at 30-40 m after 15 years, which is a cooling of about 14 \ |degree|\ C. These values are read from the figure and are accurate to about 0.5 \ |degree|\ C.

The cold region is much smaller than the region invaded by the CO2: after 15 years the temperature is back to about 45 \ |degree|\ C beyond 100-150 m, while the CO2 front is at 550-600 m. The thermal front is retarded with respect to the gas front by the heat capacity of the rock and of the immobile water, which have to be cooled before the cooling can move forward. A simple estimate is the ratio between the heat capacity of the gas (:math:`\phi\,S_g\,\rho_g c_{p,g} \approx 4\times 10^{4}` J/m\ :sup:`3`/K, for typical values of the CO2 density and heat capacity) and the total heat capacity of the rock, water and gas (about :math:`2\times 10^{6}` J/m\ :sup:`3`/K): the thermal front moves approximately 50 times more slowly in volume, which is a ratio of :math:`\sqrt{1/50} \approx 0.15` in radius, so that a CO2 front at 550-600 m corresponds to a thermal front at 80-100 m, in agreement with the figure. This is the transient effect mentioned above: the matrix attenuates the cooling, but the thermal front keeps propagating, and the minimum temperature keeps decreasing because the pressure drop that drives the JT cooling increases with time. The weak undershoot below 45 \ |degree|\ C (a few tenths of a degree) at 15 years between the thermal front and the CO2 front is the cooling due to the small pressure drop in this region.

The magnitude of the minimum temperature is larger than a simple estimate based on the JT coefficient at the initial conditions. The pressure difference between the well and the location of the minimum at 15 years is about 10 bar (72 bar at the well, 61-62 bar at 30-40 m), so that a cooling of 14 \ |degree|\ C requires an apparent JT coefficient of at least 1.4 K/bar, larger than the 0.9 K/bar of CO2 at 45 \ |degree|\ C and 50 bar. This is consistent with the strong increase of :math:`\mu_{JT}` when CO2 approaches its critical point (31 \ |degree|\ C and 73.8 bar), which is the case in the cold region, and it means that the cooling is sensitive to the thermodynamic model.

.. plot:: ../inputFiles/compositionalMultiphaseFlow/benchmarks/jouleThomsonCooling/plot_joule_thomson_profile.py plot_temperature
    :align: center
    :caption: Radial profiles of the temperature at the output times (lines: *GEOS*; dots: results of Figure 7 of Oldenburg (2007), digitized). There is no reference curve for 1 day. The minimum temperature after 15 years (about 31 \ |degree|\ C at 30-40 m) is reproduced, while *GEOS* is warmer than the paper close to the well (r < 10 m).

The comparison with the points taken from the paper (dots) shows a good overall agreement, with some differences that are worth noting:

- **Good agreement on the cold front and on the minimum temperature.** The position and the steepness of the thermal front are close to the digitized points at all times (for instance the rise from 35 to 45 \ |degree|\ C between about 10 and 30 m after 1 year, and from 35 to 40 \ |degree|\ C between about 60 and 100 m after 15 years), and the minimum temperature after 15 years, about 31 \ |degree|\ C at 30-40 m, is the same as in the paper. This minimum is the quantity that matters for the risk of hydrate or ice formation, and it remains far above the temperature of hydrate stability. The thermal front and the CO2 front also have the same relative position as in the paper (the cooled zone stays inside the CO2-rich region).
- **Warmer temperatures in the near-well region.** For :math:`r` smaller than about 10 m, the digitized points are lower than the *GEOS* curves. At early times (10 to 90 days), the paper gives 34-35 \ |degree|\ C at 1.5-3 m, about 2 \ |degree|\ C below the *GEOS* curves, and after 15 years the paper gives 35.5 \ |degree|\ C at 3 m and 33.8 \ |degree|\ C at 7 m, while *GEOS* is about 3-4 \ |degree|\ C warmer. The paper also shows a minimum that is more flat and more extended (about 31-34 \ |degree|\ C between 7 and 40 m at 15 years), whereas the *GEOS* profile decreases almost linearly with :math:`\log r` up to about 10 m and then drops more rapidly toward the minimum. *GEOS* therefore predicts a weaker cooling than the paper in the first meters around the injector.
- **Small-scale features at early times.** The 1 day and 10 day curves show small kinks at 0.2-0.8 m, at the scale of the first cells (a few centimeters to decimeters). They have no counterpart in the paper, which has no data below 1.5 m, and they are related to the resolution of the mesh and of the time step at the beginning of the injection.
- **No reference at 1 day.** The dark blue curve is only given by *GEOS*.

The origin of these differences has not been isolated, and several effects are likely to contribute:

- the thermodynamic model. The JT effect is very sensitive to the equation of state in the near-critical region reached by the CO2 in the cold zone, and the enthalpy depends on the ideal-gas heat capacities and on the departure functions. *GEOS* uses Peng-Robinson with polynomial ideal-gas heat capacities, whereas the thermodynamic properties of EOS7C may differ in these respects, and a few percent of difference on :math:`\mu_{JT}` or on :math:`c_p` close to the critical point leads to a difference of several degrees of cooling;
- the representation of the injector. In the paper, the CO2 is injected as a source term at 45 \ |degree|\ C in the first grid block, whereas *GEOS* injects it through a thermal well, which has its own pressure and temperature state and is connected to the reservoir by a perforation with a prescribed transmissibility. This changes the pressure drop and the enthalpy balance in the first cells, which can modify the temperature a few meters away, where the heat is exchanged with the rock;
- the physical processes that are not in the *GEOS* model: the capillary pressure and the molecular diffusion of the paper, and the slightly different gas relative permeability (Brooks-Corey instead of Corey);
- the accuracy of the reference: the dots are read from a figure with a logarithmic radial axis, so that the uncertainty is about 0.5 \ |degree|\ C on the temperature and 10-15 % on the radius, which is a significant part of the differences at large radius but cannot explain the 3-4 \ |degree|\ C close to the well.

.. note::
    The cooling can be compared with the upper bound :math:`\mu_{JT}\,\Delta P`. The relevant pressure difference is the one between the well and the point of interest, and not the total pressure increase at the well, since the gas is injected at the reservoir temperature and cools as it expands toward lower pressures. As discussed above, the cooling observed in the simulation corresponds to an apparent JT coefficient larger than the one at the initial conditions, because of the proximity of the critical point, even though the heat released by the matrix, which is important in the first days and progressively exhausted, attenuates the cooling.

------------------------------------------------------------------
To go further
------------------------------------------------------------------

The other cases of the paper can be obtained from this input by changing a few parameters:

- the permeability (``rockPerm``): :math:`10^{-13}` and :math:`10^{-14}` m\ :sup:`2` for Figures 5 and 6, and the summary of Figure 8;
  the injection rate is kept unchanged, so the pressure drop and the JT cooling vary with the permeability;
- the porosity (``rockPorosity``), to study the compensation between the permeability and the heat supplied by the matrix (Figure 9);
- the injection temperature and the initial reservoir temperature (``injectionTemperature``, ``initialTemperature``) for the
  cases of Figures 10 and 11 (20 \ |degree|\ C injection in a 75 \ |degree|\ C reservoir). The initial compositions must then be computed again, as explained above.

**Feedback on this example**

For any feedback on this example, please submit a `GitHub issue on the project's GitHub page <https://github.com/GEOS-DEV/GEOS/issues>`_.

------------------------------------------------------------------
References
------------------------------------------------------------------

- Oldenburg, C. M. 2007. Joule-Thomson cooling due to CO2 injection into natural gas reservoirs. *Energy Conversion and Management*, 48(6), 1808-1815. https://doi.org/10.1016/j.enconman.2007.01.010
- Lemmon, E. W., McLinden, M. O., and Friend, D. G. 1998. Thermophysical Properties of Fluid Systems. In *NIST Chemistry WebBook*, NIST Standard Reference Database, eds. Linstrom, P. J. and Mallard, W. G., Vol. 69. Gaithersburg, Maryland, USA: National Institute of Standards and Technology. https://webbook.nist.gov/chemistry/fluid/
