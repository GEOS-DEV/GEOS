.. _jouleThomsonCooling:

##########################################################################
Joule-Thomson Cooling during CO2 Injection into a Gas Reservoir
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

where :math:`c_p` and :math:`v` are the specific heat capacity and the specific volume of the gas. For CO2 at 45 °C (318.15 K) and 50 bar, :math:`\mu_{JT}` is approximately 0.9 K/bar according to the NIST reference equation of state, and it decreases when the pressure increases (approximately 0.4 K/bar at 100 bar). The product :math:`\mu_{JT}\,\Delta P` is thus the maximum cooling that can be obtained for a given pressure drop. In a reservoir, this maximum is never reached: the rock grains and the residual water, with a volumetric heat capacity much larger than that of the gas, release heat as the gas cools down. This is a transient effect: the matrix progressively loses its heat during a prolonged expansion. The cooling therefore depends not only on the pressure drop, but also on the porosity, the residual water saturation, the heat capacity of the grains and the thermal conductivity.

**Reservoir model of the paper**

The paper considers an idealized depleted gas reservoir, with the geometry and the properties given in the paper. It is represented by a one-dimensional radially symmetric system, 50 m thick and 1130 m in radius, as illustrated below.

.. code-block:: console

    injection well (r = 0)  --->  CO2-rich gas  | CH4-rich gas  --->  production (r ~ 1100 m)

The reservoir is initially filled with methane at 50 bar and 45 °C. The water is immobile (the water saturation is equal to the residual saturation :math:`S_{lr} = 0.2`), and the only mobile phase is the gas. Pure CO2 is injected at a constant rate of 3 kg/s and at the reservoir temperature. Methane is produced at the outer boundary, at 0.56 kg/s, to represent the production of four surrounding wells located at a distance of 1100 m in a five-spot pattern. The system is closed: no fluid or heat crosses the outer boundary. Because the injection temperature is equal to the initial reservoir temperature, every temperature change in the reservoir is caused by the JT effect (attenuated by the heat provided by the rock and the water), which makes this setup well suited to isolate and validate this effect.

The parameters of the problem are summarized in the following table.

.. list-table:: Parameters of the problem
    :widths: 30 70
    :header-rows: 1
    :class: align-top

    * - Parameter
      - Value
    * - Geometry
      - 1D radial, :math:`R = 1130` m, thickness 50 m

        GEOS is modeled as a 1° sector, :math:`r \in [0.1, 1130]` m.
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
      - 45 °C (318.15 K)
    * - Initial fluids
      - CH4 gas, water at :math:`S_l = S_{lr} = 0.2`

        GEOS initializes this via mass fractions: :math:`z_{CO_2} = 0`, :math:`z_{CH_4} = 0.115401`, :math:`z_{H_2O} = 0.884599`.
    * - CO2 injection
      - 3 kg/s of CO2 at 45 °C
    * - CH4 production
      - 0.56 kg/s at the outer boundary (1100 m)
    * - Simulated time
      - 15 years
    * - Equation of state
      - Peng-Robinson for the gas phase

        The paper uses EOS7C/GasEOS, while GEOS explicitly uses Soreide-Whitson for the aqueous phase.

..
  **What the paper finds for this case**

..
  For a constant injection rate, the pressure build-up near the well scales inversely with the permeability. The paper reports a cooling of about 2 °C for 100 mD, 8 °C for 10 mD, and about 14 °C for 5 mD (Figure 8) after 15 years of injection. For 5 mD (Figure 7), the cooled region extends over more than 100 m after 15 years and it remains inside the region occupied by nearly pure CO2. The minimum temperature remains around 31 °C, well above the temperatures at which CO2 hydrates are stable at these pressures (below approximately 10 °C), so that the problem is not expected to be a concern for CSEGR in this setting.

------------------------------------------------------------------------
Input base
------------------------------------------------------------------------

This benchmark test is based on the XML file located below:

.. code-block:: console

  inputFiles/compositionalMultiphaseFlow/benchmarks/jouleThomsonCooling/jouleThomsonCooling.xml

The mesh is a thin wedge of cylinder (``InternalWellbore``), 1° wide, which is refined near the injection well by the radial bias ``rBias``. The 360° of the paper are accounted for by dividing all the well rates by 360.

.. literalinclude:: ../../../../../../../inputFiles/compositionalMultiphaseFlow/benchmarks/jouleThomsonCooling/jouleThomsonCooling.xml
    :language: xml
    :start-after: <!-- SPHINX_MESH -->
    :end-before: <!-- SPHINX_MESH_END -->

The mesh is only one cell thick in the angular and vertical directions, so the flow is purely radial. The two wells are vertical and are located in the mid-plane of the wedge (:math:`\theta = 0.5°`): the injector is at :math:`r = 0.125` m, in the first cells of the mesh, and the producer is at :math:`r \approx 1097` m, close to the outer boundary at 1130 m. In the paper, the production is taken from the outermost grid block. Since there is a single cell over the thickness, the position of the perforation along the well (``distanceFromHead``) has no influence on the solution. The two wells are listed in the ``ElementRegions`` together with the reservoir region. All the regions share the same fluid, relative permeability and thermal conductivity models:

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

- ``equationsOfState`` lists one equation of state for each phase, in the order of ``phaseNames``. The gas phase uses the   Peng-Robinson (PR) cubic equation of state, as in the paper (EOS7C computes the properties of the CO2-CH4-H2O gas mixtures with PR). The aqueous phase uses the Soreide-Whitson variant of PR, which modifies the attractive term for water and uses salinity-dependent interaction parameters, in order to represent the mutual solubilities of the gases and of water. The aqueous phase stays at its residual saturation in this problem and acts mostly as a heat reservoir, so that the solubility effects are weak, but the model remains fully compositional.
- the component data (``componentMolarWeight`` in kg/mol, critical pressures and temperatures, acentric factors) are the standard values of CO2, CH4 and H2O, given in this order, which is the order of ``componentNames``. The same order is used for the ``injectionStream``, the ``globalCompFraction`` components and the binary and heat capacity coefficients. 
- the aqueous phase density and viscosity follow the brine correlations of Phillips et al., with ``salinity`` equal to 0 (fresh water).
- ``pressureCoordinates`` and ``temperatureCoordinates`` define the pressure-temperature grid on which the fluid properties are evaluated. They must bracket all the states visited during the simulation: here from 1 bar (the minimum bottom-hole pressure of the producer) to 300 bar (the maximum bottom-hole pressure of the injector), and from -10 °C to 110 °C. The states actually visited are between about 5 and 7 MPa, and between 31 and 45 °C.

.. _jouleThomsonCoolingHeatCapacity:

**Heat capacity coefficients**

The enthalpy of a phase is the sum of an ideal-gas contribution and of a departure (residual) contribution given by the equation of state, which accounts for the non-ideal behavior of the real gas. The ideal-gas part is obtained by integrating the isobaric ideal-gas heat capacity of each component, which is provided by ``componentHeatCapacityCoefficients`` (one row per component, in the order of ``componentNames``):

.. literalinclude:: ../../../../../../../inputFiles/compositionalMultiphaseFlow/benchmarks/jouleThomsonCooling/jouleThomsonCooling.xml
    :language: xml
    :start-after: <!-- SPHINX_FLUID -->
    :end-before: <!-- SPHINX_FLUID_END -->
    :lines: 11-15
    :dedent: 6

Each row contains the five coefficients :math:`a_0, \ldots, a_4` of a polynomial in the temperature difference
with respect to 25 °C:

.. math::

    c_p^{0}(T) = \sum_{i=0}^{4} a_i \, (T - 298.15)^i \qquad [\mathrm{J/mol/K}], \quad T \text{ in K}.

The first coefficient :math:`a_0` is therefore the ideal-gas molar heat capacity at 25 °C (37.14 J/mol/K for CO2, 35.71 for CH4 and 33.59 for H2O). The following table gives the resulting values at the reservoir temperature, which are in agreement with the reference ideal-gas heat capacities of the NIST equations of state within 0.01 J/mol/K over the whole temperature range of the fluid tables (-10 to 110 °C).

.. list-table:: Ideal-gas heat capacity :math:`c_p^0` at 318.15 K
    :widths: 20 20 20 20
    :header-rows: 1

    * - Component
      - CO2
      - CH4
      - H2O
    * - :math:`c_p^0` (J/mol/K)
      - 38.04
      - 36.50
      - 33.69

These coefficients are important for two reasons. First, the JT coefficient is inversely proportional to the heat capacity (see the formula above), and the heat capacity of the gas, which is the sum of the ideal-gas part and of the real-gas departure, controls how much the gas cools for a given enthalpy variation. Second, the heat capacity sets the thermal inertia of the fluids against that of the rock, which determines how fast the cooling is attenuated by the matrix.

.. _jouleThomsonBinaryCoeff:

**Binary interaction parameters**

The Peng-Robinson equation of state uses a mixing rule for the attractive term of the mixture in which the cross term of each pair of components :math:`(i,j)` is scaled by :math:`(1 - k_{ij})`. The binary interaction parameters :math:`k_{ij}` are given by ``componentBinaryCoeff`` as a symmetric matrix with a zero diagonal, in the order CO2, CH4, H2O:

.. literalinclude:: ../../../../../../../inputFiles/compositionalMultiphaseFlow/benchmarks/jouleThomsonCooling/jouleThomsonCooling.xml
    :language: xml
    :start-after: <!-- SPHINX_FLUID -->
    :end-before: <!-- SPHINX_FLUID_END -->
    :lines: 10
    :dedent: 6

so that :math:`k_{CO_2-CH_4} = 0.1` and :math:`k_{CO_2-H_2O} = k_{CH_4-H_2O} = 0.03`. The CO2-CH4 coefficient is the one that matters most for this problem, since the gas phase is a CO2-CH4 mixture in the region where the CO2 front progresses and displaces the methane. The coefficients involving water are small, which is consistent with the fact that the gas phase contains only a very small amount of water vapor at 45 °C; the interactions of water with the dissolved gases in the aqueous phase are handled by the specific correlations of the Soreide-Whitson formulation.

**Rock, relative permeability and thermal conductivity**

The rock is modeled as a slightly compressible porous medium with a constant permeability. Its heat capacity is given as
a volumetric value, :math:`\rho_R C_R = 2600 \times 1000 = 2.6 \times 10^{6}` J/m\ :sup:`3`/K, which is the product
of the grain density and the grain heat capacity of the paper. The matrix contributes :math:`(1-\phi)\rho_R C_R` per unit of bulk volume to the
energy accumulation, which is the term that attenuates the JT cooling in the paper.

.. literalinclude:: ../../../../../../../inputFiles/compositionalMultiphaseFlow/benchmarks/jouleThomsonCooling/jouleThomsonCooling.xml
    :language: xml
    :start-after: <!-- SPHINX_ROCK -->
    :end-before: <!-- SPHINX_ROCK_END -->

The relative permeabilities are Brooks-Corey power laws with the residual saturations of the paper. With a water residual saturation of 0.2
and an initial water saturation equal to 0.2, the water relative permeability is zero and the water does not move:

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
    - the aqueous phase uses the Soreide-Whitson equation of state, whereas EOS7C uses PR for the whole system;
    - the producer is a well with a prescribed surface gas rate located at 1097 m, instead of a source term in the outermost grid block.

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
  ``HydrostaticEquilibrium`` is required, unlike in the SPE11 example) since the model is one cell thick and the problem is
  controlled by the radial flow.
- **Temperature, 45 °C (318.15 K).** This is the reservoir temperature of the paper, and it is the same as the injection
  temperature. The initial thermal state is therefore in equilibrium with the injected fluid, and **any temperature
  deviation observed afterwards is the result of the JT expansion** (and of its attenuation by the heat of the rock),
  not of a thermal contrast between the injected fluid and the reservoir.
- **Composition.** The reservoir contains methane and water, with no CO2. In *GEOS*, the saturations are not primary inputs: the initial
  phase distribution is computed by a flash from the global composition, pressure and temperature. The paper
  prescribes an initial water saturation equal to the residual saturation, :math:`S_w = 0.2`, so that the water is immobile and the only mobile
  phase is gas. The global mass fractions are calibrated to obtain exactly this saturation. At 50 bar and 45 °C, the density of CH4 gas
  is :math:`\rho_g = 32.36` kg/m\ :sup:`3` and the density of water is :math:`\rho_w = 992.3` kg/m\ :sup:`3`, so that

  .. math::

      z_{CH_4} = \frac{\rho_g (1 - S_w)}{\rho_g (1 - S_w) + \rho_w S_w} = \frac{32.36 \times 0.8}{32.36 \times 0.8 + 992.3 \times 0.2} \approx 0.1154,
      \qquad z_{H_2O} = 1 - z_{CH_4} \approx 0.8846,

  which are the values of ``initZ1`` and ``initZ2`` (``initZ0`` is the CO2 mass fraction and is equal to zero). The mutual
  solubility between the phases is neglected in this estimate, and its effect on the saturation is very small.

.. note::
    These mass fractions are specific to the initial pressure and temperature. If these values are modified (for example
    to reproduce the cases of the paper with a reservoir at 75 °C), the densities change and the compositions must be
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
      - Constant-rate injection of 3 kg/s. The rate is divided by 360 because the model is a 1° sector
        (:math:`3 \times 1/360`).
    * - Injected fluid
      - ``injectionStream`` = {0.9998, 1e-4, 1e-4}, ``injectionTemperature`` = 318.15 K
      - Pure CO2 at :math:`T_{inj}` = 45 °C, equal to the reservoir temperature. The stream is
        99.98 % CO2 (in mass), with traces of CH4 and H2O, which is equivalent to pure CO2 for this problem.
    * - Injector safeguard
      - ``MaximumBHPConstraint``, 300 bar at the well head (``referenceElevation`` = 50 m)
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
      - ``MinimumBHPConstraint``, 1 bar
      - Not in the paper. It is inactive here since the reservoir pressure is above 50 bar.
    * - Crossflow
      - ``enableCrossflow="0"`` for both wells
      - The injector can only inject and the producer can only produce.

The mass rate of the producer in the paper and the surface volumetric rate used here are equivalent as long as the
produced gas is methane. A simple volumetric balance indicates that, after 15 years, the CO2 front remains at a few hundred metres
from the injector, so the producer, at 1097 m, only produces methane during the simulation, and the rate is
constant in mass. The ``writeCSV`` flag exports the time history of the rates, of the bottom-hole pressure and of the
temperature of each well, which is useful to monitor the injection overpressure.

------------------------------------------------------------------------
Events and output
------------------------------------------------------------------------

The simulation covers 15 years. The outputs are requested at the times of the curves shown in Figure 7: the initial state, then 10, 20,
30 and 90 days, 1 year and 15 years (a year is 3.15576e7 s). The output events are listed before the solver event so that
they record the state at their target time. The solver is applied with a maximum time step of 30 days and the time step
is controlled by the targets defined in the flow solver. Statistics are computed once per year, and restart files are written
every 5 years.

.. literalinclude:: ../../../../../../../inputFiles/compositionalMultiphaseFlow/benchmarks/jouleThomsonCooling/jouleThomsonCooling.xml
    :language: xml
    :start-after: <!-- SPHINX_EVENTS -->
    :end-before: <!-- SPHINX_EVENTS_END -->

.. note::
    The time tags in *GEOS* are in seconds. The ``maxTime`` of the events is slightly larger than 15 years, so that the
    last output event is reached.

------------------------------------------------------------------------
Post-processing and reporting
------------------------------------------------------------------------

The results are inspected with Paraview as in the usual *GEOS* workflow. To compare with Figure 7, the pressure, the temperature and the CO2 mass
fraction in the gas phase can be plotted along the radial direction (for instance with *Plot Over Line* along the first
cell row of the wedge) for the output times listed above. Note that the plots of the paper show the first 100 m of the system.

The main features to retrieve are:

- a cold region around the injector whose temperature decreases and whose extent increases with time. After 15 years the cooling is
  about 14 °C with respect to the initial 45 °C (minimum temperature around 31 °C), and the cooled region extends over
  more than 100 m;
- a pressure increase near the well, of about 20 bar after 15 years, which decreases quickly with the distance to the well;
- a cooled and pressurized region that remains inside the region where the gas is almost pure CO2;
- a temperature that remains at 45 °C farther from the well.

.. note::
    The cooling can be compared with the upper bound :math:`\mu_{JT}\,\Delta P`. With :math:`\mu_{JT}` of the order of 0.8 K/bar
    and a pressure increase of about 20 bar, the bound is of the order of 15 °C, close to the observed cooling.
    The cooling observed in the simulation is of the same order since the heat released by the matrix is important
    in the first days but is progressively exhausted, while the pressure keeps increasing.

.. .. _jouleThomson_T:
.. .. figure:: ./pictures/jouleThomsonCooling_T.png
..     :align: center
..     :width: 500
..     :figclass: align-center
..
..     Temperature, pressure and CO2 mass fraction profiles at the output times.

------------------------------------------------------------------
To go further
------------------------------------------------------------------

The other cases of the paper can be obtained from this input by changing a few parameters:

- the permeability (``rockPerm``): :math:`10^{-13}` and :math:`10^{-14}` m\ :sup:`2` for Figures 5 and 6, and the summary of Figure 8;
  the injection rate is kept unchanged, so the pressure drop and the JT cooling vary with the permeability;
- the porosity (``rockPorosity``), to study the compensation between the permeability and the heat supplied by the matrix (Figure 9);
- the injection temperature and the initial reservoir temperature (``injectionTemperature``, ``initT``) for the
  cases of Figures 10 and 11 (20 °C injection in a 75 °C reservoir). The initial compositions must then be computed again, as explained above.

**Feedback on this example**

For any feedback on this example, please submit a `GitHub issue on the project's GitHub page <https://github.com/GEOS-DEV/GEOS/issues>`_.
