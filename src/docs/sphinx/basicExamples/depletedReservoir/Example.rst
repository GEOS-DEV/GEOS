.. |oslash| unicode:: U+00F8

.. |degree| unicode:: U+00B0

.. _CO2InjectionDepletedReservoirExample:

##########################################################################
CO2 Injection into a Depleted Reservoir
##########################################################################


**Context**

Depleted hydrocarbon fields are among the most attractive targets for geological carbon storage. A field that has trapped gas over geological time has proven its seal, its storage volume and its injectivity are constrained by decades of production data, and the pressure drop caused by production leaves room for the injected fluid. The same history makes the problem demanding to simulate. The reservoir is under-pressured, it still contains hydrocarbons that mix with the injected CO2, and the sequence of field operations (production, water injection, shut-in, CO2 injection) imposes successive pressure and saturation transients on a heterogeneous structure.

This example sets up such a case on a dome-shaped reservoir that contains the residual hydrocarbons of a previous production. It illustrates how to combine in a single *GEOS* input

- a compositional model with four components (CH4, CO2, H2S and H2O) distributed in a gas phase and an aqueous phase, available with two alternative descriptions of the phase equilibrium: pre-tabulated equilibrium constants (K-values), or an equation of state (Peng-Robinson for the gas and S\ |oslash|\ reide-Whitson for the brine);
- an initialization by capillary-gravity equilibrium around a prescribed gas-water contact in a heterogeneous structure;
- three rock types, the top one being an inactive seal that is part of the mesh but not of the flow problem;
- a multi-stage well schedule driven by time tables: production, water injection, shut-in and injection of a CO2-CH4 mixture, with rate targets and pressure limits on every well.

The post-processing is part of the example: the VTK files written by *GEOS* are converted into one HDF5 file, and all the pictures of this page, including the map of the dome, are generated from this file when the documentation is built.

------------------------------------------------------------------------
Brief case description
------------------------------------------------------------------------

**Dome model**

The reservoir is a dome. A sand body of good quality forms the base and the core of the structure, a clay layer of poorer quality covers it, and a caprock closes the top. The clay is thick on the flanks and thins out over the crest, where the sand comes close to the caprock. The figure below shows the rock types in the vertical section that contains the wells. It is not a stored image: it is rebuilt from the cell regions of the simulation output.

.. .. plot:: plot_dome.py plot_rock_types
    :align: center
    :caption: Rock types of the dome, in the vertical x-z section that contains the wells. The dashed line is the gas-water contact, the triangles are the wells (injectors pointing up, producer pointing down). The caprock is part of the mesh but is inactive in the flow problem.

**Operations**

The model starts from a depleted state, with the gas-water contact at -4400 m. The field is then operated in four successive stages, summarized in the table below and in the first figure of the Results section.

.. list-table:: Operating schedule
    :widths: 12 18 20 20 30
    :header-rows: 1
    :class: align-top

    * - Time (yr)
      - Producer
      - Water injector
      - Gas injector
      - Stage
    * - 0 - 5
      - 15 m\ :sup:`3`/s
      - closed
      - closed
      - **Depletion**: fluid is produced from the crest region
    * - 5 - 8
      - 15 m\ :sup:`3`/s
      - 0.2 m\ :sup:`3`/s
      - closed
      - **Water flooding**: water is injected from the flank while production continues
    * - 8 - 13
      - closed
      - closed
      - closed
      - **Shut-in**: pressure and fluids re-equilibrate
    * - 13 - 18
      - closed
      - closed
      - 20 m\ :sup:`3`/s
      - **Gas injection**: 85 % CO2 and 15 % CH4 (mole fractions)
    * - 18 - 20
      - closed
      - closed
      - closed
      - **Post-injection**: migration of the plume

All the rates are volumetric rates at surface conditions (101325 Pa and 288.71 K).

The parameters of the problem are summarized in the following table.

.. list-table:: Parameters of the problem
    :widths: 30 70
    :header-rows: 1
    :class: align-top

    * - Parameter
      - Value
    * - Temperature
      - 344.15 K (71 \ |degree|\ C), uniform and constant (isothermal)
    * - Brine salinity
      - 2.3 (molality of NaCl)
    * - Components
      - CH4, CO2, H2S, H2O
    * - Phases
      - gas and water (aqueous)
    * - Rock regions
      - SAND (base and core), CLAY (middle layer), CAPROCK (inactive seal)
    * - Porosity
      - SAND 0.14375, CLAY 0.125 at the reference pressure of 30 MPa

        Rock compressibility :math:`3.4\times 10^{-10}` Pa\ :sup:`-1`, water compressibility :math:`4.1483\times 10^{-10}` Pa\ :sup:`-1`
    * - Permeability
      - SAND: :math:`k_x = k_y = 3.12\times 10^{-13}` m\ :sup:`2` (316 mD), :math:`k_z = 3.12\times 10^{-14}` m\ :sup:`2` (31.6 mD)

        CLAY: :math:`k_x = k_y = 3.12\times 10^{-14}` m\ :sup:`2` (31.6 mD), :math:`k_z = 3.12\times 10^{-15}` m\ :sup:`2` (3.16 mD)
    * - Saturation end points
      - Irreducible water saturation 0.20 (SAND) and 0.28 (CLAY)

        Critical gas saturation 0.08 (SAND) and 0.12 (CLAY)
    * - Initial state
      - Capillary-gravity equilibrium, 27 MPa at -4453.75 m, gas-water contact at -4400 m
    * - Wells
      - Water injector and gas injector at (x, y) = (125, 25) m, perforated at z = -4505 m

        Producer at (x, y) = (2056.25, 25) m, perforated at z = -4126.25 m
    * - Simulated time
      - 20 years

------------------------------------------------------------------------
Input files
------------------------------------------------------------------------

The case is split in several XML files that are assembled with the ``Included`` mechanism, so that the physical model, the mesh and the fluid description can be exchanged independently. They are located in

.. code-block:: console

  inputFiles/compositionalMultiphaseWell

.. list-table:: Input files
    :widths: 38 62
    :header-rows: 1
    :class: align-top

    * - File
      - Content
    * - ``dome_soreide_whitson_base.xml``
      - Everything except the mesh, with the equation-of-state fluid model (Peng-Robinson and S\ |oslash|\ reide-Whitson)
    * - ``dome_kvalue_base.xml``
      - Same as above, with the fluid model based on tabulated K-values (the tables are read from the ``tables`` directory)
    * - ``dome_mesh_bench.xml``, ``dome_mesh_smoke.xml``
      - Mesh file and well trajectories, for the benchmark resolution and for a quick functional check
    * - ``dome_properties.xml``
      - Porosity and permeability of every rock type
    * - ``dome_soreide_whitson_bench.xml``, ``dome_soreide_whitson_smoke.xml``, ``dome_kvalue_bench.xml``, ``dome_kvalue_smoke.xml``
      - Drivers that include one base file and one mesh file. Each one is a complete input

The results shown below are obtained with ``dome_soreide_whitson_bench.xml``. The run is started in the usual way

.. code-block:: console

  geos -i dome_soreide_whitson_bench.xml

The sections below follow the order of the input and quote the files with the equation-of-state fluid, except where the K-value variant differs. The marked parts of the XML files are included directly in the text.

------------------------------------------------------------------------
Mesh and wells
------------------------------------------------------------------------

The mesh is read from a VTK image file (``dome_bench.vti``) with ``VTKMesh``. The rock type is a cell attribute of this file, with the values 1, 2 and 3, and *GEOS* splits the mesh accordingly into the cell blocks ``1_hexahedra`` (sand), ``2_hexahedra`` (clay) and ``3_hexahedra`` (caprock). The three wells are declared in the mesh with ``InternalWell``: each one is a vertical polyline with a single perforation, and is linked by ``wellControlsName`` to the well definition of the solver section. The water injector and the gas injector are two separate wells that share the same trajectory, because they operate at different times.

.. literalinclude:: ../../../../../inputFiles/compositionalMultiphaseWell/dome_mesh_bench.xml
    :language: xml
    :start-after: <!-- SPHINX_DOME_MESH_BEGIN -->
    :end-before: <!-- SPHINX_DOME_MESH_END -->
    :dedent: 8
    :emphasize-lines: 3,10,16,23,29,36,42

The injectors are at the bottom of the left flank, in the sand and below the gas-water contact. The producer is on the right flank, above the contact. With this arrangement the injected fluid has to travel up the structure to reach the gas that remains from the depletion.

------------------------------------------------------------------------
Regions and rock properties
------------------------------------------------------------------------

**Element regions**

Each cell block is assigned to a region, with its own list of constitutive models. The caprock is only given the rock model: it has no fluid, relative permeability or capillary pressure, and is therefore excluded from the flow problem (it is also absent from the ``targetRegions`` of the solvers). It acts as a perfect seal and is only there to close the geometry. The wells have their own regions, which carry the fluid and relative permeability models used in the wellbore.

.. literalinclude:: ../../../../../inputFiles/compositionalMultiphaseWell/dome_soreide_whitson_base.xml
    :language: xml
    :start-after: <!-- SPHINX_DOME_ELEMENT_REGIONS_BEGIN -->
    :end-before: <!-- SPHINX_DOME_ELEMENT_REGIONS_END -->
    :dedent: 2
    :emphasize-lines: 5,10,15,19,23,27

**Rock**

The rock is a compressible porous medium (``PressurePorosity``) with a constant permeability. The values given in the model are default values: the porosity and the permeability of every region are overwritten cell by cell by the ``FieldSpecification`` tags of ``dome_properties.xml``, which scale the reference porosity and the three permeability components of each rock type:

.. literalinclude:: ../../../../../inputFiles/compositionalMultiphaseWell/dome_soreide_whitson_base.xml
    :language: xml
    :start-after: <!-- SPHINX_DOME_ROCK_BEGIN -->
    :end-before: <!-- SPHINX_DOME_ROCK_END -->
    :dedent: 4
    :emphasize-lines: 7,8

.. literalinclude:: ../../../../../inputFiles/compositionalMultiphaseWell/dome_properties.xml
    :language: xml
    :start-after: <!-- SPHINX_DOME_PROPERTIES_BEGIN -->
    :end-before: <!-- SPHINX_DOME_PROPERTIES_END -->
    :dedent: 8
    :emphasize-lines: 2,9,16,40

The sand is about ten times more permeable than the clay, and in both rocks the vertical permeability is ten times smaller than the horizontal one, which makes the layering a strong control of the vertical migration of the gas.

------------------------------------------------------------------------
Fluid models
------------------------------------------------------------------------

The fluid is a mixture of four components, whose order is the same in all the tables of the input: CH4, CO2, H2S and H2O. It is distributed in two phases, a ``gas`` phase and a ``water`` phase, by a two-phase flash. The brine density and viscosity follow the correlations of Phillips et al., with a salinity of 2.3. The pure-component data (molar mass, critical properties and acentric factor) and the binary interaction coefficients, non-zero only for the pairs that involve water, are common to the two models.

Two fluid descriptions are provided. They differ only by the way the phase equilibrium is computed, and the rest of the input is identical.

.. list-table:: The two fluid modeling approaches
    :widths: 22 39 39
    :header-rows: 1
    :class: align-top

    * -
      - Equation of state
      - K-values
    * - Input file
      - ``dome_soreide_whitson_*.xml``
      - ``dome_kvalue_*.xml``
    * - *GEOS* model
      - ``CompositionalTwoPhaseFluidPhillipsBrine``
      - ``CompositionalTwoPhaseKValueFluidPhillipsBrine``
    * - Phase equilibrium
      - Stability analysis and flash with Peng-Robinson (gas) and S\ |oslash|\ reide-Whitson (water): the equilibrium depends on the composition
      - Flash with K-values :math:`K_i = y_i / x_i` read from tables of pressure and temperature: the equilibrium does not depend on the composition
    * - Specific input
      - ``pressureCoordinates``, ``temperatureCoordinates`` and the tolerances of the stability and flash algorithms
      - ``kValueTables``, one table per component

**Equation-of-state model**

The ``equationsOfState`` attribute lists one equation of state per phase, in the order of ``phaseNames``. The gas phase uses the Peng-Robinson (PR) cubic equation of state. The aqueous phase uses the S\ |oslash|\ reide-Whitson variant of PR, which modifies the attractive term for water and uses salinity-dependent interaction parameters, so that the solubility of the gases in the brine and the vaporization of water into the gas are represented. This matters here, since the CO2 dissolves into the brine around the plume. The tolerances ``stabilityTolerance``, ``flashTolerance`` and ``stabilityThreshold`` control the phase stability test and the flash.

.. literalinclude:: ../../../../../inputFiles/compositionalMultiphaseWell/dome_soreide_whitson_base.xml
    :language: xml
    :start-after: <!-- SPHINX_DOME_FLUID_BEGIN -->
    :end-before: <!-- SPHINX_DOME_FLUID_END -->
    :dedent: 4
    :emphasize-lines: 4,5,15,16,18,19,20,21

The fluid properties are evaluated on a pressure-temperature grid, defined by ``pressureCoordinates`` (13 values, from 1 bar to 600 bar) and ``temperatureCoordinates`` (5 values, from 283.15 K to 354.15 K). The grid must bracket all the states visited during the simulation. Here it covers the reservoir temperature (344.15 K), the surface conditions used to convert the well rates (101325 Pa and 288.71 K) and the whole range of the well pressure limits (1 to 50 MPa).

**K-value model**

In the K-value variant, the fluid block is replaced by ``CompositionalTwoPhaseKValueFluidPhillipsBrine``. The component data and the equation-of-state names are kept, but the pressure and temperature grid and the tolerances are replaced by ``kValueTables``, which gives, for each component, the name of a table that returns the equilibrium constant as a function of pressure and temperature:

.. literalinclude:: ../../../../../inputFiles/compositionalMultiphaseWell/dome_kvalue_base.xml
    :language: xml
    :start-after: <!-- SPHINX_DOME_FLUID_BEGIN -->
    :end-before: <!-- SPHINX_DOME_FLUID_END -->
    :dedent: 4
    :emphasize-lines: 4,15,17

The tables are ``TableFunction`` objects whose coordinates and values are read from text files:

.. literalinclude:: ../../../../../inputFiles/compositionalMultiphaseWell/dome_kvalue_base.xml
    :language: xml
    :start-after: <!-- SPHINX_DOME_KVALUE_TABLES_BEGIN -->
    :end-before: <!-- SPHINX_DOME_KVALUE_TABLES_END -->
    :dedent: 4
    :emphasize-lines: 5,6,12,13,19,20,26,27

A K-value flash is simpler than a flash based on an equation of state. Its drawback is that the :math:`K_i` do not depend on the composition. It is a good choice when the compositions stay in a range for which the tables are representative, and the comparison of the two models on the same case is one of the purposes of this example.

**Relative permeability and capillary pressure**

The relative permeabilities and the gas-water capillary pressures are tables of the water saturation (or of the gas saturation for the gas relative permeability), with one set per rock type. The clay uses a J-function (``JFunctionCapillaryPressure``), which scales a single dimensionless curve with the permeability and the porosity of each cell and with the surface tension, so that the capillary pressure adapts to the local rock quality.

.. literalinclude:: ../../../../../inputFiles/compositionalMultiphaseWell/dome_soreide_whitson_base.xml
    :language: xml
    :start-after: <!-- SPHINX_DOME_RELPERM_BEGIN -->
    :end-before: <!-- SPHINX_DOME_RELPERM_END -->
    :dedent: 4
    :emphasize-lines: 2,4,7,9

.. literalinclude:: ../../../../../inputFiles/compositionalMultiphaseWell/dome_soreide_whitson_base.xml
    :language: xml
    :start-after: <!-- SPHINX_DOME_CAPPRES_BEGIN -->
    :end-before: <!-- SPHINX_DOME_CAPPRES_END -->
    :dedent: 4
    :emphasize-lines: 2,4,7,9,10

The curves are shown below. They show the end points given in the parameter table. The capillary pressure is of the order of 1 to 2 MPa at the irreducible water saturation, which is large compared with the gravity pressure difference between gas and water over tens of meters. This explains why a thick transition zone develops above the contact in the initial state, and why the gas needs a significant pressure to enter the clay.

.. plot:: plot_dome.py plot_saturation_functions
    :align: center
    :caption: Relative permeabilities (left) and gas-water capillary pressures (right) of the sand (solid lines) and of the clay (dashed lines), as functions of the water saturation. The capillary pressure of the clay is obtained from its J-function with the surface tension of the model and with the clay porosity and horizontal permeability.

------------------------------------------------------------------------
Flow solver
------------------------------------------------------------------------

The flow is solved by ``CompositionalMultiphaseFVM``, with a two-point flux approximation (TPFA) discretization, on the sand and clay regions only:

.. literalinclude:: ../../../../../inputFiles/compositionalMultiphaseWell/dome_soreide_whitson_base.xml
    :language: xml
    :start-after: <!-- SPHINX_DOME_FLOW_SOLVER_BEGIN -->
    :end-before: <!-- SPHINX_DOME_FLOW_SOLVER_END -->
    :dedent: 4
    :emphasize-lines: 5,6,7,10

The ``temperature`` attribute fixes the temperature of the reservoir at 344.15 K: the energy balance is not solved. ``useMass="0"`` selects a molar formulation, so that the compositions and the injection streams are expressed as mole fractions. ``targetPhaseVolFractionChangeInTimeStep`` is relaxed to 1.0, so that the time step is selected essentially from the convergence of the Newton method.

The reservoir and the wells are solved together, in a fully implicit way, by ``CompositionalMultiphaseReservoir``, which also holds the nonlinear and linear solver parameters:

.. literalinclude:: ../../../../../inputFiles/compositionalMultiphaseWell/dome_soreide_whitson_base.xml
    :language: xml
    :start-after: <!-- SPHINX_DOME_COUPLED_SOLVER_BEGIN -->
    :end-before: <!-- SPHINX_DOME_COUPLED_SOLVER_END -->
    :dedent: 4
    :emphasize-lines: 3,4,5,9,11,12,15,18,19,22

The time step is increased when the Newton method converges in less than 40 % of ``newtonMaxIter`` iterations (``timeStepIncreaseIterLimit``) and is decreased when it needs more than 90 % of them (``timeStepDecreaseIterLimit``). The time step can be cut up to 100 times in case of failure. The linear system is solved with a flexible GMRES method preconditioned by a multigrid reduction strategy (``mgr``), which is designed for the coupling of the reservoir and well equations. The linear tolerance is adapted along the Newton iterations (``krylovAdaptiveTol``): it is loose at the beginning, and tightened when the nonlinear solution gets closer.

------------------------------------------------------------------------
Initial conditions
------------------------------------------------------------------------

.. _CO2InjectionDepletedReservoirInitialState:

**Initial state**

The initial state is not uniform: the pressure, the composition and the phase saturations are set by a ``HydrostaticEquilibrium`` tag, which computes the capillary-gravity equilibrium of the whole structure from a datum and from a gas-water contact:

.. literalinclude:: ../../../../../inputFiles/compositionalMultiphaseWell/dome_soreide_whitson_base.xml
    :language: xml
    :start-after: <!-- SPHINX_DOME_EQUILIBRIUM_BEGIN -->
    :end-before: <!-- SPHINX_DOME_EQUILIBRIUM_END -->
    :dedent: 4
    :emphasize-lines: 4,5,7,9,10

- **Datum.** The pressure is 27 MPa at the elevation -4453.75 m, which is below the contact. This value is far below the pressure of about 45 MPa that a column of brine in hydrostatic equilibrium would give at this depth, as expected in a depleted field.
- **Contact.** ``phaseContacts`` places the gas-water contact at -4400 m. The vertical pressure of each phase is integrated from the datum with the densities of the fluid, and the saturation in each cell results from the capillary pressure at its elevation above the contact. The initial gas saturation therefore increases progressively from the contact up to the crest, with a transition zone whose thickness is controlled by the capillary pressure curves of each rock type. This is why sand and clay, which have different curves, do not have the same saturation at the same elevation.
- **Composition.** The overall composition is given as a function of elevation by one table per component, linearly interpolated between -5000 m and -4000 m. The lower part of the column is water-rich and the upper part is gas-rich, with a methane fraction that increases with elevation and a small amount of H2S. These tables represent the hydrocarbons that remain after the production.
- **Temperature.** A constant table, equal to the temperature of the flow solver.

.. literalinclude:: ../../../../../inputFiles/compositionalMultiphaseWell/dome_soreide_whitson_base.xml
    :language: xml
    :start-after: <!-- SPHINX_DOME_INIT_TABLES_BEGIN -->
    :end-before: <!-- SPHINX_DOME_INIT_TABLES_END -->
    :dedent: 4
    :emphasize-lines: 2,4,7,9,12,14,17,19,22,24

For instance, the tables give at the contact (-4400 m) the overall mole fractions 0.52 for CH4, 0.03 for H2S and 0.45 for H2O, with no CO2. The number of increments in the table of hydrostatic pressure is controlled by ``elevationIncrementInHydrostaticPressureTable``, here 1 m.

**Boundary conditions**

No boundary condition is imposed. All the outer boundaries of the domain are closed (no-flow), so that the reservoir is a closed system in which the fluid can only enter or leave through the wells. The pressure changes are therefore driven by the compressibility of the fluids and of the rock.

.. _CO2InjectionDepletedReservoirWellControls:

**Well controls**

The three wells are defined in the ``WellManager``. Each well has a main control (``control``), a table that gives the target rate as a function of time, and a limit on the bottom-hole pressure (BHP) that protects the well and the formation. In all of them, ``useSurfaceConditions="1"`` means that the rates are volumetric rates at 101325 Pa and 288.71 K.

The water injector injects pure water:

.. literalinclude:: ../../../../../inputFiles/compositionalMultiphaseWell/dome_soreide_whitson_base.xml
    :language: xml
    :start-after: <!-- SPHINX_DOME_WELL_WATER_INJECTOR_BEGIN -->
    :end-before: <!-- SPHINX_DOME_WELL_WATER_INJECTOR_END -->
    :dedent: 6
    :emphasize-lines: 7,8,15,19,20

The gas injector injects a mixture of 85 % CO2 and 15 % CH4. The stream is given in the order of the components, CH4, CO2, H2S and H2O:

.. literalinclude:: ../../../../../inputFiles/compositionalMultiphaseWell/dome_soreide_whitson_base.xml
    :language: xml
    :start-after: <!-- SPHINX_DOME_WELL_GAS_INJECTOR_BEGIN -->
    :end-before: <!-- SPHINX_DOME_WELL_GAS_INJECTOR_END -->
    :dedent: 6
    :emphasize-lines: 7,8,15,19,20

The producer is controlled by the volumetric rate of the gas phase, with a lower limit on its BHP:

.. literalinclude:: ../../../../../inputFiles/compositionalMultiphaseWell/dome_soreide_whitson_base.xml
    :language: xml
    :start-after: <!-- SPHINX_DOME_WELL_PRODUCER_BEGIN -->
    :end-before: <!-- SPHINX_DOME_WELL_PRODUCER_END -->
    :dedent: 6
    :emphasize-lines: 7,8,15,16,19

.. list-table:: Well controls
    :widths: 20 35 45
    :header-rows: 1
    :class: align-top

    * - Well
      - *GEOS* input
      - Meaning
    * - Water injector
      - ``control="totalVolRate"`` with ``InjectionVolumeRateConstraint``, ``injectionStream`` = {0, 0, 0, 1}
      - Injection of water at the target rate of the table ``WATER.INJECTOR.RATE``, at the reservoir temperature. The BHP must remain below 50 MPa
    * - Gas injector
      - ``control="totalVolRate"`` with ``InjectionVolumeRateConstraint``, ``injectionStream`` = {0.15, 0.85, 0, 0}
      - Injection of the CO2-CH4 mixture at the target rate of the table ``GAS.INJECTOR.RATE``, at the reservoir temperature. The BHP must remain below 50 MPa
    * - Producer
      - ``control="phaseVolRate"`` for the ``gas`` phase with ``ProductionPhaseVolumeRateConstraint``
      - Production at the target gas rate of the table ``PRODUCER.RATE``. The BHP must remain above 1 MPa
    * - Crossflow
      - ``enableCrossflow="0"`` for the injectors, ``"1"`` for the producer
      - The injectors can only inject. The producer is allowed to let fluid flow back into the formation if the pressure of a perforation is lower than that of the wellbore

The BHP constraints are safeguards. When the reservoir cannot deliver or accept the target rate within the pressure limit of the well, the rate cannot be honoured, and the pressure limit becomes the active constraint. The pressure at the well cells shown in the Results section is a useful indicator of the proximity to these limits.

**Schedule**

The target rates are given by three tables of the time (in seconds), which are piecewise constant (``interpolation="lower"``): the value at a given time is the value of the last coordinate that is not greater than this time. The first coordinate of each table, -31557600 s (one year before the start), only defines the value of the rate before the first event, which is zero.

.. literalinclude:: ../../../../../inputFiles/compositionalMultiphaseWell/dome_soreide_whitson_base.xml
    :language: xml
    :start-after: <!-- SPHINX_DOME_SCHEDULE_TABLES_BEGIN -->
    :end-before: <!-- SPHINX_DOME_SCHEDULE_TABLES_END -->
    :dedent: 4
    :emphasize-lines: 5,6,12,13,19,20

The coordinates are multiples of one year (31557600 s). Reading the tables in years gives the schedule of the first table of this page:

- the producer table, ``{-1, 0, 5, 8}`` years with the values ``{0, 15, 15, 0}``, gives a production of 15 m\ :sup:`3`/s between 0 and 8 years. The coordinate at 5 years does not change the value; it is the date of the start of the water injection;
- the water injector table, ``{-1, 5, 8, 100}`` years with the values ``{0, 0.2, 0, 0}``, gives an injection of 0.2 m\ :sup:`3`/s between 5 and 8 years;
- the gas injector table, ``{-1, 13, 18}`` years with the values ``{0, 20, 0}``, gives an injection of 20 m\ :sup:`3`/s between 13 and 18 years.

Since the producer remains open during the water flooding, the water injection is accompanied by production, and the period between 8 and 13 years, during which all the wells are closed, separates the water flooding from the gas injection.

------------------------------------------------------------------------
Events and output
------------------------------------------------------------------------

The simulation covers 20 years. The ``maxTime`` of the events is slightly larger than 20 years (632016000 s), so that the output at 20 years is reached. The events are of three kinds.

.. literalinclude:: ../../../../../inputFiles/compositionalMultiphaseWell/dome_soreide_whitson_base.xml
    :language: xml
    :start-after: <!-- SPHINX_DOME_EVENTS_BEGIN -->
    :end-before: <!-- SPHINX_DOME_EVENTS_END -->
    :dedent: 2
    :emphasize-lines: 5,8,18,23,25,26,27,30,32,33,34,37,39,40

- **Solver applications.** The coupled solver is applied by three periodic events, which differ by the maximum time step (``maxEventDt``). The maximum time step is 100 days (8640000 s) before 5 years and after 5 years and 10 days. Between 5 years and 5 years and 10 days, it is limited to 1 day (86400 s), in order to handle the sudden start of the water injection. In all the cases, the actual time step is chosen by the solver from the convergence of the Newton method, and it is truncated in order to reach the output times.
- **Outputs.** The VTK output is written every half year (15778800 s), and ``targetExactTimestep="1"`` forces the time steps to end exactly at the output times. Restart files are written every 50 cycles.
- **Statistics.** The ``CompositionalMultiphaseStatistics`` task is executed every year (31557600 s). With ``computeRegionStatistics="1"`` and ``writeCSV="1"``, it writes the average pressure and the mass of every phase and component in each region. These files give the mass balance of the injected CO2, which is not available from the maps of the Results section.

The output and the statistics tasks are defined in the following blocks:

.. literalinclude:: ../../../../../inputFiles/compositionalMultiphaseWell/dome_soreide_whitson_base.xml
    :language: xml
    :start-after: <!-- SPHINX_DOME_OUTPUTS_BEGIN -->
    :end-before: <!-- SPHINX_DOME_OUTPUTS_END -->
    :dedent: 2
    :emphasize-lines: 3,6

.. literalinclude:: ../../../../../inputFiles/compositionalMultiphaseWell/dome_soreide_whitson_base.xml
    :language: xml
    :start-after: <!-- SPHINX_DOME_STATISTICS_BEGIN -->
    :end-before: <!-- SPHINX_DOME_STATISTICS_END -->
    :dedent: 2
    :emphasize-lines: 3,6,8

------------------------------------------------------------------------
Post-processing and reporting
------------------------------------------------------------------------

The VTK output of *GEOS* is a ``.pvd`` index file that points to one ``.vtm`` file per output time, each one being split in one ``.vtu`` file per region and per MPI rank. This format is convenient for ParaView but not for scripts, and the documentation does not store pictures. Two scripts, located next to the input files, bridge the gap:

.. code-block:: console

  inputFiles/compositionalMultiphaseWell/dome_vtk_to_hdf5.py
  inputFiles/compositionalMultiphaseWell/plot_dome.py

``dome_vtk_to_hdf5.py`` converts the whole VTK output into one HDF5 file, independent of the number of MPI ranks. It reads the ``.pvd`` file, merges the ranks of each region, orders the cells in the same way at every time, and stores the following datasets (it needs the python packages ``numpy``, ``h5py`` and ``vtk``):

.. list-table:: Content of the HDF5 file
    :widths: 30 70
    :header-rows: 1
    :class: align-top

    * - Dataset
      - Content
    * - ``times_s``, ``times_years``
      - Output times
    * - ``mesh/center``, ``mesh/size``
      - Coordinates of the cell centers and dimensions of the cells (x, y, z)
    * - ``mesh/region_id``
      - Rock type of each cell: 0 for SAND, 1 for CLAY and 2 for CAPROCK. The names are in the attribute ``region_names``
    * - ``fields/pressure``
      - Pressure (Pa), array of dimensions (times, cells)
    * - ``fields/phaseVolumeFraction_gas``, ``fields/phaseVolumeFraction_water``
      - Phase saturations
    * - ``fields/globalCompFraction_CH4``, ``..._CO2``, ``..._H2S``, ``..._H2O``
      - Overall mole fractions of the components

The fields are undefined (NaN) in the caprock, which is not part of the flow problem. The conversion is run once the simulation is complete:

.. code-block:: console

  python dome_vtk_to_hdf5.py VTK.OUTPUT.pvd -o dome_results.hdf5

The option ``--list-fields`` lists the arrays of the VTK files, ``--fields`` chooses the arrays to export and ``--times-years`` keeps only the outputs closest to the given times, which reduces the size of the file.

``plot_dome.py`` defines one function per figure of this page. The documentation calls them with the ``plot`` directive of Sphinx, so that the figures are regenerated from the HDF5 file at each build. The file ``dome_results.hdf5`` must be available next to the script, or its path given in the environment variable ``DOME_HDF5``. The script reads the position of the wells, the gas-water contact and the schedule from the input files, so the figures remain consistent with the input if it is modified. It can also be run directly to write all the figures as PNG files (``python plot_dome.py --outdir figures``).

------------------------------------------------------------------------
Results
------------------------------------------------------------------------

.. The statements of this section describe the mechanisms expected from the input and
   should be checked against the figures of a new run, in particular after any change of
   the schedule, of the mesh or of the fluid model.

**Layout of the results**

The dome is long and thin: its horizontal extent is several times its height. A conventional grid of snapshots, with the variables as columns and the times as rows, would give panels that are too small to read. The results are therefore presented in a different way. Each figure shows one variable, with the key times stacked vertically, in the same horizontal axis and with a single color scale. A given color has the same meaning in all the panels of a figure, and the evolution of any location can be followed by reading the figure from top to bottom. The key times are the ends of the stages of the schedule.

.. list-table:: Key times
    :widths: 20 55 25
    :header-rows: 1
    :class: align-top

    * - Time (yr)
      - Event
      - VTK output
    * - 0
      - Initial state
      - 000000
    * - 5
      - End of the depletion, start of the water injection
      - 000108
    * - 8
      - End of the water flooding, closure of the producer
      - 000368
    * - 13
      - End of the shut-in, start of the gas injection
      - 000418
    * - 15.5
      - 2.5 years after the start of the gas injection
      - 000560
    * - 18
      - End of the gas injection
      - 000656
    * - 20
      - End of the simulation
      - 000721

On the maps, the dashed line is the gas-water contact, the triangles are the wells, the thin lines are the interface between sand and clay, and the gray area is the inactive caprock. The discussion follows the stages of the schedule. It starts with the evolution of the pressure over the whole simulation, which relates every stage to the operations.

**Operations and pressure history**

The first figure puts the schedule (top) next to the evolution of the pressure (bottom). The black curve is the average pressure of the reservoir, weighted by the volume of the cells, and the colored curves are the pressures in the cell of the injectors and in the cell of the producer. The colored bands mark the stages, which are deduced from the rate tables of the input.

.. plot:: plot_dome.py plot_history
    :align: center
    :caption: Target rates of the wells (top; the water rate has its own axis on the right) and pressure history (bottom): volume-averaged pressure of the reservoir, and pressure in the cell of the injectors and in the cell of the producer. The bands mark the stages of the schedule.

The production lowers the average pressure, and the pressure in the producer cell is below the average pressure. The injection has the opposite effect, and the pressure in the injector cell is above the average pressure. During the shut-in the pressure differences between the wells equilibrate, and the average pressure stays almost constant, since the reservoir is closed. At the end, the average pressure results from the balance of the produced and injected volumes. The pressures at the well cells give an indication of how close the wells are to their pressure limits (1 MPa for the producer and 50 MPa for the injectors): a pressure that stays at a limit would show that the target rate is not honoured.

**Pressure**

The pressure maps give the absolute pressure at the key times. The initial map is a hydrostatic distribution: the pressure increases with depth, with a gradient that changes at the gas-water contact because the gas is lighter than the water. In absolute values, the largest part of the color range is due to this vertical gradient, and the changes caused by the operations are only a fraction of it.

.. plot:: plot_dome.py plot_pressure
    :align: center
    :caption: Pressure (MPa) at the key times. The pressure increases with depth, and the gradient changes at the gas-water contact (dashed line).

The second figure removes the hydrostatic distribution and shows the change of pressure with respect to the initial state. It is the most direct picture of the effect of the operations: blue (decrease) is created by the production, and red (increase) is created by the injections.

.. plot:: plot_dome.py plot_pressure_change
    :align: center
    :caption: Change of pressure with respect to the initial state (MPa) at the key times. Blue is a decrease and red an increase of pressure, on a color scale that is symmetric about zero.

- **Depletion (0 - 5 years).** The production lowers the pressure first around the producer. The disturbance then extends toward the whole gas-filled part of the structure, and into the water leg. It travels easily in the sand and is slowed down by the clay, which has a lower permeability.
- **Water flooding (5 - 8 years).** The pressure rises around the injector, at the bottom of the left flank. As the producer is still open, a pressure gradient develops between the two wells, and the reservoir is on average repressurized less than the injector region.
- **Shut-in (8 - 13 years).** The pressure field becomes smoother, as the high and low pressure regions equilibrate. The map at 13 years is the starting point of the gas injection.
- **Gas injection (13 - 18 years).** The pressure rises around the injector, and the average pressure of the reservoir increases with the injected volume. The pressure build-up at the well is the quantity to compare with the limit of the BHP.
- **Post-injection (18 - 20 years).** The overpressure of the injector relaxes, while the plume continues to move.

**Gas saturation**

The gas saturation shows where the gas is, and its evolution is the key to understanding the other maps.

.. plot:: plot_dome.py plot_gas_saturation
    :align: center
    :caption: Gas saturation at the key times. The gas is initially located above the gas-water contact (dashed line), with a transition zone that depends on the rock type.

At the initial time the gas occupies the upper part of the dome, above the contact. The saturation increases from the contact to the crest, with a transition zone governed by the capillary pressure curves. As the sand and the clay have different curves, the saturation changes across their interface.

- During the depletion and the water flooding the injectors are below the contact, so the water that they inject mainly pushes the contact up and repressurizes the structure. The changes of gas saturation remain limited to the vicinity of the contact.
- During the gas injection the gas is injected below the contact, in the water leg. As it is lighter than the water, it rises, with a plume that is formed at the injector and moves up the structure, under the clay, toward the crest, and spreads in the sand. The plume is fed during the five years of injection and its volume grows continuously.
- After the end of the injection the plume keeps migrating and spreading under the effect of gravity, even if the pressure at the well relaxes. The inactive caprock closes the top of the domain. The clay, which has a lower permeability and a higher capillary pressure, slows down the vertical migration of the gas, and the gas tends to accumulate under it.

**Composition: CO2 and CH4**

The maps of the overall mole fractions show the plume that results from the injected mixture, and its mixing with the hydrocarbons of the field. The overall composition includes the fluid dissolved in the water. It is not the composition of the gas phase.

The CO2 fraction is zero in the whole domain until the start of the gas injection, as the initial fluid does not contain CO2. The figure therefore only shows the times after this start.

.. plot:: plot_dome.py plot_co2_fraction
    :align: center
    :caption: Overall mole fraction of CO2 during and after the gas injection. The color scale is the same in all the panels.

The CO2 maps follow the plume of the gas saturation, and give a more precise description of it: the front of CO2 is behind the front of the gas, because the gas that is first displaced is the native gas, mixed with CH4. In the cells that remain water-saturated, the CO2 fraction is small but not zero. It shows the CO2 that dissolves in the brine, which is a part of the storage that the saturation does not show.

The CH4 maps complete the description. The methane fraction of the injected mixture is lower than that of the native gas, so the plume appears as a region where the methane fraction is lowered: the native methane is displaced and diluted by the injected CO2. The CH4 maps also show the changes before the injection, since the overall methane fraction in the gas leg is changed by the production and by the water flooding.

.. plot:: plot_dome.py plot_ch4_fraction
    :align: center
    :caption: Overall mole fraction of CH4 at the key times. The native methane of the gas leg is displaced and diluted where the injected CO2 arrives.

**Conclusion**

The example shows how the stages of a depleted-field storage project are driven in *GEOS* by time tables: production, water flooding, shut-in and CO2 injection, on a heterogeneous structure initialized by capillary-gravity equilibrium. The injection of the CO2-CH4 mixture below the contact leads to a plume that rises and spreads under the clay and the caprock, while the pressure of the field is first lowered by the production, then partially restored by the injections. The two fluid descriptions can be compared on the same case by running ``dome_kvalue_bench.xml`` instead of ``dome_soreide_whitson_bench.xml``, converting its output, and building the pictures again.

------------------------------------------------------------------
To go further
------------------------------------------------------------------

The case can be modified in a few places:

- the schedule: the three tables ``PRODUCER.RATE``, ``WATER.INJECTOR.RATE`` and ``GAS.INJECTOR.RATE`` are the only input that defines the operations. For instance, replacing the coordinates ``{ -31557600, 0, 157788000, 252460800 }`` and the values ``{ 0, 15, 15, 0 }`` of the producer by ``{ -31557600, 0, 157788000 }`` and ``{ 0, 15, 0 }`` closes the producer when the water injection starts. The key times and the stages of the figures should then be reviewed;
- the injected mixture (``injectionStream`` of the gas injector), for instance with a small fraction of H2S, since the third component is already part of the model;
- the salinity (``salinity``) and the equilibrium tables, to study the solubility of CO2 in the brine;
- the mesh (``dome_mesh_smoke.xml`` instead of ``dome_mesh_bench.xml``) for a quick check of the setup.

**Feedback on this example**

For any feedback on this example, please submit a `GitHub issue on the project's GitHub page <https://github.com/GEOS-DEV/GEOS/issues>`_.

------------------------------------------------------------------
References
------------------------------------------------------------------

- Peng, D.-Y., and Robinson, D. B. 1976. A new two-constant equation of state. *Industrial & Engineering Chemistry Fundamentals*, 15(1), 59-64.
- S\ |oslash|\ reide, I., and Whitson, C. H. 1992. Peng-Robinson predictions for hydrocarbons, CO2, N2, and H2S with pure water and NaCl brine. *Fluid Phase Equilibria*, 77, 217-240.
- Phillips, S. L., Igbene, A., Fair, J. A., Ozbek, H., and Tavana, M. 1981. A technical databook for geothermal energy utilization. Lawrence Berkeley Laboratory Report LBL-12810.
