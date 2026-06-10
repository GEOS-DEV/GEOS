.. _ExampleSPE11b:

########################################################################
Assessment of CO2 Storage residual and dissolution trapping mechanism
########################################################################


**Context**

We consider the benchmark proposed in `(Nordbotten et al., 2024) <https://onepetro.org/SJ/article-pdf/29/05/2507/4090996/spe-218015-pa.pdf/1>`__
to showcase the effect of the convective mixing in addition to the usually considered residual trapping in a multi-facies resevoir.
Also, by using a thermal formulation, the effect of injecting cold CO2 can be observed.

This example can serve as a guideline to set-up the input XML deck to reproduce the published results on `SPE11 CSP <https://github.com/Simulation-Benchmarks/11thSPE-CSP>`__.
The input decks for other simulators can be found at `decks <https://github.com/Simulation-Benchmarks/11thSPE-CSP/tree/main/input_decks>`__.

.. note::
    We refer to the `official comparative website <https://moyner.github.io/SPE11-plot-test-deploy/>`__ to see all
    submitted results.

------------------------------------------------------------------------
Brief case description
------------------------------------------------------------------------

As the detailed description is available in `(Nordbotten et al.,2024) <https://onepetro.org/SJ/article-pdf/29/05/2507/4090996/spe-218015-pa.pdf/1>`__,
we will only briefly review the data set that will be used to in the following sections and let the interested reader read
the full version.

The spe11-b is a reservoir-like 2D case rescaled from the `FluidFlower <https://link.springer.com/article/10.1007/s11242-023-01977-7>`__ experiment.
It consists in a 8.4 km large and 1.2 km deep reservoir which top depth is at 2km from the surface.
A global geothermal gradient of 25 °C/km is imposed from the bottom surface at 70 °C.
The reservoir is composed of 7 facies with different properties:

.. _spe11b_facies:
.. figure:: ./pictures/spe11b_presentation.png
    :align: center
    :width: 500
    :figclass: align-center

    Schematic representation of SPE11B case and its reporting boxes. Image extracted from `arxiv version <https://arxiv.org/abs/2507.15861>`__

Here blue, red and orange boxes are representing the prescribed reporting boxes, respectively denoted box A, B and C from the description.
These locations are used to observe the initial accumulation in the first anticline, the accumulation at the top of anticlines due to fluid
migration through heterogeneous structures, and the development of convective mixing finger patterns.

------------------------------------------------------------------------
Input base
------------------------------------------------------------------------

This benchmark test is based on the XML file located below:

.. code-block:: console

  inputFiles/compositionalMultiphaseFlow/benchmarks/SPE11/b/spe11b_vti_source_00840x00120.xml

which includes

.. code-block:: console

  inputFiles/compositionalMultiphaseFlow/benchmarks/SPE11/b/spe11b_vti_source_base.xml

for all parameters not related to the mesh discretization. The full simulation involves 1000 years of initial thermal equilibration
required to let the system stabilize under the geothermal gradient before starting the injection schedule.

The well schedule consists of 50 years of injection for the bottom injector and 25 years with a 25 year delayed start for the top injector.
It is then followed by a 950 years of migration, dissolution and convection.

As we can see in the snippet below, `SourceFlux` was choosen to model the injection (rather than wellbores). It goes with a
`FieldSpecification` setting for the temperature at the same place to define the imposed injected *CO2* temperature. This also includes to object
defined elsewhere as a set ``thermalSources1`` (respectively ``thermalSources2``) and a `Function` that varies over time to represent the flux schedule.
These will be discussed in more details below.

.. literalinclude:: ../../../../../../../inputFiles/compositionalMultiphaseFlow/benchmarks/SPE11/b/spe11b_vti_source_base.xml
    :language: xml
    :start-after: <!-- SPHINX_SOURCES -->
    :end-before: <!-- SPHINX_SOURCES_END -->

------------------------------------------------------------------------
Constitutives and includes
------------------------------------------------------------------------

.. _KRPCSection:

The capillary and residual trapping in such a multi-facies reservoir occur due to the heterogeneous properties. For instance, a large entry pressure
jump exists between anticline facies \#5 and the upper layers. The box A is chosen such that it can monitor the expected trapped *CO2* dissolution.
The dissolved *CO2* makes the brine heavier which results in its sinking and starts the convective mixing instabilites.

Below we show the sets of relative permeabilities (left) and capillary pressure (right) produced from their tabulated values. The facies \#5 is highlighted
as it plays an important role in the overall trapping.

.. plot:: ../inputFiles/compositionalMultiphaseFlow/benchmarks/SPE11/b/include/plotAllKrPc.py

The tables import in the simulation deck is done as the following showing the relative permeability and capillary pressure tables for facies \#5:

.. literalinclude:: ../../../../../../../inputFiles/compositionalMultiphaseFlow/benchmarks/SPE11/b/include/kr.xml
    :language: xml
    :start-after: <!-- SPHINX_KR_TABLE_FACIES_5_START -->
    :end-before: <!-- SPHINX_KR_TABLE_FACIES_5_END -->

.. literalinclude:: ../../../../../../../inputFiles/compositionalMultiphaseFlow/benchmarks/SPE11/b/include/kr.xml
    :language: xml
    :start-after: <!-- SPHINX_PC_TABLE_FACIES_5_START -->
    :end-before: <!-- SPHINX_PC_TABLE_FACIES_5_END -->

Then, relative permeabilities and capillary pressure are defined under the `Constitutive` tag. For example, the relative permeability for the facies \#1 
is defined as

.. literalinclude:: ../../../../../../../inputFiles/compositionalMultiphaseFlow/benchmarks/SPE11/b/spe11b_vti_source_base.xml
    :language: xml
    :start-after: <!-- SPHINX_KR_CONST -->
    :end-before: <!-- SPHINX_KR_CONST_END -->

Other properties such as permeabilities, diffusivity, and rock thermal conductivity are modeled as constant per facies
(thermal conductivity definition is shown as an example):

.. literalinclude:: ../../../../../../../inputFiles/compositionalMultiphaseFlow/benchmarks/SPE11/b/spe11b_vti_source_base.xml
    :language: xml
    :start-after: <!-- SPHINX_THC_CONST -->
    :end-before: <!-- SPHINX_THC_CONST_END -->

In addition, the scaling of the `Constitutive` values is defined to represent the facies:

.. literalinclude:: ../../../../../../../inputFiles/compositionalMultiphaseFlow/benchmarks/SPE11/b/include/properties_vti.xml
    :language: xml
    :start-after: <!-- SPHINX_THC_FIELD -->
    :end-before: <!-- SPHINX_THC_FIELD_END -->

.. note::

    Their `FieldSpecification` scaling definition is linked to the `Constitutive` definition through the generated field name, for instance *thermalCond_rockThermalConductivity*.
    The prefix being extracted from the `Constitutive` `name` attributes.

The permeabilities and porosities are defined similarly to the rock thermal conductivities using the scaling approach.
The diffusivity being non-heterogeneous in the specifications of the case, it is solely defined in its `Constitutive` model.

Eventually, one has to pick a model for both brine and *CO2* densities and viscosities with respect to varying pressure, temperature and
composition. It is done using the `CO2BrinePhillipsThermalFluid` tag. It lists files that include the related parameters. For instance, *CO2* density
is modeled through *Span Wagner* model, while the viscosities are obtained via *Fenghour* model (see :ref:`CO2-EOS` for more details).

.. literalinclude:: ../../../../../../../inputFiles/compositionalMultiphaseFlow/benchmarks/SPE11/b/spe11b_vti_source_base.xml
    :language: xml
    :start-after: <!-- SPHINX_EOS_CONST -->
    :end-before: <!-- SPHINX_EOS_CONST_END -->

The model conformity to NIST values is assessed in the following plots:

.. plot:: ../inputFiles/compositionalMultiphaseFlow/benchmarks/SPE11/b/tables/plotNIST.py

Fluid model comparison with NIST data. The dotted lines represent *GEOS* values while the solid-and-stars are the values from NIST database.

The plots demonstrate the validity of the selected models with respect to the provided NIST tables over the range of pressures and temperatures.

------------------------------------------------------------------------
Flow solver
------------------------------------------------------------------------

The `CompositionalMultiphaseFVM` solver uses the finite volume discretization on the `targetRegions`. The `discretization` will be TPFA


.. literalinclude:: ../../../../../../../inputFiles/compositionalMultiphaseFlow/benchmarks/SPE11/b/spe11b_vti_source_base.xml
    :language: xml
    :start-after: <!-- SPHINX_TPFA -->
    :end-before: <!-- SPHINX_TPFA_END -->

and here is linked-by-name to the `TwoPointFluxApproximation` XML block. The *Phase Potential Upwind (PPU)* scheme is used.

.. literalinclude:: ../../../../../../../inputFiles/compositionalMultiphaseFlow/benchmarks/SPE11/b/spe11b_vti_source_base.xml
    :language: xml
    :start-after: <!-- SPHINX_FLOW_SOL -->
    :end-before: <!-- SPHINX_FLOW_SOL_END -->


Note the important `isThermal` and `useMass` flags are enabled. The `targetPhaseVolFractionChangeInTimeStep` parameter controls the time step selection
and the `maxCompFractionChange` parameter controls the solution chopping strategy in the solver. The `NonlinearSolverParameters` contains the nonlinear tolerance,
*line search* parameters, and the parameters for the time step increase and decrease. The `LinearSolverParameters` contains parameters for the linear solver:
the standard pick here is *Flexible GMRES* iterative solve with *Multi Grid Reduction* preconditioner.

------------------------------------------------------------------------
Initial and boundary conditions
------------------------------------------------------------------------
.. _BCSection:

In the following, we will provide more details about the initialization and describe how to set both the injections and the left and right buffers used to mimic large
connected aquifers.

As mentioned above, the injection description refers to a cell set to be applied on and a `Function` defining incoming flux over time.

The cell sets are defined using boxes in the discretization-specific file (also root file):

.. literalinclude:: ../../../../../../../inputFiles/compositionalMultiphaseFlow/benchmarks/SPE11/b/spe11b_vti_source_00840x00120.xml
    :language: xml
    :start-after: <!-- SPHINX_GEOM -->
    :end-before: <!-- SPHINX_GEOM_END -->

The lower-interpolated 1D function over time serve as a scaler for the incoming flux:

.. literalinclude:: ../../../../../../../inputFiles/compositionalMultiphaseFlow/benchmarks/SPE11/b/spe11b_vti_source_base.xml
    :language: xml
    :start-after: <!-- SPHINX_FUNCTIONS -->
    :end-before: <!-- SPHINX_FUNCTIONS_END -->

.. note::
    Note that the units here are in *mass* as the `useMass=1` attribute is set in `Solvers`.

Using `Functions` we also imposed the initial geothermal gradient as well as initial compositions.
Then we need to define the `HydrostaticEquilibrium` to be computed, as shown below:

.. literalinclude:: ../../../../../../../inputFiles/compositionalMultiphaseFlow/benchmarks/SPE11/b/spe11b_vti_source_base.xml
    :language: xml
    :start-after: <!-- SPHINX_EQUIL -->
    :end-before: <!-- SPHINX_EQUIL_END -->

Finally, the the equilibration stage is defined in the `Events` by setting the simulation to start at *-1000 years*:

.. literalinclude:: ../../../../../../../inputFiles/compositionalMultiphaseFlow/benchmarks/SPE11/b/spe11b_vti_source_00840x00120.xml
    :language: xml
    :start-after: <!-- SPHINX_EVENTS -->
    :end-before: <!-- SPHINX_EVENTS_END -->

.. note::
    The time tags in *GEOS* are in seconds. Here a split between equilibration, start of the injection and post-injection is chosen as
    the maximum time steps in these stages are quite different. Specifically, separate ``solverApplication1`` definition is ensuring stability as the injection starts.

The boundary conditions for the domain are imposed as follows. The top and bottom temperatures are set using the `FieldSpecification` on `top` and `bottom` geometrical sets defined above:

.. literalinclude:: ../../../../../../../inputFiles/compositionalMultiphaseFlow/benchmarks/SPE11/b/include/dirichlet_boundary_vti.xml
    :language: xml
    :start-after: <!-- SPHINX_TOP_BOTTOM -->
    :end-before: <!-- SPHINX_TOP_BOTTOM_END -->

*Fictive* aquifers are used as buffers to damp the pressure buildup
(see `(Nordbotten et al.,2024) <https://onepetro.org/SJ/article-pdf/29/05/2507/4090996/spe-218015-pa.pdf/1>`__, for details about these buffers definition).
They are set similar to the sources definition but instead of using the boxes for cell sets,
we use tags created during the construction of the mesh, namely ``12_hexahedra`` to ``15_hexahedra``. It is then easy to
rescale the volumes in those regions:

.. literalinclude:: ../../../../../../../inputFiles/compositionalMultiphaseFlow/benchmarks/SPE11/b/include/properties_vti.xml
    :language: xml
    :start-after: <!-- SPHINX_BUFFER_FIELDS -->
    :end-before: <!-- SPHINX_BUFFER_FIELDS_END -->

The scaling factor depends on the discretization and hence that is defined in the discretization-specific file:

.. literalinclude:: ../../../../../../../inputFiles/compositionalMultiphaseFlow/benchmarks/SPE11/b/spe11b_vti_source_00840x00120.xml
    :language: xml
    :start-after: <!-- SPHINX_BUFFER_FN -->
    :end-before: <!-- SPHINX_BUFFER_FN_END -->


------------------------------------------------------------------------
Post-processing and reporting
------------------------------------------------------------------------

First, inspection of the results is done using Paraview as in the usual GEOS workflow. Here, we report temperature, saturation and dissolved CO2 fraction for 500 years and 1000 years:

.. _spe11b_T_500:
.. figure:: ./pictures/spe11b_thermal_T.0100.png
    :align: center
    :figwidth: 45%
    :figclass: align-center

    Thermal state after 500 years

.. _spe11b_T_1000:
.. figure:: ./pictures/spe11b_thermal_T.0200.png
    :align: center
    :figwidth: 45%
    :figclass: align-center

    Thermal state after 1000 years

.. .. _spe11b_s_500:
.. .. figure:: ./pictures/spe11b_thermal_sat.0100.png
..     :align: left
..     :width: 250
..     :figclass: align-left
..
.. .. _spe11b_s_1000:
.. .. figure:: ./pictures/spe11b_thermal_sat.0200.png
..     :align: right
..     :width: 250
..     :figclass: align-right
..
.. .. _spe11b_x_500:
.. .. figure:: ./pictures/spe11b_thermal_xcp.0100.png
..     :align: left
..     :width: 250
..     :figclass: align-left
..
.. .. _spe11b_x_1000:
.. .. figure:: ./pictures/spe11b_thermal_xcp.0200.png
..     :align: right
..     :width: 250
..     :figclass: align-right


The saturation distribution is showing that, for this mesh resolution, almost all gaseous CO2 is trapped and then dissolved after 1000 years
of the injection scenario. The reporting of the dissolved fraction clearly shows evidence of on-set of convective mixing
fingers, with heavier CO2-saturated brine sinking. In particular, the dissolved CO2 sinks and accumulates in the lower left part of the domain, 
perpendicular to the first injector. This behavior is similar for the temperature distributions.


GEOS run can be post-treated leveraging `Paraview <https://www.paraview.org/>`_ (as mentioned in other tutorials) or `pyvtk <https://pypi.org/project/PyVTK/>`_. Here is an example of the latter.
It is built around two main driver `write()` and `plot()` functions which respectively write the integrated mass balance in
the reporting boxes A and B and re-read the generated file and plot it to png. It requires usual dependencies such as `re`, `os`, `pandas`,
`numpy`, `scipy (for interpolation) <https://docs.scipy.org/doc/scipy/tutorial/interpolate.html>`_ and `pyplot`

.. literalinclude:: ./example_script.py
    :language: python
    :start-after: #SPHINX_IMPORT_BEGIN
    :end-before: #SPHINX_IMPORT_END

in addition to `pyvtk <https://pypi.org/project/PyVTK/>`_. The main drivers hereafter,

.. literalinclude:: ./example_script.py
    :language: python
    :start-after: #SPHINX_MAIN_BEGIN
    :end-before: #SPHINX_MAIN_END

are using simple data structure to gather info and to link GEOS' name to name used in formulas,

.. literalinclude:: ./example_script.py
    :language: python
    :start-after: #SPHINX_DATA_BEGIN
    :end-before: #SPHINX_DATA_END

It uses a string interpreted function `process_keys(formula, fielddict)` to build the observable from the fields loaded from GEOS :

.. literalinclude:: ./example_script.py
    :language: python
    :start-after: #SPHINX_UTILS_BEGIN
    :end-before: #SPHINX_UTILS_END

and a tool set to interpolate and integrate them over the correct boxes,

.. literalinclude:: ./example_script.py
    :language: python
    :start-after: #SPHINX_INTERPOLATE_BEGIN
    :end-before: #SPHINX_INTERPOLATE_END

Eventually, the per-time driver, load usefull specific fields from the correct time in the multi-block VTK format GEOS is using,

.. literalinclude:: ./example_script.py
    :language: python
    :start-after: #SPHINX_PTIME_BEGIN
    :end-before: #SPHINX_PTIME_END

and this function is mapped over the schedule of the simulation.


------------------------------------------------------------------
To go further
------------------------------------------------------------------

.. The more complex 3D case is presented in  :ref:`Placeholder`.

**Feedback on this example**

For any feedback on this example, please submit a `GitHub issue on the project's GitHub page <https://github.com/GEOS-DEV/GEOS/issues>`_.




