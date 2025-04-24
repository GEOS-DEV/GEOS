.. _ExampleSPE11b:

########################################################################
Assessment of CO2 Storage residual and dissolution trapping mechanism
########################################################################


**Context**

In this article, we consider the benchmark proposed in `(Nordbotten et al.,2024) <https://norceresearch.brage.unit.no/norceresearch-xmlui/bitstream/handle/11250/3137806/spe-218015-pa.pdf?sequence=1>`__
to showcase the effect of the convective mixing in addition to the usually considered residual trapping in a multi-facies resevoir.
Using a thermal formulation, the effect of injecting cold CO2 can be observed.

This example can serve as a guideline to set-up the input XML deck to reproduce the published results on `SPE11 CSP <https://github.com/Simulation-Benchmarks/11thSPE-CSP>`__.
The reader can also find `decks <https://github.com/Simulation-Benchmarks/11thSPE-CSP/tree/main/input_decks>`__ for
other open-source simulators and can try to benchmark a subset of them,
even if the use of other simulators is not in the scope of this article.

.. note::
    Interested reader can refer to the `official comparative website <https://moyner.github.io/SPE11-plot-test-deploy/>`__ to see all
    submitted results

------------------------------------------------------------------------
Brief case description
------------------------------------------------------------------------

As the detailed description is available in `(Nordbotten et al.,2024) <https://norceresearch.brage.unit.no/norceresearch-xmlui/bitstream/handle/11250/3137806/spe-218015-pa.pdf?sequence=1>`__,
we will only briefly mentioned data that will be referred to in the following sections and let the interested reader read
the full version.

The spe11-b is a reservoir-like 2D case rescaled from the `FluidFlower <https://arxiv.org/pdf/2302.10986>`__ experiment.
It consists in a 8.4 km large and 1.2 km deep reservoir which top depth is at 2km from free surface.
A global geothermal gradient of 25K/km is imposed from the bottom surface at 70 Celsius.
As its precursor it is composed of 7 facies with different properties.

.. _spe11b_facies:
.. figure:: ./pictures/spe11b_presentation.png
    :align: center
    :width: 500
    :figclass: align-center

Here blue, red and orange boxes are materialization of the prescribed reporting boxes, respectively denoted box A , B and C in the description.
They are places for observing first anticline accumulation, top anticlines accumulation through heterogeneous structure's dripping and
convective mixing finger structures.

------------------------------------------------------------------------
Input base
------------------------------------------------------------------------

This benchmark test is based on the XML file located below:

.. code-block:: console

  inputFiles/compositionalMultiphaseFlow/benchmarks/SPE11/b/spe11b_vti_source_00840x00120.xml

it includes

.. code-block:: console

  inputFiles/compositionalMultiphaseFlow/benchmarks/SPE11/b/spe11b_vti_source_base.xml

for all its non-discretization-related parameters. The full simulation involves a 1000 years of thermal equilibration.
Hydrostatic equilibration indeed does not includes thermal effects and it is then required to let the system equilibrate
under the geothermal gradient before starting the injection schedule.

The injection schedule is of 50 years from the bottom injector and 25 years with a 25 year delayed start from the top injector.
It is then followed by a 950 years of migration, dissolution and convection.

.. literalinclude:: ../../../../../../../inputFiles/compositionalMultiphaseFlow/benchmarks/SPE11/b/spe11b_vti_source_base.xml
    :language: xml
    :start-after: <!-- SPHINX_SOURCES -->
    :end-before: <!-- SPHINX_SOURCES_END -->

As we can see in the snippet above, `SourceFlux` is the modeling choice for incoming fluxes (rather than wellbores). It goes with a
`FieldSpecification` setting the temperature at the same place to mimic imposed injected *CO2* temperature. This includes to object
left to be defined a set ``thermalSources1`` (respectively ``thermalSources2``) and a `Function` that can vary flux over time. These will
be discussed in the :ref:`BCSection`.


------------------------------------------------------------------------
Constitutives and includes
------------------------------------------------------------------------

.. _KRPCSection:

The capillary and residual trapping in such a multi-facies reservoir occurs due to the heterogeneous properties. For instance, jumping from
anticline facies \#5 to the upper layers a large entry pressure has to be reached. The box A is the design such that it can monitor the
expected trapped *CO2* dissolving and making the surrounding brine heavier. This then will result in sinking brine and start convective
mixing instabilites.

.. plot:: ../inputFiles/compositionalMultiphaseFlow/benchmarks/SPE11/b/include/plotAllKrPc.py

Above are reported the sets of relative permeabilities (left) and capillary pressure (right) from their tabulated values. The facies \#5 is highlighted
as playing an important role in the overall trapping. The import in the simulation deck is done through the ``include/`` folder

.. literalinclude:: ../../../../../../../inputFiles/compositionalMultiphaseFlow/benchmarks/SPE11/b/include/kr.xml
    :language: xml
    :start-after: <!-- SPHINX_KR_TABLE -->
    :end-before: <!-- SPHINX_KR_TABLE_END -->

and through the definition of relative permeabilities and capillary pressure under the `Constitutive` tag. As an example, definition of
 relperm for the facies \#1 is reported,


.. literalinclude:: ../../../../../../../inputFiles/compositionalMultiphaseFlow/benchmarks/SPE11/b/spe11b_vti_source_base.xml
    :language: xml
    :start-after: <!-- SPHINX_KR_CONST -->
    :end-before: <!-- SPHINX_KR_CONST_END -->

There is also properties such as permeabilities, diffusivity and rock thermal conductivity that are modeled as constant per facies
and then needs to be defined as such. They then will be defined as

.. literalinclude:: ../../../../../../../inputFiles/compositionalMultiphaseFlow/benchmarks/SPE11/b/include/properties_vti.xml
    :language: xml
    :start-after: <!-- SPHINX_THC_FIELD -->
    :end-before: <!-- SPHINX_THC_FIELD_END -->

scaling their `Constitutive` 's values that are defined,

.. literalinclude:: ../../../../../../../inputFiles/compositionalMultiphaseFlow/benchmarks/SPE11/b/spe11b_vti_source_base.xml
    :language: xml
    :start-after: <!-- SPHINX_THC_CONST -->
    :end-before: <!-- SPHINX_THC_CONST_END -->

.. note::

    Their `FieldSpecification` scaling definition is linked to their `Constitutive` 's through the generated field name, for instance *thermalCond_rockThermalConductivity*.
    The prefix being extracted from the `Constitutive` 's `name` attributes.

 .. note::

     The permeabilities and porosities being also heterogeneous per facies, we let the reader inspect the files in ``include/`` and make themselves aware
     of the places where fields are imported and defined similarly to the rock thermal conductivities.

The diffusivity being non heterogeneous in the specifications of the case, it is solely defined in its `Constitutive` model.

Eventually, one has to pick a model for both brine and *CO2* densities and viscosities with respect to varying pressure, temperature and
composition.

.. literalinclude:: ../../../../../../../inputFiles/compositionalMultiphaseFlow/benchmarks/SPE11/b/spe11b_vti_source_base.xml
    :language: xml
    :start-after: <!-- SPHINX_EOS_CONST -->
    :end-before: <!-- SPHINX_EOS_CONST_END -->

It is done using the `CO2BrinePhillipsThermalFluid` tag. It lists files that includes table generator's parameters. For instance, *CO2*'s density
is modeled through *Span Wagner* model, while the viscosities are obtained via *Fenghour* model. For more detail, read :ref:`CO2-EOS`.
Their conformity to NIST values is assessed in the following plot,

.. plot:: ../inputFiles/compositionalMultiphaseFlow/benchmarks/SPE11/b/tables/plotNIST.py

 Above is the one to one comparison with respect to pressure range used for a variety of pressures.

------------------------------------------------------------------------
Flow solver
------------------------------------------------------------------------

Now that the physics in our case is set, with its global and spatialized parameters defined, let's dive into `FlowSolver`
tunning. It is how we connect together these definition, up to a initial and boundary conditions (discussed in the next section).

------------------------------------------------------------------------
Initial and boundary conditions
------------------------------------------------------------------------
.. _BCSection:

In the following, we will cover how to set both the injections and the left and right buffers used to mimic large
connected aquifers. This will also cover the initialization.


As mentioned above, the injection description lacks a `Geometry` cell set to be applied on and a `Function` varying incoming flux over time.
The first point is defined in the discretization-specific file (also root file) and the latter is found in the base file.

.. literalinclude:: ../../../../../../../inputFiles/compositionalMultiphaseFlow/benchmarks/SPE11/b/spe11b_vti_source_00840x00120.xml
    :language: xml
    :start-after: <!-- SPHINX_GEOM -->
    :end-before: <!-- SPHINX_GEOM_END -->

.. note::
    Note that the units here are in *mass* as the `useMass=1` attribute is set in `Solver` used.
    Default will be volumetric.


The lower-interpolated 1D function over time serve as a scaler for the incoming flux.

.. literalinclude:: ../../../../../../../inputFiles/compositionalMultiphaseFlow/benchmarks/SPE11/b/spe11b_vti_source_base.xml
    :language: xml
    :start-after: <!-- SPHINX_FUNCTIONS -->
    :end-before: <!-- SPHINX_FUNCTIONS_END -->

From the `Functions` we can also imposed the initial geo-thermal gradient as well as initial compositions.
Then we need an `HydrostaticEquilibrium` to be computed from the inital pressure and composition, as shown below:

.. literalinclude:: ../../../../../../../inputFiles/compositionalMultiphaseFlow/benchmarks/SPE11/b/spe11b_vti_source_base.xml
    :language: xml
    :start-after: <!-- SPHINX_EQUIL -->
    :end-before: <!-- SPHINX_EQUIL_END -->

Then, the only thing left to fully set the initialization is for the simulation to start at *-1000 years* as it is done in the `Events`

.. literalinclude:: ../../../../../../../inputFiles/compositionalMultiphaseFlow/benchmarks/SPE11/b/spe11b_vti_source_base.xml
    :language: xml
    :start-after: <!-- SPHINX_EVENTS -->
    :end-before: <!-- SPHINX_EVENTS_END -->

.. note::
    The time tags in *GEOS* are in seconds. Here a split between equilibration, start of the injection and injection-post-injection is chosen as
    the maximal dt from experience in these time range is quite different. ``solverApplication1`` is then ensuring stability as the injection starts.

We are then left with the imposition of domain boundary conditions. The top and bottom temperatures are imposed with the help of the subfile
included as seen in :ref:`KRPCSection`. Content of the inclusion is FieldSpecification on geometrical sets defined earlier in this section.

.. literalinclude:: ../../../../../../../inputFiles/compositionalMultiphaseFlow/benchmarks/SPE11/b/include/dirichlet_boundary_vti.xml
    :language: xml
    :start-after: <!-- SPHINX_TOP_BOTTOM -->
    :end-before: <!-- SPHINX_TOP_BOTTOM_END -->

The last point to be tackled is how to set the *fictive* aquifers that are used to damped pressure build up that would occur otherwise
in such a 2D constrained domain. It is done as for the sources to the exception that, instead of capturing the cellset thanks to a box,
it has already been tagged as a region in the construction of the mesh, namely ``12_hexahedra`` to ``15_hexahedra``. It is then easy to
rescale the volume in those regions,

.. literalinclude:: ../../../../../../../inputFiles/compositionalMultiphaseFlow/benchmarks/SPE11/b/include/properties_vti.xml
    :language: xml
    :start-after: <!-- SPHINX_BUFFER_FIELDS -->
    :end-before: <!-- SPHINX_BUFFER_FIELDS_END -->

that then is scaled in the discretization-specific file, as the scaling factor depends on the minimal distance that the discretization is
able to capture, hence the discretization itself. (see.  `(Nordbotten et al.,2024) <https://norceresearch.brage.unit.no/norceresearch-xmlui/bitstream/handle/11250/3137806/spe-218015-pa.pdf?sequence=1>`__, for details about this buffer definition).

.. literalinclude:: ../../../../../../../inputFiles/compositionalMultiphaseFlow/benchmarks/SPE11/b/spe11b_vti_source_00840x00120.xml
    :language: xml
    :start-after: <!-- SPHINX_BUFFER_FN -->
    :end-before: <!-- SPHINX_BUFFER_FN_END -->

We are now set to launch the simulation.

------------------------------------------------------------------------
Post-treating and prescribed reports
------------------------------------------------------------------------

First inspection of the results is done using Paraview as in the usual GEOS workflow. A short time after injection, the pressure
build up has been damped by the boundary *buffers* and is only inspected as a dynamic response.

Hereafter are reported respectively, temperature, saturation and dissolved CO2 fraction for 500 years and 1000 years.

.. _spe11b_T_500:
.. figure:: ./pictures/spe11b_thermal_T.0100.png
    :align: left
    :width: 250
    :figclass: align-left
.. _spe11b_T_1000:
.. figure:: ./pictures/spe11b_thermal_T.0200.png
    :align: right
    :width: 250
    :figclass: align-right
.. _spe11b_s_500:
.. figure:: ./pictures/spe11b_thermal_sat.0100.png
    :align: left
    :width: 250
    :figclass: align-left
.. _spe11b_s_1000:
.. figure:: ./pictures/spe11b_thermal_sat.0200.png
    :align: right
    :width: 250
    :figclass: align-right
.. _spe11b_x_500:
.. figure:: ./pictures/spe11b_thermal_xcp.0100.png
    :align: left
    :width: 250
    :figclass: align-left
.. _spe11b_x_1000:
.. figure:: ./pictures/spe11b_thermal_xcp.0200.png
    :align: right
    :width: 250
    :figclass: align-right

Firstly, saturation is showing that, for this mesh resolution, almost all gaseous CO2 is trapped then dissolved after 1000 years
of the injection scenario. Then the reporting of the dissolved fraction clearly shows evidence of on-set of convective mixing
fingers, with heavier CO2-saturated brine sinking. Reader can in particular focus on the lower left part of the domain, perpendicular
to the first injector. Here dissolved CO2 will sink and accumulate. This observation is repeated looking at the temperature maps report.


Here we will draft base of reporting scripts ...

------------------------------------------------------------------
To go further
------------------------------------------------------------------

.. The more complex 3D case is presented in  :ref:`Placeholder`.

**Feedback on this example**

For any feedback on this example, please submit a `GitHub issue on the project's GitHub page <https://github.com/GEOS-DEV/GEOS/issues>`_.




