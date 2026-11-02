.. _FieldSpecification:

Initial and Boundary Conditions
==================================

Initial conditions, boundary conditions and source terms are all declared in the
``<FieldSpecifications>`` block. Each entry says *which value* is applied to *which field*,
on *which mesh objects*, and *when*.

Generic field specification
---------------------------

The basic element is ``FieldSpecification``. The example below, from
``inputFiles/singlePhaseFlow/3D_10x10x10_compressible_base.xml``, sets the initial pressure
in all cells, then holds the pressure fixed in the cells of the ``source`` and ``sink`` sets:

.. literalinclude:: ../../../../inputFiles/singlePhaseFlow/3D_10x10x10_compressible_base.xml
  :language: xml
  :start-after: <!-- SPHINX_TUT_INT_HEX_FIELDS -->
  :end-before: <!-- SPHINX_TUT_INT_HEX_FIELDS_END -->

Its main attributes are:

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Attribute
     - Meaning
   * - ``objectPath``
     - Where the field lives: ``ElementRegions/<region>/<cellBlock>`` (or ``ElementRegions/<region>``,
       ``ElementRegions``) for cell fields, ``nodeManager`` for nodal fields (e.g. displacement),
       ``faceManager`` for face fields.
   * - ``fieldName``
     - Name of the field, e.g. ``pressure``, ``temperature``, ``globalCompFraction``,
       ``totalDisplacement``. The fields registered on each mesh object are listed in the
       :doc:`/docs/sphinx/datastructure/CompleteXMLSchema`.
   * - ``component``
     - For vector or tensor fields, the component to set (0, 1, 2 for x, y, z; for compositions,
       the component index in the order defined in the fluid model). ``-1`` (default) for scalar fields.
   * - ``setNames``
     - Required. The sets of mesh objects the value is applied to (see `Sets`_ below).
   * - ``scale``
     - The value. If ``functionName`` is given, the applied value is ``scale`` × the function value.
   * - ``functionName``
     - Optional function from ``<Functions>`` (see :ref:`FunctionManager`). A function of time makes the
       condition time-dependent; a function of space (e.g. of ``elementCenter``) makes it heterogeneous.
   * - ``initialCondition``
     - ``1``: apply once, at initialization (an initial condition). ``0`` (default): apply at every time step
       (a boundary condition).
   * - ``beginTime``, ``endTime``
     - Time interval during which a boundary condition is active.

``Dirichlet`` is an alias of ``FieldSpecification`` with the same attributes.

Specialized conditions
----------------------

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Element
     - Use
   * - ``SourceFlux``
     - Mass (or molar) rate injected into or produced from the cells of the target sets.
       The rate given by ``scale`` (× ``functionName``) is the total for the set, split equally among its cells.
       No ``fieldName`` is needed. See `Units and sign conventions`_ for the sign.
   * - ``Traction``
     - Traction applied to faces (``objectPath="faceManager"``) in mechanics. ``tractionType`` is ``vector``
       (``scale`` × ``direction``), ``normal`` (``scale`` × outward face normal, i.e. a pressure) or
       ``stress`` (``inputStress`` · face normal).
   * - ``HydrostaticEquilibrium``
     - Initializes pressure (and phase distribution) in hydrostatic equilibrium. See :ref:`EquilibriumInitialCondition`.
   * - ``Aquifer``
     - Analytical (Carter-Tracy) aquifer attached to reservoir boundary faces. See :ref:`AquiferBoundaryCondition`.
   * - ``PML``
     - Perfectly matched absorbing layer for wave propagation solvers.

Sets
----

``setNames`` refers to named sets of nodes, faces and cells. They come from three sources:

* **Built-in sets.** ``all`` contains every object. ``InternalMesh`` also creates the node sets
  ``xneg``, ``xpos``, ``yneg``, ``ypos``, ``zneg``, ``zpos`` on the six faces of the domain.
* **Geometric objects.** Each object in ``<Geometry>`` (``Box``, ``Cylinder``, ...) creates a set with the
  same name. It contains the nodes located inside the object, and the faces and cells *all* of whose nodes are inside.
  Make boxes slightly larger than the region they should capture, as in the example above.
* **Imported sets.** With ``VTKMesh``, node sets stored in the file can be imported with ``nodesetNames``,
  and each cell block / region of the file can be targeted through ``objectPath``.

Units and sign conventions
--------------------------

GEOS uses SI units throughout (Pa, m, s, kg, K); see the units table in :ref:`InputFiles`.

* **Stress and strain: compression is negative.** A compressive in-situ stress or a confining pressure
  applied as a stress must be entered with a negative value.
* **Traction:** with ``tractionType="normal"`` the traction is ``scale`` times the *outward* face normal,
  so a negative ``scale`` pushes on the boundary (compression), as a fluid pressure would.
* **SourceFlux: negative values inject, positive values produce.** For example, ``scale="-0.03"`` injects
  0.03 kg/s (or mol/s if the flow solver uses moles) into the target set.

.. toctree::
   :maxdepth: 1

   EquilibriumInitialCondition

   AquiferBoundaryCondition
