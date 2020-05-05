.. _TutorialSinglePhaseFlowWithInternalMesh:

#####################################################
First steps with GEOSX: the single-phase flow solver
#####################################################

**Context**

In this tutorial, we use a single-phase flow solver (see :ref:`SinglePhaseFlow`)
from GEOSX to solve for pressure propagation on a simple discretized 10x10x10 cube mesh.
A pressure source term will be set on one side of the cube, along a face, and
a sink term will be set on the opposite face of the cube.

**Objectives**

At the end of this tutorial you will know:

  - the basic structure of XML input files used by GEOSX,
  - how to run GEOSX on a simple case requiring no external input files,
  - the basic syntax of a solver block for single-phase problems,
  - how to control output and visualize results.


**Input file**

This tutorial uses no external input files and everything required is
contained within a single GEOSX input file.
The xml input file for this test case is located at:

.. code-block:: console

  src/coreComponents/physicsSolvers/fluidFlow/integratedTests/singlePhaseFlow/3D_10x10x10_compressible.xml


------------------------------------
GEOSX input files
------------------------------------

GEOSX runs by reading user input information from one or several XML files.
For instance, if everything we need to run is contained in a file called ``my_input.xml``,
GEOSX runs this file by executing:

.. code-block:: console

  /your/path/to/GEOSX -i /your/path/to/my_input.xml
  
The ``-i`` flag indicates the path to the XML input file.


XML files store information in a tree-like structure using nested blocks of information called *elements*.
In GEOSX, the root of this tree structure is an element called *Problem*. It defines the problem we wish to solve.
All elements in an XML file have names (commonly called *tags*) and properties (commonly called *attributes*).
A typical GEOSX input file contains the following XML tags:


 #. :ref:`Solver <Solver_tag_single_phase_internal_mesh>`
 #. :ref:`Mesh <Mesh_tag_single_phase_internal_mesh>`
 #. :ref:`Geometry <Geometry_tag_single_phase_internal_mesh>`
 #. :ref:`Events <Events_tag_single_phase_internal_mesh>`
 #. :ref:`NumericalMethods <NumericalMethods_tag_single_phase_internal_mesh>`
 #. :ref:`ElementRegions <ElementRegions_tag_single_phase_internal_mesh>`
 #. :ref:`Constitutive <Constitutive_tag_single_phase_internal_mesh>`
 #. :ref:`FieldSpecifications <FieldSpecifications_tag_single_phase_internal_mesh>`
 #. :ref:`Functions<Functions_tag_single_phase_internal_mesh>`
 #. :ref:`Outputs <Outputs_tag_single_phase_internal_mesh>`


In addition to the data required to solve the problem,
it is a best practice to start an XML files with an optional file (called a *schema*)
that specifies the writing conventions used in the XML file.

.. literalinclude:: ../../../../coreComponents/physicsSolvers/fluidFlow/integratedTests/singlePhaseFlow/3D_10x10x10_compressible.xml
  :language: xml
  :start-after: <!-- SPHINX_TUT_INT_HEX_PROBLEM_OPEN -->
  :end-before: <!-- SPHINX_TUT_INT_HEX_PROBLEM_OPEN_END -->

The attributes ``xmlns:xsi`` and ``xsi:noNamespaceSchemaLocation`` are used to define the file format and schema.
While optional, they may be used to configure various xml validation tools.


.. _Solver_tag_single_phase_internal_mesh:

Defining a solver
-----------------

GEOSX is a multi-physics tool. To find solution to different categories of physical problem
(diffusion, convection, deformation, etc.), GEOSX uses combinations of numerical solvers.
The XML **Solvers** tag is used to list and parameterize these solvers.
In GEOSX, different combinations of solvers can be applied
in different regions of the mesh at different moments of the simulation.


Here, to keep things simple, we use one type of solver in the entire domain and
for the entire duration of the simulation.
The solver we are specifying here is a single-phase flow solver.
In GEOSX, such a solver is identified by a **SinglePhaseFVM** element (part of a family of cell-centered single-phase finite volume methods).


The XML block used to define this single-phase finite volume solver is shown here:

.. literalinclude:: ../../../../coreComponents/physicsSolvers/fluidFlow/integratedTests/singlePhaseFlow/3D_10x10x10_compressible.xml
  :language: xml
  :start-after: <!-- SPHINX_TUT_INT_HEX_SOLVERS -->
  :end-before: <!-- SPHINX_TUT_INT_HEX_SOLVERS_END -->


Each type of solver has a specific set of parameters that are required and
some parameters that are optional. Optional values are usually set with sensible default values.

To start, we see that our solver is registered with a user-chosen name
(here ``SinglePhaseFlow``, but it could be anything).
This is a common practice in GEOSX: users need to give names to objects they define.
This is analogous to instantiating a class in C++ code.
From now on, this unique name will be used as the *handle* to this specific flow solver instance.
This handle will also be used inside the code
to point to this specific instance of a SinglePhaseFVM solver.


Then, we set a solver-specific level of on-console logging (``logLevel`` set to 1 here).
Higher values will lead to more console output and/or intermediate results saved to files.


For solvers of the ``SinglePhaseFVM`` family, we must specify a discretization scheme.
Here, we use a Two-Point Flux Approximation (TPFA) finite volume discretization scheme, a typical discretization scheme
used in cartesian cell-centered finite volume methods.

We have also specified a collection of fluids, rocks, and
target regions of the mesh on which this solver will be applied (``Region2``).
Curly brackets are used in GEOSX inputs to indicate collections of values (sets or lists).
The curly brackets used here are necessary, even if the collection contains a single value.

Finally, note that other XML elements can be nested inside the ``Solvers`` element.
Here, we use specific XML elements to set values for numerical tolerances
(the solver has converged when numerical residuals are smaller than these tolerances)
and for the maximum number of iterations allowed to reach convergence.




.. _Mesh_tag_single_phase_internal_mesh:

Specifying a computational mesh
----------------------------------

The single-phase flow solver in GEOSX uses cell-centered finite volume computations.
We must thus a define a mesh (or grid) to perform numerical calculations on.
The **Mesh** element allows users to specify this support.

There are two approaches to specifying meshes in GEOSX: internal or external.
The external approach consists of importing mesh files created outside of GEOSX, such as a
corner-point grids or generic unstructured grids.
This external approach is generally used when using real data and
geological models with complex shapes and structures.


The internal approach uses a convenient tool in GEOSX called the internal mesh generator.
The internal mesh generator creates simple geometric grids directly inside GEOSX
from a small number of parameters. It does not require any external file information.

In this tutorial, to keep things as simple and self-contained as possible,
we use the internal mesh generator. We parameterize it with the **InternalMesh** element.

.. literalinclude:: ../../../../coreComponents/physicsSolvers/fluidFlow/integratedTests/singlePhaseFlow/3D_10x10x10_compressible.xml
  :language: xml
  :start-after: <!-- SPHINX_TUT_INT_HEX_MESH -->
  :end-before: <!-- SPHINX_TUT_INT_HEX_MESH_END -->


We create a mesh registered as ``mesh1``.
Just like for solvers, this registration name is chosen by the user.

Then, we specify the collection of elements *types* that this mesh contains.
Tetrahedra, hexahedra, wedges, prisms are examples of element types.
If a mesh contains different types of elements (a hybrid mesh),
we should indicate this here by listing all unique types of elements in curly brackets.
Keeping things simple, our element collection has only one type of element: a ``C3D8`` type.
This value is a code taken from the usual finite
element nomenclature. It represents a general purpose linear 8-node brick element (linear hexahedron).


Last, we specify the spatial arrangement of the mesh elements.
The mesh defined here goes from coordinate x=0 to x=10 in the x-direction, with ``nx=10`` subdivisions along this segment.
The same is true for the y-dimension and the z-dimension.
We therefore have a cube of 10x10x10 elements with a bounding box defined by corner coordinates (0,0,0) and (10,10,10).


.. image:: cube_mesh_10x10x10.png


.. _Geometry_tag_single_phase_internal_mesh:

Geometry tag
-----------------

The **Geometry** tag is useful to define specific parts of a mesh and assign properties to them.
Here, for instance, we use two **Box** elements to specify where our source and sink pressure terms are located.
We want the source to be all elements along the x=0 face of the domain, and the sink to be all the elements at x=10.
Later in the file, we will assign a high pressure to the source box, and a low pressure to the sink box.

Note that for an element to be considered **inside** a geometric region, it needs to have all its vertices inside the region.
This explains why we need to extend the geometry limits to 0.01 beyond the coordinates of the elements, to be sure to catch all vertices.

.. literalinclude:: ../../../../coreComponents/physicsSolvers/fluidFlow/integratedTests/singlePhaseFlow/3D_10x10x10_compressible.xml
  :language: xml
  :start-after: <!-- SPHINX_TUT_INT_HEX_GEOMETRY -->
  :end-before: <!-- SPHINX_TUT_INT_HEX_GEOMETRY_END -->

There are several methods to achieve similar conditions (Dirichlet boundary condition on faces, etc.).
The **Box** defined here is one of the simplest approach.
Just like meshes and solvers, boxes are named objects and will be registered and used using their handle (``source`` and ``sink``).
We can refer to their handle later in the input file when assigning property values to them.


.. image:: cube_initial.png


.. _Events_tag_single_phase_internal_mesh:

Specifying events
------------------------

In GEOSX, we call **Events** anything that happens at a set time, or a set frequency (**PeriodicEvents**).
Events are very important and useful elements in GEOSX,
and a dedicated section just for events is necessary to gives them the treatment they deserve.


But for now, we focus on three simple types of events: the time at which we wish the simulation to end (``maxTime``),
the times at which we want the solver to perform computations,
and the times we wish to have simulation output values reported.


In GEOSX, all times are specified in **seconds**, so here ``maxTime=5000.0`` means that the simulation will run from time 0 to time 5,000 seconds.


If we focus on the two periodic events, we see :

 #. A periodic solver application: this event is registered here as ``solverApplications`` (user-defined name). With the attribute ``forceDt=20``, it forces the solver to compute results at every 20 second time intervals. We know what this event does by looking at its ``target`` attribute: here, from time 0 to ``maxTime`` and with a forced time step of 20 seconds, we instruct GEOSX to call the solver registered as ``SinglePhaseFlow``. Note the hierarchical structure of the target formulation, using '/' to indicate a specific named instance (``SinglePhaseFlow``) of an element (``Solvers``). Also note that if the solver needs to take smaller time steps than 20 seconds (for numerical convergence, for instance) it is allowed to do so. But it will have to compute results for every 20 seconds increments between time zero and ``maxTime`` regardless of possible intermediate time steps required.
 #. An output event: this event is used for reporting purposes and forces GEOSX to write out results at specific frequencies. Here, we need to see results at every 100 seconds increments. The ``targetExactTimestep=1`` flag is used to instruct GEOSX that this output event must be always be done jointly with a full application of solvers at the output time, even if the solvers were not synchronized with the outputs. In other words, with this flag set to 1, an output event will force an application of solvers, possibly in addition to the periodic events requested directly by solvers.

.. literalinclude:: ../../../../coreComponents/physicsSolvers/fluidFlow/integratedTests/singlePhaseFlow/3D_10x10x10_compressible.xml
  :language: xml
  :start-after: <!-- SPHINX_TUT_INT_HEX_EVENTS -->
  :end-before: <!-- SPHINX_TUT_INT_HEX_EVENTS_END -->


.. _NumericalMethods_tag_single_phase_internal_mesh:

Defining Numerical Methods
----------------------------------

GEOSX comes with a number of useful numerical methods.
Here, for instance, in the Solvers elements, we have specified that we use a two-point flux approximation
as discretization scheme for the finite volume single-phase solver.
To use this scheme, we need to supply more details in the **NumericalMethods** element.

.. literalinclude:: ../../../../coreComponents/physicsSolvers/fluidFlow/integratedTests/singlePhaseFlow/3D_10x10x10_compressible.xml
  :language: xml
  :start-after: <!-- SPHINX_TUT_INT_HEX_NUM_METHODS -->
  :end-before: <!-- SPHINX_TUT_INT_HEX_NUM_METHODS_END -->

Briefly, the ``fieldName`` attribute specifies which property will be used for flux computations,
the``boundaryFieldName`` attribute specifies that for Dirichlet boundary conditions,
the pressure at the element face value is used.
Last, the ``coefficientName`` attribute is used for the stencil transmissibility computations.

Note that in GEOSX, we are distinguishing solvers from numerical methods
and their parameterization are independent. We can thus solve have
multiple solvers using the same numerical scheme but with different tolerances, for instance.


.. _ElementRegions_tag_single_phase_internal_mesh:

Defining regions in the mesh
-----------------------------------

Regions are important in GEOSX to specify the material properties of elements.
The **ElementRegions** element is used here to list all the regions used in the simulation.
Here we use only a single region to represent the entire domain (named ``Region2``),
with a collection of elements containing only the ``cb1`` blocks defined in the mesh section.
We must also specify the material contained in that region (here, two materials are used: ``water`` and ``rock``; their properties will be defined next).


.. literalinclude:: ../../../../coreComponents/physicsSolvers/fluidFlow/integratedTests/singlePhaseFlow/3D_10x10x10_compressible.xml
  :language: xml
  :start-after: <!-- SPHINX_TUT_INT_HEX_ELEM_REGIONS -->
  :end-before: <!-- SPHINX_TUT_INT_HEX_ELEM_REGIONS_END -->


.. _Constitutive_tag_single_phase_internal_mesh:

Defining material properties with constitutive laws
---------------------------------------------------------------------

The **Constitutive** element allows to list all elements contained in the domain.
Here, the physical properties of the elements defined as ``water`` and ``rock`` are specified under this tag,
each with a specific type (a ``CompressibleSinglePhaseFluid`` for the water, and a ``PoreVolumeCompressibleSolid`` for the rock).

.. literalinclude:: ../../../../coreComponents/physicsSolvers/fluidFlow/integratedTests/singlePhaseFlow/3D_10x10x10_compressible.xml
  :language: xml
  :start-after: <!-- SPHINX_TUT_INT_HEX_CONSTITUTIVE -->
  :end-before: <!-- SPHINX_TUT_INT_HEX_CONSTITUTIVE_END -->

.. _FieldSpecifications_tag_single_phase_internal_mesh:

Defining properties with the FieldSpecifications
---------------------------------------------------------------------

Here, fields such as porosity, permeability, source and sink terms or initial field values are specified.
Our test case exhibits an anisotropic homogeneous permeability which components are so that:
  - permeability in the x-direction: ``permx``, constant value of 1.0e-12 m\ :sup:`2` (1 Darcy), and is considered the 0\ :sup:`th` component of the ``permeability`` vector,
  - permeability in the y-direction: ``permy``, constant value of 1.0e-12 m\ :sup:`2` (1 Darcy),
  - a lower permeability in the z-direction: ``permz``, constant value of 1.0e-15 m\ :sup:`2` (1 mD)

The ``setNames`` node value specifies the geometric zone where the value should be applied.
These directional permeabilities are followed by all the other field initializations. Please note the change in ``component`` node value as we are dealing with a permeability diagonal tensor.
The other fields to be specified are a constant homogeneous reference porosity for the whole domain, initial pressure, and source and sink term pressures.

.. literalinclude:: ../../../../coreComponents/physicsSolvers/fluidFlow/integratedTests/singlePhaseFlow/3D_10x10x10_compressible.xml
  :language: xml
  :start-after: <!-- SPHINX_TUT_INT_HEX_FIELDS -->
  :end-before: <!-- SPHINX_TUT_INT_HEX_FIELDS_END -->

.. _Functions_tag_single_phase_internal_mesh:

Here we leave ``Functions`` element empty. The description of the Functions elements is detailed in another tutorial.

.. literalinclude:: ../../../../coreComponents/physicsSolvers/fluidFlow/integratedTests/singlePhaseFlow/3D_10x10x10_compressible.xml
  :language: xml
  :start-after: <!-- SPHINX_TUT_INT_HEX_BLANKS -->
  :end-before: <!-- SPHINX_TUT_INT_HEX_BLANKS_END -->

.. _Outputs_tag_single_phase_internal_mesh:

Specifying the output formats
----------------------------------

In order to get the results from simulation written to visualization and post-processing file(s), we specify the output path:

.. literalinclude:: ../../../../coreComponents/physicsSolvers/fluidFlow/integratedTests/singlePhaseFlow/3D_10x10x10_compressible.xml
  :language: xml
  :start-after: <!-- SPHINX_TUT_INT_HEX_OUTPUTS -->
  :end-before: <!-- SPHINX_TUT_INT_HEX_OUTPUTS_END -->

And this concludes our XML file:

.. literalinclude:: ../../../../coreComponents/physicsSolvers/fluidFlow/integratedTests/singlePhaseFlow/3D_10x10x10_compressible.xml
  :language: xml
  :start-after: <!-- SPHINX_TUT_INT_HEX_PROBLEM_CLOSE -->
  :end-before: <!-- SPHINX_TUT_INT_HEX_PROBLEM_CLOSE_END -->

------------------------------------
Runnning GEOSX
------------------------------------

The command to run GEOSX is

``path/to/geosx -i path/to/this/xml_file.xml``

Note that all paths for files included in the XML file are relative to this XML file. While running GEOSX, it will log status info in the console output.

For internal mesh generation,

.. code-block:: sh

  GEOS must be configured to use Python to use parameters, symbolic math, etc. in input files
  Adding Solver of type SinglePhaseFlow, named SinglePhaseFlow
  Adding Mesh: InternalMesh, mesh1
  Adding Geometric Object: Box, source
  Adding Geometric Object: Box, sink
  Adding Event: PeriodicEvent, solverApplications
  Adding Event: PeriodicEvent, outputs
  Adding Output: Silo, siloOutput
  Adding Object CellElementRegion named Region2 from ObjectManager::Catalog.


The time iteration are then logged until the end of the simulation

.. code-block:: sh

  Running simulation
  Time: 0s, dt:20s, Cycle: 0
  Attempt: 0, Newton: 0, R = 5.6703
  Attempt: 0, Newton: 1, R = 0.000207606
  Attempt: 0, Newton: 2, R = 9.87966e-11
  Time: 20s, dt:20s, Cycle: 1
  Attempt: 0, Newton: 0, R = 0.0680544
  Attempt: 0, Newton: 1, R = 5.30163e-05
  Attempt: 0, Newton: 2, R = 5.0784e-12
  Time: 4960s, dt:20s, Cycle: 248
  Attempt: 0, Newton: 0, R = 9.33817e-07
  Time: 4980s, dt:20s, Cycle: 249
  Attempt: 0, Newton: 0, R = 9.33817e-07
  Cleaning up events

  init time = 0.043643s, run time = 4.0304s

------------------------------------
Visualization of results
------------------------------------


All results are written in a format compatible with `VisIt
<https://wci.llnl.gov/simulation/computer-codes/visit/>`_.

For instance, here are reported diagonal pressure profile from sink to source blocks with the time being increased (on the left) and the 3D plot of the transient pressure gradient to the linear solution (on the right)

.. image:: Plots.png
   :width: 400px

.. image:: IntHexMovie.mpg
   :width: 500px



All results are written in a format compatible with `VisIt
<https://wci.llnl.gov/simulation/computer-codes/visit/>`_.

For more details on the single-phase flow solvers, please see :ref:`SinglePhaseFlow`.
