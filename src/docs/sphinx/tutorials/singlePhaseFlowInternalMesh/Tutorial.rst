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

  - the structure of the XML input files used by GEOSX,
  - how to run GEOSX on a simple case requiring no external input files,
  - the basic syntax of a solver block for single-phase problems,
  - how to export and visualize results.


**Input file**

This tutorial uses no external input files and everything required is contained within the GEOSX input file.
The xml input file for this test case is located at:

.. code-block:: console

  src/coreComponents/physicsSolvers/fluidFlow/integratedTests/singlePhaseFlow/3D_10x10x10_compressible.xml


------------------------------------
GEOSX input files
------------------------------------

GEOSX runs by reading user input information from one or multiple XML files.
For instance, if everything we need to run is contained in a file called ``my_input.xml``,
GEOSX runs this file by executing:

.. code-block:: console

  /your/path/to/GEOSX -i my_input.xml
  
The ``-i`` flag indicates the path to the XML input file.


XML files store information in a tree-like structure using nested blocks of information called *elements*.
In GEOSX, the root of the tree structure is an element called *Problem* (it defines the problem we wish to solve).
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
 #. :ref:`Functions and Partition <Functions_tag_single_phase_internal_mesh>`
 #. :ref:`Outputs <Outputs_tag_single_phase_internal_mesh>`


In addition to the data required to solve the problem we wish to address,
it is a best practice to start an XML files with a convention (called a *schema*)
to specify the types of writing conventions used in the XML file.

.. literalinclude:: ../../../../coreComponents/physicsSolvers/fluidFlow/integratedTests/singlePhaseFlow/3D_10x10x10_compressible.xml
  :language: xml
  :start-after: <!-- SPHINX_TUT_INT_HEX_PROBLEM_OPEN -->
  :end-before: <!-- SPHINX_TUT_INT_HEX_PROBLEM_OPEN_END -->

The attributes ``xmlns:xsi`` and ``xsi:noNamespaceSchemaLocation`` are used to define the file format and schema.
While optional, they may be used to configure various xml validation tools.


.. _Solver_tag_single_phase_internal_mesh:

Defining a solver
-----------------

GEOSX is a multi-physics tool. The solution to each physical problem
(diffusion, convection, deformation, etc.) is found using one or many numerical solvers and
the **Solvers** tag is used to list and parameterize them.
Several solvers can be defined in the input file,
and different combinations of solvers can be applied
in different regions of the mesh at different moments of the simulation.


Here, to keep things simple, we use one type of solver in the entire domain and
for the entire duration of the simulation. This solver is going to perform a single-phase flow solve.
In GEOSX, such a solver is identified by a **SinglePhaseFVM** element (part of a family of cell-centered single-phase finite volume methods).


The XML block used to define this single-phase finite volume solver is shown here:

.. literalinclude:: ../../../../coreComponents/physicsSolvers/fluidFlow/integratedTests/singlePhaseFlow/3D_10x10x10_compressible.xml
  :language: xml
  :start-after: <!-- SPHINX_TUT_INT_HEX_SOLVERS -->
  :end-before: <!-- SPHINX_TUT_INT_HEX_SOLVERS_END -->

Each type of solver has a specific set of parameters that are required and
some that are optional (usually with sensible default values).
Here, for instance, we see that our solver is registered with a user-chosen name (``SinglePhaseFlow``).
From now on, this unique name will be used in the parameter file and in the code to point to this instance of a SinglePhaseFVM solver.
We also set a solver-specific level of on-console logging (set to 1 here).

For solvers of the ``SinglePhaseFVM`` family, we must specify a discretization scheme.
Here, we use a Two-Point Flux Approximation (TPFA) finite volume discretization scheme, a typical discretization scheme
used in cartesian cell-centered finite volume methods.

We have also specified a collection of fluids, rocks, and
target regions of the mesh on which this solver will be applied (``Region2``).
The curly brackets used here are necessary, even if the collection contains a single value.

Finally, note that other XML elements can be nested inside the ``Solvers`` element.
Here, for instance, we specify values for numerical tolerances
(the solver has converged when numerical residuals are smaller than these tolerances)
and for the maximum number of iterations allowed to reach convergence.




.. _Mesh_tag_single_phase_internal_mesh:

Specifying a computational Mesh
----------------------------------

The single-phase flow solver in GEOSX uses cell-centered finite volume computations.
We must thus a define a mesh (or grid) as a numerical support to work on.
The **Mesh** element allows users to specify this support.

There are two approaches to specifying meshes in GEOSX: internal or external.
The external approach consists of importing mesh files created outside of GEOSX, such as a
corner-point grids or generic unstructured grids.
The internal approach consists of using a tool in GEOSX called the internal mesh generator.
The internal mesh generator is useful to create grids on-the-fly
from a small number of geometric parameters, without external file imports.

In this tutorial, to keep things as simple and self-contained as possible,
we use GEOSX's internal mesh generator, parameterized in the **InternalMesh** element.

.. literalinclude:: ../../../../coreComponents/physicsSolvers/fluidFlow/integratedTests/singlePhaseFlow/3D_10x10x10_compressible.xml
  :language: xml
  :start-after: <!-- SPHINX_TUT_INT_HEX_MESH -->
  :end-before: <!-- SPHINX_TUT_INT_HEX_MESH_END -->


Here, we create a mesh registered as ``mesh1`` (this registration name is chosen by the user),
consisting of elements of type ``C3D8`` (this is a code from the usual finite
element nomenclature for a general purpose linear 8-node brick element).
Elements in the x-direction go from x=0 to x=10. The x-dimension is divided into nx=10 elements.
The same is true for the y-dimension and the z-dimension.

We therefore have a cube of 10x10x10 elements with a bounding box defined by corner coordinates (0,0,0) and (10,10,10).


.. image:: cube_mesh_10x10x10.png


.. _Geometry_tag_single_phase_internal_mesh:

Geometry tag
-----------------

The **Geometry** tag is useful to point to specific parts of a mesh and assign properties to them.
Here, for instance, we use two **Box** elements to specify where our source and sink pressure terms are located.
We want the source to be all elements along the x=0 face of the domain, and the sink to be all the elements at x=10.

Note that for an element to be considered **inside** a geometric region, it needs to have all vertices inside the region.
This explains why we need to extend the geometry limits to 0.01 beyond the minimum and maximum coordinates, to be sure to encompass the entire elements.

.. literalinclude:: ../../../../coreComponents/physicsSolvers/fluidFlow/integratedTests/singlePhaseFlow/3D_10x10x10_compressible.xml
  :language: xml
  :start-after: <!-- SPHINX_TUT_INT_HEX_GEOMETRY -->
  :end-before: <!-- SPHINX_TUT_INT_HEX_GEOMETRY_END -->

There are several methods to achieve similar conditions (Dirichlet boundary condition on faces, etc.).
The **Box** defined here is one of the simplest approach.
Boxes defined here are named objects, and will be registered and used using their names (``source`` and ``sink``).


.. image:: cube_initial.png


.. _Events_tag_single_phase_internal_mesh:

Events tag
---------------
The Event tag includes the final time of our simulation under ``maxTime`` node. Under *PeriodicEvent* embededd tags, we can set:

 #. which solver has to be called (among the child tag defined under the above mentinoned *Solver* tag) with its initial time step defined as the ``forceDt`` node value.
 #. under which ``timeFrequency`` will we need to output results (targeting the settings defined under some child tag of the below explained *Output* tag).

.. literalinclude:: ../../../../coreComponents/physicsSolvers/fluidFlow/integratedTests/singlePhaseFlow/3D_10x10x10_compressible.xml
  :language: xml
  :start-after: <!-- SPHINX_TUT_INT_HEX_EVENTS -->
  :end-before: <!-- SPHINX_TUT_INT_HEX_EVENTS_END -->

.. _NumericalMethods_tag_single_phase_internal_mesh:

NumericalMethods tag
------------------------

The two-point flux approximation, which was first introduced under the *Solver>SinglePhaseFlow* child tag as the value of ``discretization`` node, is defined here.

.. literalinclude:: ../../../../coreComponents/physicsSolvers/fluidFlow/integratedTests/singlePhaseFlow/3D_10x10x10_compressible.xml
  :language: xml
  :start-after: <!-- SPHINX_TUT_INT_HEX_NUM_METHODS -->
  :end-before: <!-- SPHINX_TUT_INT_HEX_NUM_METHODS_END -->

Here the ``boundaryFieldName`` node specifies that for Dirichlet boundary conditions the face located value is considered. The ``coefficientName`` node refers to the field which has to be considered in the stencil computation.

.. _ElementRegions_tag_single_phase_internal_mesh:

Element Regions tag
---------------------

This block defines regions.
Here, the entire field is one region called ``Domain``,
and contains ``water`` and ``rock`` only.

.. literalinclude:: ../../../../coreComponents/physicsSolvers/fluidFlow/integratedTests/singlePhaseFlow/3D_10x10x10_compressible.xml
  :language: xml
  :start-after: <!-- SPHINX_TUT_INT_HEX_ELEM_REGIONS -->
  :end-before: <!-- SPHINX_TUT_INT_HEX_ELEM_REGIONS_END -->

.. _Constitutive_tag_single_phase_internal_mesh:

Constitutive tag
---------------------

The physical properties of ``water`` and ``rock`` elements can be found and set under this tag.

.. literalinclude:: ../../../../coreComponents/physicsSolvers/fluidFlow/integratedTests/singlePhaseFlow/3D_10x10x10_compressible.xml
  :language: xml
  :start-after: <!-- SPHINX_TUT_INT_HEX_CONSTITUTIVE -->
  :end-before: <!-- SPHINX_TUT_INT_HEX_CONSTITUTIVE_END -->

.. _FieldSpecifications_tag_single_phase_internal_mesh:

FieldSpecifications tag
---------------------------
Here, fields such as porosity, permeability, source and sink terms or initial field values are specified. Our test case exhibits an anisotropic homogeneous permeability which components are so that:
  - permeability in the x-direction: ``permx``, constant value of 1.0e-12 m\ :sup:`2` (100 mD), and is considered the 0\ :sup:`th` component of the ``permeability`` vector,
  - permeability in the y-direction: ``permy``, constant value of 1.0e-12 m\ :sup:`2` (100 mD),
  - a lower permeability in the z-direction: ``permz``, constant value of 1.0e-15 m\ :sup:`2` (10 mD)

The ``setNames`` node value specifies the geometric zone where the value should be applied.
These directional permeabilities are followed by all the other field initializations. Please note the change in ``component`` node value as we are dealing with a permeability diagonal tensor.
The other fields to be specified are a constant homogeneous reference porosity for the whole domain, initial pressure, and source and sink term pressures.

.. literalinclude:: ../../../../coreComponents/physicsSolvers/fluidFlow/integratedTests/singlePhaseFlow/3D_10x10x10_compressible.xml
  :language: xml
  :start-after: <!-- SPHINX_TUT_INT_HEX_FIELDS -->
  :end-before: <!-- SPHINX_TUT_INT_HEX_FIELDS_END -->

.. _Functions_tag_single_phase_internal_mesh:

Here we leave ``Functions`` and ``Partition`` tags unspecified as the description of their use are detailed in other tutorials.

.. literalinclude:: ../../../../coreComponents/physicsSolvers/fluidFlow/integratedTests/singlePhaseFlow/3D_10x10x10_compressible.xml
  :language: xml
  :start-after: <!-- SPHINX_TUT_INT_HEX_BLANKS -->
  :end-before: <!-- SPHINX_TUT_INT_HEX_BLANKS_END -->

.. _Outputs_tag_single_phase_internal_mesh:

Outputs tag
----------------
In order to get the results from simulation written to file, we specify the output path:

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
