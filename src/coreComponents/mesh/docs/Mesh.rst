============
Meshes
============

The purpose of this document is to explain how users and developers interact with mesh data.
This section describes how meshes are handled and stored in GEOSX.

There are two possible methods for generating a mesh:
either by using GEOSX's internal mesh generator (for Cartesian meshes only),
or by importing meshes from various common mesh file formats.
This latter options allows one to work with more complex geometries,
such as unstructured meshes comprised of a variety of element types (polyhedral elements).

************************
Mesh Data Structure
************************

GEOSX uses a hierarchical class structure to store the mesh.
To illustrate it, we show a model with two regions (Top and Bottom).

.. image:: ../../../coreComponents/mesh/docs/model.png

This model can be meshed with different types of polyhedra.
Here, the model is meshed with pyramids,
tetrahedra, hexahedra and wedges.

.. image:: ../../../coreComponents/mesh/docs/mesh_multi.png

This mesh can then be split across different computational nodes.
The resulting split consists of several ``DomainPartition`` of the mesh (here, two).
These partitions are not necessarily identical to the two regions (Top and Bottom).

.. image:: ../../../coreComponents/mesh/docs/mesh_domain.png

Each ``DomainPartition`` is handled by the ``DomainPartition`` class in GEOSX.
Each MPI process has its own instance of a ``DomainPartition``.

A ``DomainPartition`` can contain several ``MeshBody``.

A ``MeshBody`` is a partition of the mesh assigned to a specific set of physical laws.
For instance, we could imagine a ``MeshBody`` where we only want to solve for flow 
and not for elasticity, while in another area, we would have a ``MeshBody`` where both solvers are used.

A ``MeshBody`` can contain several ``MeshLevel``.
``MeshLevel`` objects are used to define levels in multi level computations.


For clarity here, we assume that there is only one ``MeshBody`` and one ``MeshLevel``.

A ``MeshLevel`` contains several managers that handle the mesh data structure.

- ``NodeManager`` handles the nodes,
- ``EdgeManager`` handles the edges,
- ``FaceManager`` handles the facets,
- ``ElementRegionManager`` handles the different regions and the polyhedra that compose them.

NodeManager
===============

The ``NodeManager`` contains information on all nodes (vertices)
of the ``MeshLevel`` and of the
``DomainPartition`` it belongs to.
Its size is equal to the number of nodes in this ``DomainPartition``/``MeshLevel``.

EdgeManager
===============

In GEOSX, an edge is a segment between two nodes.
The ``EdgeManager`` contains information on all edges
of the ``MeshLevel`` and of the ``DomainPartition`` it belongs to.
Its size is equal to the number of edges in this ``DomainPartition``/``MeshLevel``.

The following picture shows the edges of the ``DomainPartition 1`` in the mesh:

.. image:: ../../../coreComponents/mesh/docs/edges_domain1.png

FaceManager
===============

In GEOSX, a face is the interface between two polyhedra. 
The ``FaceManagers`` contains informations on all faces of the ``MeshLevel`` and of the
``DomainPartitions`` it belongs to. 
Its size is equal to the number of faces in this ``DomainPartition``/``MeshLevel``.

ElementRegionManager
========================

The ``ElementRegionManager`` handles all polyhedral elements of the ``DomainPartition``/``MeshLevel`` it belongs to.
An ``ElementRegion`` in GEOSX is thus a collection of polyhedral elements managed by an ``ElementRegionManager``. 
In the example above, the ``ElementRegionManager`` of one ``DomainPartition``/``MeshLevel`` manages two ``ElementRegion`` instances: one corresponding to the Bottom region, and one corresponding to the Top region.

The element geometry information is stored in the ``CellElementSubRegion``.
An ``ElementRegion`` can contain several ``CellElementSubRegion``.
There is one such ``CellElementSubRegion`` for each element type.
In our example, there are four distinct elements types (hexahedra, tetrahedra, wedges and pyramids).
As a consequence, our ``ElementRegion`` will contain four different ``CellElementSubRegion`` instances:
one for all hexahedra, one for all tetrahedra, one for all wedges, one for all pyramids.

Ghosting structure
==================

To ease the communication between ``DomainPartition`` objects across MPI processes,
GEOSX computes ghost elements.
Ghost elements provide an overlap between two adjacents ``DomainPartition`` objects.
From now on, we will distinguish between *owned* elements (that belong to the domain)
and *ghost* elements (than belong to the neighboring domain).

.. image:: ../../../coreComponents/mesh/docs/split.png

.. warning::
  Asking for the size of the ``NodeManager``, ``EdgeManager``, ``FaceManager`` or a ``CellElementSubRegion``
  will return the number of owned elements plus the number of ghost elements.

The complete mesh data structure is shown in the next picture.

.. image:: ../../../coreComponents/mesh/docs/diag.png


************************
Internal Mesh Generation
************************

The Internal Mesh Generator allows one to quickly build simple cartesian grids and divide
them into several regions.
The following is an example XML ``<Mesh>`` block:

.. code-block:: xml

  <Mesh>
    <InternalMesh name="mesh"
                  elementTypes="C3D8"
                  xCoords="0, 1"
                  yCoords="0, 1"
                  zCoords="0, 2, 6"
                  nx="1"
                  ny="1"
                  nz="2, 4"
                  cellBlockNames="cb1 cb2"/>
  </Mesh>

- ``name`` the name of the mesh body
- ``elementTypes`` the type of the elements that will be generated.
- ``xCoord`` List of ``x`` coordinates of the boundaries of the ``CellBlocks``
- ``yCoord`` List of ``y`` coordinates of the boundaries of the ``CellBlocks``
- ``zCoord`` List of ``z`` coordinates of the boundaries of the ``CellBlocks``
- ``nx`` List containing the number of cells in ``x`` direction within the ``CellBlocks``
- ``ny`` List containing the number of cells in ``y`` direction within the ``CellBlocks``
- ``nz`` List containing the number of cells in ``z`` direction within the ``CellBlocks``
- ``cellBlockNames`` List containing the names of the ``CellBlocks``

The previous sample of XML file will generate a vertical beam with two ``CellBlocks``
(one in red and one in blue in the following picture)

.. image:: ../../../coreComponents/mesh/docs/beam.png


**************************
Using an External Mesh
**************************

Supported Formats
=================

GEOSX provides features to run simulations on unstructured meshes.
It uses PAMELA_ to read the external meshes and its API to write
it into the GEOSX mesh data structure.

The supported mesh format are:

- The GMSH_ file format (.msh v2).
- The MEDIT_ file format (.mesh)
- The ECLIPSE file formats (.egrid, .grdecl)

The supported mesh elements for volume elements consist of the following:

- 4 nodes tetrahedra,
- 5 nodes pyramids,
- 6 nodes wedges,
- 8 nodes hexahedra,

The mesh can be divided in several regions.
These regions are intended to support different physics
or to define different constitutive properties.

- For the GMSH file format, the regions are defined using the `elementary geometrical tags`_
  provided by GMSH.
- For the MEDIT file format, the regions are defined using the tag of the element.
- For the ECLIPSE file formats, the regions have to be first defined using the ECLIPSE software.

Importing the Mesh
==================

Several blocks are involved to import an external mesh into GEOSX, defined in the XML input file.
These are the ``<Mesh>`` block and the ``<ElementRegions>`` block.

The mesh block has the following syntax:

.. code-block:: xml

  <Mesh>
    <PAMELAMeshGenerator name="MyMeshName"
                         file="/path/to/the/mesh/file.msh"/>
  </Mesh>

We advise users to use absolute path to the mesh file.

GEOSX uses ``ElementRegions`` to support different physics
or to define different constitutive properties.
An ``ElementRegion`` is defined as a set of ``CellBlocks``.
A ``CellBlock`` is an ensemble of elements with the same element geometry.

.. image:: mesh.svg

In the example presented above, the mesh is is composed of two regions (*Region 0* and *Region 1*).
Each region contains 3 ``CellBlocks``.

The ``ElementRegions`` are defined as below :

.. code-block:: xml

  <ElementRegions>
    <ElementRegion name="Top" cellBlocks="0_HEX 0_WEDGE 0_TETRA" materialList="water rock"/>
    <ElementRegion name="Bot" cellBlocks="1_HEX 1_WEDGE 1_TETRA" materialList="water rock"/>
  </ElementRegions>

You have to use the following syntax to declare your ``CellBlocks`` :

.. code-block:: none

  indexOfTheRegionWithinTheMesh_typeOfTheElement

The keywords for the element types are :

- TETRA
- WEDGE
- PYR
- HEX

.. _PAMELA: https://github.com/GEOSX/PAMELA
.. _GMSH: http://gmsh.info
.. _MEDIT: https://people.sc.fsu.edu/~jburkardt/data/medit/medit.html
.. _`elementary geometrical tags`: http://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format-version-2
