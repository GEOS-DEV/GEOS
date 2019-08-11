============
Meshes
============

This section describes how meshes are handled.
We first briefly describe how meshes are stored within GEOSX, 
in order to clarify how the user can interact with mesh data.

There are then two options for generating a mesh.  GEOSX can internally
generate simple (Cartesian) meshes.  For more complex geometries, we support
a number of mesh import file formats.  This latter options allows one to work
with unstructured mesh data and a variety of element types.

************************
Mesh Data Structure
************************

TODO: explain elementRegions and cellBlocks, and then refer to more extensive developer guide documentation.

************************
Internal Mesh Generation
************************

The Internal Mesh Generator allows one to quickly build simple cartesian grids and divide
them into several regions.  The following is an example XML ``<mesh>`` block:

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
- ``ny`` List containing the number of cells in ``x`` direction within the ``CellBlocks``
- ``nz`` List containing the number of cells in ``x`` direction within the ``CellBlocks``
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

The supported mesh elements are, for volume elements:

- 4 nodes tetrahedra
- 5 nodes pyramids
- 6 nodes wedges
- 8 nodes hexahedra

The mesh can be divided in several regions.
These regions are intended
to support different physics or to define different constitutive properties.

- For the GMSH file format, the regions are defined using the `elementary geometrical tags`_
  provided by GMSH
- For the MEDIT file format, the regions are defined using the tag of the element
- For the ECLIPSE file formats, the regions have to be first defined using the ECLIPSE software

Importing the Mesh
==================

Several blocks are involved to import an external mesh into GEOSX, defined in the XML input file.
These are the ``<Mesh>`` block and the ``<ElementRegions>`` block.

The mesh block has the following syntax.

.. code-block:: xml

  <Mesh>
    <PAMELAMeshGenerator name="MyMeshName"
                         file="/path/to/the/mesh/file.msh"/>
  </Mesh>

We strongly recommand to use absolute path to the mesh file.

GEOSX uses ``ElementRegions`` to support different physics, or to define different constitutive properties.
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
