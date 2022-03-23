.. _Meshes:

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
Internal Mesh Generation
************************

Basic Example
=================

The Internal Mesh Generator allows one to quickly build simple cartesian grids and divide
them into several regions.  The following attributes are supported in the input block for InternalMesh:

.. include:: /coreComponents/schema/docs/InternalMesh.rst


The following is an example XML ``<mesh>`` block, which will generate a vertical beam with two ``CellBlocks`` (one in red and one in blue in the following picture).

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

.. image:: ../../../coreComponents/mesh/docs/beam.png


.. _Mesh_bias:

Mesh Bias
===========

The internal mesh generator is capable of producing meshes with element sizes that vary smoothly over space.
This is achieved by specifying ``xBias``, ``yBias``, and/or ``zBias`` fields.
(Note: if present, the length of these must match ``nx``, ``ny``, and ``nz``, respectively, and each individual value must be in the range (-1, 1).)

For a given element block, the average element size will be

.. math::
   dx_{average}[i] = \frac{xCoords[i+1]-xCoords[i]}{nx[i]},

the element on the left-most side of the block will have size

.. math::
   dx_{left}[i] = (1 + xBias[i]) \cdot dx_{average}[i],

and the element on the right-most side will have size

.. math::
   dx_{right}[i] = (1 - xBias[i]) \cdot dx_{average}[i].


The following are the two most common scenarios that occur while designing a mesh with bias:

1. The size of the block and the element size on an adjacent region are known.  Assuming that we are to the left of the target block, the appropriate bias would be:

.. math::
   xBias[i] = 1 - \frac{nx[i] \cdot dx_{left}[i+1]}{xCoords[i+1]-xCoords[i]}

2. The bias of the block and the element size on an adjacent region are known.  Again, assuming that we are to the left of the target block, the appropriate size for the block would be:

.. math::
   xCoords[i+1]-xCoords[i] = \frac{nx[i] \cdot dx_{left}[i+1]}{1 - xBias[i]}


The following is an example of a mesh block along each dimension, and an image showing the corresponding mesh.  Note that there is a core region of elements with zero bias, and that the transitions between element blocks are smooth.

.. literalinclude:: ../../../../inputFiles/solidMechanics/sedov_with_bias.xml
  :language: xml
  :start-after: <!-- SPHINX_MESH_BIAS -->
  :end-before: <!-- SPHINX_MESH_BIAS_END -->

.. image:: ../../../coreComponents/mesh/docs/mesh_with_bias.png


Advanced Cell Block Specification
==================================
It's possible to generate more complex ``CellBlock`` using the ``InternalMeshGenerator``.
For instance, the staircase example is a model which is often used in GEOSX as an integrated
test. It defines ``CellBlocks`` in the three directions to generate a staircase-like model
with the following code.

.. code-block:: xml

  <Mesh>
    <InternalMesh name="mesh1"
                  elementTypes="{C3D8}"
                  xCoords="{0, 5, 10}"
                  yCoords="{0, 5, 10}"
                  zCoords="{0, 2.5, 5, 7.5, 10}"
                  nx="{5, 5}"
                  ny="{5, 5}"
                  nz="{3, 3, 3, 3}"
                  cellBlockNames="{b00,b01,b02,b03,b04,b05,b06,b07,b08,b09,b10,b11,b12,b13,b14,b15}"/>
  </Mesh>

  <ElementRegions>
     <CellElementRegion name="Channel"
                    cellBlocks="{b08,b00,b01,b05,b06,b14,b15,b11}"
                    materialList="{fluid1, rock, relperm}"/>
     <CellElementRegion name="Barrier"
                    cellBlocks="{b04,b12,b13,b09,b10,b02,b03,b07}"
                    materialList="{}"/>
  </ElementRegions>

Thus, the generated mesh will be :

.. figure:: ../../../coreComponents/mesh/docs/staircase.svg
   :align: center
   :width: 500

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
- The ECLIPSE file formats (.egrid, .grdecl)

The supported mesh elements for volume elements consist of the following:

- 4 nodes tetrahedra,
- 5 nodes pyramids,
- 6 nodes wedges,
- 8 nodes hexahedra,

The mesh can be divided in several regions.
These regions are intended to support different physics
or to define different constitutive properties.

- For the GMSH file format, the regions are defined using the `physical entity names`_
  provided by GMSH.
- For the ECLIPSE file formats, the regions have to be first defined using the ECLIPSE software.

.. _ImportingExternalMesh:

Importing the Mesh
==================

Importing regions
*****************

Several blocks are involved to import an external mesh into GEOSX, defined in the XML input file.
These are the ``<Mesh>`` block and the ``<ElementRegions>`` block.

The mesh block has the following syntax:

.. code-block:: xml

  <Mesh>
    <PAMELAMesh name="MyMeshName"
                file="/path/to/the/mesh/file.msh"/>
  </Mesh>

We advise users to use absolute path to the mesh file.

GEOSX uses ``ElementRegions`` to support different physics
or to define different constitutive properties.
An ``ElementRegion`` is defined as a set of ``CellBlocks``.
A ``CellBlock`` is an ensemble of elements with the same element geometry.

.. image:: mesh.svg

In the example presented above, the mesh is is composed of two regions (*Top* and *Bot*).
Each region contains 3 ``CellBlocks``.

The ``ElementRegions`` are defined as below :

.. code-block:: xml

  <ElementRegions>
    <ElementRegion name="Top" cellBlocks="Top_HEX Top_WEDGE Top_TETRA" materialList="water rock"/>
    <ElementRegion name="Bot" cellBlocks="Bot_HEX Bot_WEDGE Bot_TETRA" materialList="water rock"/>
  </ElementRegions>

You have to use the following syntax to declare your ``CellBlocks`` :

.. code-block:: none

  nameOfTheRegionWithinTheMesh_typeOfTheElement

The keywords for the element types are :

- TETRA
- WEDGE
- PYR
- HEX

If the regions are not named in the file (it happens with all the eclipse grids and several GMSH mesh
files), the name of the region is ``DEFAULT``, e.g:

.. code-block:: xml

  <ElementRegions>
    <ElementRegion name="Default" cellBlocks="DEFAULT_HEX" materialList="water rock"/>
  </ElementRegions>

Using the gmsh file format, regions can be easily named
as a preprocessed step using the gmsh software of directly editing the file following the syntax
defined in the documentation_.

An example of a gmsh file with all the physical regions defined is used in :ref:`TutorialFieldCase`.

Importing surfaces
******************

Surfaces are imported throught point sets in GEOSX. This feature is supported using only the gmsh file format.
In the same way than the regions, the surfaces of interests can be defined using the `physical entity names`_.
The surfaces are automatically import in GEOSX if they exist in the gmsh file.
Within GEOSX, the point set will have the same name than the one given in the file. This name can be used
again to impose boundary condition. For instance, if a surface is named "Bottom" and the user wants to
impose a Dirichlet boundary condition of 0 on it, it can be easily done using this syntax.

.. code-block:: xml

  <FieldSpecification
    name="zconstraint"
    objectPath="nodeManager"
    fieldName="Velocity"
    component="2"
    scale="0.0"
    setNames="{ Bottom }"/>

The name of the surface of interest appears under the keyword ``setNames``. Again, an example of a gmsh file
with the surfaces fully defined is available within :ref:`TutorialFieldCase`.

.. _PAMELA: https://github.com/GEOSX/PAMELA
.. _GMSH: http://gmsh.info
.. _documentation: https://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format-version-2-_0028Legacy_0029
.. _`physical entity names`: https://gmsh.info/doc/texinfo/gmsh.html#Elementary-entities-vs-physical-groups
