
GEOS Mesh Tools
--------------------------

The `geosx_mesh_tools` python package includes tools for converting meshes from common formats (abaqus, etc.) to those that can be read by GEOS (gmsh, vtk).
See :ref:`PythonToolsSetup` for details on setup instructions, and :ref:`ExternalMeshUsage` for a detailed description of how to use external meshes in GEOS.
The available console scripts for this package and its API are described below.


convert_abaqus
^^^^^^^^^^^^^^

Compile an xml file with advanced features into a single file that can be read by GEOS.

.. argparse::
   :module: geosx_mesh_tools.main
   :func: build_abaqus_converter_input_parser
   :prog: convert_abaqus


.. note::
    For vtk format meshes, the user also needs to determine the region ID numbers and names of nodesets to import into GEOS.
    The following shows how these could look in an input XML file for a mesh with three regions (*REGIONA*, *REGIONB*, and *REGIONC*) and six nodesets (*xneg*, *xpos*, *yneg*, *ypos*, *zneg*, and *zpos*):


.. code-block:: xml

    <Problem>
      <Mesh>
        <VTKMesh
          name="external_mesh"
          file="mesh.vtu"
          regionAttribute="REGIONA-REGIONB-REGIONC"
          nodesetNames="{ xneg, xpos, yneg, ypos, zneg, zpos }"/>
      </Mesh>

      <ElementRegions>
        <CellElementRegion
          name="ALL"
          cellBlocks="{ 0_tetrahedra, 1_tetrahedra, 2_tetrahedra }"
          materialList="{ water, porousRock }"
          meshBody="external_mesh"/>
      </ElementRegions>
    </Problem>


API
^^^

.. automodule:: geosx_mesh_tools.abaqus_converter
    :members:
