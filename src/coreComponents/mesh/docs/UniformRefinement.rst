.. _UniformMeshRefinement:

Uniform refinement of a VTK mesh
===============================

Set ``uniformRefinement`` on ``VTKMesh`` to refine a coarse mesh before GEOS
builds connectivity, ghosts, numerical methods and degrees of freedom::

  <Mesh>
    <VTKMesh
      name="mesh"
      file="coarse.vtu"
      scatterMethod="rcb"
      partitionRefinement="1"
      uniformRefinement="2" />
  </Mesh>

The level count is a nonnegative integer. Two levels subdivide the first level's
children again. The default is zero, which preserves the existing import path.
Negative and fractional values are input errors.

``partitionRefinement`` controls the existing coarse graph partitioning.
``uniformRefinement`` runs after scatter, optional graph partitioning and
redistribution of surfaces and fractures. Each volume child stays on its
parent's rank. Fine connectivity is constructed locally and shared point IDs are
reconciled with the relevant ranks. Refinement does not gather the fine mesh or
repartition its volume children.

Supported cells
---------------

.. list-table:: Children in one level
   :header-rows: 1
   :widths: 35 65

   * - Coarse volume
     - Children
   * - Tetrahedron
     - 8 tetrahedra
   * - Hexahedron or voxel
     - 8 hexahedra
   * - Wedge
     - 8 wedges
   * - Pyramid
     - 6 pyramids and 4 tetrahedra
   * - Polygonal prism with N sides, N = 5 through 11
     - 2N hexahedra

The recognized shapes can also use VTK polyhedron encoding. Cells must have
valid supported topology and positive, admissible geometry. Arbitrary polyhedra,
high-order cells and curved-surface projection require other meshing methods.
Main-mesh line and vertex cells are diagnosed as unsupported at positive levels;
separate well generators keep their existing role.

The coarse volume mesh must be conforming. A quad opposite two triangles is a
nonmatching interface, rather than a common quad face. Shared triangle faces
produce four triangles; quads produce four quads; polygonal prism caps produce
N quads. The surface and volume templates use the same edge and face points.
Coincident points with distinct IDs retain their distinct topology.

Choose numerical methods that support the generated child types. In particular,
pyramid refinement introduces tetrahedra and prism refinement introduces hexes.

Regions and fields
------------------

Existing ``cellBlocks`` selections refer to the original coarse source blocks.
Their descendants stay in the selected material region. For example, selecting
``1_pyramids`` includes ``1_pyramids__refined_tetrahedra``; those tetrahedra remain
separate from cells descended from an original ``1_tetrahedra`` block. Attribute
and wildcard selections retain their coarse meaning. All descendants must still
belong to exactly one ``CellElementRegion``.

Ordinary cell values are intensive and copy from parent to child. Total quantities
require an explicit extensive transfer policy, which uses measured child-volume
or child-area fractions. The internal transfer-policy API supplies that policy;
there is no XML attribute for selecting extensive arrays in this interface.
Floating point data on points uses the geometry's affine interpolation. Integral
point labels require equal supporting values; differing labels cause an input
error. Names, types, components and active VTK attribute roles are preserved.

For unsigned node-set masks named by ``nodesetNames``, a new point joins the set
when every defining support corner is a member. Explicit surface marker labels
copy to each surface child. Geometric sets such as ``Box`` evaluate the final
GEOS mesh in the normal initialization phase.

Original main point IDs survive refinement. New points receive exact integral
IDs and shared copies agree. Active child cell IDs are allocated in new ranges;
root and immediate-parent lineage is stored separately. IDs, ghost flags,
extraction provenance and ``collocated_nodes`` use specialized handling, rather
than ordinary numerical interpolation.

Fractures and restart
---------------------

Named face blocks refine their actual incident volume-face traces. New
``collocated_nodes`` buckets reference the new main-mesh points on each actual
side, including distinct coincident sides. An endpoint bucket is not combined
with unrelated endpoint IDs to invent new edges. Auxiliary point IDs have a
separate namespace per block; auxiliary cell ranges are disjoint across blocks
before the normal importer applies its main-mesh offset.

Use the original coarse input and the same level count when restarting. GEOS
restores the saved fine state through its normal checkpoint path. A changed
``uniformRefinement`` count is rejected, and already refined input containing
refinement metadata is diagnosed before another positive refinement pass.

Compatibility and scaling
-------------------------

Positive refinement currently rejects ``structuredIndexAttribute``: a parent's
logical IJK cannot describe all its children. Coordinate transforms must have
finite, nonzero scales and preserve orientation. Refinement does not relax the
coarse partitioner's restrictions, including nonempty coarse volume partitions
and supported fracture/partitioner combinations.

Work and storage grow quickly with the level count. Hex and tetrahedron cell
counts multiply by eight each level; pyramids have mixed descendants. Checked
counts and ID ranges reject representational overflow. Sharing discovery and
full-face validation operate on coarse topology once. Later levels retain shared
interface metadata and use cached neighbors. The existing root memory for coarse
input reading remains part of the baseline.

See :ref:`uniform-refinement-scaling` for component measurements and the
one-, two- and four-node Dane test workflow. Component timings do not measure
full solver initialization or fracture/ghost/DOF costs.
