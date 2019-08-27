###############################################################################
DoF Manager
###############################################################################

This will contains a description of the DoF manager in GEOSX.

Brief description
========================

The main aim of the Degrees of freedom (DoF) Manager class is to handle all
degrees of freedom associated to elements, faces and nodes.
Moreover, it creates an unique map between local and global DoFs.

Key concepts are locations and connectors.
Locations, that can be elements, faces or nodes, represent where the DoF is
assigned.
For example, a DoF for pressure in a two-point flux approximation will be on a
cell (i.e. element), while a displacement DoF for structural equations will be
on a node.
The counterparts of locations are connectors, that are the geometrical entities
that link together different DoFs are create the sparsity pattern.
Connectors can be elements, faces, edges, nodes or none.
Using the same example as before, connectors will be faces and cells,
respectively.
The case of a mass matrix, where every element is linked only to itself, is an
example when there are no connectors, i.e. these have to be set to none.

DoFs are handled through sparse matrices.
Inside the node, a simple compress row storage (CRS) matrix type is implemented
in the class itself and used to collect connectivity-location information.
Among nodes, the ``ParallelMatrix`` data structure, i.e. a distributed CSR, is
used.
CSR matrices contain information about the sparsity pattern only, thus binary
matrices are used.

DoFs located on a mesh object are owned by the same rank that owns the object
in parallel mesh partitioning.
DoFs are assigned global indices as follows:

1. Contiguous across MPI ranks
2. Within each rank, sequentially across present fields
3. Within each field, sequentially across owned locations
4. Within each location, according to DoF component number

DoF Manager allocates a separate "Dof Index" array for each field on the mesh.
It is an array of global indices, where each value represents the first DoF
index for that field and location (or equivalently, the row and column offset
of that location's equations and variables for the field in the global matrix).
For example, if index array a field with 3 components contains the value N,
global DoF numbers for that location will be N, N+1, N+2.
DoF on ghosted locations have the same indices as on the owning rank.
The array is stored under a generated key, which can be queried from the DoF
manager, and is typically used in system assembly.

This last assumption (sequential numbering of components) cannot be lifted
as it will necessarily increase memory footprint and bandwidth.
Likewise, the first principle (rank-wise ordering) is unlikely to change
as it is a requirement of many linear solver packages to have contiguously
numbered matrix rows on each rank.
However, other numbering schemes (e.g. matrix bandwidth reducing permutations)
may be implemented in future for points 2 and 3.

Methods
========================

The main methods of ``DoF Manager`` are:

* ``setMesh``: sets which portion of the mesh the DoF manager instance is
  referring to.
  ``domain`` identifies the global mesh, while ``meshLevelIndex`` and
  ``meshBodyIndex`` determine the specific level and body, respectively.

 .. code-block:: c

  void setMesh(DomainPartition * const domain,
               localIndex const meshLevelIndex = 0,
               localIndex const meshBodyIndex = 0);

* ``addField``: creates a new set of DoF, labeled ``field``, with specific
  ``location`` and ``connectivity``.
  Default number of ``components`` is ``1``, like for pressure in flux.
  Default ``regions`` is the empty string, meaning all domain.

 .. code-block:: c

  void addField(string const & field,
                Location const location,
                Connectivity const connectivity,
                localIndex const components,
                string_array const & regions);

* ``addCoupling``: creates a coupling between two fields (``rowField`` and
  ``colField``) according to a given ``connectivity`` in the regions defined
  by ``regions``.
  Of course, both fields (row and column) must have already been defined on
  the regions where is required the coupling among them.
  Default value for ``regions`` is the whole intersection between the regions
  where the first and the second fields are defined.
  This method can create also the coupling between ``colField`` and
  ``rowField``, i.e. the transpose of the rectangular sparsity pattern.
  This behaviour is avoided by default, but can be activated setting
  ``symmetric`` to ``true``.

 .. code-block:: c

  void addCoupling(string const & rowField,
                   string const & colField,
                   Connectivity const connectivity,
                   string_array const & regions,
                   bool const symmetric);

* ``close``: finish populating field and coupling information and apply DoF
  numbering

 .. code-block:: c

  void close();

* ``clear``: removes all fields, releases memory and re-opens the DofManager

 .. code-block:: c

  void clear();

* ``getSparsityPattern``: gets the sparsity ``pattern`` for the given
  ``rowField`` and ``colField``.
  Default case is the complete sparsity pattern, for all DoFs.

 .. code-block:: c

  void getSparsityPattern( SparsityPattern & pattern,
                           string const & rowField = "",
                           string const & colField = "") const;

Minor methods are:

* ``numGlobalDofs``: returns the total number of DoFs across all processors for
  the specified ``field``.
  Whenever no ``field`` is specified, it returns the total number of DoFs for
  all fields.

 .. code-block:: c

  globalIndex numGlobalDofs( string const & field = "" ) const;

* ``numLocalDofs``: returns the number of DoFs on this process for the
  specified ``field``.
  As for ``numGlobalDofs``, when ``field`` is lacking, all fields are
  considered.

 .. code-block:: c

  localIndex numLocalDofs( string const & field = "" ) const;

* ``printConnectivityLocationPattern``: prints the connectivity-location
  pattern for ``field``.
  Unless a filename is provided, the default behaviour is to print on screen.

 .. code-block:: c

  void printConnectivityLocationPattern( string const & field,
                                         string const & fileName = "" ) const;


* ``printParallelMatrix``: prints a ``ParallelMatrix`` on the file named
  ``filename`` using MatrixMarket format (MTX file).

 .. code-block:: c

  void printParallelMatrix( ParallelMatrix const & matrix,
                            string const & filename ) const;

Example
=======

Here we show how the sparsity pattern is computed for a simple 2D quadrilateral
mesh with 6 elements.
Unknowns are pressure, located on the element center, and displacements (*x*
and *y* components), located on the nodes.
For fluxes, a two-point flux approximation (TPFA) is used.
The representation of the sparsity pattern of the :math:`\mathsf{C_L}` matrix
(connectors/locations) for the simple mesh, shown in :numref:`meshFig`, is
reported in :numref:`CLFig`.
It can be notices that the two unknowns for the displacements *x* and *y* are
grouped together.
Elements are the connectivity for DoF on nodes (Finite Element Method for
displacements) and on elements (pressures).
Faces are the connectivity for DoF on elements (Finite Volume Method for
pressure), being the flux computation based on the pressure on the two adjacent
elements.

.. _meshFig:
.. figure:: /coreComponents/linearAlgebraInterface/docs/images/mesh2D.svg
   :align: center
   :width: 250
   :figclass: align-center

   Small 2D quadrilateral mesh used for this examples.
   Nodes are label with black numbers, elements with light gray numbers and
   faces with italic dark gray numbers.

.. _CLFig:
.. figure:: /coreComponents/linearAlgebraInterface/docs/images/CL.svg
   :align: center
   :width: 500
   :figclass: align-center

   Sparsity pattern of the binary matrix connections/locations.

The global sparsity pattern, shown in :numref:`patternFig`, is obtained through
the symbolic multiplication of the transpose of the matrix :math:`\mathsf{C_L}`
and the matrix itself, i.e. :math:`\mathsf{P = C_L^T C_L}`.

.. _patternFig:
.. figure:: /coreComponents/linearAlgebraInterface/docs/images/pattern.svg
   :align: center
   :width: 400
   :figclass: align-center

   Sparsity pattern of the global matrix, where red and green entries are
   related to the displacement field and to the pressure field, respectively.
   Blue entries represent coupling blocks.

Real mesh and patterns
======================

Now we build the pattern of the Jacobian matrix for a simple 3D mesh, shown in
:numref:`meshCubeFig`. Fields are:

- displacement (location: node, connectivity: element) defined on the blue,
  orange and red regions;
- pressure (location: element, connectivity: face) defined on the green, orange
  and red regions;
- mass matrix (location: element, connectivity: element) defined on the green
  region only.

Moreover, following coupling are imposed:

- displacement-pressure (connectivity: element) on the orange region only;
- pressure-mass matrix and transpose (connectivity: element) everywhere it is
  possibile.

.. _meshCubeFig:
.. figure:: /coreComponents/linearAlgebraInterface/docs/images/meshCube3D.svg
   :align: center
   :width: 400
   :figclass: align-center

   Real mesh used to compute the Jacobian pattern.

:numref:`globalPatterFig` shows the global pattern with the field-based
ordering of unknowns.
Different colors mean different fields.
Red unkwnons are associated with displacement, yellow ones with pressure and
blue ones with mass matrix.
Orange means the coupling among displacement and pressure, while green is the
symmetric coupling among pressure and mass matrix.

.. _globalPatterFig:
.. figure:: /coreComponents/linearAlgebraInterface/docs/images/global.svg
   :align: center
   :width: 400
   :figclass: align-center

   Global pattern with field-based ordering.
   Red is associated with displacement unknowns, yellow with pressure ones and
   blue with those of mass matrix field.
   Orange means the coupling among displacement and pressure, while green is
   the symmetric coupling among pressure and mass matrix.

:numref:`permutedPatterFig` shows the global pattern with the MPI rank-based
ordering of unknowns.
In this case, just two processes are used.
Again, different colors indicate different ranks.

.. _permutedPatterFig:
.. figure:: /coreComponents/linearAlgebraInterface/docs/images/permutedGlobal.svg
   :align: center
   :width: 400
   :figclass: align-center

   Global pattern with MPI rank-based ordering.
   Red unkwnons are owned by rank 0 and green ones by rank 1.
   Blue indicates the coupling among the two processes.
