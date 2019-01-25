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
Connectors can be elements, faces, node or none.
Using the same example as before, connectors will be faces and cells,
respectively.
The case of a mass matrix, where every element is linked only to itself, is an
example when there are no connectors, i.e. these have to be set to none.

DoFs are handled through sparse matrices (see CSR description ...).
CSR matrices contain information about the sparsity pattern only, thus binary
matrices are used.

Methods
========================

The main functionalities of ``DoF Manager`` are:

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

 .. code-block:: c

  void addField(string const & field,
                Location const location,
                Connectivity const connectivity,
                localIndex const components = 1,
                string_array const & regions = string_array() );

* ``addCoupling``: creates a coupling between two fields (``rowField`` and
  ``colField``) according to a given ``connectivity``.
  It can creates also the coupling between ``colField`` and ``rowField``, i.e.
  it transposes the rectangular sparsity pattern.
  This behaviour is avoided by default, but can be activated setting
  ``symmetric`` to ``true``.

 .. code-block:: c

  void addCoupling(string const & rowField,
                   string const & colField,
                   Connectivity const connectivity,
                   bool const symmetric = false);

Minor methods are:

* ``numGlobalDofs``: returns the total number of DoFs across all processors for
  the specified ``field``.
  Whenever no ``field`` is specified, it returns the total number of DoFs for
  all fields.

 .. code-block:: c

  globalIndex numGlobalDofs( string const & field = "all" ) const;

* ``numLocalDofs``: returns the number of DoFs on this process for the
  specified ``field``.
  As for ``numGlobalDofs``, when ``field`` is lacking, all fields are
  considered.

 .. code-block:: c

  localIndex numLocalDofs( string const & field = "all" ) const;

* ``getSparsityPattern``: gets the sparsity ``pattern`` for the given
  ``rowField`` and ``colField``.
  Default case is the complete sparsity pattern, for all DoFs.

 .. code-block:: c

  void getSparsityPattern( SparsityPattern & pattern,
                           string const & rowField = "all",
                           string const & colField = "all") const;

* ``getIndices``: gets global ``indices`` for DoFs with a given local ``index``
  and linked by a specific ``connectivity``.
  When ``field`` is not set, all DoFs assigned to the local ``index`` are
  considered.
  In case of DoF located on elements, there is the need to set ``region`` and
  ``subregion``.

 .. code-block:: c

  void getIndices( globalIndex_array & indices,
                   Connectivity const connectivity,
                   localIndex const region,
                   localIndex const subregion,
                   localIndex const index,
                   string const & field = "all") const;

  void getIndices( globalIndex_array & indices,
                   Connectivity const connectivity,
                   localIndex const index,
                   string const & field = "all") const;

Example
=======

Here we show how the sparsity pattern is computed for a simple 2D quadrilateral
mesh with 6 elements.
Unknowns are pressure, located on the element center, and displacements (*x*
and *y* components), located on the nodes.
For fluxes, a two-point flux approximation (TPFA) is used.
The representation of the sparsity pattern of the :math:`\mathsf{C_L}` matrix
(connectors/locations) for the simple mesh, shown in :numref:`mesh`, is
reported in :numref:`CL`.
It can be notices that the two unknowns for the displacements *x* and *y* are
grouped together.
Elements are the connectivity for DoF on nodes (Finite Element Method for
displacements) and on elements (pressures).
Faces are the connectivity for DoF on elements (Finite Volume Method for
pressure), being the flux computation based on the pressure on the two adjacent
elements.

.. _mesh:
.. figure:: ../../coreComponents/linearAlgebraInterface/docs/images/mesh.svg
   :align: center
   :width: 250
   :figclass: align-center

   Small 2D quadrilateral mesh used for this examples.
   Nodes are label with black numbers, elements with light gray numbers and
   faces with italic dark gray numbers.

.. _CL:
.. figure:: ../../coreComponents/linearAlgebraInterface/docs/images/CL.svg
   :align: center
   :width: 500
   :figclass: align-center

   Sparsity pattern of the binary matrix connections/locations.

The global sparsity pattern, shown in :numref:`pattern`, is obtained through
the symbolic multiplication of the transpose of the matrix :math:`\mathsf{C_L}`
and the matrix itself, i.e. :math:`\mathsf{P = C_L^T C_L}`.

.. _pattern:
.. figure:: ../../coreComponents/linearAlgebraInterface/docs/images/pattern.svg
   :align: center
   :width: 400
   :figclass: align-center

   Sparsity pattern of the global matrix, where red and green entries are
   related to the displacement field and to the pressure field, respectively.
   Blue entries represent coupling blocks.
