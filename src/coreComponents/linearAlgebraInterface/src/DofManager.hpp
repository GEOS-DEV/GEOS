/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file DofManager.hpp
 */

#ifndef SRC_CORECOMPONENTS_LINEARALGEBRAINTERFACE_SRC_DOFMANAGER_HPP_
#define SRC_CORECOMPONENTS_LINEARALGEBRAINTERFACE_SRC_DOFMANAGER_HPP_

#include "common/DataTypes.hpp"
#include "managers/DomainPartition.hpp"
#include "MPI_Communications/CommunicationTools.hpp"
#include "MPI_Communications/NeighborCommunicator.hpp"
#include "mesh/MeshLevel.hpp"
#include "dataRepository/ManagedGroup.hpp"
#include "dataRepository/MappedVector.hpp"
#include "TrilinosInterface.hpp"

namespace geosx
{

/**
 * Defines a simple GEOSX sparsity pattern.
 * It is intended to be a lightweight way to communicate between the
 * DoFManager and the LAI implementation, without relying on a specific
 * LAI choice.
 *
 * TODO: replace with proper CrsArray
 */
struct SparsityPattern
{
  // Set default values (empty matrix)
  SparsityPattern() :
      nRows( 0 ),
      nCols( 0 )
  {
  }
  localIndex nRows; //<! number of rows
  localIndex nCols; //<! number of columns
  localIndex_array rowLengths; //<! row lengths, size numLocalRows
  globalIndex_array colIndices; //<! packed column indices, size numLocalNonZeros
};

using ParallelMatrix = typename TrilinosInterface::ParallelMatrix;

/**
 * The DoFManager is responsible for allocating global dofs, constructing
 * sparsity patterns, and generally simplifying the interaction between
 * PhysicsSolvers and linear algebra operations.
 */
class DofManager
{
public:

  /**
   * Constructor
   */
  DofManager();

  /**
   * Destructor
   */
  ~DofManager() = default;

  /**
   * Enumeration of geometric objects for support location.  Note that this enum
   * is nearly identical to Connectivity, but we keep both for code readability
   * in function calls.
   */
  enum class Location
  {
    Elem, Face, Node
  };

  /**
   * Enumeration of geometric objects for connectivity type.  Note that this enum
   * is nearly identical to Location, but we keep both for code readability
   * in function calls.
   */
  enum class Connectivity
  {
    Elem, Face, Node, None
  };

  /**
   * Assign a mesh
   */
  void setMesh( DomainPartition * const domain,
                localIndex const meshLevelIndex = 0,
                localIndex const meshBodyIndex = 0 );

  /**
   * Add fields.
   * The user can add a field with a support location, connectivity type, string key, number of scalar components,
   * and a list of element regions over which the field is active.  If the region list is empty, it is assumed
   * the field exists on all regions.
   *
   * The connectivity type is used to infer the sparsity pattern that connects degrees of freedom.
   * If LC denotes a boolean connectivity graph between support locations L and connectors C, the desired sparsity
   * pattern will be computed as LC*CL.  For example, for a TPFA discretization we have dofs located at cell centers,
   * and connected through adjacent faces.  In this example, LC is the cell-to-face connectivity, and LC*CL is the
   * desired TPFA sparsity pattern.  More generally,
   *
   * Example 1 = ("displacement",NODE,ELEM,3) for a Q1 finite-element interpolation for elasticity
   * Example 2 = ("pressure",ELEM,FACE,1) for a scalar TPFA-type approximation
   * Example 3 = ("pressure",ELEM,NODE,1) for a scalar MPFA-type approximation
   * Example 4 = ("mass",ELEM,NONE,1) for a diagonal-only sparsity pattern (no connectivitys)
   *
   * When the number of components is greater than one, we always assume they are tightly coupled to one another
   * and form a dense block.  The sparsity pattern LC*CL is then interpreted as the super-node pattern, containing
   * dense sub-blocks.
   */

  /**
   * Just an interface to allow only three parameters
   */
  void addField( string const & field,
                 Location const location,
                 Connectivity const connectivity );

  /**
   * Just another interface to allow four parameters (no regions)
   */
  void addField( string const & field,
                 Location const location,
                 Connectivity const connectivity,
                 localIndex const components );

  /**
   * Just another interface to allow four parameters (no components)
   */
  void addField( string const & field,
                 Location const location,
                 Connectivity const connectivity,
                 string_array const & regions );

  /**
   * The real function, allowing the creation of self-connected blocks
   */
  void addField( string const & field,
                 Location const location,
                 Connectivity const connectivity,
                 localIndex const components,
                 string_array const & regions );

  /**
   * Add coupling between two fields.
   * The connectivity argument defines how the two fields couple. If the first field has support location A,
   * the second field has support location B, and the connecting object is C, the sparsity pattern will be
   * defined as (AC)(CB).  The final argument indicates if the coupling is symmetric, in the sense that there
   * is a two-way coupling between the fields.  Without this argument, a nonzero block will be added to the
   * system matrix for block AB, but block BA will remain zero (one-way coupling).
   *
   * Example 1 = ("node_field","elem_field", ELEM, true) couples all dofs sharing a common element (two-way coupling)
   * Example 2 = ("node_field_1","node_field_2", NODE, true) couples all dofs sharing a common node (two-way coupling)
   * Example 3 = ("node_field_1","face_field", NODE, false) couples nodal dofs to adjacent faces (one-way coupling)
   */

  /**
   * Just an interface to allow only three parameters
   */
  void addCoupling( string const & rowField,
                    string const & colField,
                    Connectivity const connectivity ) const;

  /**
   * Just another interface to allow four parameters (no symmetry)
   */
  void addCoupling( string const & rowField,
                    string const & colField,
                    Connectivity const connectivity,
                    string_array const & regions ) const;

  /**
   * Just another interface to allow four parameters (no regions)
   */
  void addCoupling( string const & rowField,
                    string const & colField,
                    Connectivity const connectivity,
                    bool const symmetric ) const;

  /**
   * The real function, allowing the creation of coupling blocks
   */
  void addCoupling( string const & rowField,
                    string const & colField,
                    Connectivity const connectivity,
                    string_array const & regions,
                    bool const symmetric ) const;

  /**
   * Return global number of dofs across all processors. If field argument is empty, return monolithic size.
   */
  globalIndex numGlobalDofs( string const & field = "" ) const;

  /**
   * Return local number of dofs on this processor.  If field argument is empty, return monolithic size
   */
  localIndex numLocalDofs( string const & field = "" ) const;

  /**
   * Get a sparsity pattern.  Without additional arguments, this function routines the sparsity pattern for the
   * monolithic matrix.  Sub-patterns can be extracted, however, using row and column field keys.
   */
  void getSparsityPattern( ParallelMatrix & locLocDistr,
                           string const & rowField = "",
                           string const & colField = "" ) const;

  /**
   * Get global indices for dofs connected by the connector type.  We have two versions, since cells need
   * three indices while faces and nodes only need two.  This keeps the interface the same, but we will only
   * implement appropriate combinations.
   *
   * Example 1 = getIndices(indices,ELEM,er,esr,ei,"pressure") = get pressure indices connected to this cell
   * Example 2 = getIndices(indices,FACE,fi,"pressure") = get pressure indices connected to this face
   * Example 3 = getIndices(indices,NODE,ni,"pressure") = get pressure indices connected to this node
   */
  void getIndices( globalIndex_array & indices,
                   Connectivity const connectivity,
                   localIndex const region,
                   localIndex const subregion,
                   localIndex const index,
                   string const & field = "" ) const;

  /**
   * Get global indices for dofs connected by the connector type.  We have two versions, since cells need
   * three indices while faces and nodes only need two.  This keeps the interface the same, but we will only
   * implement appropriate combinations.
   *
   * Example 1 = getIndices(indices,ELEM,er,esr,ei,"pressure") = get pressure indices connected to this cell
   * Example 2 = getIndices(indices,FACE,fi,"pressure") = get pressure indices connected to this face
   * Example 3 = getIndices(indices,NODE,ni,"pressure") = get pressure indices connected to this node
   */
  void getIndices( globalIndex_array & indices,
                   Connectivity const connectivity,
                   localIndex const index,
                   string const & field = "" ) const;

  /**
   * Print the global connectivity matrix
   */
  void printConnectivityMatrix() const;

  /**
   * Print the connectivity-location pattern for a specific field
   */
  void printConnectivityLocationPattern( string const & field, string const & fileName = "" ) const;

  /**
   * Print the given parallel matrix in Matrix Market format (MTX file)
   */
  void printParallelMatrix( ParallelMatrix const & matrix,
                            string const & filename ) const;

  /**
   * Print a CSR pattern on file or on screen
   */
  void printSparsityPattern( SparsityPattern const & pattern, string const & fileName = "" ) const;

  /**
   * Release internal storage
   */
  void cleanUp() const;

private:
  /**
   *  Limit on max number of fields
   */
  localIndex const static MAX_NUM_FIELDS = 10;

  /**
   * Pointer to domain manager
   */
  DomainPartition * m_domain = nullptr;

  /**
   * Pointer to corresponding MeshLevel
   */
  MeshLevel * m_meshLevel = nullptr;

  /**
   * Field description
   */
  struct FieldDescription
  {
    string name; //!< field name
    array1d<string> regionNames; //!< active element regions
    array1d<ElementRegion*> regionPtrs; //!< saved pointers to active regions
    Location location; //!< support location
    localIndex numComponents; //!< number of vector components
    string key; //!< string key for index array
    string docstring; //!< documentation string
    localIndex numLocalNodes; //!< number of local nodes
    localIndex numLocalRows; //!< number of local rows
    globalIndex numGlobalRows; //!< number of ghost rows
    globalIndex firstLocalRow; //!< first row on this processor (without field offset)
    globalIndex fieldOffset; //!< global row offset for multi-field problems
    globalIndex firstLocalConnectivity; //!< first connector on this processor
    globalIndex numGlobalConnectivity; //!< number of connector for this field
    ParallelMatrix* connLocPattern; //!< pattern for the connectivity-location matrix
  };

  /**
   * Array of field descriptions
   */
  array1d<FieldDescription> m_fields;

  /**
   * Table of connectivities within and between fields
   */
  array2d<Connectivity> m_connectivity;

  /**
   * Definifion for entries of sparse matrices collection
   */
  typedef std::pair<ParallelMatrix*, ParallelMatrix*> matrixPair;

  /**
   * Table of sparsity patterns within and between fields
   */
  array2d<matrixPair> m_sparsityPattern;

  /**
   * Number of MPI ranks
   */
  int mpiSize;

  /**
   * This mpi rank
   */
  int mpiRank;

  /**
   * Check if string key is already being used
   */
  bool keyInUse( string const & key ) const;

  /**
   * Get field index from string key
   */
  localIndex fieldIndex( string const & key ) const;

  /**
   * Create index array
   */
  void createIndexArray_NodeOrFaceVersion( FieldDescription & field,
                                           localIndex_array const & activeRegionsInput = localIndex_array() ) const;

  /**
   * Create element index array
   */
  void createIndexArray_ElemVersion( FieldDescription & field ) const;

  /**
   * Create sparsity pattern for a field with itself (diagonal entries in the
   * connectivity matrix)
   */
  void addDiagSparsityPattern( SparsityPattern & connLocPatt,
                               localIndex const & fieldIdx,
                               Connectivity const connectivity,
                               localIndex_array const & activeRegionsInput = localIndex_array() ) const;

  /**
   * Create sparsity pattern for two fields (extra-diagonal entries in the
   * connectivity matrix)
   */
  void addExtraDiagSparsityPattern( ParallelMatrix *& rowConnLocPattDistr,
                                    ParallelMatrix *& colConnLocPattDistr,
                                    localIndex const & rowFieldIndex,
                                    localIndex const & colFieldIndex,
                                    localIndex_array const & rowActiveRegions,
                                    localIndex_array const & colActiveRegions,
                                    Connectivity const connectivity ) const;

  /**
   * Definifion for entries of sparse matrix in COO format
   */
  typedef std::pair<localIndex, globalIndex> indexPair;

  /**
   * Structure used to create CSR matrix from COO format
   */
  struct pairComparison
  {
    inline bool operator()( const indexPair& lhs, const indexPair& rhs ) const
                            {
      if( lhs.first < rhs.first )
        return true;
      else if( lhs.first == rhs.first )
        return lhs.second < rhs.second;
      else
        return false;
    }
    ;
  };

  /**
   * Convert a sparse matrix in COO format in the CSR version
   */
  void vectorOfPairsToCSR( array1d<indexPair> const & pairs,
                           localIndex const nRows,
                           localIndex const nCols,
                           SparsityPattern & pattern ) const;

  /**
   * Just a local INT_MAX
   */
  globalIndex const globalIndexMax = std::numeric_limits<globalIndex>::max();

  /**
   * Map a global row index to local row index
   */
  localIndex ParallelMatrixGetLocalRowID( EpetraMatrix const &A, globalIndex const index ) const;

  /**
   * Performe a matrix matrix product with Parallel Matrix
   */
  void MatrixMatrixMultiply( EpetraMatrix const &A,
                             bool const transA,
                             EpetraMatrix const &B,
                             bool const transB,
                             EpetraMatrix &C,
                             bool const call_FillComplete = true ) const;

};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_LINEARALGEBRAINTERFACE_SRC_DOFMANAGER_HPP_ */
