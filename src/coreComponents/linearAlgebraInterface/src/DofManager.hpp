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

#include <numeric>
#include "common/DataTypes.hpp"
#include "managers/DomainPartition.hpp"
#include "MPI_Communications/CommunicationTools.hpp"
#include "MPI_Communications/NeighborCommunicator.hpp"
#include "mesh/MeshLevel.hpp"
#include "dataRepository/ManagedGroup.hpp"
#include "dataRepository/MappedVector.hpp"
#include "InterfaceTypes.hpp"

namespace geosx
{

/**
 * @struct Dof_SparsityPattern
 *
 * @brief Defines a simple GEOSX sparsity pattern.
 * It is intended to be a lightweight way to communicate between the
 * DoFManager and the LAI implementation, without relying on a specific
 * LAI choice.
 *
 * TODO: replace with proper CrsArray
 */
struct Dof_SparsityPattern
{
  // Set default values (empty matrix)
  Dof_SparsityPattern() :
      nRows( 0 ),
      nCols( 0 )
  {
  }
  localIndex nRows; //!< number of rows
  localIndex nCols; //!< number of columns
  localIndex_array rowLengths; //!< row lengths, size numLocalRows
  globalIndex_array colIndices; //!< packed column indices, size numLocalNonZeros
  localIndex_array nnzEntries; //!< packed values (of type localIndex), size numLocalNonZeros
};

/**
 * @class DofManager
 * @brief The DoFManager is responsible for allocating global dofs, constructing
 * sparsity patterns, and generally simplifying the interaction between
 * PhysicsSolvers and linear algebra operations.
 */
class DofManager
{
public:

  /**
   * @brief Constructor.
   *
   * @param [in] verbosity Optional localIndex setting the verbosity level.
   *                       - 0: nothing (default)
   *                       - >0: minimal info
   */
  DofManager( localIndex const verbosity = 0 );

  /**
   * @brief Destructor.
   */
  ~DofManager()
  {
    cleanUp();
  }

  /**
   * @brief Enumeration of geometric objects for support location. Note that this
   * enum is nearly identical to Connectivity, but we keep both for code readability
   * in function calls.
   */
  enum class Location
  {
    Elem, //!< location is element (like pressure in finite volumes)
    Face, //!< location is face (like flux in mixed finite elements)
    Node, //!< location is node (like displacements in finite elements)
    USER_DEFINED //!< user defined location (for input connectivity pattern)
  };

  /**
   * @brief Enumeration of geometric objects for connectivity type. Note that this
   * enum is nearly identical to Location, but we keep both for code readability
   * in function calls.
   */
  enum class Connectivity
  {
    Elem, //!< connectivity is element (like in finite elements)
    Face, //!< connectivity is face (like in finite volumes TPFA)
    Node, //!< connectivity is node (like in finite volumes MPFA)
    None, //!< there is no connectivity (self connected field, like a mass matrix)
    USER_DEFINED //!< user defined connectivity (for input connectivity pattern)
  };

  /**
   * @brief Define the possible status of a location: undefined or belonging to ghost
   * locations of this MPI process.
   */
  enum class LocationStatus : globalIndex
  {
    notAssigned = -1, //!< location is not yet assigned
    notMyGhostLocation = -2 //!< location is not part of local ghost locations
  };

  /**
   * @brief Assign a mesh.
   *
   * @param [in] domain DomainPartition the input domain.
   * @param [in] meshLevelIndex Optional localIndex the mesh level.
   * @param [in] meshBodyIndex Optional localIndex the body level.
   */
  void setMesh( DomainPartition * const domain,
                localIndex const meshLevelIndex = 0,
                localIndex const meshBodyIndex = 0 );

  /**
   * @brief Just an interface to allow only three parameters.
   *
   * @param [in] field string the name of the field.
   * @param [in] location Location where it is defined.
   * @param [in] connectivity Connectivity through what it is connected.
   */
  void addField( string const & field,
                 Location const location,
                 Connectivity const connectivity );

  /**
   * @brief Just another interface to allow four parameters (no regions, default is everywhere).
   *
   * @param [in] field string the name of the field.
   * @param [in] location Location where it is defined.
   * @param [in] connectivity Connectivity through what it is connected.
   * @param [in] components localIndex number of components (for vector fields).
   */
  void addField( string const & field,
                 Location const location,
                 Connectivity const connectivity,
                 localIndex const components );

  /**
   * @brief Just another interface to allow four parameters (no components, default is 1).
   *
   * @param [in] field string the name of the field.
   * @param [in] location Location where it is defined.
   * @param [in] connectivity Connectivity through what it is connected.
   * @param [in] regions string_array where this field is defined.
   */
  void addField( string const & field,
                 Location const location,
                 Connectivity const connectivity,
                 string_array const & regions );

  /**
   * @brief The user can add a field with a support location, connectivity type, string key, number of scalar
   * components, and a list of element regions over which the field is active. If the region list is empty,
   * it is assumed the field exists on all regions.
   *
   * The connectivity type is used to infer the sparsity pattern that connects degrees of freedom.
   * If LC denotes a boolean connectivity graph between support locations L and connectors C, the desired sparsity
   * pattern will be computed as LC*CL. For example, for a TPFA discretization we have dofs located at cell centers,
   * and connected through adjacent faces. In this example, LC is the cell-to-face connectivity, and LC*CL is the
   * desired TPFA sparsity pattern. More generally,
   *
   * - Example 1 = ("displacement",NODE,ELEM,3) for a Q1 finite-element interpolation for elasticity
   * - Example 2 = ("pressure",ELEM,FACE,1) for a scalar TPFA-type approximation
   * - Example 3 = ("pressure",ELEM,NODE,1) for a scalar MPFA-type approximation
   * - Example 4 = ("mass",ELEM,NONE,1) for a diagonal-only sparsity pattern (no connectivitys)
   *
   * When the number of components is greater than one, we always assume they are tightly coupled to one another
   * and form a dense block. The sparsity pattern LC*CL is then interpreted as the super-node pattern, containing
   * dense sub-blocks.
   *
   * @param [in] field string the name of the field.
   * @param [in] location Location where it is defined.
   * @param [in] connectivity Connectivity through what it is connected.
   * @param [in] components localIndex number of components (for vector fields).
   * @param [in] regions string_array where this field is defined.
   */
  void addField( string const & field,
                 Location const location,
                 Connectivity const connectivity,
                 localIndex const components,
                 string_array const & regions );

  /**
   * @brief Interface to allow only two parameters.
   *
   * @param [in] field string the name of the field.
   * @param [in] connLocInput ParallelMatrix input LC pattern.
   */
  void addField( string const & field,
                 ParallelMatrix const & connLocInput );

  /**
   * @brief Just another interface to allow three parameters (no connectivity, default is USER_DEFINED).
   *
   * @param [in] field string the name of the field.
   * @param [in] connLocInput ParallelMatrix input LC pattern.
   * @param [in] components localIndex number of components (for vector fields).
   */
  void addField( string const & field,
                 ParallelMatrix const & connLocInput,
                 localIndex const components );

  /**
   * @brief Just another interface to allow three parameters (no components, default is 1).
   *
   * @param [in] field string the name of the field.
   * @param [in] connLocInput ParallelMatrix input LC pattern.
   * @param [in] connectivity Connectivity through what it is connected.
   */
  void addField( string const & field,
                 ParallelMatrix const & connLocInput,
                 Connectivity const connectivity );

  /**
   * @brief addField with an input pattern.
   * Allow the usage of a predefined location-connection pattern (user-defined).
   *
   * @param [in] field string the name of the field.
   * @param [in] connLocInput ParallelMatrix input LC pattern.
   * @param [in] components localIndex number of components (for vector fields).
   * @param [in] connectivity Connectivity through what it is connected.
   */
  void addField( string const & field,
                 ParallelMatrix const & connLocInput,
                 localIndex const components,
                 Connectivity const connectivity );

  /**
   * @brief Just an interface to allow only three parameters.
   *
   * @param [in] rowField string the name of the row field.
   * @param [in] colField string the name of the col field.
   * @param [in] connectivity Connectivity through what they are connected.
   */
  void addCoupling( string const & rowField,
                    string const & colField,
                    Connectivity const connectivity );

  /**
   * @brief Just another interface to allow four parameters (no symmetry, default is true).
   *
   * @param [in] rowField string the name of the row field.
   * @param [in] colField string the name of the col field.
   * @param [in] connectivity Connectivity through what they are connected.
   * @param [in] regions string_array where this coupling is defined.
   */
  void addCoupling( string const & rowField,
                    string const & colField,
                    Connectivity const connectivity,
                    string_array const & regions );

  /**
   * @brief Just another interface to allow four parameters (no regions, default is everywhere).
   *
   * @param [in] rowField string the name of the row field.
   * @param [in] colField string the name of the col field.
   * @param [in] connectivity Connectivity through what they are connected.
   * @param [in] symmetric bool is it symmetric, i.e., both row-col and col-row?
   */
  void addCoupling( string const & rowField,
                    string const & colField,
                    Connectivity const connectivity,
                    bool const symmetric );

  /**
   * @brief Add coupling between two fields.
   * The connectivity argument defines how the two fields couple. If the first field has support location A,
   * the second field has support location B, and the connecting object is C, the sparsity pattern will be
   * defined as (AC)(CB). The final argument indicates if the coupling is symmetric, in the sense that there
   * is a two-way coupling between the fields. Without this argument, a nonzero block will be added to the
   * system matrix for block AB, but block BA will remain zero (one-way coupling).
   *
   * - Example 1 = ("node_field","elem_field", ELEM, true) couples all dofs sharing a common element (two-way coupling)
   * - Example 2 = ("node_field_1","node_field_2", NODE, true) couples all dofs sharing a common node (two-way coupling)
   * - Example 3 = ("node_field_1","face_field", NODE, false) couples nodal dofs to adjacent faces (one-way coupling)
   *
   * @param [in] rowField string the name of the row field.
   * @param [in] colField string the name of the col field.
   * @param [in] connectivity Connectivity through what they are connected.
   * @param [in] regions string_array where this coupling is defined.
   * @param [in] symmetric bool is it symmetric, i.e., both row-col and col-row?
   */
  void addCoupling( string const & rowField,
                    string const & colField,
                    Connectivity const connectivity,
                    string_array const & regions,
                    bool const symmetric );

  /**
   * @brief Return the key used to record the field in the DofManager.
   *
   * @param [in] field string the name of the field.
   */
  string getKey( string const & field ) const;

  /**
   * @brief Return global number of dofs across all processors. If field argument is empty, return
   * monolithic size.
   *
   * @param [in] field Optional string the name of the field.
   */
  globalIndex numGlobalDofs( string const & field = "" ) const;

  /**
   * @brief Return local number of dofs on this processor. If field argument is empty, return
   * monolithic size.
   *
   * @param [in] field Optional string the name of the field.
   */
  localIndex numLocalDofs( string const & field = "" ) const;

  /**
   * @brief Return the sum of local dofs across all previous processors w.r.t. to the calling one for
   * the specified field.
   *
   * @param [in] field Optional string the name of the field.
   */
  localIndex offsetLocalDofs( string const & field = "" ) const;

  /**
   * @brief Set a sparsity pattern. Without additional arguments, this function provides the sparsity
   * pattern for the monolithic matrix. Sub-patterns can be extracted, however, using row and column
   * field keys.
   *
   * @param [out] locLocDistr ParallelMatrix the location-location sparsity pattern (LC*CL).
   * @param [in]  rowField Optional string the name of the row field.
   * @param [in]  colField Optional string the name of the col field.
   */
  void setSparsityPattern( ParallelMatrix & locLocDistr,
                           string const & rowField = "",
                           string const & colField = "" ) const;

  /**
   * @brief Set a sparsity pattern. Low level version.
   *
   * @param [out] locLocDistr ParallelMatrix the location-location sparsity pattern (LC*CL).
   * @param [in]  rowFieldIndex localIndex row field index (-1 means all fields).
   * @param [in]  colFieldIndex localIndex col field index (-1 means all fields).
   */
  void setSparsityPattern( ParallelMatrix & locLocDistr,
                           localIndex const rowFieldIndex,
                           localIndex const colFieldIndex ) const;

  /**
   * @brief Allocate a vector. Without additional arguments, this function provides the a vector
   * consistent with the sparsity pattern for the monolithic matrix. Sub-vectors can be extracted,
   * however, using row and column field keys.
   *
   * @param [out] vector ParallelVector the output vector.
   * @param [in]  rowField Optional string the name of the row field.
   * @param [in]  colField Optional string the name of the col field.
   */
  void setVector( ParallelVector & vector,
                  string const & rowField = "",
                  string const & colField = "" ) const;

  /**
   * @brief Allocate a vector. Low level version.
   *
   * @param [out] vector ParallelVector the output vector.
   * @param [in]  rowFieldIndex localIndex row field index (-1 means all fields).
   * @param [in]  colFieldIndex localIndex col field index (-1 means all fields).
   */
  void setVector( ParallelVector & vector,
                  localIndex const rowFieldIndex,
                  localIndex const colFieldIndex ) const;

  /**
   * @brief Get global indices for dofs connected by the connector type. We have two versions, since cells
   * need three indices while faces and nodes only need two. This keeps the interface the same, but we will
   * only implement appropriate combinations.
   *
   * - Example 1 = getIndices(indices,ELEM,er,esr,ei,"pressure") = get pressure indices connected to this cell
   * - Example 2 = getIndices(indices,FACE,fi,"pressure") = get pressure indices connected to this face
   * - Example 3 = getIndices(indices,NODE,ni,"pressure") = get pressure indices connected to this node
   *
   * @param [out] indices globalIndex_array array gathering requested DoF indices.
   * @param [in]  connectivity Connectivity connectivity used to retrieve DoF indices.
   * @param [in]  region localIndex region of the location.
   * @param [in]  subregion localIndex subregion of the location.
   * @param [in]  index localIndex index of the location.
   * @param [in]  field Optional string field name for which DoFs are required.
   */
  void getIndices( globalIndex_array & indices,
                   Connectivity const connectivity,
                   localIndex const region,
                   localIndex const subregion,
                   localIndex const index,
                   string const & field = "" ) const;

  /**
   * @brief Get global indices for dofs connected by the connector type. We have two versions, since cells
   * need three indices while faces and nodes only need two. This keeps the interface the same, but we will
   * only implement appropriate combinations.
   *
   * - Example 1 = getIndices(indices,ELEM,er,esr,ei,"pressure") = get pressure indices connected to this cell
   * - Example 2 = getIndices(indices,FACE,fi,"pressure") = get pressure indices connected to this face
   * - Example 3 = getIndices(indices,NODE,ni,"pressure") = get pressure indices connected to this node
   *
   * @param [out] indices globalIndex_array array gathering requested DoF indices.
   * @param [in]  connectivity Connectivity connectivity used to retrieve DoF indices.
   * @param [in]  index localIndex index of the location.
   * @param [in]  field Optional string field name for which DoFs are required.
   */
  void getIndices( globalIndex_array & indices,
                   Connectivity const connectivity,
                   localIndex const index,
                   string const & field = "" ) const;

  /**
   * @brief Create the permutation that collects together all DoFs of each MPI process.
   *
   * @param [out] permutation ParallelMatrix the permutation matrix.
   */
  void createPermutation( ParallelMatrix & permutation ) const;

  /**
   * @brief Permute the GLOBAL sparsity pattern (location-location). Low level interface.
   *
   * @param [in]  locLocDistr ParallelMatrix the field-based sparsity pattern.
   * @param [in]  permutation ParallelMatrix the permutation matrix.
   * @param [out] permutedMatrix ParallelMatrix the permuted sparsity pattern.
   */
  void permuteSparsityPattern( ParallelMatrix const & locLocDistr,
                               ParallelMatrix const & permutation,
                               ParallelMatrix & permutedMatrix ) const;

  /**
   * @brief Copy values from DOFs to nodes.
   *
   * @param [in]    vector ParallelVector the input vector.
   * @param [in]    field string the name of the field.
   * @param [inout] manager dataRepository::ManagedGroup manager where values are stored.
   */
  void copyVectorToField( ParallelVector const & vector,
                          string const & field,
                          dataRepository::ManagedGroup * const manager ) const;

  /**
   * @brief Copy values from nodes to DOFs.
   *
   * @param [out] vector ParallelVector vector where values are stored.
   * @param [in]  field string the name of the field.
   * @param [in]  manager dataRepository::ManagedGroup the input manager.
   */
  void copyFieldToVector( ParallelVector & vector,
                          string const & field,
                          dataRepository::ManagedGroup const * const manager ) const;

  /**
   * @brief Print the global connectivity matrix.
   */
  void printConnectivityMatrix() const;

  /**
   * @brief Print the connectivity-location pattern for a specific field.
   *
   * @param [in] field string the name of the field.
   * @param [in] fileName Optional string the name of the output file (if empty, fileName is
   *                      formed based on field name).
   */
  void printConnectivityLocationPattern( string const & field, string const & fileName = "" ) const;

  /**
   * @brief Print a CSR pattern on file or on screen.
   *
   * @param [in] pattern Dof_SparsityPattern local sparsity pattern.
   * @param [in] fileName Optional string the name of the output file (empty means stdout).
   */
  void printSparsityPattern( Dof_SparsityPattern const & pattern, string const & fileName = "" ) const;

  /**
   * @brief Get m_doubleSync.
   */
  inline bool needDoubleSync() const {
    return m_doubleSync;
  }

private:
  /**
   * Verbosity level
   */
  localIndex m_verbosity = 0;

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
   * To fix the case when a processor handles just one layer of cells
   */
  bool m_doubleSync = false;

  /**
   * Field description
   */
  struct FieldDescription
  {
    FieldDescription()
    {}
    
    ~FieldDescription()
    {}

    FieldDescription( FieldDescription const & src ) :
      name( src.name ),
      regionNames( src.regionNames ),
      regionPtrs( src.regionPtrs ),
      location( src.location ),
      numComponents( src.numComponents ),
      key( src.key ),
      docstring( src.docstring ),
      numLocalNodes( src.numLocalNodes ),
      numLocalRows( src.numLocalRows ),
      numLocalConnectivity( src.numLocalConnectivity ),
      numGlobalRows( src.numGlobalRows ),
      firstLocalRow( src.firstLocalRow ),
      fieldOffset( src.fieldOffset ),
      firstLocalConnectivity( src.firstLocalConnectivity ),
      connLocPattern( src.connLocPattern )
    {}

    FieldDescription( FieldDescription && src ) :
      name( std::move( src.name ) ),
      regionNames( std::move( src.regionNames ) ),
      regionPtrs( std::move( src.regionPtrs ) ),
      location( std::move( src.location ) ),
      numComponents( std::move( src.numComponents ) ),
      key( std::move( src.key ) ),
      docstring( std::move( src.docstring ) ),
      numLocalNodes( std::move( src.numLocalNodes ) ),
      numLocalRows( std::move( src.numLocalRows ) ),
      numLocalConnectivity( std::move( src.numLocalConnectivity ) ),
      numGlobalRows( std::move( src.numGlobalRows ) ),
      firstLocalRow( std::move( src.firstLocalRow ) ),
      fieldOffset( std::move( src.fieldOffset ) ),
      firstLocalConnectivity( std::move( src.firstLocalConnectivity ) ),
      connLocPattern( std::move( src.connLocPattern ) )
    {}

    string name; //!< field name
    array1d<string> regionNames; //!< active element regions
    array1d<ElementRegion*> regionPtrs; //!< saved pointers to active regions
    Location location; //!< support location
    localIndex numComponents; //!< number of vector components
    string key; //!< string key for index array
    string docstring; //!< documentation string
    localIndex numLocalNodes; //!< number of local nodes
    localIndex numLocalRows; //!< number of local rows
    localIndex numLocalConnectivity; //!< number of local connectors
    globalIndex numGlobalRows; //!< number of ghost rows
    globalIndex firstLocalRow; //!< first row on this processor (without field offset)
    globalIndex fieldOffset; //!< global row offset for multi-field problems
    globalIndex firstLocalConnectivity; //!< first connector on this processor
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
   * Definition for entries of sparse matrices collection
   */
  struct matrixPair
  {
    ParallelMatrix* first;
    ParallelMatrix* second;
  };

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
   * Initialize data structure for connectivity and sparsity pattern
   */
  void initializeDataStructure();

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
                                           localIndex_array const & activeRegionsInput = localIndex_array() );

  /**
   * Create element index array
   */
  void createIndexArray_ElemVersion( FieldDescription & field ) const;

  /**
   * Create sparsity pattern for a field with itself (diagonal entries in the
   * connectivity matrix)
   */
  void addDiagSparsityPattern( Dof_SparsityPattern & connLocPatt,
                               localIndex const & fieldIdx,
                               Connectivity const connectivity,
                               localIndex_array const & activeRegionsInput = localIndex_array() );

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
                                    Connectivity const connectivity );

  /**
   * Definifion for entries of sparse matrix in COO format
   */
  typedef std::tuple<localIndex, globalIndex, localIndex> indexPair;

  /**
   * Compare structure used to create CSR matrix from COO format
   */
  struct pairComparison
  {
    inline bool operator()( const indexPair& lhs, const indexPair& rhs ) const
    {
      if( std::get<0>( lhs ) < std::get<0>( rhs ) )
        return true;
      else if( std::get<0>( lhs ) == std::get<0>( rhs ) )
        return std::get<1>( lhs ) < std::get<1>( rhs );
      else
        return false;
    }
  };

  /**
   * Compare second element of a pair
   */
  struct pairSecondComparison
  {
    inline bool operator()( const indexPair& lhs, const indexPair& rhs ) const
    {
      return ( std::get<1>( lhs ) < std::get<1>( rhs ) );
    }
  };

  /**
   * Convert a sparse matrix in COO format in the CSR version
   */
  void vectorOfPairsToCSR( array1d<indexPair> const & pairs,
                           localIndex const nRows,
                           localIndex const nCols,
                           Dof_SparsityPattern & pattern ) const;

  /**
   * Release internal storage
   */
  void cleanUp();
};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_LINEARALGEBRAINTERFACE_SRC_DOFMANAGER_HPP_ */
