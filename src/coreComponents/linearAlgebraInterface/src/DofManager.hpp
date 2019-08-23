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
#include "InterfaceTypes.hpp"

namespace geosx
{

class DomainPartition;
class MeshLevel;
class ObjectManagerBase;

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
   * @brief Enumeration of geometric objects for support location. Note that this
   * enum is nearly identical to Connectivity, but we keep both for code readability
   * in function calls.
   */
  enum class Location
  {
    Elem, //!< location is element (like pressure in finite volumes)
    Face, //!< location is face (like flux in mixed finite elements)
    Edge, //!< location is edge (like flux between fracture elements)
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
    Edge, //!< connectivity is edge (like fracture element connectors)
    Node, //!< connectivity is node (like in finite volumes MPFA)
    None, //!< there is no connectivity (self connected field, like a mass matrix)
    USER_DEFINED //!< user defined connectivity (for input connectivity pattern)
  };

  /**
   * Field description
   */
  struct FieldDescription
  {
    FieldDescription()
      : numComponents( 1 )
    {}

    string name; //!< field name
    array1d<string> regionNames; //!< active element regions
    Location location; //!< support location
    localIndex numComponents; //!< number of vector components
    string key; //!< string key for index array
    string docstring; //!< documentation string
    localIndex numLocalNodes; //!< number of local nodes
    localIndex numLocalRows; //!< number of local rows
    globalIndex numGlobalRows; //!< number of global rows
    globalIndex firstLocalRow; //!< field's first row on current processor (in its block, not considering other fields)
    globalIndex fieldOffset; //!< global offset of field's DOFs on current processor for multi-field problems
  };

  /**
   * @brief Constructor.
   *
   * @param [in] name a unique name for this DoF manager
   * @param [in] verbosity Optional localIndex setting the verbosity level.
   *                       - 0: nothing (default)
   *                       - >0: minimal info
   */
  DofManager( string name, localIndex const verbosity = 0 );

  /**
   * @brief Destructor.
   */
  ~DofManager() = default;

  /**
   * Remove all fields
   */
  void clear();

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
   * @param [in] fieldName string the name of the field.
   * @param [in] location Location where it is defined.
   * @param [in] connectivity Connectivity through what it is connected.
   */
  void addField( string const & fieldName,
                 Location const location,
                 Connectivity const connectivity );

  /**
   * @brief Just another interface to allow four parameters (no regions, default is everywhere).
   *
   * @param [in] fieldName string the name of the field.
   * @param [in] location Location where it is defined.
   * @param [in] connectivity Connectivity through what it is connected.
   * @param [in] components localIndex number of components (for vector fields).
   */
  void addField( string const & fieldName,
                 Location const location,
                 Connectivity const connectivity,
                 localIndex const components );

  /**
   * @brief Just another interface to allow four parameters (no components, default is 1).
   *
   * @param [in] fieldName string the name of the field.
   * @param [in] location Location where it is defined.
   * @param [in] connectivity Connectivity through what it is connected.
   * @param [in] regions string_array where this field is defined.
   */
  void addField( string const & fieldName,
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
  void addField( string const & fieldName,
                 Location const location,
                 Connectivity const connectivity,
                 localIndex const components,
                 string_array const & regions );

  /**
   * @brief Interface to allow only two parameters.
   *
   * @param [in] fieldName string the name of the field.
   * @param [in] connLocInput ParallelMatrix input LC pattern.
   */
  void addField( string const & fieldName,
                 ParallelMatrix const & connLocInput );

  /**
   * @brief Just another interface to allow three parameters (no connectivity, default is USER_DEFINED).
   *
   * @param [in] fieldName string the name of the field.
   * @param [in] connLocInput ParallelMatrix input LC pattern.
   * @param [in] components localIndex number of components (for vector fields).
   */
  void addField( string const & fieldName,
                 ParallelMatrix const & connLocInput,
                 localIndex const components );

  /**
   * @brief Just another interface to allow three parameters (no components, default is 1).
   *
   * @param [in] fieldName string the name of the field.
   * @param [in] connLocInput ParallelMatrix input LC pattern.
   * @param [in] connectivity Connectivity through what it is connected.
   */
  void addField( string const & fieldName,
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
  void addField( string const & fieldName,
                 ParallelMatrix const & connLocInput,
                 localIndex const components,
                 Connectivity const connectivity );

  /**
   * @brief Just an interface to allow only three parameters.
   *
   * @param [in] rowFieldName string the name of the row field.
   * @param [in] colFieldName string the name of the col field.
   * @param [in] connectivity Connectivity through what they are connected.
   */
  void addCoupling( string const & rowFieldName,
                    string const & colFieldName,
                    Connectivity const connectivity );

  /**
   * @brief Just another interface to allow four parameters (no symmetry, default is true).
   *
   * @param [in] rowFieldName string the name of the row field.
   * @param [in] colFieldName string the name of the col field.
   * @param [in] connectivity Connectivity through what they are connected.
   * @param [in] regions string_array where this coupling is defined.
   */
  void addCoupling( string const & rowFieldName,
                    string const & colFieldName,
                    Connectivity const connectivity,
                    string_array const & regions );

  /**
   * @brief Just another interface to allow four parameters (no regions, default is everywhere).
   *
   * @param [in] rowFieldName string the name of the row field.
   * @param [in] colFieldName string the name of the col field.
   * @param [in] connectivity Connectivity through what they are connected.
   * @param [in] symmetric bool is it symmetric, i.e., both row-col and col-row?
   */
  void addCoupling( string const & rowFieldName,
                    string const & colFieldName,
                    Connectivity const connectivity,
                    bool const symmetric );

  /**
   * @brief Finish populating fields and apply appropriate dof renumbering
   *
   * This function must be called after all field and coupling information has been added.
   * It adjusts DoF index arrays to account for presence of other fields (in a global monolithic fashion).
   *
   * @note After DofManager has been closed, new fields and coupling cannot be added, until
   *       @ref clear or @ref setMesh is called.
   *
   * @note After close() is called, the meaning of FieldDescription::fieldOffset changes from
   *       "global offset of field's block in a global field-wise ordered (block) system" to
   *       "global offset of field's block on current processor in a rank-wise ordered system".
   *       This meaning is consistent with its use throughout. For example, this is the row/col
   *       global offset used to insert the field's sparsity block into a global coupled system.
   */
  void close();

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
   * @param [in] rowFieldName string the name of the row field.
   * @param [in] colFieldName string the name of the col field.
   * @param [in] connectivity Connectivity through what they are connected.
   * @param [in] regions string_array where this coupling is defined.
   * @param [in] symmetric bool is it symmetric, i.e., both row-col and col-row?
   */
  void addCoupling( string const & rowFieldName,
                    string const & colFieldName,
                    Connectivity const connectivity,
                    string_array const & regions,
                    bool const symmetric );

  /**
   * @brief Return the key used to record the field in the DofManager.
   *
   * @param [in] fieldName string the name of the field.
   */
  string getKey( string const & fieldName ) const;

  /**
   * @brief Return global number of dofs across all processors. If field argument is empty, return
   * monolithic size.
   *
   * @param [in] fieldName Optional string the name of the field.
   */
  globalIndex numGlobalDofs( string const & fieldName = "" ) const;

  /**
   * @brief Return local number of dofs on this processor. If field argument is empty, return
   * monolithic size.
   *
   * @param [in] fieldName Optional string the name of the field.
   */
  localIndex numLocalDofs( string const & fieldName = "" ) const;

  /**
   * @brief Return the sum of local dofs across all previous processors w.r.t. to the calling one for
   * the specified field.
   *
   * @param [in] fieldName Optional string the name of the field.
   */
  localIndex offsetLocalDofs( string const & fieldName = "" ) const;

  /**
   * @brief Set a sparsity pattern. Without additional arguments, this function provides the sparsity
   * pattern for the monolithic matrix. Sub-patterns can be extracted, however, using row and column
   * field keys.
   *
   * @param [out] locLocDistr ParallelMatrix the location-location sparsity pattern (LC*CL).
   * @param [in]  rowFieldName Optional string the name of the row field.
   * @param [in]  colFieldName Optional string the name of the col field.
   */
  void setSparsityPattern( ParallelMatrix & locLocDistr,
                           string const & rowFieldName = "",
                           string const & colFieldName = "" ) const;

  /**
   * @brief Set a sparsity pattern. Low level version.
   *
   * @param [out] matrix ParallelMatrix the location-location sparsity pattern (LC*CL).
   * @param [in]  rowFieldIndex localIndex row field index (-1 means all fields).
   * @param [in]  colFieldIndex localIndex col field index (-1 means all fields).
   */
  void setSparsityPattern( ParallelMatrix & matrix,
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
                  string const & fieldName = "" ) const;

  /**
   * @brief Allocate a vector. Low level version.
   *
   * @param [out] vector ParallelVector the output vector.
   * @param [in]  rowFieldIndex localIndex row field index (-1 means all fields).
   * @param [in]  colFieldIndex localIndex col field index (-1 means all fields).
   */
  void setVector( ParallelVector & vector,
                  localIndex const fieldIndex ) const;

  /**
   * @brief Copy values from DOFs to nodes.
   *
   * @tparam FIELD_OP operation to perform (see FieldSpecificationOps.hpp)
   * @tparam POLICY execution policy for the kernel
   * @param vector source LA vector
   * @param srcFieldName name of the source field (as defined in DofManager)
   * @param scalingFactor a factor to scale vector values by
   * @param manager mesh object manager that contains the target field (subregion for elements)
   * @param dstFieldName name of the destination field (view wrapper key on the manager)
   * @param loCompIndex index of starting DoF component (for partial copy)
   * @param hiCompIndex index past the ending DoF component (for partial copy)
   *
   * @note [@p loCompIndex , @p hiCompIndex) form a half-open interval.
   *       Negative value of @p hiCompIndex means use full number of field components
   */
  void copyVectorToField( ParallelVector const & vector,
                          string const & srcFieldName,
                          real64 const scalingFactor,
                          ObjectManagerBase * const manager,
                          string const & dstFieldName,
                          localIndex const loCompIndex = 0,
                          localIndex const hiCompIndex = -1 ) const;

  /**
   * @brief Add values from DOFs to nodes.
   *
   * @tparam FIELD_OP operation to perform (see FieldSpecificationOps.hpp)
   * @tparam POLICY execution policy for the kernel
   * @param vector source LA vector
   * @param srcFieldName name of the source field (as defined in DofManager)
   * @param scalingFactor a factor to scale vector values by
   * @param manager mesh object manager that contains the target field (subregion for elements)
   * @param dstFieldName name of the destination field (view wrapper key on the manager)
   * @param loCompIndex index of starting DoF component (for partial copy)
   * @param hiCompIndex index past the ending DoF component (for partial copy)
   *
   * @note [@p loCompIndex , @p hiCompIndex) form a half-open interval.
   *       Negative value of @p hiCompIndex means use full number of field components
   */
  void addVectorToField( ParallelVector const & vector,
                         string const & srcFieldName,
                         real64 const scalingFactor,
                         ObjectManagerBase * const manager,
                         string const & dstFieldName,
                         localIndex const loCompIndex = 0,
                         localIndex const hiCompIndex = -1 ) const;

  /**
   * @brief Copy values from nodes to DOFs.
   *
   * @tparam FIELD_OP operation to perform (see FieldSpecificationOps.hpp)
   * @tparam POLICY execution policy for the kernel
   * @param manager mesh object manager that contains the target field (subregion for elements)
   * @param srcFieldName name of the source field (view wrapper key on the manager)
   * @param scalingFactor a factor to scale vector values by
   * @param vector target LA vector
   * @param dstFieldName name of the destination field (as defined in DofManager)
   * @param loCompIndex index of starting DoF component (for partial copy)
   * @param hiCompIndex index past the ending DoF component (for partial copy)
   *
   * @note [@p loCompIndex , @p hiCompIndex) form a half-open interval.
   *       Negative value of @p hiCompIndex means use full number of field components
   */
  void copyFieldToVector( ObjectManagerBase const * const manager,
                          string const & srcFieldName,
                          real64 const scalingFactor,
                          ParallelVector & vector,
                          string const & dstFieldName,
                          localIndex const loCompIndex = 0,
                          localIndex const hiCompIndex = -1 ) const;

  /**
   * @brief Add values from nodes to DOFs.
   *
   * @tparam FIELD_OP operation to perform (see FieldSpecificationOps.hpp)
   * @tparam POLICY execution policy for the kernel
   * @param manager mesh object manager that contains the target field (subregion for elements)
   * @param srcFieldName name of the source field (view wrapper key on the manager)
   * @param scalingFactor a factor to scale vector values by
   * @param vector target LA vector
   * @param dstFieldName name of the destination field (as defined in DofManager)
   * @param loCompIndex index of starting DoF component (for partial copy)
   * @param hiCompIndex index past the ending DoF component (for partial copy)
   *
   * @note [@p loCompIndex , @p hiCompIndex) form a half-open interval.
   *       Negative value of @p hiCompIndex means use full number of field components
   */
  void addFieldToVector( ObjectManagerBase const * const manager,
                         string const & srcFieldName,
                         real64 const scalingFactor,
                         ParallelVector & vector,
                         string const & dstFieldName,
                         localIndex const loCompIndex = 0,
                         localIndex const hiCompIndex = -1 ) const;

  /**
   * @brief Print the global connectivity matrix.
   */
  void printConnectivityMatrix( std::ostream & os = std::cout ) const;

  /**
   * @brief Print the connectivity-location pattern for a specific field.
   *
   * @param [in] fieldName string the name of the field.
   * @param [in] fileName Optional string the name of the output file (if empty, fileName is
   *                      formed based on field name).
   */
  void printConnectivityLocationPattern( string const & fieldName, string const & fileName = "" ) const;

private:

  /**
   * @brief Initialize data structure for connectivity and sparsity pattern
   */
  void initializeDataStructure();

  /**
   * @brief Check if string key is already being used
   */
  bool keyInUse( string const & key ) const;

  /**
   * @brief Get field index from string key
   */
  localIndex getFieldIndex( string const & key ) const;

  /**
   * @brief Create index array for the field
   */
  template< typename ... SUBREGIONTYPES >
  void createIndexArray( FieldDescription & field );

  /**
   * @brief Remove an index array for the field
   */
  template< typename ... SUBREGIONTYPES >
  void removeIndexArray( FieldDescription const & field );

  /**
   * @brief Create a connector-location sparsity pattern for a field
   */
  void makeConnLocPattern( FieldDescription const & fieldDesc,
                           Connectivity const connectivity,
                           array1d <string> const & regions,
                           ParallelMatrix & connLocPattern );

  /**
   * @brief Populate the sparsity pattern for a coupling block between given fields.
   * @param locLocDistr the sparsity to be filled
   * @param rowFieldIndex index of row field (must be non-negative)
   * @param colFieldIndex index of col field (must be non-negative)
   *
   * This private function is used as a building block by higher-level SetSparsityPattern()
   */
  void setSparsityPatternOneBlock( ParallelMatrix & locLocDistr,
                                   localIndex const rowFieldIndex,
                                   localIndex const colFieldIndex ) const;

  /**
   * @brief Generic implementation for @ref copyVectorToField and @ref addVectorToField
   * @tparam FIELD_OP operation to perform (see FieldSpecificationOps.hpp)
   * @tparam POLICY execution policy for the kernel
   * @param vector source LA vector
   * @param srcFieldName name of the source field (as defined in DofManager)
   * @param scalingFactor a factor to scale vector values by
   * @param manager mesh object manager that contains the target field (subregion for elements)
   * @param dstFieldName name of the destination field (view wrapper key on the manager)
   * @param loCompIndex index of starting DoF component (for partial copy)
   * @param hiCompIndex index past the ending DoF component (for partial copy)
   *
   * @note [@p loCompIndex , @p hiCompIndex) form a half-open interval.
   *       Negative value of @p hiCompIndex means use full number of field components
   */
  template< typename FIELD_OP, typename POLICY >
  void vectorToField( ParallelVector const & vector,
                      string const & srcFieldName,
                      real64 const scalingFactor,
                      ObjectManagerBase * const manager,
                      string const & dstFieldName,
                      localIndex const loCompIndex,
                      localIndex const hiCompIndex ) const;

  /**
   * @brief Generic implementation for @ref copyFieldToVector and @ref addFieldToVector
   * @tparam FIELD_OP operation to perform (see FieldSpecificationOps.hpp)
   * @tparam POLICY execution policy for the kernel
   * @param manager mesh object manager that contains the target field (subregion for elements)
   * @param srcFieldName name of the source field (view wrapper key on the manager)
   * @param scalingFactor a factor to scale vector values by
   * @param vector target LA vector
   * @param dstFieldName name of the destination field (as defined in DofManager)
   * @param loCompIndex index of starting DoF component (for partial copy)
   * @param hiCompIndex index past the ending DoF component (for partial copy)
   *
   * @note [@p loCompIndex , @p hiCompIndex) form a half-open interval.
   *       Negative value of @p hiCompIndex means use full number of field components
   */
  template< typename FIELD_OP, typename POLICY >
  void fieldToVector( ObjectManagerBase const * const manager,
                      string const & srcFieldName,
                      real64 const scalingFactor,
                      ParallelVector & vector,
                      string const & dstFieldName,
                      localIndex const loCompIndex,
                      localIndex const hiCompIndex ) const;

  /**
   * Name of the manager
   */
  string m_name;

  /**
   * Verbosity level
   */
  localIndex m_verbosity;

  /**
   *  Limit on max number of fields
   */
  static localIndex constexpr MAX_NUM_FIELDS = 10;

  /**
   * Pointer to domain manager
   */
  DomainPartition * m_domain = nullptr;

  /**
   * Pointer to corresponding MeshLevel
   */
  MeshLevel * m_mesh = nullptr;

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
    std::unique_ptr<ParallelMatrix> first;
    std::unique_ptr<ParallelMatrix> second;
  };

  /**
   * Table of sparsity patterns within and between fields
   */
  array2d<matrixPair> m_sparsityPattern;

  /**
   * Indicates that the manager is closed for adding new fields
   */
  bool m_closed;
};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_LINEARALGEBRAINTERFACE_SRC_DOFMANAGER_HPP_ */
