/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file DofManager.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_MULTISCALE_DOFMANAGER_HPP
#define GEOSX_LINEARALGEBRA_MULTISCALE_DOFMANAGER_HPP

#include "linearAlgebra/DofManager.hpp"

namespace geos
{
namespace multiscale
{

class MeshObjectManager;

/**
 * @brief Degree-of-freedom manager that works with multiscale mesh levels.
 *
 * Simplified implementation and capabilities compared to the geos::DofManager.
 */
class DofManager
{
public:

  /**
   * @brief Constructor.
   *
   * @param [in] name a unique name for this DoF manager
   */
  explicit DofManager( string name );

  /**
   * @brief Deleted copy constructor.
   */
  DofManager( DofManager const & ) = delete;

  /**
   * @brief Move constructor.
   */
  DofManager( DofManager && ) = default;

  /**
   * @brief Deleted copy assignment.
   * @return
   */
  DofManager & operator=( DofManager const & ) = delete;

  /**
   * @brief Defaulted move assignment.
   * @return
   */
  DofManager & operator=( DofManager && ) = default;

  /**
   * @brief Remove all previously added fields.
   */
  void clear();

  /**
   * @brief Assign a domain.
   * @param domain the target domain
   */
  void setDomain( DomainPartition & domain );

  /**
   * @brief Add a field.
   * @param fieldName name of the field
   * @param components number of components
   * @param manager the mesh manager the field is defined on
   */
  void addField( string const & fieldName,
                 integer const components,
                 MeshObjectManager & manager );

  /**
   * @brief Finish populating fields and apply appropriate dof renumbering.
   */
  void reorderByRank();

  /**
   * @brief @return the number of registered fields.
   */
  std::size_t numFields() const
  {
    return m_fields.size();
  }

  /**
   * @brief Check if string key is already being used.
   * @param name field key to check
   * @return flag true if exists
   */
  bool fieldExists( string const & name ) const;

  /**
   * @brief Return the key used to record the field in the DofManager.
   * @param [in] fieldName string the name of the field.
   * @return string indicating name of the field.
   */
  string const & key( string const & fieldName ) const;

  /**
   * @brief @return number of global dofs of a given field.
   * @param fieldName the name of the field
   */
  globalIndex numGlobalDofs( string const & fieldName ) const;

  /**
   * @brief @return total number of dofs across all fields and processors.
   */
  globalIndex numGlobalDofs() const;

  /**
   * @brief @return number of local dofs of a given field.
   * @param fieldName name of the field
   */
  localIndex numLocalDofs( string const & fieldName ) const;

  /**
   * @brief @return total number of dofs across all fields on current processor.
   */
  localIndex numLocalDofs() const;

  /**
   * @brief @return local offset of field's block on current processor in the system matrix.
   * @param [in] fieldName name of the field.
   */
  localIndex localOffset( string const & fieldName ) const;

  /**
   * @brief @return rank offset of given field, i.e. total number of dofs of this field on previous processors.
   * @param fieldName name of the field
   */
  globalIndex rankOffset( string const & fieldName ) const;

  /**
   * @brief @return rank offset of current processor, i.e. total number of dofs on previous processors.
   */
  globalIndex rankOffset() const;

  /**
   * @brief @return number of components in a given field.
   * @param fieldName name of the field
   */
  integer numComponents( string const & fieldName ) const;

  /**
   * @brief @return number of dof components across all fields.
   */
  integer numComponents() const;

  /**
   * @brief @return global offset of field's block on current processor in the system matrix.
   * @param [in] fieldName name of the field.
   */
  globalIndex globalOffset( string const & fieldName ) const;

  /**
   * @brief @return reference to the mesh object manager for field @p fieldName.
   * @param fieldName name name of the field.
   */
  MeshObjectManager const & manager( string const & fieldName ) const;

  /**
   * @brief @return the physical domain
   */
  DomainPartition & domain() { return *m_domain; }

  /**
   * @brief Create a matrix that restricts vectors and matrices to a subset of DOFs
   * @tparam MATRIX type of matrix used for restrictor
   * @param fieldName a fields to select
   * @param comm the MPI communicator to use in the operator
   * @param transpose if @p true, the transpose (prolongation) operator will be created
   * @param restrictor resulting operator
   *
   * @note Can only be called after reorderByRank(), since global DOF indexing is required
   *       for the restrictor to make sense.
   */
  template< typename MATRIX >
  void makeRestrictor( string const & fieldName,
                       MPI_Comm const & comm,
                       bool transpose,
                       MATRIX & restrictor ) const;

private:

  /**
   * Field description
   */
  struct FieldDescription
  {
    string name;                   ///< field name
    string key;                    ///< string key for index array
    MeshObjectManager * manager{}; ///< Pointer to mesh manager
    integer numComponents = 1;     ///< number of vector components
    localIndex numLocalDof = 0;    ///< number of local rows
    localIndex localOffset = 0;    ///< local offset of field on current processor (number of local DoFs preceding it)
    globalIndex numGlobalDof = 0;  ///< number of global rows
    globalIndex blockOffset = 0;   ///< offset of this field's block in a block-wise ordered system
    globalIndex rankOffset = 0;    ///< field's first DoF on current processor (within its block, ignoring other fields)
    globalIndex globalOffset = 0;  ///< global offset of field's DOFs on current processor for multi-field problems
  };

  /**
   * @brief Get field description from string key
   */
  FieldDescription const & getField( string const & name ) const;

  /**
   * @brief Create index array for the field
   * @param field the field descriptor
   */
  void createIndexArray( FieldDescription const & field );

  /// Name of the manager (unique, for unique identification of index array keys)
  string m_name;

  /// Pointer to domain
  DomainPartition * m_domain = nullptr;

  /// Array of field descriptions
  std::vector< FieldDescription > m_fields;

  /// Flag indicating that DOFs have been reordered rank-wise.
  bool m_reordered = false;
};

} // geosx
} // multiscale

#endif //GEOSX_LINEARALGEBRA_MULTISCALE_DOFMANAGER_HPP
