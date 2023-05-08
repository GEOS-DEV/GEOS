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
 * @file MeshLevel.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_MULTISCALE_MULTISCALEMESHLEVEL_HPP
#define GEOSX_LINEARALGEBRA_MULTISCALE_MULTISCALEMESHLEVEL_HPP

#include "common/DataTypes.hpp"
#include "linearAlgebra/DofManager.hpp"
#include "MeshObjectManager.hpp"
#include "linearAlgebra/utilities/LinearSolverParameters.hpp"

namespace geos
{

class DomainPartition;
class MeshLevel;

namespace multiscale
{

/**
 * @brief Multiscale mesh level.
 *
 * Mesh representation differs from the main code in several ways:
 * 1. Only stores part of the mesh on which corresponding physics problem is defined.
 * 2. Regions/subregions are erased, all elements numbered contiguously both locally and globally.
 * 3. Only elements and nodes are stored (no faces/edges).
 * 4. Element shapes are not defined (as coarse level elements are agglomerated).
 *
 */
class MeshLevel
{
public:

  /**
   * @brief Constructor.
   * @param name the level name
   */
  explicit MeshLevel( string const & name );

  /**
   * @brief @return level's cell manager
   */
  MeshObjectManager & cellManager() { return m_cellManager; }

  /**
   * @brief @return level's cell manager
   */
  MeshObjectManager const & cellManager() const { return m_cellManager; }

  /**
   * @brief @return level's node manager
   */
  MeshObjectManager & nodeManager() { return m_nodeManager; }

  /**
   * @brief @return level's node manager
   */
  MeshObjectManager const & nodeManager() const { return m_nodeManager; }

  /**
   * @brief Construct the finest level mesh representation at this level.
   * @param domain the physical domain
   * @param support the support domain (list of bodies/regions) for the physical field
   */
  void buildFineMesh( DomainPartition & domain,
                      Span< geos::DofManager::FieldSupport const > const support );

  /**
   * @brief Construct a coarse mesh representation at this level.
   * @param fineMesh previous (fine) level mesh
   * @param params coarsening parameters
   * @param boundaryNodeSets list of global domain boundary node set names
   */
  void buildCoarseMesh( multiscale::MeshLevel & fineMesh,
                        LinearSolverParameters::Multiscale::Coarsening const & params,
                        array1d< string > const & boundaryNodeSets );

  /**
   * @brief Push cell data to finer levels
   * @param fieldNames the list of field names to write
   * @param depth counter used to track recursion depth (leave default value of 0)
   */
  void writeCellData( std::vector< string > const & fieldNames, int depth = 0 ) const;

  /**
   * @brief Push node data to finer levels
   * @param fieldNames the list of field names to write
   * @param depth counter used to track recursion depth (leave default value of 0)
   */
  void writeNodeData( std::vector< string > const & fieldNames, int depth = 0 ) const;

  /**
   * @brief @return the level name
   */
  string const & name() const { return m_name; }

  /**
   * @brief @return pointer to the physical domain object
   */
  DomainPartition * domain() const { return m_domain; }

  /**
   * @brief @return pointer to parent (fine) mesh if one exists
   */
  multiscale::MeshLevel * fineMesh() const { return m_fineMesh; }

private:

  void writeCellDataFine( std::vector< string > const & fieldNames, int depth ) const;
  void writeCellDataCoarse( std::vector< string > const & fieldNames, int depth ) const;
  void writeNodeDataFine( std::vector< string > const & fieldNames, int depth ) const;
  void writeNodeDataCoarse( std::vector< string > const & fieldNames, int depth ) const;

  string m_name; ///< Unique name prefix

  conduit::Node m_rootNode;     ///< not used per se, but needed to satisfy Group constructor
  dataRepository::Group m_root; ///< needed for ObjectManagerBase constructor

  MeshObjectManager m_cellManager; ///< Cell manager
  MeshObjectManager m_nodeManager; ///< Node manager

  DomainPartition * m_domain{};                ///< Pointer to domain required to access communicators
  geos::DofManager::FieldSupport m_support{}; ///< List of regions in the source GEOSX mesh
  multiscale::MeshLevel * m_fineMesh{};        ///< Pointer to parent fine mesh (saved on each coarse level)
};

} // namespace multiscale
} // namespace geos

#endif //GEOSX_LINEARALGEBRA_MULTISCALE_MULTISCALEMESHLEVEL_HPP
