/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file PartitionerManager.hpp
 * @brief Manager for domain partitioning strategies
 */

#ifndef GEOS_PARTITIONER_PARTITIONERMANAGER_HPP_
#define GEOS_PARTITIONER_PARTITIONERMANAGER_HPP_

#include "dataRepository/Group.hpp"
#include "DomainPartitioner.hpp"

namespace geos
{

/**
 * @class PartitionerManager
 * @brief Manager for domain partitioning strategies
 *
 * PartitionerManager is a container for a single partitioner instance.
 *
 * The manager supports multiple partitioner types:
 * - Geometric partitioners (CartesianPartitioner, ParticleCartesianPartitioner)
 * - Mesh partitioners (CellGraphPartitioner, LayeredGraphPartitioner)
 *
 * Only one partitioner can be active at a time.
 *
 * XML Structure:
 * @code{.xml}
 * <Partitioner>
 *   <CellGraphPartitioner name="cellPartitioner"
 *                         engine="parmetis"
 *                         numRefinements="5" />
 * </Partitioner>
 * @endcode
 *
 * or
 *
 * @code{.xml}
 * <Partitioner>
 *   <CartesianPartitioner name="cartPartitioner"
 *                         partitionCounts="{ 4, 4, 2 }" />
 * </Partitioner>
 * @endcode
 *
 * Usage:
 * @code
 * PartitionerManager & partMgr = problemManager.getGroup<PartitionerManager>("Partitioner");
 *
 * if( partMgr.hasPartitioner() )
 * {
 *   DomainPartitioner & part = partMgr.getPartitioner();
 *
 *   // Query partitioner type
 *   if( auto* meshPart = dynamic_cast<MeshPartitioner*>(&part) )
 *   {
 *     // Use mesh partitioner
 *   }
 *   else if( auto* geomPart = dynamic_cast<GeometricPartitioner*>(&part) )
 *   {
 *     // Use geometric partitioner
 *   }
 * }
 * @endcode
 */
class PartitionerManager : public dataRepository::Group
{
public:

  /**
   * @brief Constructor
   * @param name The name of this manager (typically "Partitioner")
   * @param parent The parent group
   */
  PartitionerManager( string const & name,
                      dataRepository::Group * const parent );

  /**
   * @brief Destructor
   */
  ~PartitionerManager() override;

  /**
   * @brief Create a child partitioner from XML
   *
   * Called by the XML parser when it encounters a partitioner element.
   * Uses the DomainPartitioner catalog to create the appropriate concrete type.
   *
   * @param childKey The catalog key (partitioner type, e.g., "CellGraphPartitioner")
   * @param childName The instance name (e.g., "cellPartitioner")
   * @return Pointer to the created partitioner
   */
  virtual Group * createChild( string const & childKey,
                               string const & childName ) override;

  /**
   * @brief Expand catalogs for schema generation
   *
   * During schema generation, this creates one instance of each registered
   * partitioner type to document available options.
   */
  virtual void expandObjectCatalogs() override;

  /**
   * @brief Check if a partitioner has been defined
   * @return true if a partitioner exists
   */
  bool hasPartitioner() const;

  /**
   * @brief Get the active partitioner
   *
   * Returns the single partitioner defined in the XML.
   *
   * @return Reference to the partitioner
   */
  DomainPartitioner & getPartitioner();

  /**
   * @brief Get the active partitioner (const version)
   *
   * @return Const reference to the partitioner
   */
  DomainPartitioner const & getPartitioner() const;
};

} // namespace geos

#endif // GEOS_PARTITIONER_PARTITIONERMANAGER_HPP_
