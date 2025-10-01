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
 */

#ifndef GEOS_PARTITIONER_PARTITIONERMANAGER_HPP_
#define GEOS_PARTITIONER_PARTITIONERMANAGER_HPP_

#include "PartitionerBase.hpp"
#include "dataRepository/Group.hpp"


namespace geos
{

/**
 * @class PartitionerManager
 * @brief This class manages the partitioner objects in GEOSX (Cartesian, ParMetis, PTScotch, etc.)
 */
class PartitionerManager : public dataRepository::Group
{
public:

  /**
   * @brief Constructor for the PartitionerManager object.
   * @param[in] name the name of the PartitionerManager object in the repository
   * @param[in] parent the parent group of the PartitionerManager object being constructed
   */
  PartitionerManager( string const & name,
                      Group * const parent );

  virtual ~PartitionerManager() override;


  /**
   * @brief Create a new sub-mesh.
   * @param[in] childKey the key of the new object in the ObjectCatalog
   * @param[in] childName the name of the new object in the collection of sub-meshes
   * @return A pointer to the Group node in the dataRepository of the new object created
   */
  virtual Group * createChild( string const & childKey, string const & childName ) override;

  /// This function is used to expand any catalogs in the data structure
  virtual void expandObjectCatalogs() override;

  /**
   * @brief Check if a partitioner has been defined
   * @return True if a partitioner exists, false otherwise
   */
  bool hasPartitioner() const;

  /**
   * @brief Get the active partitioner (assumes only one partitioner is defined in XML)
   * @return Reference to the active partitioner
   */
  PartitionerBase & getPartitioner();

  /**
   * @brief Get the active partitioner (const version)
   * @return Const reference to the active partitioner
   */
  PartitionerBase const & getPartitioner() const;

private:

  /**
   * @brief Deleted default constructor of the PartitionerManager
   */
  PartitionerManager() = delete;

};

} // namespace geos

#endif // GEOS_PARTITIONER_PARTITIONERMANAGER_HPP_
