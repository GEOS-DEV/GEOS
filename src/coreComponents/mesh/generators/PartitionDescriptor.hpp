/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file ParrtitionDescriptor.hpp
 */

#ifndef GEOSX_MESH_PARTITIONDESCRIPTOR_H_
#define GEOSX_MESH_PARTITIONDESCRIPTOR_H_

#include <set>
#include "mesh/mpiCommunications/SpatialPartition.hpp"

namespace geosx
{

/**
 * @class PartitionDescriptor
 * @brief Simple utility to retrieve partition information in case of Metis or Spatial partition.
 */
class PartitionDescriptor
{
public:
  /**
   * @brief indicate if the partition is described using a Metis Neighbor list.
   */
  bool hasMetisNeighborList() const { return m_hasMetisNeighborList; }
  /**
   * @brief Sets the boolean that indicates if the partition is described using a Metis Neighbor list.
   */
  void setHasMetisNeighborList( bool value ) { m_hasMetisNeighborList = value; }
  /**
   * @brief Gets a reference to the list of metis neighbor list.
   */
  std::set< int > const & getMetisNeighborList() const { return m_metisNeighborList; }

  /**
   * @brief Sets the list of metis neighbor list.
   */
  void setMetisNeighborList( std::set< int > const & metisNeighborList ) { m_metisNeighborList = metisNeighborList; }
  template< class InputIt >
  void setMetisNeighborList( InputIt begin, InputIt end ) { m_metisNeighborList.insert( begin, end ); }

  /**
   * @brief indicate if the partition is described using a spatial partition.
   */
  bool hasSpatialParition() const { return !m_hasMetisNeighborList; }
  /**
   * @brief Sets the boolean that indicates if the partition is described using a Metis Neighbor list.
   */
  void setHasSpatialPartition( bool value ) { m_hasMetisNeighborList = !value; }
  /**
   * @brief Returns a reference to the spatialPartition
   */
  SpatialPartition const & getSpatialPartition() const { return m_spatialPartition; }

  /**
   * @brief Sets the spatialPartition
   */
  void setSpatialPartition( SpatialPartition const & spatialPartition ) { m_spatialPartition = spatialPartition; }
private:
  bool m_hasMetisNeighborList;         //< Indicate if we use metis neighbor list or spatial partition to describe the partition
  std::set< int > m_metisNeighborList;   //< The list of neighbors computed wwith metis
  SpatialPartition m_spatialPartition; //< The spatial partition
};

}
#endif /* GEOSX_MESH_PARTITIONDESCRIPTOR_H_ */
