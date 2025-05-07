
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
 * @file CartesianPartitioner.hpp
 */

#ifndef GEOS_PARTITIONER_CARTESIANPARTITIONER_HPP_
#define GEOS_PARTITIONER_CARTESIANPARTITIONER_HPP_

#include "PartitionerBase.hpp"

#include "common/DataTypes.hpp"

#include <mpi.h>



constexpr int nsdof = 3;
namespace geos
{



class CartesianPartitioner : public PartitionerBase {
public:
    void partition() override;
    std::vector<int> getNeighbors() const override;
    int getColor() const override;



 /**
   * @brief Boolean like array of length 3 (space dimensions).
   *
   * 1 means periodic.
   */
  array1d< int > m_Periodic;
  /// ijk partition indexes
  //array1d< int > m_coords;

  /**
   * @brief Recursively builds neighbors if an MPI cartesian topology is used 
   * @param idim Dimension index in the cartesian.
   * @param cartcomm Communicator with cartesian structure.
   * @param ncoords Cartesian coordinates of a process (assumed to be of length 3).
   *
   */
  void addNeighbors( const unsigned int idim,
                     MPI_Comm & cartcomm,
                     int * ncoords );





};


} // namespace geos

#endif // GEOS_PARTITIONER_CARTESIANPARTITIONER_HPP_
