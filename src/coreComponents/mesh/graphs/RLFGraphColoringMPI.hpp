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
 * @file RLFGraphColoringMPI.hpp
 */

#ifndef GEOS_GRAPH_RLFGRAPHCOLORINGMPI_HPP_
#define GEOS_GRAPH_RLFGRAPHCOLORINGMPI_HPP_

#include <vector>
#include <unordered_set>
#include <cstddef>
#include "common/DataTypes.hpp"
#include "GraphTools.hpp"
#include "GraphColoringBase.hpp"

namespace geos
{
namespace graph
{


class RLFGraphColoringMPI : public GraphColoringBase
{
public:
  RLFGraphColoringMPI( MPI_Comm comm = MPI_COMM_WORLD );
  ~RLFGraphColoringMPI();

  size_t getNumberOfColors( const std::vector< int > & colors ) const;
  size_t getNumberOfColors( const int color ) const;
  bool isColoringValid( const std::vector< camp::idx_t > & localAdjncy, const int localColor ) const;


  std::vector< int > colorGraph( const std::vector< camp::idx_t > & localXadj, const std::vector< camp::idx_t > & localAdjncy ) override;

  // Simplified version assuming one node per rank
  int colorGraph( const std::vector< camp::idx_t > & localAdjncy ) override;
};

} // namespace graph
} // namespace geos

#endif // GEOS_GRAPH_RLFGRAPHCOLORINGMPI_HPP_
