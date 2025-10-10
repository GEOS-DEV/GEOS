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
 * @file CellGraphPartitioner.hpp
 */

#ifndef GEOS_MESH_MPICOMMUNICATIONS_CELLGRAPHPARTITIONER_HPP
#define GEOS_MESH_MPICOMMUNICATIONS_CELLGRAPHPARTITIONER_HPP

#include "MeshPartitioner.hpp"

namespace geos
{

/**
 * @class CellGraphPartitioner
 * @brief Partition mesh using cell dual graph
 *
 * Uses pure graph-based partitioning:
 * 1. Build dual graph (cells connected if they share face)
 * 2. Partition graph using ParMETIS/PTScotch
 */
class CellGraphPartitioner : public MeshPartitioner
{
public:

  CellGraphPartitioner( string const & name, Group * const parent );

  virtual ~CellGraphPartitioner() override;

  static string catalogName() { return "CellGraphPartitioner"; }

  virtual string getInfoString() const override;

protected:

  /**
   * @brief Compute partition assignment using cell graph
   *
   * @param mesh Input mesh
   * @param comm MPI communicator
   * @return Partition assignment for each element (3D + 2D)
   */
  virtual array1d< int64_t > computePartitioning( vtk::AllMeshes & mesh,
                                                  MPI_Comm const comm ) override;
};

} // namespace geos

#endif
