/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file CPMeshComm.hpp
 */

#ifndef GEOSX_MESHUTILITIES_CPMESH_CPMESHCOMM_HPP_
#define GEOSX_MESHUTILITIES_CPMESH_CPMESHCOMM_HPP_

#include "common/DataTypes.hpp"
#include "dataRepository/ObjectCatalog.hpp"

namespace geosx
{

namespace CPMesh
{

class CPMeshComm
{

public:

  CPMeshComm( string const & name );

  /// Default virtual destructor
  virtual ~CPMeshComm() = default;

  /// Default copy constructor
  CPMeshComm( CPMeshComm const & ) = default;

  /// Default move constructor
  CPMeshComm( CPMeshComm && ) = default;

  /// Deleted copy assignment operator
  CPMeshComm & operator=( CPMeshComm const & ) = delete;

  /// Deleted move assignment operator
  CPMeshComm & operator=( CPMeshComm && ) = delete;

  using CatalogInterface = dataRepository::CatalogInterface< CPMeshComm,
                                                             string const & >;
  static typename CatalogInterface::CatalogType & getCatalog()
  {
    static CatalogInterface::CatalogType catalog;
    return catalog;
  }

  static string catalogName() { return "CPMeshComm"; }

  virtual string getCatalogName() const { return catalogName(); }

  /**
   * @brief Rank 0 partitions the domain, and the communication pattern is set up
   * @details This involves multiple steps:
   *  1) Rank 0 broadcasts the full domain boundaries (nX, nY, nZ)
   *  2) Rank 0 partitions the domain in boxes using a simple approach
   *  3) Rank 0 scatters the partition boundaries
   *  4) Each rank uses the structure of the mesh to figure out the boundary vertices it owns, and whom it is communicating to
   */
  void setupMPIPartition( CPMeshData & cPMeshData ) const;

  /**
   * @brief Each rank gathers the number of cells/vertices owned by the other ranks
   */
  localIndex gatherLocalCountForEachRank( localIndex const countOnMyRank ) const;

  /**
   * @brief The global indices of the vertices located on the boundaries of each MPI partition are synchronized
   * @details Each rank does the following:
   *  - It receives the global indices of the boundary vertices owned by other ranks
   *  - It sends the global indices of the boundary vertices that it owns
   */
  void synchronizeBoundaryVertices( CPMeshData & m_cPMeshData ) const;

private:

  /**
   * @brief Rank 0 (that has read the number of cells) broadcasts full domain boundaries (nX, nY, nZ)
   */
  void broadcastDomainBoundaries( localIndex const myRank,
                                  CPMeshData & cPMeshData ) const;

  /**
   * @brief Rank 0 partitions the full domain and scatters the MPI partition boundaries
   */
  void scatterPartitionBoundaries( localIndex const myRank,
                                   localIndex const size,
                                   CPMeshData & cPMeshData ) const;

  /**
   * @brief Each rank loops over the boundaries of its MPI partition, and collects the vertices it owns
   * @details The determination of vertex ownership is done using the ijk logic. For each pillar, we can
   * explicitly compute the (i,j,k) indices of the surrounding cells. The pillar is assigned to the rank
   * (owning one of the surrounding cells) with the smallest index.
   */
  void fillOwnershipArrays( CPMeshData & cPMeshData ) const;

  /// Name of the mesh
  string m_meshName;

};

} // end namespace CPMesh

} // end namespace geosx

#endif //GEOSX_MESHUTILITIES_CPMESH_CPMESHCOMM_HPP_
