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
#include "meshUtilities/CPMesh/CPMeshData.hpp"

namespace geosx
{

namespace CPMesh
{

/**
 * @class CPMeshComm
 * @brief This class is in charge of MPI communications during the construction of the mesh
 */
class CPMeshComm
{

public:

  /**
   * @brief Constructor.
   * @param name the name of the class
   */
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

  /// using alias for templated Catalog CPMeshBuilder type
  using CatalogInterface = dataRepository::CatalogInterface< CPMeshComm,
                                                             string const & >;

  /**
   * @brief Getter for the singleton Catalog object
   * @return a static reference to the Catalog object
   */
  static typename CatalogInterface::CatalogType & getCatalog()
  {
    static CatalogInterface::CatalogType catalog;
    return catalog;
  }

  /**
   * @brief Define the catalog name for this class
   * @return the catalog name
   */
  static string catalogName() { return "CPMeshComm"; }

  /**
   * @brief Const getter for catalog name of this class
   * @return the catalog name
   */
  virtual string getCatalogName() const { return catalogName(); }

  /**
   * @brief Const getter for the list of my neighbor ranks
   * @return a set containing the indices of neighbors
   */
  std::set< int > const & neighborsList() const { return m_neighborsList; }

  /**
   * @brief Rank 0 partitions the domain, and the communication pattern is set up
   * @param[out] meshDims the dimensions of the problem / the MPI partition
   * @details This involves multiple steps:
   *  1) Rank 0 broadcasts the full domain boundaries (nX, nY, nZ)
   *  2) Rank 0 partitions the domain in boxes using a simple approach
   *  3) Rank 0 scatters the partition boundaries
   *  4) Rank 0 uses the structure of the mesh to figure out the communication pattern
   */
  void setupMPIPartition( CPMeshDimensions & meshDims );

private:

  /**
   * @brief Rank 0 (that has read the number of cells) broadcasts full domain dimensions (nX, nY, nZ)
   * @param[inout] meshDims the struct describing the local and global mesh dimensions
   */
  void broadcastDomainDimensions( CPMeshDimensions & meshDims ) const;

  /**
   * @brief Rank 0 partitions the full domain and scatters the MPI partition boundaries
   * @param[inout] meshDims the struct describing the local and global mesh dimensions
   */
  void scatterPartitionBoundaries( CPMeshDimensions & meshDims );


  /// Index of the main rank
  localIndex m_mainRank;

  /// Index of my rank
  localIndex m_myRank;

  /// Number of ranks
  localIndex m_size;

  /// Indices of my neighbors
  std::set< int > m_neighborsList;

  /// Name of the mesh
  string m_meshName;

};

} // end namespace CPMesh

} // end namespace geosx

#endif //GEOSX_MESHUTILITIES_CPMESH_CPMESHCOMM_HPP_
