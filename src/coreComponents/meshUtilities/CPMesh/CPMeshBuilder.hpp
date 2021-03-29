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
 * @file CPMeshBuilder.hpp
 */

#ifndef GEOSX_MESHUTILITIES_CPMESH_CPMESHBUILDER_HPP_
#define GEOSX_MESHUTILITIES_CPMESH_CPMESHBUILDER_HPP_

#include "dataRepository/ObjectCatalog.hpp"
#include "meshUtilities/CPMesh/CPMeshParser.hpp"
#include "meshUtilities/CPMesh/CPMeshData.hpp"
#include "meshUtilities/CPMesh/CPMeshComm.hpp"

namespace geosx
{

namespace CPMesh
{

class CPMeshBuilder
{

public:

  CPMeshBuilder( string const & name );

  /// Default virtual destructor
  virtual ~CPMeshBuilder() = default;

  /// Default copy constructor
  CPMeshBuilder( CPMeshBuilder const & ) = default;

  /// Default move constructor
  CPMeshBuilder( CPMeshBuilder && ) = default;

  /// Deleted copy assignment operator
  CPMeshBuilder & operator=( CPMeshBuilder const & ) = delete;

  /// Deleted move assignment operator
  CPMeshBuilder & operator=( CPMeshBuilder && ) = delete;

  using CatalogInterface = dataRepository::CatalogInterface< CPMeshBuilder,
                                                             string const & >;
  static typename CatalogInterface::CatalogType & getCatalog()
  {
    static CatalogInterface::CatalogType catalog;
    return catalog;
  }

  static string catalogName() { return "CPMeshBuilder"; }

  virtual string getCatalogName() const { return catalogName(); }

  /**
   * @brief Driver function of the class: reads a corner-point mesh in the Eclipse format and constructs a conforming mesh
   * @param filePath the GRDECL file that contains the description of the corner-point mesh
   */
  void buildMesh( Path const & filePath );

private:

  /**
   * @brief Once the file is read, process the Eclipse keywords and fill the datastructure describing the conforming mesh
   * @details The post-processing proceeds in the following steps:
   *   1) For each corner-point cell, compute the position of the eight corner-point vertices
   *   2) Filter out duplicate vertices
   *   3) Assign global cell indices
   *   4) Assign global vertex indices
   *   5) Construct faces
   */
  void postProcessMesh();

  /**
   * @brief Loop over corner-point cells to compute the position of the eight corner-point vertices in each cell
   */
  void buildHexahedralMesh();

  /**
   * @brief Using ZCORN and COORD, compute the position of the eight corner-point vertices in cell (i,j,k)
   * @param i index in the x-direction
   * @param j index in the y-direction
   * @param k index in the z-direction
   * @param nXLocal number of cells in the x-direction on my rank
   * @param nYLocal number of cells in the y-direction on my rank
   * @param coord local content of the COORD keyword
   * @param zcorn local content of the ZCORN keyword
   * @param xPos x-position of the eight corner-point vertices
   * @param yPos y-position of the eight corner-point vertices
   * @param zPos z-position of the eight corner-point vertices
   */
  bool processActiveHexahedron( localIndex const i,
                                localIndex const j,
                                localIndex const k,
                                localIndex const nXLocal,
                                localIndex const nYLocal,
                                array1d< real64 > const & coord,
                                array1d< real64 > const & zcorn,
                                array1d< real64 > & xPos,
                                array1d< real64 > & yPos,
                                array1d< real64 > & zPos ) const;

  /**
   * @brief Loop over the corner-point vertices and filter out duplicates. Assign a local index to each unique vertex
   */
  void filterVertices();

  /**
   * @brief Assign global cell indices
   * @details This is done in two steps:
   *  1) Gather the number of cells owned by each rank
   *  2) Using the offset computed at step 1), each rank computes the global index of each local cell
   */
  void assignGlobalActiveCellIndices();

  /**
   * @brief Assign global cell indices
   * @details This is done in three steps:
   *  1) Gather the number of vertices owned by each rank
   *  2) Using the offset computed at step 1), each rank computes the global index of each local vertex
   *  3) Synchronize global indices on the boundaries of the MPI partitions
   */
  void assignGlobalVertexIndices();

  /**
   * @brief For the non-conforming case, build faces (not implemented yet)
   */
  void buildFaces();

  /// Object in charge of parsing the GRDECL file
  CPMeshParser m_cPMeshParser;

  /// Object holding the data local to each rank
  CPMeshData m_cPMeshData;

  /// Object holding the communication info
  CPMeshComm m_cPMeshComm;

  /// Name of the mesh
  string m_meshName;

};

} // end namespace CPMesh

} // end namespace geosx

#endif //GEOSX_MESHUTILITIES_CPMESH_CPMESHBUILDER_HPP_
