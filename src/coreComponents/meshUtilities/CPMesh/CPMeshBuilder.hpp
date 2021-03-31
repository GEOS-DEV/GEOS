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

/**
 * @class CPMeshBuilder
 * @brief This class is in charge of constructing a conforming mesh from an Eclipse corner-point mesh
 */
class CPMeshBuilder
{

public:

  /**
   * @brief Constructor.
   * @param name the name of the class
   */
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

/// using alias for templated Catalog CPMeshBuilder type
  using CatalogInterface = dataRepository::CatalogInterface< CPMeshBuilder,
                                                             string const & >;

  /**
   * @brief Getter for the singleton Catalog object
   * @return a static reference to the Catalog object
   *
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
  static string catalogName() { return "CPMeshBuilder"; }

  /**
   * @brief Const getter for catalog name of this class
   * @return the catalog name
   */
  virtual string getCatalogName() const { return catalogName(); }

  /**
   * @brief Driver function of the class: reads a corner-point mesh in the Eclipse format and constructs a conforming mesh
   * @param filePath the GRDECL file that contains the description of the corner-point mesh
   */
  void buildMesh( Path const & filePath );

  // vertices

  /**
   * @brief Const getter for the (unique) vertex positions
   * @return the 2D array of vertex positions
   */
  arrayView2d< real64 const > vertices() const
  { return m_meshVertices.m_vertices.toViewConst(); }

  /**
   * @brief Const getter for the map from corner-point vertex to (unique) vertex
   * @return the map from corner-point vertex to (unique) vertex
   */
  arrayView1d< localIndex const > cPVertexToVertex() const
  { return m_meshVertices.m_cPVertexToVertex.toViewConst(); }

  /**
   * @brief Const getter for the map from local corner-point vertex index to global corner-point vertex index
   * @return the map from local corner-point vertex index to global corner-point vertex index
   */
  arrayView1d< globalIndex const > cPVertexToGlobalCPVertex() const
  { return m_meshVertices.m_cPVertexToGlobalCPVertex.toViewConst(); }

  /**
   * @brief Const getter for the map from local vertex index to global vertex index
   * @return the map from local vertex index to global corner-point vertex index
   */
  arrayView1d< globalIndex const > vertexToGlobalVertex() const
  { return m_meshVertices.m_vertexToGlobalVertex.toViewConst(); }

  // cells

  /**
   * @brief Const getter for the map from active cell inside partition to active cell
   * @return the map from active cell inside partition to active cell
   */
  arrayView1d< localIndex const > activeCellInsidePartitionToActiveCell() const
  { return m_meshCells.m_activeCellInsidePartitionToActiveCell.toViewConst(); }

  /**
   * @brief Const getter for the map from active cell to cell
   * @return the map from active cell to cell
   */
  arrayView1d< localIndex const > activeCellToCell() const
  { return m_meshCells.m_activeCellToCell.toViewConst(); }

  /**
   * @brief Const getter for the map from cell to its first corner-point vertex
   * @return the map from cell to its first corner-point vertex
   */
  arrayView1d< localIndex const > cellToCPVertices() const
  { return m_meshCells.m_cellToCPVertices.toViewConst(); }

  /**
   * @brief Const getter for the map from local active cell index to global cell index
   * @return the map from from local active cell index to global cell index
   */
  arrayView1d< globalIndex const > activeCellToGlobalCell() const
  { return m_meshCells.m_activeCellToGlobalCell.toViewConst(); }

  // properties

  /**
   * @brief Const getter for the porosity field (includes inactive cells!)
   * @return an array of porosity values
   */
  arrayView1d< real64 const > porosityField() const { return m_meshParser.poro().toViewConst(); }

  /**
   * @brief Const getter for the porosity field (includes inactive cells!)
   * @return an array of permeability values
   */
  arrayView2d< real64 const > permeabilityField() const { return m_meshParser.perm().toViewConst(); }

  // MPI information

  /**
   * @brief Const getter for the topological mesh information
   * @return the topological mesh information encapsulated in an object of type CPMeshData
   */
  std::set< int > const & neighborsList() const { return m_meshComm.neighborsList(); }


private:

  /**
   * @brief Once the file is read, process the Eclipse keywords and fill the datastructure describing the conforming mesh
   * @details The post-processing proceeds in the following steps:
   *   1) For each corner-point cell, compute the position of the eight corner-point vertices
   *   2) Filter out duplicate vertices
   *   3) Construct faces (no implemented yet)
   */
  void postProcessMesh();

  /**
   * @brief Loop over corner-point cells to compute the position of the eight corner-point vertices in each cell
   */
  void buildCornerPointCells();

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
                                localIndex const iMinOverlap,
                                localIndex const iMaxOverlap,
                                localIndex const jMinOverlap,
                                localIndex const jMaxOverlap,
                                array1d< real64 > const & coord,
                                array1d< real64 > const & zcorn,
                                array1d< real64 > & xPos,
                                array1d< real64 > & yPos,
                                array1d< real64 > & zPos,
                                array1d< bool > & cPVertexIsInside ) const;

  /**
   * @brief Loop over the corner-point vertices and filter out duplicates.
   */
  void filterVertices();

  /**
   * @brief For the non-conforming case, build faces (not implemented yet)
   */
  void buildFaces();


  /// Object holding the mesh (and MPI partition) dimensions
  CPMeshDimensions m_meshDims;

  /// Object storing vertex data
  CPMeshVertices m_meshVertices;

  /// Object storing face data
  CPMeshFaces m_meshFaces;

  /// Object storing cell data
  CPMeshCells m_meshCells;

  /// Object in charge of parsing the GRDECL file
  CPMeshParser m_meshParser;

  /// Object holding the communication info
  CPMeshComm m_meshComm;

  /// Name of the mesh
  string m_meshName;

};

} // end namespace CPMesh

} // end namespace geosx

#endif //GEOSX_MESHUTILITIES_CPMESH_CPMESHBUILDER_HPP_
