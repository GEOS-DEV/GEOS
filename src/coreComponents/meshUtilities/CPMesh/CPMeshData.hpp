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
 * @file CPMeshData.hpp
 */


#ifndef GEOSX_MESHUTILITIES_CPMESH_CPMESHDATA_HPP_
#define GEOSX_MESHUTILITIES_CPMESH_CPMESHDATA_HPP_

#include "common/DataTypes.hpp"
#include "dataRepository/ObjectCatalog.hpp"

namespace geosx
{

namespace CPMesh
{

class CPMeshData
{

public:

  CPMeshData( string const & name );

  /// Default virtual destructor
  virtual ~CPMeshData() = default;

  /// Default copy constructor
  CPMeshData( CPMeshData const & ) = default;

  /// Default move constructor
  CPMeshData( CPMeshData && ) = default;

  /// Deleted copy assignment operator
  CPMeshData & operator=( CPMeshData const & ) = delete;

  /// Deleted move assignment operator
  CPMeshData & operator=( CPMeshData && ) = delete;

  using CatalogInterface = dataRepository::CatalogInterface< CPMeshData,
                                                             string const & >;
  static typename CatalogInterface::CatalogType & getCatalog()
  {
    static CatalogInterface::CatalogType catalog;
    return catalog;
  }

  static string catalogName() { return "CPMeshData"; }

  virtual string getCatalogName() const { return catalogName(); }

  void defineDomainBoundaries( localIndex const nX, localIndex const nY, localIndex const nZ );

  void definePartitionBoundaries( localIndex const iMin, localIndex const jMin,
                                  localIndex const iMax, localIndex const jMax );

  void definePartitionOverlaps( localIndex const iMinOverlap, localIndex const jMinOverlap,
                                localIndex const iMaxOverlap, localIndex const jMaxOverlap );

  // For now, I put all the arrays in this class for convenience
  // Later, I may split them into different classes (ownership vectors to CPMeshComm, keywords vectors to CPMeshParsers)
  // Also, we may not need all these getters, they are really cumbersome (use public members instead)

  localIndex nX() const { return m_nX; }
  localIndex nY() const { return m_nY; }
  localIndex nZ() const { return m_nZ; }

  localIndex nXLocal() const { return (m_iMax+m_iMaxOverlap)-(m_iMin-m_iMinOverlap)+1; }
  localIndex nYLocal() const { return (m_jMax+m_jMaxOverlap)-(m_jMin-m_jMinOverlap)+1; }
  localIndex nZLocal() const { return m_nZ; }

  localIndex iMinLocal() const { return m_iMin-m_iMinOverlap; }
  localIndex jMinLocal() const { return m_jMin-m_jMinOverlap; }
  localIndex iMaxLocal() const { return m_iMax+m_iMaxOverlap; }
  localIndex jMaxLocal() const { return m_jMax+m_jMaxOverlap; }

  localIndex iMinOwned() const { return m_iMin; }
  localIndex jMinOwned() const { return m_jMin; }
  localIndex iMaxOwned() const { return m_iMax; }
  localIndex jMaxOwned() const { return m_jMax; }



  array1d< real64 > & coord() { return m_coord; }
  array1d< real64 > const & coord() const { return m_coord; }

  array1d< real64 > & zcorn() { return m_zcorn; }
  array1d< real64 > const & zcorn() const { return m_zcorn; }

  array1d< localIndex > & actnum() { return m_actnum; }
  array1d< localIndex > const & actnum() const { return m_actnum; }



  localIndex & nOwnedVertices() { return m_nOwnedVertices; }
  localIndex const & nOwnedVertices() const { return m_nOwnedVertices; }

  localIndex & nLocalVertices() { return m_nLocalVertices; }
  localIndex nLocalVertices() const { return m_nLocalVertices; }

  array2d< real64 > & cPVertices() { return m_cPVertices; }
  array2d< real64 > const & cPVertices() const { return m_cPVertices; }

  array1d< bool > & cPVertexIsOwned() { return m_cPVertexIsOwned; }
  array1d< bool > const & cPVertexIsOwned() const { return m_cPVertexIsOwned; }

  array2d< real64 > & vertices() { return m_vertices; }
  array2d< real64 > const & vertices() const { return m_vertices; }

  array1d< bool > & vertexIsOwned() { return m_vertexIsOwned; }
  array1d< bool > const & vertexIsOwned() const { return m_vertexIsOwned; }

  array1d< localIndex > & localCPVertexToLocalVertex() { return m_localCPVertexToLocalVertex; }
  array1d< localIndex > const & localCPVertexToLocalVertex() const { return m_localCPVertexToLocalVertex; }

  array1d< globalIndex > & localVertexToGlobalVertex() { return m_localVertexToGlobalVertex; }
  array1d< globalIndex > const & localVertexLocalToGlobalVertex() const { return m_localVertexToGlobalVertex; }



  localIndex & nOwnedActiveCells() { return m_nOwnedActiveCells; }
  localIndex nOwnedActiveCells() const { return m_nOwnedActiveCells; }

  localIndex & nLocalActiveCells() { return m_nLocalActiveCells; }
  localIndex nLocalActiveCells() const { return m_nLocalActiveCells; }

  array1d< bool > & localCellIsOwned() { return m_localCellIsOwned; }
  array1d< bool > const & localCellIsOwned() const { return m_localCellIsOwned; }

  array1d< localIndex > & localActiveCellToLocalCell() { return m_localActiveCellToLocalCell; }
  array1d< localIndex > const & localActiveCellToLocalCell() const { return m_localActiveCellToLocalCell; }

  array1d< localIndex > & localCellToLocalCPVertices() { return m_localCellToLocalCPVertices; }
  array1d< localIndex > const & localCellToLocalCPVertices() const { return m_localCellToLocalCPVertices; }

  array1d< globalIndex > & localActiveCellToGlobalActiveCell() { return m_localActiveCellToGlobalActiveCell; }
  array1d< globalIndex > const & localActiveCellToGlobalActiveCell() const { return m_localActiveCellToGlobalActiveCell; }

private:

  /// global information

  /// total number of cells in the x direction
  localIndex m_nX;

  /// total number of cells in the y direction
  localIndex m_nY;

  /// total number of cells in the z direction
  localIndex m_nZ;



  /// local information: definition of the MPI partition

  /// beginning of the MPI partition in the x direction
  localIndex m_iMin;

  /// beginning of the MPI partition in the y direction
  localIndex m_jMin;

  /// end of the MPI partition in the x direction
  localIndex m_iMax;

  /// end of the MPI partition in the y direction
  localIndex m_jMax;

  /// overlap in the x-min direction
  localIndex m_iMinOverlap;

  /// overlap in the y-min direction
  localIndex m_jMinOverlap;

  /// overlap in the x-max direction
  localIndex m_iMaxOverlap;

  /// overlap in the y-max direction
  localIndex m_jMaxOverlap;



  /// local information: corner-point information

  /// local content of the COORD keyword in the GRDECL file
  array1d< real64 > m_coord;

  /// local content of the ZCORN keyword in the GRDECL file
  array1d< real64 > m_zcorn;

  /// local content of the ACTNUM keyword in the GRDECL file
  array1d< localIndex > m_actnum;



  /// node data

  /// number of vertices owned by my rank
  localIndex m_nOwnedVertices;

  /// number of vertices owned by my rank
  localIndex m_nLocalVertices;

  /// original (duplicated) CPG vertices
  array2d< real64 > m_cPVertices;

  /// true if a cPVertex is owned by my rank, false otherwise
  array1d< bool > m_cPVertexIsOwned;

  /// vertices obtained by filtering out duplicates in m_cPVertices
  array2d< real64 > m_vertices;

  /// true if a vertex is owned by my rank, false otherwise
  array1d< bool > m_vertexIsOwned;

  /// map from cPVertex local index to unique (filtered) vertex
  array1d< localIndex > m_localCPVertexToLocalVertex;

  /// map from vertex local index to vertex global index
  array1d< globalIndex > m_localVertexToGlobalVertex;



  /// cell data

  /// number of owned active and valid cells
  localIndex m_nOwnedActiveCells;

  /// number of local active and valid cells
  localIndex m_nLocalActiveCells;

  /// true if a cell is owned by my rank, false otherwise
  array1d< bool > m_localCellIsOwned;

  /// map from active cell local index to cell local index
  array1d< localIndex > m_localActiveCellToLocalCell;

  /// map from local index to original CPG nodes
  array1d< localIndex > m_localCellToLocalCPVertices;

  /// map from active cell local index to active cell global index
  array1d< globalIndex > m_localActiveCellToGlobalActiveCell;



  /// name of the mesh
  string m_meshName;

};

} // end namespace CPMesh

} // end namespace geosx

#endif //GEOSX_MESHUTILITIES_CPMESH_CPMESHDATA_HPP_
