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
 * @file CPMeshParser.hpp
 */

#ifndef GEOSX_MESHUTILITIES_CPMESH_CPMESHPARSER_HPP_
#define GEOSX_MESHUTILITIES_CPMESH_CPMESHPARSER_HPP_

#include "codingUtilities/StringUtilities.hpp"
#include "common/DataTypes.hpp"
#include "dataRepository/ObjectCatalog.hpp"
#include "meshUtilities/CPMesh/CPMeshData.hpp"

namespace geosx
{

namespace CPMesh
{

class CPMeshParser
{

public:

  CPMeshParser( string const & name ) { GEOSX_UNUSED_VAR( name ); }

  /// Default virtual destructor
  virtual ~CPMeshParser() = default;

  /// Default copy constructor
  CPMeshParser( CPMeshParser const & ) = default;

  /// Default move constructor
  CPMeshParser( CPMeshParser && ) = default;

  /// Deleted copy assignment operator
  CPMeshParser & operator=( CPMeshParser const & ) = delete;

  /// Deleted move assignment operator
  CPMeshParser & operator=( CPMeshParser && ) = delete;

  using CatalogInterface = dataRepository::CatalogInterface< CPMeshParser,
                                                             string const & >;
  static typename CatalogInterface::CatalogType & getCatalog()
  {
    static CatalogInterface::CatalogType catalog;
    return catalog;
  }

  static string catalogName() { return "CPMeshParser"; }

  virtual string getCatalogName() const { return catalogName(); }

  /**
   * @brief Rank 0 reads the keyword DIMENS or SPECGRID to get the total number of cells (nX, nY, nZ)
   */
  void readNumberOfCells( Path const & filePath,
                          localIndex & nX,
                          localIndex & nY,
                          localIndex & nZ ) const;

  /**
   * @brief Each rank reads its part of the mesh
   * @details This involves reading the main keywords (COORD, ZCORN, ACTNUM), plus the property keywords
   */
  void readMesh( Path const & filePath,
                 CPMeshData & cpMeshData ) const;

private:

  /**
   * @brief Utility (temporary) functions that reads all the data below a keyword
   * @details This is a carry-over from PAMELA, but this function should go if we want to load only local data
   * (i.e., data corresponding to the MPI domain owned by each rank)
   */
  std::string extractDataBelowKeyword( std::istringstream & stringBlock ) const;

  /**
   * @brief Read the local content of the COORD keyword (i.e, data corresponding to the MPI domain owned by each rank)
   * @details This function will have to be rewritten
   * @param meshFile the content of the mesh file
   * @param cPMeshData the datastructure describing CP mesh
   */
  void readLocalCOORD( std::istringstream & meshFile,
                       CPMeshData & cPMeshData ) const;

  /**
   * @brief Read the local content of the ZCORN keyword (i.e, data corresponding to the MPI domain owned by each rank)
   * @details This function will have to be rewritten
   * @param meshFile the content of the mesh file
   * @param cPMeshData the datastructure describing CP mesh
   */
  void readLocalZCORN( std::istringstream & meshFile,
                       CPMeshData & cPMeshData ) const;

  /**
   * @brief Read the local content of the ACTNUM keyword (i.e, data corresponding to the MPI domain owned by each rank)
   * @details This function will have to be rewritten
   * @param meshFile the content of the mesh file
   * @param cPMeshData the datastructure describing CP mesh
   */
  void readLocalACTNUM( std::istringstream & meshFile,
                        CPMeshData & cPMeshData ) const;

  /**
   * @brief If ACTNUM is absent from the file, fill the actnum vector with ones (all cells are active)
   * @param cPMeshData the datastructure describing CP mesh
   */
  void fillLocalACTNUM( CPMeshData & cPMeshData ) const;

  /// name of the mesh
  string m_meshName;

};

} // end namespace CPMesh

} // end namespace geosx

#endif //GEOSX_MESHUTILITIES_CPMESH_CPMESHPARSER_HPP_
