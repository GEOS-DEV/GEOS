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
 * @file OutputUtilities.hpp
 */

#ifndef GEOSX_MESHUTILITIES_CPMESH_UTILITIES_OUTPUTUTILITIES_HPP_
#define GEOSX_MESHUTILITIES_CPMESH_UTILITIES_OUTPUTUTILITIES_HPP_

#include "meshUtilities/CPMesh/CPMeshData.hpp"

namespace geosx
{

namespace CPMesh
{

namespace OutputUtilities
{

/**
 * @brief Debug function that outputs a vtk file representing the new mesh
 * @param[in] meshVertices the struct holding vertex information
 * @param[in] meshFaces the struct holding face information
 * @param[in] meshCells the struct holding cell information
 */
void outputDebugVTKFile( CPMeshVertices const & meshVertices,
                         CPMeshFaces const & meshFaces,
                         CPMeshCells const & meshCells )
{
  GEOSX_UNUSED_VAR( meshFaces ); // not needed for now

  std::ofstream myfile;
  myfile.open ( "debug.vtk" );
  myfile << "# vtk DataFile Version 3.0\n";
  myfile << "debug mesh\n";
  myfile << "ASCII\n";
  myfile << "DATASET UNSTRUCTURED_GRID\n";

  array2d< real64 > const & vertices = meshVertices.m_vertices;
  myfile << "POINTS " << vertices.size() << " float\n";
  for( localIndex iVertex = 0; iVertex < vertices.size( 0 ); ++iVertex )
  {
    myfile << vertices( iVertex, 0 ) << " " << vertices( iVertex, 1 ) << " " << vertices( iVertex, 2 ) << "\n";
  }

  array1d< localIndex > const & activeCellToCell = meshCells.m_activeCellToCell;
  array1d< localIndex > const & cellToCPVertices = meshCells.m_cellToCPVertices;
  array1d< localIndex > const & cPVertexToVertex = meshVertices.m_cPVertexToVertex;
  localIndex const nActiveCells = activeCellToCell.size();

  myfile << "CELLS " << nActiveCells << " " << 9*nActiveCells << "\n";
  for( localIndex iActiveCell = 0; iActiveCell < nActiveCells; ++iActiveCell )
  {
    localIndex const iLocalCell = activeCellToCell( iActiveCell );
    localIndex const iFirstCPVertex = cellToCPVertices( iLocalCell );
    myfile << "8 ";
    for( localIndex pos = 0; pos < 8; ++pos )
    {
      myfile << cPVertexToVertex( iFirstCPVertex + pos ) << " ";
    }
    myfile << "\n";
  }

  myfile << "CELL_TYPES " << nActiveCells << std::endl;
  for( localIndex iActiveCell = 0; iActiveCell < nActiveCells; ++iActiveCell )
  {
    myfile << "12\n";
  }
  myfile.close();
}


} // end namespace OutputUtilities

} // end namespace CPMesh

} // end namespace geosx

#endif //GEOSX_MESHUTILITIES_CPMESH_OUTPUTUTILITIES_HPP_
