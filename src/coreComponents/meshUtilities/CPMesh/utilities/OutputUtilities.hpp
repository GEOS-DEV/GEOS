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

void outputDebugVTKFile( CPMeshData const & cPMeshData )
{
  std::ofstream myfile;
  myfile.open ( "debug.vtk" );
  myfile << "# vtk DataFile Version 3.0\n";
  myfile << "debug mesh\n";
  myfile << "ASCII\n";
  myfile << "DATASET UNSTRUCTURED_GRID\n";

  array2d< real64 > const & vertices = cPMeshData.vertices();
  myfile << "POINTS " << vertices.size() << " float\n";
  for( localIndex iVertex = 0; iVertex < vertices.size( 0 ); ++iVertex )
  {
    myfile << vertices( iVertex, 0 ) << " " << vertices( iVertex, 1 ) << " " << vertices( iVertex, 2 ) << "\n";
  }

  localIndex const nLocalActiveCells = cPMeshData.nLocalActiveCells();
  array1d< localIndex > const & localActiveCellToLocalCell = cPMeshData.localActiveCellToLocalCell();
  array1d< localIndex > const & localCellToLocalCPVertices = cPMeshData.localCellToLocalCPVertices();
  array1d< localIndex > const & localCPVertexToLocalVertex = cPMeshData.localCPVertexToLocalVertex();

  myfile << "CELLS " << nLocalActiveCells << " " << 9*nLocalActiveCells << "\n";
  for( localIndex iLocalActiveCell = 0; iLocalActiveCell < nLocalActiveCells; ++iLocalActiveCell )
  {
    localIndex const iLocalCell = localActiveCellToLocalCell( iLocalActiveCell );
    localIndex const iFirstCPVertex = localCellToLocalCPVertices( iLocalCell );
    myfile << "8 ";
    for( localIndex pos = 0; pos < 8; ++pos )
    {
      myfile << localCPVertexToLocalVertex( iFirstCPVertex + pos ) << " ";
    }
    myfile << "\n";
  }

  myfile << "CELL_TYPES " << nLocalActiveCells << std::endl;
  for( localIndex iLocalActiveCell = 0; iLocalActiveCell < nLocalActiveCells; ++iLocalActiveCell )
  {
    myfile << "12\n";
  }
  myfile.close();
}


} // end namespace OutputUtilities

} // end namespace CPMesh

} // end namespace geosx

#endif //GEOSX_MESHUTILITIES_CPMESH_OUTPUTUTILITIES_HPP_
