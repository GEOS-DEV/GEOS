/*	
 * ------------------------------------------------------------------------------------------------------------	
 * SPDX-License-Identifier: LGPL-2.1-only	
 *	
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC	
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University	
 * Copyright (c) 2018-2020 Total, S.A	
 * Copyright (c) 2020-     GEOSX Contributors	
 * All right reserved	
 *	
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.	
 * ------------------------------------------------------------------------------------------------------------	
 */

#include "common/DataTypes.hpp"

namespace geosx
{

void getFaceNodes( string const & elementType,
                   localIndex const iElement,
                   localIndex const iFace,
                   array2d< localIndex, cells::NODE_MAP_PERMUTATION > const & elementToNodes,
                   array1d <localIndex> & nodeIndices )
{
  if( elementType == "C3D8" )
  {
    nodeIndices.resize( 4 );
    if( iFace == 0 )
    {
      nodeIndices[0] = elementToNodes[iElement][0];
      nodeIndices[1] = elementToNodes[iElement][1];
      nodeIndices[2] = elementToNodes[iElement][5];
      nodeIndices[3] = elementToNodes[iElement][4];
    }
    else if( iFace == 1 )
    {
      nodeIndices[0] = elementToNodes[iElement][0];
      nodeIndices[1] = elementToNodes[iElement][2];
      nodeIndices[2] = elementToNodes[iElement][3];
      nodeIndices[3] = elementToNodes[iElement][1];
    }
    else if( iFace == 2 )
    {
      nodeIndices[0] = elementToNodes[iElement][0];
      nodeIndices[1] = elementToNodes[iElement][4];
      nodeIndices[2] = elementToNodes[iElement][6];
      nodeIndices[3] = elementToNodes[iElement][2];

    }
    else if( iFace == 3 )
    {
      nodeIndices[0] = elementToNodes[iElement][1];
      nodeIndices[1] = elementToNodes[iElement][3];
      nodeIndices[2] = elementToNodes[iElement][7];
      nodeIndices[3] = elementToNodes[iElement][5];
    }
    else if( iFace == 4 )
    {
      nodeIndices[0] = elementToNodes[iElement][3];
      nodeIndices[1] = elementToNodes[iElement][2];
      nodeIndices[2] = elementToNodes[iElement][6];
      nodeIndices[3] = elementToNodes[iElement][7];
    }
    else if( iFace == 5 )
    {
      nodeIndices[0] = elementToNodes[iElement][4];
      nodeIndices[1] = elementToNodes[iElement][5];
      nodeIndices[2] = elementToNodes[iElement][7];
      nodeIndices[3] = elementToNodes[iElement][6];
    }
  }
  else if( elementType == "C3D6" )
  {
    if( iFace == 0 )
    {
      nodeIndices.resize( 4 );
      nodeIndices[0] = elementToNodes[iElement][0];
      nodeIndices[1] = elementToNodes[iElement][1];
      nodeIndices[2] = elementToNodes[iElement][5];
      nodeIndices[3] = elementToNodes[iElement][4];
    }
    else if( iFace == 1 )
    {
      nodeIndices.resize( 4 );
      nodeIndices[0] = elementToNodes[iElement][0];
      nodeIndices[1] = elementToNodes[iElement][2];
      nodeIndices[2] = elementToNodes[iElement][3];
      nodeIndices[3] = elementToNodes[iElement][1];
    }
    else if( iFace == 2 )
    {
      nodeIndices.resize( 3 );
      nodeIndices[0] = elementToNodes[iElement][0];
      nodeIndices[1] = elementToNodes[iElement][2];
      nodeIndices[2] = elementToNodes[iElement][4];
    }
    else if( iFace == 3 )
    {
      nodeIndices.resize( 3 );
      nodeIndices[0] = elementToNodes[iElement][1];
      nodeIndices[1] = elementToNodes[iElement][3];
      nodeIndices[2] = elementToNodes[iElement][5];
    }
    else if( iFace == 4 )
    {
      nodeIndices.resize( 4 );
      nodeIndices[0] = elementToNodes[iElement][2];
      nodeIndices[1] = elementToNodes[iElement][3];
      nodeIndices[2] = elementToNodes[iElement][5];
      nodeIndices[3] = elementToNodes[iElement][4];
    }
  }
  else if( elementType == "C3D4" )
  {
    nodeIndices.resize( 3 );
    if( iFace == 0 )
    {
      nodeIndices[0] = elementToNodes[iElement][0];
      nodeIndices[1] = elementToNodes[iElement][2];
      nodeIndices[2] = elementToNodes[iElement][1];
    }
    else if( iFace == 1 )
    {
      nodeIndices[0] = elementToNodes[iElement][0];
      nodeIndices[1] = elementToNodes[iElement][1];
      nodeIndices[2] = elementToNodes[iElement][3];
    }
    else if( iFace == 2 )
    {
      nodeIndices[0] = elementToNodes[iElement][0];
      nodeIndices[1] = elementToNodes[iElement][3];
      nodeIndices[2] = elementToNodes[iElement][2];
    }
    else if( iFace == 3 )
    {
      nodeIndices[0] = elementToNodes[iElement][1];
      nodeIndices[1] = elementToNodes[iElement][2];
      nodeIndices[2] = elementToNodes[iElement][3];
    }
  }
  else if( elementType == "C3D5" )
  {
    if( iFace == 0 )
    {
      nodeIndices.resize( 4 );
      nodeIndices[0] = elementToNodes[iElement][0];
      nodeIndices[1] = elementToNodes[iElement][1];
      nodeIndices[2] = elementToNodes[iElement][2];
      nodeIndices[3] = elementToNodes[iElement][3];
    }
    else if( iFace == 1 )
    {
      nodeIndices.resize( 3 );
      nodeIndices[0] = elementToNodes[iElement][0];
      nodeIndices[1] = elementToNodes[iElement][1];
      nodeIndices[2] = elementToNodes[iElement][4];
    }
    else if( iFace == 2 )
    {
      nodeIndices.resize( 3 );
      nodeIndices[0] = elementToNodes[iElement][1];
      nodeIndices[1] = elementToNodes[iElement][2];
      nodeIndices[2] = elementToNodes[iElement][4];
    }
    else if( iFace == 3 )
    {
      nodeIndices.resize( 3 );
      nodeIndices[0] = elementToNodes[iElement][2];
      nodeIndices[1] = elementToNodes[iElement][3];
      nodeIndices[2] = elementToNodes[iElement][4];
    }
    else if( iFace == 4 )
    {
      nodeIndices.resize( 3 );
      nodeIndices[0] = elementToNodes[iElement][3];
      nodeIndices[1] = elementToNodes[iElement][0];
      nodeIndices[2] = elementToNodes[iElement][4];
    }
  }
  else
  {
    GEOSX_ERROR( "Error. Don't know what kind of element this is and cannot build faces." );
  }
}

}
