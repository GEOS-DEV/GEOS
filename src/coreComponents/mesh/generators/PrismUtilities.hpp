/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2020-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#ifndef GEOSX_MESH_GENERATORS_PRISMUTILITIES_HPP_
#define GEOSX_MESH_GENERATORS_PRISMUTILITIES_HPP_

#include "common/DataTypes.hpp"
#include "common/Span.hpp"

namespace geosx
{

template< localIndex N >
localIndex getFaceNodesPrism( localIndex const faceNum,
                              arraySlice1d< localIndex const, cells::NODE_MAP_USD - 1 > const & elemNodes,
                              Span< localIndex > const faceNodes )
{
  if( faceNum == 0 )
  {
    GEOSX_ERROR_IF_LT( faceNodes.size(), 4 );
    faceNodes[0] = elemNodes[0];
    faceNodes[1] = elemNodes[1];
    faceNodes[2] = elemNodes[N+1];
    faceNodes[3] = elemNodes[N];
    return 4;
  }
  else if( faceNum == 1 )
  {
    GEOSX_ERROR_IF_LT( faceNodes.size(), N );
    faceNodes[0] = elemNodes[0];
    for( localIndex i = 1; i <  N; ++i )
    {
      faceNodes[i] = elemNodes[N-i];
    }
    return N;
  }
  else if( faceNum == 2 )
  {
    GEOSX_ERROR_IF_LT( faceNodes.size(), 4 );
    faceNodes[0] = elemNodes[0];
    faceNodes[1] = elemNodes[N];
    faceNodes[2] = elemNodes[N*2-1];
    faceNodes[3] = elemNodes[N-1];
    return 4;
  }
  else if( faceNum >= 3 && faceNum <= N )
  {
    GEOSX_ERROR_IF_LT( faceNodes.size(), 4 );
    faceNodes[0] = elemNodes[faceNum-2];
    faceNodes[1] = elemNodes[faceNum-1];
    faceNodes[2] = elemNodes[N+faceNum-1];
    faceNodes[3] = elemNodes[N+faceNum-2];
    return 4;
  }
  else if( faceNum == N + 1 )
  {
    GEOSX_ERROR_IF_LT( faceNodes.size(), N );
    for( localIndex i = 0; i <  N; ++i )
    {
      faceNodes[i] = elemNodes[i+N];
    }
    return N;
  }
  else
  {
    GEOSX_ERROR( "Invalid local face index for element of type Wedge: " << faceNum );
    return 0;
  }

}

}

#endif // GEOSX_MESH_GENERATORS_PRISMUTILITIES_HPP_
