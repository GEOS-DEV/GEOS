/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

// Source includes
#include "codingUtilities/UnitTestUtilities.hpp"
#include "common/logger/Logger.hpp"
#include "finiteVolume/mimeticInnerProducts/MimeticInnerProductBase.hpp"
#include "finiteVolume/mimeticInnerProducts/TPFAInnerProduct.hpp"
#include "finiteVolume/mimeticInnerProducts/QuasiTPFAInnerProduct.hpp"
#include "finiteVolume/mimeticInnerProducts/SimpleInnerProduct.hpp"
#include "finiteVolume/mimeticInnerProducts/BdVLMInnerProduct.hpp"
#include "mainInterface/initialization.hpp"
#include "mesh/FaceManager.hpp"
#include "mesh/utilities/ComputationalGeometry.hpp"
#include <iostream>
#include <iomanip>

// TPL includes
#include <gtest/gtest.h>

using namespace geos;
using namespace geos::mimeticInnerProduct;
using namespace geos::computationalGeometry;
using namespace geos::testing;

struct InnerProductType
{
  static constexpr integer TPFA = 0;
  static constexpr integer QUASI_TPFA = 1;
  static constexpr integer QUASI_TPFA_WITH_MULTIPLIERS = 2;
  static constexpr integer SIMPLE = 3;
  static constexpr integer SIMPLE_WITH_MULTIPLIERS = 4;
  static constexpr integer BDVLM = 5;
  static constexpr integer BDVLM_WITH_MULTIPLIERS = 6;
};

static constexpr real64 consistency_tol = 1e-15;

void compareTransmissibilityMatrices( arraySlice2d< real64 const > const & transMatrix,
                                      arraySlice2d< real64 const > const & transMatrixRef )
{
  for( localIndex ifaceLoc = 0; ifaceLoc < transMatrix.size( 0 ); ++ifaceLoc )
  {
    for( localIndex jfaceLoc = 0; jfaceLoc < transMatrix.size( 1 ); ++jfaceLoc )
    {
      checkRelativeError( transMatrix( ifaceLoc, jfaceLoc ),
                          transMatrixRef( ifaceLoc, jfaceLoc ),
                          1e-15 );
    }
  }
}

void compareInnerProductMatrices( arraySlice2d< real64 const > const & M,
                                  arraySlice2d< real64 const > const & Mref )
{
  for( localIndex ifaceLoc = 0; ifaceLoc < M.size( 0 ); ++ifaceLoc )
  {
    for( localIndex jfaceLoc = 0; jfaceLoc < M.size( 1 ); ++jfaceLoc )
    {
      checkRelativeError( M( ifaceLoc, jfaceLoc ),
                          Mref( ifaceLoc, jfaceLoc ),
                          1e-15 );
    }
  }
}

void computeVolumeAndCenter( array2d< real64, nodes::REFERENCE_POSITION_PERM > const & nodePosition,
                             array1d< localIndex > const & toNodes,
                             real64 ( & elemCenter )[3],
                             real64 & elemVolume )
{
  localIndex const numNodes = toNodes.size();

  GEOS_ERROR_IF( numNodes != 8 && numNodes != 4,
                 "This number of nodes is not supported in the test yet" );

  LvArray::tensorOps::fill< 3 >( elemCenter, 0.0 );

  if( numNodes == 8 )
  {
    real64 Xlocal[8][3];
    for( integer a = 0; a < numNodes; ++a )
    {
      Xlocal[a][0] = nodePosition( toNodes( a ), 0 );
      Xlocal[a][1] = nodePosition( toNodes( a ), 1 );
      Xlocal[a][2] = nodePosition( toNodes( a ), 2 );
      LvArray::tensorOps::add< 3 >( elemCenter, Xlocal[a] );
    }
    elemVolume = computationalGeometry::hexahedronVolume( Xlocal );
  }
  else if( numNodes == 4 )
  {
    real64 Xlocal[4][3];
    for( integer a = 0; a < numNodes; ++a )
    {
      Xlocal[a][0] = nodePosition( toNodes( a ), 0 );
      Xlocal[a][1] = nodePosition( toNodes( a ), 1 );
      Xlocal[a][2] = nodePosition( toNodes( a ), 2 );
      LvArray::tensorOps::add< 3 >( elemCenter, Xlocal[a] );
    }
    elemVolume = computationalGeometry::tetrahedronVolume( Xlocal );
  }
  LvArray::tensorOps::scale< 3 >( elemCenter, 1.0 / numNodes );
}

void makeHexa( array2d< real64, nodes::REFERENCE_POSITION_PERM > & nodePosition,
               FaceManager::NodeMapType & faceToNodes,
               array1d< localIndex > & elemToFaces,
               real64 ( & elemCenter )[3],
               real64 & elemVolume,
               real64 ( & elemPerm )[3],
               real64 & lengthTolerance,
               integer const ipType,
               arraySlice2d< real64 > const & transMatrixRef )
{
  localIndex constexpr numNodes = 8;
  localIndex constexpr numFaces = 6;
  localIndex constexpr numNodesPerFace = 4;

  lengthTolerance = 1.73205e-8;

  elemPerm[ 0 ] = 1e-12;
  elemPerm[ 1 ] = 2e-12;
  elemPerm[ 2 ] = 3e-12;

  // elem-to-faces map
  elemToFaces.resize( numFaces );
  elemToFaces( 0 ) = 1;
  elemToFaces( 1 ) = 2;
  elemToFaces( 2 ) = 0;
  elemToFaces( 3 ) = 5;
  elemToFaces( 4 ) = 4;
  elemToFaces( 5 ) = 3;

  // face-to-nodes map
  faceToNodes.resize( numFaces );
  for( localIndex iface = 0; iface < numFaces; ++iface )
  {
    faceToNodes.resizeArray( iface, numNodesPerFace );
  }
  faceToNodes( 0, 0 ) = 0;
  faceToNodes( 0, 1 ) = 1;
  faceToNodes( 0, 2 ) = 3;
  faceToNodes( 0, 3 ) = 2;
  faceToNodes( 1, 0 ) = 0;
  faceToNodes( 1, 1 ) = 4;
  faceToNodes( 1, 2 ) = 5;
  faceToNodes( 1, 3 ) = 1;
  faceToNodes( 2, 0 ) = 0;
  faceToNodes( 2, 1 ) = 2;
  faceToNodes( 2, 2 ) = 6;
  faceToNodes( 2, 3 ) = 4;
  faceToNodes( 3, 0 ) = 1;
  faceToNodes( 3, 1 ) = 5;
  faceToNodes( 3, 2 ) = 7;
  faceToNodes( 3, 3 ) = 3;
  faceToNodes( 4, 0 ) = 6;
  faceToNodes( 4, 1 ) = 2;
  faceToNodes( 4, 2 ) = 3;
  faceToNodes( 4, 3 ) = 7;
  faceToNodes( 5, 0 ) = 4;
  faceToNodes( 5, 1 ) = 6;
  faceToNodes( 5, 2 ) = 7;
  faceToNodes( 5, 3 ) = 5;

  // node position
  nodePosition.resize( numNodes, 3 );
  nodePosition( 0, 0 ) = 0;
  nodePosition( 0, 1 ) = 0;
  nodePosition( 0, 2 ) = 0;
  nodePosition( 1, 0 ) = 0.75;
  nodePosition( 1, 1 ) = 0;
  nodePosition( 1, 2 ) = 1;
  nodePosition( 2, 0 ) = 0;
  nodePosition( 2, 1 ) = 1;
  nodePosition( 2, 2 ) = 0;
  nodePosition( 3, 0 ) = 0.75;
  nodePosition( 3, 1 ) = 1;
  nodePosition( 3, 2 ) = 1;
  nodePosition( 4, 0 ) = 1;
  nodePosition( 4, 1 ) = 0;
  nodePosition( 4, 2 ) = 0;
  nodePosition( 5, 0 ) = 1.75;
  nodePosition( 5, 1 ) = 0;
  nodePosition( 5, 2 ) = 1;
  nodePosition( 6, 0 ) = 1;
  nodePosition( 6, 1 ) = 1;
  nodePosition( 6, 2 ) = 0;
  nodePosition( 7, 0 ) = 1.75;
  nodePosition( 7, 1 ) = 1;
  nodePosition( 7, 2 ) = 1;

  array1d< localIndex > toNodes;
  toNodes.resize( numNodes );
  toNodes( 0 ) = 0;
  toNodes( 1 ) = 4;
  toNodes( 2 ) = 2;
  toNodes( 3 ) = 6;
  toNodes( 4 ) = 1;
  toNodes( 5 ) = 5;
  toNodes( 6 ) = 3;
  toNodes( 7 ) = 7;

  computeVolumeAndCenter( nodePosition,
                          toNodes,
                          elemCenter,
                          elemVolume );

  // reference transmissibility matrix computed with MRST
  if( ipType == InnerProductType::TPFA )
  {
    transMatrixRef( 0, 0 ) = 4e-12;
    transMatrixRef( 1, 1 ) = 3.84e-12;
    transMatrixRef( 2, 2 ) = 2e-12;
    transMatrixRef( 3, 3 ) = 2e-12;
    transMatrixRef( 4, 4 ) = 4e-12;
    transMatrixRef( 5, 5 ) = 3.84e-12;
  }
  else if( ipType == InnerProductType::QUASI_TPFA )
  {
    transMatrixRef( 0, 0 ) = 4e-12;

    transMatrixRef( 1, 1 ) = 6e-12;
    transMatrixRef( 1, 2 ) = -2.250e-12;
    transMatrixRef( 1, 3 ) = 2.250e-12;

    transMatrixRef( 2, 1 ) = -2.225e-12;
    transMatrixRef( 2, 2 ) = 5.375e-12;
    transMatrixRef( 2, 5 ) = 2.225e-12;

    transMatrixRef( 3, 1 ) = 2.25e-12;
    transMatrixRef( 3, 3 ) = 5.375e-12;
    transMatrixRef( 3, 5 ) = -2.25e-12;

    transMatrixRef( 4, 4 ) = 4e-12;

    transMatrixRef( 5, 2 ) = 2.25e-12;
    transMatrixRef( 5, 5 ) = 6e-12;
    transMatrixRef( 5, 3 ) = -2.25e-12;
  }
  else if( ipType == InnerProductType::QUASI_TPFA_WITH_MULTIPLIERS )
  {
    transMatrixRef( 0, 0 ) = 4.00e-12;

    transMatrixRef( 1, 1 ) = 4.817e-12;
    transMatrixRef( 1, 2 ) = -1.829e-12;
    transMatrixRef( 1, 3 ) = 0.094e-12;
    transMatrixRef( 1, 5 ) = 0.851e-12;

    transMatrixRef( 2, 1 ) = -1.829e-12;
    transMatrixRef( 2, 2 ) = 3.991e-12;
    transMatrixRef( 2, 3 ) = 0.008e-12;
    transMatrixRef( 2, 5 ) = 1.315e-12;

    transMatrixRef( 3, 1 ) = 0.094e-12;
    transMatrixRef( 3, 2 ) = 0.008e-12;
    transMatrixRef( 3, 3 ) = 0.213e-12;
    transMatrixRef( 3, 5 ) = -0.068e-12;

    transMatrixRef( 4, 4 ) = 4e-12;

    transMatrixRef( 5, 1 ) = 0.851e-12;
    transMatrixRef( 5, 2 ) = 1.315e-12;
    transMatrixRef( 5, 5 ) = 3.703e-12;
    transMatrixRef( 5, 3 ) = -0.068e-12;
  }
  else if( ipType == InnerProductType::SIMPLE )
  {
    transMatrixRef( 0, 0 ) = 8e-12;
    transMatrixRef( 0, 4 ) = 4e-12;

    transMatrixRef( 1, 1 ) = 9e-12;
    transMatrixRef( 1, 2 ) = -2.250e-12;
    transMatrixRef( 1, 3 ) = 2.250e-12;
    transMatrixRef( 1, 5 ) = 3.e-12;

    transMatrixRef( 2, 1 ) = -2.225e-12;
    transMatrixRef( 2, 2 ) = 12.06e-12;
    transMatrixRef( 2, 3 ) = 6.69e-12;
    transMatrixRef( 2, 5 ) = 2.25e-12;

    transMatrixRef( 3, 1 ) = 2.25e-12;
    transMatrixRef( 3, 2 ) = 6.69e-12,
    transMatrixRef( 3, 3 ) = 12.06e-12;
    transMatrixRef( 3, 5 ) = -2.25e-12;

    transMatrixRef( 4, 4 ) = 8e-12;
    transMatrixRef( 4, 0 ) = 4e-12;

    transMatrixRef( 5, 1 ) = 3e-12;
    transMatrixRef( 5, 2 ) = 2.25e-12;
    transMatrixRef( 5, 5 ) = 9e-12;
    transMatrixRef( 5, 3 ) = -2.25e-12;
  }
  else if( ipType == InnerProductType::SIMPLE_WITH_MULTIPLIERS )
  {
    transMatrixRef( 0, 0 ) = 8e-12;
    transMatrixRef( 0, 4 ) = 4e-12;

    transMatrixRef( 1, 1 ) = 7.494e-12;
    transMatrixRef( 1, 2 ) = -2.757e-12;
    transMatrixRef( 1, 3 ) = 0.066e-12;
    transMatrixRef( 1, 5 ) = 2.530e-12;

    transMatrixRef( 2, 1 ) = -2.757e-12;
    transMatrixRef( 2, 2 ) = 5.499e-12;
    transMatrixRef( 2, 3 ) = 0.088e-12;
    transMatrixRef( 2, 5 ) = 1.548e-12;

    transMatrixRef( 3, 1 ) = 0.066e-12;
    transMatrixRef( 3, 2 ) = 0.088e-12;
    transMatrixRef( 3, 3 ) = 0.218e-12;
    transMatrixRef( 3, 5 ) = -0.037e-12;

    transMatrixRef( 4, 4 ) = 8e-12;
    transMatrixRef( 4, 0 ) = 4e-12;

    transMatrixRef( 5, 1 ) = 2.53e-12;
    transMatrixRef( 5, 2 ) = 1.548e-12;
    transMatrixRef( 5, 5 ) = 5.317e-12;
    transMatrixRef( 5, 3 ) = -0.037e-12;

  }
  else if( ipType == InnerProductType::BDVLM )
  {
    transMatrixRef( 0, 0 ) =  4.240e-12;
    transMatrixRef( 0, 4 ) =  0.240e-12;

    transMatrixRef( 1, 1 ) =  5.240e-12;
    transMatrixRef( 1, 2 ) = -2.250e-12;
    transMatrixRef( 1, 3 ) =  2.250e-12;
    transMatrixRef( 1, 5 ) = -0.760e-12;

    transMatrixRef( 2, 1 ) = -2.225e-12;
    transMatrixRef( 2, 2 ) =  6.187e-12;
    transMatrixRef( 2, 3 ) =  0.812e-12;
    transMatrixRef( 2, 5 ) =  2.25e-12;

    transMatrixRef( 3, 1 ) =  2.250e-12;
    transMatrixRef( 3, 2 ) =  0.812e-12;
    transMatrixRef( 3, 3 ) =  6.187e-12;
    transMatrixRef( 3, 5 ) = -2.225e-12;

    transMatrixRef( 4, 4 ) =  4.240e-12;
    transMatrixRef( 4, 0 ) =  0.240e-12;

    transMatrixRef( 5, 1 ) = -0.760e-12;
    transMatrixRef( 5, 2 ) =  2.25e-12;
    transMatrixRef( 5, 5 ) =  5.240e-12;
    transMatrixRef( 5, 3 ) = -2.225e-12;
  }
  else if( ipType == InnerProductType::BDVLM_WITH_MULTIPLIERS )
  {
    transMatrixRef( 0, 0 ) =  4.240e-12;
    transMatrixRef( 0, 4 ) =  0.240e-12;

    transMatrixRef( 1, 1 ) =  4.179e-12;
    transMatrixRef( 1, 2 ) = -1.924e-12;
    transMatrixRef( 1, 3 ) =  0.082e-12;
    transMatrixRef( 1, 5 ) =  0.233e-12;

    transMatrixRef( 2, 1 ) = -1.924e-12;
    transMatrixRef( 2, 2 ) =  4.364e-12;
    transMatrixRef( 2, 3 ) =  0.029e-12;
    transMatrixRef( 2, 5 ) =  1.489e-12;

    transMatrixRef( 3, 1 ) =  0.082e-12;
    transMatrixRef( 3, 2 ) =  0.029e-12;
    transMatrixRef( 3, 3 ) =  0.214e-12;
    transMatrixRef( 3, 5 ) = -0.064e-12;

    transMatrixRef( 4, 4 ) =  4.240e-12;
    transMatrixRef( 4, 0 ) =  0.240e-12;

    transMatrixRef( 5, 1 ) =  0.233e-12;
    transMatrixRef( 5, 2 ) =  1.489e-12;
    transMatrixRef( 5, 5 ) =  3.288e-12;
    transMatrixRef( 5, 3 ) = -0.064e-12;
  }
}


void makeTetra( array2d< real64, nodes::REFERENCE_POSITION_PERM > & nodePosition,
                FaceManager::NodeMapType & faceToNodes,
                array1d< localIndex > & elemToFaces,
                real64 ( & elemCenter )[3],
                real64 & elemVolume,
                real64 ( & elemPerm )[3],
                real64 & lengthTolerance,
                integer const ipType,
                arraySlice2d< real64 > const & transMatrixRef )
{
  localIndex constexpr numNodes = 4;
  localIndex constexpr numFaces = 4;
  localIndex constexpr numNodesPerFace = 3;

  lengthTolerance = 1.73205e-8;

  elemPerm[ 0 ] = 1e-12;
  elemPerm[ 1 ] = 2e-12;
  elemPerm[ 2 ] = 3e-12;

  // elem-to-faces map
  elemToFaces.resize( numFaces );
  elemToFaces( 0 ) = 0;
  elemToFaces( 1 ) = 1;
  elemToFaces( 2 ) = 2;
  elemToFaces( 3 ) = 3;

  // face-to-nodes map
  faceToNodes.resize( numFaces );
  for( localIndex iface = 0; iface < numFaces; ++iface )
  {
    faceToNodes.resizeArray( iface, numNodesPerFace );
  }
  faceToNodes( 0, 0 ) = 0;
  faceToNodes( 0, 1 ) = 1;
  faceToNodes( 0, 2 ) = 2;
  faceToNodes( 1, 0 ) = 0;
  faceToNodes( 1, 1 ) = 1;
  faceToNodes( 1, 2 ) = 3;
  faceToNodes( 2, 0 ) = 1;
  faceToNodes( 2, 1 ) = 2;
  faceToNodes( 2, 2 ) = 3;
  faceToNodes( 3, 0 ) = 2;
  faceToNodes( 3, 1 ) = 0;
  faceToNodes( 3, 2 ) = 3;

  // node position
  nodePosition.resize( numNodes, 3 );
  nodePosition( 0, 0 ) = 0;
  nodePosition( 0, 1 ) = 0;
  nodePosition( 0, 2 ) = 0;
  nodePosition( 1, 0 ) = 1;
  nodePosition( 1, 1 ) = 0;
  nodePosition( 1, 2 ) = 0;
  nodePosition( 2, 0 ) = 0;
  nodePosition( 2, 1 ) = 1;
  nodePosition( 2, 2 ) = 0;
  nodePosition( 3, 0 ) = 0;
  nodePosition( 3, 1 ) = 0;
  nodePosition( 3, 2 ) = 3;

  array1d< localIndex > toNodes;
  toNodes.resize( numNodes );
  toNodes( 0 ) = 0;
  toNodes( 1 ) = 1;
  toNodes( 2 ) = 2;
  toNodes( 3 ) = 3;

  computeVolumeAndCenter( nodePosition,
                          toNodes,
                          elemCenter,
                          elemVolume );

  // reference transmissibility matrix computed with MRST

  if( ipType == InnerProductType::TPFA )
  {
    transMatrixRef( 0, 0 ) = 1.952e-12;
    transMatrixRef( 1, 1 ) = 5.684e-12;
    transMatrixRef( 2, 2 ) = 9.818e-12;
    transMatrixRef( 3, 3 ) = 2.842e-12;
  }
  else if( ipType == InnerProductType::QUASI_TPFA )
  {
    transMatrixRef( 0, 0 ) =  5.25e-12;
    transMatrixRef( 0, 1 ) =  3.75e-12;
    transMatrixRef( 0, 2 ) =  2.25e-12;
    transMatrixRef( 0, 3 ) =  3.75e-12;

    transMatrixRef( 1, 0 ) =  3.75e-12;
    transMatrixRef( 1, 1 ) =  12.75e-12;
    transMatrixRef( 1, 2 ) = -5.25e-12;
    transMatrixRef( 1, 3 ) =  3.75e-12;

    transMatrixRef( 2, 0 ) =  2.25e-12;
    transMatrixRef( 2, 1 ) = -5.25e-12;
    transMatrixRef( 2, 2 ) =  18.75e-12;
    transMatrixRef( 2, 3 ) = -0.75e-12;

    transMatrixRef( 3, 0 ) =  3.75e-12;
    transMatrixRef( 3, 1 ) =  3.75e-12;
    transMatrixRef( 3, 2 ) = -0.75e-12;
    transMatrixRef( 3, 3 ) =  8.25e-12;
  }
  else if( ipType == InnerProductType::SIMPLE )
  {
    transMatrixRef( 0, 0 ) =  6.21e-12;
    transMatrixRef( 0, 1 ) =  4.71e-12;
    transMatrixRef( 0, 2 ) =  3.21e-12;
    transMatrixRef( 0, 3 ) =  4.71e-12;

    transMatrixRef( 1, 0 ) =  4.71e-12;
    transMatrixRef( 1, 1 ) =  13.71e-12;
    transMatrixRef( 1, 2 ) = -4.29e-12;
    transMatrixRef( 1, 3 ) =  4.71e-12;

    transMatrixRef( 2, 0 ) =  3.21e-12;
    transMatrixRef( 2, 1 ) = -4.29e-12;
    transMatrixRef( 2, 2 ) =  19.71e-12;
    transMatrixRef( 2, 3 ) =  0.21e-12;

    transMatrixRef( 3, 0 ) =  4.71e-12;
    transMatrixRef( 3, 1 ) =  4.71e-12;
    transMatrixRef( 3, 2 ) =  0.21e-12;
    transMatrixRef( 3, 3 ) =  9.21e-12;
  }
  else if( ipType == InnerProductType::BDVLM )
  {
    transMatrixRef( 0, 0 ) =  2.99e-12;
    transMatrixRef( 0, 1 ) =  1.49e-12;
    transMatrixRef( 0, 2 ) = -0.01e-12;
    transMatrixRef( 0, 3 ) =  1.49e-12;

    transMatrixRef( 1, 0 ) =  1.49e-12;
    transMatrixRef( 1, 1 ) = 10.49e-12;
    transMatrixRef( 1, 2 ) = -7.51e-12;
    transMatrixRef( 1, 3 ) =  1.49e-12;

    transMatrixRef( 2, 0 ) = -0.01e-12;
    transMatrixRef( 2, 1 ) = -7.51e-12;
    transMatrixRef( 2, 2 ) = 16.49e-12;
    transMatrixRef( 2, 3 ) = -3.01e-12;

    transMatrixRef( 3, 0 ) =  1.49e-12;
    transMatrixRef( 3, 1 ) =  1.49e-12;
    transMatrixRef( 3, 2 ) = -3.01e-12;
    transMatrixRef( 3, 3 ) =  5.99e-12;

  }
}

void makeCube( array2d< real64, nodes::REFERENCE_POSITION_PERM > & nodePosition,
               FaceManager::NodeMapType & faceToNodes,
               array1d< localIndex > & elemToFaces,
               real64 ( & elemCenter )[3],
               real64 & elemVolume,
               real64 ( & elemPerm )[3],
               real64 & lengthTolerance,
               integer const ipType,
               arraySlice2d< real64 > const & transMatrixRef )
{
  localIndex constexpr numNodes = 8;
  localIndex constexpr numFaces = 6;
  localIndex constexpr numNodesPerFace = 4;

  lengthTolerance = 1e-12;

  // permeability
  elemPerm[0] = 1.0;
  elemPerm[1] = 1.0;
  elemPerm[2] = 1.0;

  // elem-to-faces
  elemToFaces.resize( numFaces );
  for( localIndex i = 0; i < numFaces; ++i )
    elemToFaces( i ) = i;

  // face-to-nodes
  faceToNodes.resize( numFaces );
  for( localIndex f = 0; f < numFaces; ++f )
    faceToNodes.resizeArray( f, numNodesPerFace );

  // vertices
  nodePosition.resize( numNodes, 3 );

  nodePosition( 0, 0 ) = 0; nodePosition( 0, 1 ) = 0; nodePosition( 0, 2 ) = 0;
  nodePosition( 1, 0 ) = 1; nodePosition( 1, 1 ) = 0; nodePosition( 1, 2 ) = 0;
  nodePosition( 2, 0 ) = 1; nodePosition( 2, 1 ) = 1; nodePosition( 2, 2 ) = 0;
  nodePosition( 3, 0 ) = 0; nodePosition( 3, 1 ) = 1; nodePosition( 3, 2 ) = 0;
  nodePosition( 4, 0 ) = 0; nodePosition( 4, 1 ) = 0; nodePosition( 4, 2 ) = 1;
  nodePosition( 5, 0 ) = 1; nodePosition( 5, 1 ) = 0; nodePosition( 5, 2 ) = 1;
  nodePosition( 6, 0 ) = 1; nodePosition( 6, 1 ) = 1; nodePosition( 6, 2 ) = 1;
  nodePosition( 7, 0 ) = 0; nodePosition( 7, 1 ) = 1; nodePosition( 7, 2 ) = 1;

  // x = 0
  faceToNodes( 0, 0 ) = 0; faceToNodes( 0, 1 ) = 3; faceToNodes( 0, 2 ) = 7; faceToNodes( 0, 3 ) = 4;

  // x = 1
  faceToNodes( 1, 0 ) = 1; faceToNodes( 1, 1 ) = 2; faceToNodes( 1, 2 ) = 6; faceToNodes( 1, 3 ) = 5;

  // y = 0
  faceToNodes( 2, 0 ) = 0; faceToNodes( 2, 1 ) = 1; faceToNodes( 2, 2 ) = 5; faceToNodes( 2, 3 ) = 4;

  // y = 1
  faceToNodes( 3, 0 ) = 3; faceToNodes( 3, 1 ) = 2; faceToNodes( 3, 2 ) = 6; faceToNodes( 3, 3 ) = 7;

  // z = 0
  faceToNodes( 4, 0 ) = 0; faceToNodes( 4, 1 ) = 1; faceToNodes( 4, 2 ) = 2; faceToNodes( 4, 3 ) = 3;

  // z = 1
  faceToNodes( 5, 0 ) = 4; faceToNodes( 5, 1 ) = 5; faceToNodes( 5, 2 ) = 6; faceToNodes( 5, 3 ) = 7;

  // center & volume
  array1d< localIndex > toNodes;
  toNodes.resize( numNodes );
  for( localIndex i = 0; i < numNodes; ++i )
    toNodes( i ) = i;

  computeVolumeAndCenter( nodePosition,
                          toNodes,
                          elemCenter,
                          elemVolume );

  // reference matrix (MRST / MATLAB)
  if( ipType == InnerProductType::TPFA ||
      ipType == InnerProductType::QUASI_TPFA )
  {
    // diagonal
    transMatrixRef( 0, 0 ) = 0.5;
    transMatrixRef( 1, 1 ) = 0.5;
    transMatrixRef( 2, 2 ) = 0.5;
    transMatrixRef( 3, 3 ) = 0.5;
    transMatrixRef( 4, 4 ) = 0.5;
    transMatrixRef( 5, 5 ) = 0.5;
  }
  else if( ipType == InnerProductType::SIMPLE ||
           ipType == InnerProductType::BDVLM )
  {
    // SIMPLE = TPFA + stabilization
    transMatrixRef( 0, 0 ) = 0.333333333333;
    transMatrixRef( 1, 1 ) = 0.333333333333;
    transMatrixRef( 2, 2 ) = 0.333333333333;
    transMatrixRef( 3, 3 ) = 0.333333333333;
    transMatrixRef( 4, 4 ) = 0.333333333333;
    transMatrixRef( 5, 5 ) = 0.333333333333;

    transMatrixRef( 0, 1 ) = -0.166666666666;
    transMatrixRef( 1, 0 ) = -0.166666666666;
    transMatrixRef( 2, 3 ) = -0.166666666666;
    transMatrixRef( 3, 2 ) = -0.166666666666;
    transMatrixRef( 4, 5 ) = -0.166666666666;
    transMatrixRef( 5, 4 ) = -0.166666666666;
  }
}

template< localIndex NF, typename ARRAY_VIEW_T >
static void runConsistencyTest( array2d< real64, nodes::REFERENCE_POSITION_PERM > const & nodePosition,
                                FaceManager::NodeMapType const & faceToNodes,
                                real64 const elemPerm[3],
                                ARRAY_VIEW_T const & transMatrix )
{
  real64 N[NF][3], C[NF][3], TC[NF][3], K[3][3];
  real64 faceCenter[3], faceNormal[3];

  // full tensor K
  real64 perm[3] = { elemPerm[0], elemPerm[1], elemPerm[2] };
  MimeticInnerProductHelpers::makeFullTensor( perm, K );

  // compute face normals
  for( localIndex iface = 0; iface < NF; ++iface )
  {
    computationalGeometry::centroid_3DPolygon(
      faceToNodes[iface], nodePosition.toViewConst(),
      faceCenter, faceNormal );

    for( int d = 0; d < 3; ++d )
    {
      N[iface][d] = faceNormal[d];
    }
  }

  // C = N * K
  LvArray::tensorOps::Rij_eq_AikBkj< NF, 3, 3 >( C, N, K );

  // TC = T * C
  LvArray::tensorOps::Rij_eq_AikBkj< NF, 3, NF >( TC, transMatrix, C );

  // define diffMat and measure its norm
  real64 diffMat[NF][3];
  real64 diffNorm = 0.0;

  for( localIndex i = 0; i < NF; ++i )
  {
    for( localIndex j = 0; j < 3; ++j )
    {
      diffMat[i][j] = C[i][j] - TC[i][j];
    }
  }

  for( std::ptrdiff_t i = 0; i < NF; ++i )
  {
    diffNorm += LvArray::tensorOps::l2NormSquared< 3 >( diffMat[i] );
  }
  diffNorm = std::sqrt( diffNorm );

  // set tol = 1e-11 because of absolute permeability employed makeHexa is of order of 1e-12
  EXPECT_LT( diffNorm, 1e-11 );
}

TEST( testMimeticInnerProducts, TPFA_hexa )
{
  localIndex constexpr NF = 6;

  array2d< real64, nodes::REFERENCE_POSITION_PERM > nodePosition;
  FaceManager::NodeMapType faceToNodes;
  array1d< localIndex > elemToFaces;
  real64 elemCenter[3];
  real64 elemPerm[3];
  real64 elemVolume = 0;
  real64 lengthTolerance = 0;
  stackArray2d< real64, NF *NF > transMatrixRef( NF, NF );

  makeHexa( nodePosition,
            faceToNodes,
            elemToFaces,
            elemCenter,
            elemVolume,
            elemPerm,
            lengthTolerance,
            InnerProductType::TPFA,
            transMatrixRef );

  stackArray2d< real64, NF *NF > transMatrix( NF, NF );
  array1d< real64 > transMultiplier( NF );
  transMultiplier.setValues< parallelHostPolicy >( 1.0 );

  stackArray1d< real64, 3 > center( 3 );
  center[0] = elemCenter[0];
  center[1] = elemCenter[1];
  center[2] = elemCenter[2];
  real64 const perm[ 3 ] = { elemPerm[0], elemPerm[1], elemPerm[2] };

  TPFAInnerProduct::compute< NF >( nodePosition.toViewConst(),
                                   transMultiplier.toViewConst(),
                                   faceToNodes.toViewConst(),
                                   elemToFaces.toSliceConst(),
                                   center,
                                   elemVolume,
                                   perm,
                                   lengthTolerance,
                                   transMatrix.toSlice() );

  compareTransmissibilityMatrices( transMatrix, transMatrixRef );

  runConsistencyTest< NF >( nodePosition,
                            faceToNodes,
                            elemPerm,
                            transMatrix.toViewConst() );
}


TEST( testMimeticInnerProducts, QTPFA_hexa )
{
  localIndex constexpr NF = 6;

  array2d< real64, nodes::REFERENCE_POSITION_PERM > nodePosition;
  FaceManager::NodeMapType faceToNodes;
  array1d< localIndex > elemToFaces;
  real64 elemCenter[3];
  real64 elemPerm[3];
  real64 elemVolume = 0;
  real64 lengthTolerance = 0;
  stackArray2d< real64, NF *NF > transMatrixRef( NF, NF );

  makeHexa( nodePosition,
            faceToNodes,
            elemToFaces,
            elemCenter,
            elemVolume,
            elemPerm,
            lengthTolerance,
            InnerProductType::QUASI_TPFA_WITH_MULTIPLIERS,
            transMatrixRef );

  stackArray2d< real64, NF *NF > transMatrix( NF, NF );
  array1d< real64 > transMultiplier( NF );
  transMultiplier.setValues< parallelHostPolicy >( 1.0 );
  transMultiplier[0] = 0.9;
  transMultiplier[5] = 0.1;
  transMultiplier[3] = 0.8;

  stackArray1d< real64, 3 > center( 3 );
  center[0] = elemCenter[0];
  center[1] = elemCenter[1];
  center[2] = elemCenter[2];
  real64 const perm[ 3 ] = { elemPerm[0], elemPerm[1], elemPerm[2] };

  QuasiTPFAInnerProduct::compute< NF >( nodePosition.toViewConst(),
                                        transMultiplier.toViewConst(),
                                        faceToNodes.toViewConst(),
                                        elemToFaces.toSliceConst(),
                                        center,
                                        elemVolume,
                                        perm,
                                        lengthTolerance,
                                        transMatrix.toSlice() );

  compareTransmissibilityMatrices( transMatrix, transMatrixRef );

  runConsistencyTest< NF >( nodePosition,
                            faceToNodes,
                            elemPerm,
                            transMatrix.toViewConst() );
}

TEST( testMimeticInnerProducts, Simple_hexa )
{
  localIndex constexpr NF = 6;

  array2d< real64, nodes::REFERENCE_POSITION_PERM > nodePosition;
  FaceManager::NodeMapType faceToNodes;
  array1d< localIndex > elemToFaces;
  real64 elemCenter[3] = { 0.0 };
  real64 elemPerm[3] = { 0.0 };
  real64 elemVolume = 0;
  real64 lengthTolerance = 0;
  stackArray2d< real64, NF *NF > transMatrixRef( NF, NF );

  makeHexa( nodePosition,
            faceToNodes,
            elemToFaces,
            elemCenter,
            elemVolume,
            elemPerm,
            lengthTolerance,
            InnerProductType::SIMPLE_WITH_MULTIPLIERS,
            transMatrixRef );

  stackArray2d< real64, NF *NF > transMatrix( NF, NF );
  array1d< real64 > transMultiplier( NF );
  transMultiplier.setValues< parallelHostPolicy >( 1.0 );
  transMultiplier[0] = 0.9;
  transMultiplier[5] = 0.1;
  transMultiplier[3] = 0.8;


  stackArray1d< real64, 3 > center( 3 );
  center[0] = elemCenter[0];
  center[1] = elemCenter[1];
  center[2] = elemCenter[2];
  real64 const perm[ 3 ] = { elemPerm[0], elemPerm[1], elemPerm[2] };

  SimpleInnerProduct::compute< NF >( nodePosition.toViewConst(),
                                     transMultiplier.toViewConst(),
                                     faceToNodes.toViewConst(),
                                     elemToFaces.toSliceConst(),
                                     center,
                                     elemVolume,
                                     perm,
                                     lengthTolerance,
                                     transMatrix.toSlice() );

  compareTransmissibilityMatrices( transMatrix, transMatrixRef );

  runConsistencyTest< NF >( nodePosition,
                            faceToNodes,
                            elemPerm,
                            transMatrix.toViewConst() );
}

TEST( testMimeticInnerProducts, BdVLM_hexa )
{
  localIndex constexpr NF = 6;

  array2d< real64, nodes::REFERENCE_POSITION_PERM > nodePosition;
  FaceManager::NodeMapType faceToNodes;
  array1d< localIndex > elemToFaces;
  real64 elemCenter[3] = { 0.0 };
  real64 elemPerm[3] = { 0.0 };
  real64 elemVolume = 0;
  real64 lengthTolerance = 0;
  stackArray2d< real64, NF *NF > transMatrixRef( NF, NF );

  makeHexa( nodePosition,
            faceToNodes,
            elemToFaces,
            elemCenter,
            elemVolume,
            elemPerm,
            lengthTolerance,
            InnerProductType::BDVLM_WITH_MULTIPLIERS,
            transMatrixRef );

  stackArray2d< real64, NF *NF > transMatrix( NF, NF );
  array1d< real64 > transMultiplier( NF );
  transMultiplier.setValues< parallelHostPolicy >( 1.0 );
  transMultiplier[0] = 0.9;
  transMultiplier[5] = 0.1;
  transMultiplier[3] = 0.8;


  stackArray1d< real64, 3 > center( 3 );
  center[0] = elemCenter[0];
  center[1] = elemCenter[1];
  center[2] = elemCenter[2];
  real64 const perm[ 3 ] = { elemPerm[0], elemPerm[1], elemPerm[2] };

  BdVLMInnerProduct::compute< NF >( nodePosition.toViewConst(),
                                    transMultiplier.toViewConst(),
                                    faceToNodes.toViewConst(),
                                    elemToFaces.toSliceConst(),
                                    center,
                                    elemVolume,
                                    perm,
                                    lengthTolerance,
                                    transMatrix.toSlice() );

  compareTransmissibilityMatrices( transMatrix, transMatrixRef );

  runConsistencyTest< NF >( nodePosition,
                            faceToNodes,
                            elemPerm,
                            transMatrix.toViewConst() );
}


TEST( testMimeticInnerProducts, TPFA_tetra )
{
  localIndex constexpr NF = 4;

  array2d< real64, nodes::REFERENCE_POSITION_PERM > nodePosition;
  FaceManager::NodeMapType faceToNodes;
  array1d< localIndex > elemToFaces;
  real64 elemCenter[3];
  real64 elemPerm[3];
  real64 elemVolume = 0;
  real64 lengthTolerance = 0;
  stackArray2d< real64, NF *NF > transMatrixRef( NF, NF );

  makeTetra( nodePosition,
             faceToNodes,
             elemToFaces,
             elemCenter,
             elemVolume,
             elemPerm,
             lengthTolerance,
             InnerProductType::TPFA,
             transMatrixRef );

  stackArray2d< real64, NF *NF > transMatrix( NF, NF );
  array1d< real64 > transMultiplier( NF );
  transMultiplier.setValues< parallelHostPolicy >( 1.0 );

  stackArray1d< real64, 3 > center( 3 );
  center[0] = elemCenter[0];
  center[1] = elemCenter[1];
  center[2] = elemCenter[2];
  real64 const perm[ 3 ] = { elemPerm[0], elemPerm[1], elemPerm[2] };

  TPFAInnerProduct::compute< NF >( nodePosition.toViewConst(),
                                   transMultiplier.toViewConst(),
                                   faceToNodes.toViewConst(),
                                   elemToFaces.toSliceConst(),
                                   center,
                                   elemVolume,
                                   perm,
                                   lengthTolerance,
                                   transMatrix.toSlice() );

  compareTransmissibilityMatrices( transMatrix, transMatrixRef );

  runConsistencyTest< NF >( nodePosition,
                            faceToNodes,
                            elemPerm,
                            transMatrix.toViewConst() );
}


TEST( testMimeticInnerProducts, QTPFA_tetra )
{
  localIndex constexpr NF = 4;

  array2d< real64, nodes::REFERENCE_POSITION_PERM > nodePosition;
  FaceManager::NodeMapType faceToNodes;
  array1d< localIndex > elemToFaces;
  real64 elemCenter[3];
  real64 elemPerm[3];
  real64 elemVolume = 0;
  real64 lengthTolerance = 0;
  stackArray2d< real64, NF *NF > transMatrixRef( NF, NF );

  makeTetra( nodePosition,
             faceToNodes,
             elemToFaces,
             elemCenter,
             elemVolume,
             elemPerm,
             lengthTolerance,
             InnerProductType::QUASI_TPFA,
             transMatrixRef );

  stackArray2d< real64, NF *NF > transMatrix( NF, NF );
  array1d< real64 > transMultiplier( NF );
  transMultiplier.setValues< parallelHostPolicy >( 1.0 );

  stackArray1d< real64, 3 > center( 3 );
  center[0] = elemCenter[0];
  center[1] = elemCenter[1];
  center[2] = elemCenter[2];
  real64 const perm[ 3 ] = { elemPerm[0], elemPerm[1], elemPerm[2] };

  QuasiTPFAInnerProduct::compute< NF >( nodePosition.toViewConst(),
                                        transMultiplier.toViewConst(),
                                        faceToNodes.toViewConst(),
                                        elemToFaces.toSliceConst(),
                                        center,
                                        elemVolume,
                                        perm,
                                        lengthTolerance,
                                        transMatrix.toSlice() );

  compareTransmissibilityMatrices( transMatrix, transMatrixRef );

  runConsistencyTest< NF >( nodePosition,
                            faceToNodes,
                            elemPerm,
                            transMatrix.toViewConst() );
}

TEST( testMimeticInnerProducts, Simple_tetra )
{
  localIndex constexpr NF = 4;

  array2d< real64, nodes::REFERENCE_POSITION_PERM > nodePosition;
  FaceManager::NodeMapType faceToNodes;
  array1d< localIndex > elemToFaces;
  real64 elemCenter[3] = { 0.0 };
  real64 elemPerm[3] = { 0.0 };
  real64 elemVolume = 0;
  real64 lengthTolerance = 0;
  stackArray2d< real64, NF *NF > transMatrixRef( NF, NF );

  makeTetra( nodePosition,
             faceToNodes,
             elemToFaces,
             elemCenter,
             elemVolume,
             elemPerm,
             lengthTolerance,
             InnerProductType::SIMPLE,
             transMatrixRef );

  stackArray2d< real64, NF *NF > transMatrix( NF, NF );
  array1d< real64 > transMultiplier( NF );
  transMultiplier.setValues< parallelHostPolicy >( 1.0 );

  stackArray1d< real64, 3 > center( 3 );
  center[0] = elemCenter[0];
  center[1] = elemCenter[1];
  center[2] = elemCenter[2];
  real64 const perm[ 3 ] = { elemPerm[0], elemPerm[1], elemPerm[2] };

  SimpleInnerProduct::compute< NF >( nodePosition.toViewConst(),
                                     transMultiplier.toViewConst(),
                                     faceToNodes.toViewConst(),
                                     elemToFaces.toSliceConst(),
                                     center,
                                     elemVolume,
                                     perm,
                                     lengthTolerance,
                                     transMatrix.toSlice() );

  compareTransmissibilityMatrices( transMatrix, transMatrixRef );

  runConsistencyTest< NF >( nodePosition,
                            faceToNodes,
                            elemPerm,
                            transMatrix.toViewConst() );
}

TEST( testMimeticInnerProducts, BdVLMtetra )
{
  localIndex constexpr NF = 4;

  array2d< real64, nodes::REFERENCE_POSITION_PERM > nodePosition;
  FaceManager::NodeMapType faceToNodes;
  array1d< localIndex > elemToFaces;
  real64 elemCenter[3] = { 0.0 };
  real64 elemPerm[3] = { 0.0 };
  real64 elemVolume = 0;
  real64 lengthTolerance = 0;
  stackArray2d< real64, NF *NF > transMatrixRef( NF, NF );

  makeTetra( nodePosition,
             faceToNodes,
             elemToFaces,
             elemCenter,
             elemVolume,
             elemPerm,
             lengthTolerance,
             InnerProductType::BDVLM,
             transMatrixRef );

  stackArray2d< real64, NF *NF > transMatrix( NF, NF );
  array1d< real64 > transMultiplier( NF );
  transMultiplier.setValues< parallelHostPolicy >( 1.0 );

  stackArray1d< real64, 3 > center( 3 );
  center[0] = elemCenter[0];
  center[1] = elemCenter[1];
  center[2] = elemCenter[2];
  real64 const perm[ 3 ] = { elemPerm[0], elemPerm[1], elemPerm[2] };

  BdVLMInnerProduct::compute< NF >( nodePosition.toViewConst(),
                                    transMultiplier.toViewConst(),
                                    faceToNodes.toViewConst(),
                                    elemToFaces.toSliceConst(),
                                    center,
                                    elemVolume,
                                    perm,
                                    lengthTolerance,
                                    transMatrix.toSlice() );

  compareTransmissibilityMatrices( transMatrix, transMatrixRef );

  runConsistencyTest< NF >( nodePosition,
                            faceToNodes,
                            elemPerm,
                            transMatrix.toViewConst() );
}

TEST( testMimeticInnerProducts, TPFAM_cube )
{
  localIndex constexpr NF = 6;

  array2d< real64, nodes::REFERENCE_POSITION_PERM > nodePosition;
  FaceManager::NodeMapType faceToNodes;
  array1d< localIndex > elemToFaces;
  real64 elemCenter[3] = { 0.0 };
  real64 elemPerm[3] = { 0.0 };
  real64 elemVolume = 0;
  real64 lengthTolerance = 0;

  stackArray2d< real64, NF * NF > Mref( NF, NF );
  Mref.setValues< parallelHostPolicy >( 0.0 );

  makeCube( nodePosition,
            faceToNodes,
            elemToFaces,
            elemCenter,
            elemVolume,
            elemPerm,
            lengthTolerance,
            InnerProductType::TPFA,
            Mref.toSlice() );

  stackArray2d< real64, NF * NF > M( NF, NF );

  stackArray1d< real64, 3 > center( 3 );
  center[0] = elemCenter[0];
  center[1] = elemCenter[1];
  center[2] = elemCenter[2];

  real64 const perm[ 3 ] = { elemPerm[0], elemPerm[1], elemPerm[2] };

  TPFAInnerProduct::computeM< NF >( nodePosition.toViewConst(),
                                    faceToNodes.toViewConst(),
                                    elemToFaces.toSliceConst(),
                                    center,
                                    1,
                                    perm,
                                    lengthTolerance,
                                    M.toSlice() );

  compareInnerProductMatrices( M.toSliceConst(),
                               Mref.toSliceConst() );
}

TEST( testMimeticInnerProducts, QuasiTPFAM_cube )
{
  localIndex constexpr NF = 6;

  array2d< real64, nodes::REFERENCE_POSITION_PERM > nodePosition;
  FaceManager::NodeMapType faceToNodes;
  array1d< localIndex > elemToFaces;
  real64 elemCenter[3] = { 0.0 };
  real64 elemPerm[3] = { 0.0 };
  real64 elemVolume = 0;
  real64 lengthTolerance = 0;

  stackArray2d< real64, NF * NF > Mref( NF, NF );
  Mref.setValues< parallelHostPolicy >( 0.0 );

  makeCube( nodePosition,
            faceToNodes,
            elemToFaces,
            elemCenter,
            elemVolume,
            elemPerm,
            lengthTolerance,
            InnerProductType::QUASI_TPFA,
            Mref.toSlice() );

  stackArray2d< real64, NF * NF > M( NF, NF );

  stackArray1d< real64, 3 > center( 3 );
  center[0] = elemCenter[0];
  center[1] = elemCenter[1];
  center[2] = elemCenter[2];

  real64 const perm[ 3 ] = { elemPerm[0], elemPerm[1], elemPerm[2] };

  QuasiTPFAInnerProduct::computeM< NF >( nodePosition.toViewConst(),
                                         faceToNodes.toViewConst(),
                                         elemToFaces.toSliceConst(),
                                         center,
                                         1,
                                         perm,
                                         lengthTolerance,
                                         M.toSlice() );

  compareInnerProductMatrices( M.toSliceConst(),
                               Mref.toSliceConst() );
}

TEST( testMimeticInnerProducts, SimpleM_cube )
{
  localIndex constexpr NF = 6;

  array2d< real64, nodes::REFERENCE_POSITION_PERM > nodePosition;
  FaceManager::NodeMapType faceToNodes;
  array1d< localIndex > elemToFaces;
  real64 elemCenter[3] = { 0.0 };
  real64 elemPerm[3] = { 0.0 };
  real64 elemVolume = 0;
  real64 lengthTolerance = 0;

  stackArray2d< real64, NF * NF > Mref( NF, NF );
  Mref.setValues< parallelHostPolicy >( 0.0 );

  makeCube( nodePosition,
            faceToNodes,
            elemToFaces,
            elemCenter,
            elemVolume,
            elemPerm,
            lengthTolerance,
            InnerProductType::SIMPLE,
            Mref.toSlice() );

  stackArray2d< real64, NF * NF > M( NF, NF );

  stackArray1d< real64, 3 > center( 3 );
  center[0] = elemCenter[0];
  center[1] = elemCenter[1];
  center[2] = elemCenter[2];

  real64 const perm[ 3 ] = { elemPerm[0], elemPerm[1], elemPerm[2] };

  SimpleInnerProduct::computeM< NF >( nodePosition.toViewConst(),
                                      faceToNodes.toViewConst(),
                                      elemToFaces.toSliceConst(),
                                      center,
                                      1,
                                      perm,
                                      lengthTolerance,
                                      M.toSlice() );

  compareInnerProductMatrices( M.toSliceConst(),
                               Mref.toSliceConst() );
}

TEST( testMimeticInnerProducts, BdVLMM_cube )
{
  localIndex constexpr NF = 6;

  array2d< real64, nodes::REFERENCE_POSITION_PERM > nodePosition;
  FaceManager::NodeMapType faceToNodes;
  array1d< localIndex > elemToFaces;
  real64 elemCenter[3] = { 0.0 };
  real64 elemPerm[3] = { 0.0 };
  real64 elemVolume = 0;
  real64 lengthTolerance = 0;

  stackArray2d< real64, NF * NF > M( NF, NF );
  M.setValues< parallelHostPolicy >( 0.0 );

  stackArray2d< real64, NF * NF > Mref( NF, NF );
  Mref.setValues< parallelHostPolicy >( 0.0 );

  makeCube( nodePosition,
            faceToNodes,
            elemToFaces,
            elemCenter,
            elemVolume,
            elemPerm,
            lengthTolerance,
            InnerProductType::BDVLM,
            Mref.toSlice() );

  stackArray1d< real64, 3 > center( 3 );
  center[0] = elemCenter[0];
  center[1] = elemCenter[1];
  center[2] = elemCenter[2];

  real64 const perm[ 3 ] = { elemPerm[0], elemPerm[1], elemPerm[2] };

  BdVLMInnerProduct::computeM< NF >( nodePosition.toViewConst(),
                                     faceToNodes.toViewConst(),
                                     elemToFaces.toSliceConst(),
                                     center,
                                     elemVolume,
                                     perm,
                                     lengthTolerance,
                                     M.toSlice() );

  compareInnerProductMatrices( M.toSliceConst(),
                               Mref.toSliceConst() );
}


//======================== Linear Pressure Recovery Test =============================
// Three cases: without distortion (unit cube cell), with distortion: (1) planar (2) nonplanar
enum class DistortionMode { None, Planar, NonPlanar };

static inline void makeUnitCube( double x_start,
                                 array2d< real64, nodes::REFERENCE_POSITION_PERM > & node,
                                 FaceManager::NodeMapType & faceTonodes,
                                 array1d< localIndex > & elemTofaces,
                                 real64 (& center)[3], real64 & vol )
{
  node.resize( 8, 3 );
  // vertices on z = 0 plane
  node( 0, 0 ) = x_start + 0; node( 0, 1 ) = 0; node( 0, 2 ) = 0;
  node( 1, 0 ) = x_start + 1; node( 1, 1 ) = 0; node( 1, 2 ) = 0;
  node( 2, 0 ) = x_start + 0; node( 2, 1 ) = 1; node( 2, 2 ) = 0;
  node( 3, 0 ) = x_start + 1; node( 3, 1 ) = 1; node( 3, 2 ) = 0;

  // vertices on z = 1 plane
  node( 4, 0 ) = x_start + 0; node( 4, 1 ) = 0; node( 4, 2 ) = 1;
  node( 5, 0 ) = x_start + 1; node( 5, 1 ) = 0; node( 5, 2 ) = 1;
  node( 6, 0 ) = x_start + 0; node( 6, 1 ) = 1; node( 6, 2 ) = 1;
  node( 7, 0 ) = x_start + 1; node( 7, 1 ) = 1; node( 7, 2 ) = 1;

  faceTonodes.resize( 6 );
  for( int f = 0; f < 6; ++f )
  {
    faceTonodes.resizeArray( f, 4 );
  }

  // define each face by its 4 node indices (consistent ordering within the cell)
  faceTonodes( 0, 0 ) = 0; faceTonodes( 0, 1 ) = 1; faceTonodes( 0, 2 ) = 3; faceTonodes( 0, 3 ) = 2; // z=0
  faceTonodes( 1, 0 ) = 0; faceTonodes( 1, 1 ) = 4; faceTonodes( 1, 2 ) = 5; faceTonodes( 1, 3 ) = 1; // y=0
  faceTonodes( 2, 0 ) = 0; faceTonodes( 2, 1 ) = 2; faceTonodes( 2, 2 ) = 6; faceTonodes( 2, 3 ) = 4; // x=0
  faceTonodes( 3, 0 ) = 1; faceTonodes( 3, 1 ) = 5; faceTonodes( 3, 2 ) = 7; faceTonodes( 3, 3 ) = 3; // x=1
  faceTonodes( 4, 0 ) = 6; faceTonodes( 4, 1 ) = 2; faceTonodes( 4, 2 ) = 3; faceTonodes( 4, 3 ) = 7; // y=1
  faceTonodes( 5, 0 ) = 4; faceTonodes( 5, 1 ) = 6; faceTonodes( 5, 2 ) = 7; faceTonodes( 5, 3 ) = 5; // z=1

  elemTofaces.resize( 6 );
  for( int i = 0; i < 6; ++i )
  {
    elemTofaces( i ) = i;
  }

  // compute cell center and volume
  array1d< localIndex > order;
  order.resize( 8 );
  order( 0 ) = 0;
  order( 1 ) = 4;
  order( 2 ) = 2;
  order( 3 ) = 6;
  order( 4 ) = 1;
  order( 5 ) = 5;
  order( 6 ) = 3;
  order( 7 ) = 7;

  real64 X[8][3]; center[0] = center[1] = center[2] = 0;
  for( int a = 0; a < 8; ++a )
  {
    X[a][0] = node( order( a ), 0 );
    X[a][1] = node( order( a ), 1 );
    X[a][2] = node( order( a ), 2 );
    center[0] += X[a][0];
    center[1] += X[a][1];
    center[2] += X[a][2];
  }

  vol = hexahedronVolume( X );
  center[0] /= 8.0;
  center[1] /= 8.0;
  center[2] /= 8.0;
}

// create planar case for distortion test
static inline void makeDistortedPlanar( array2d< real64, nodes::REFERENCE_POSITION_PERM > & nodeL,
                                        array2d< real64, nodes::REFERENCE_POSITION_PERM > & nodeR,
                                        FaceManager::NodeMapType const & faceL,
                                        FaceManager::NodeMapType const & faceR,
                                        real64 eps )
{
  // indicate shared face
  localIndex const fL = 3, fR = 2;
  // L face(3) node order: {1,5,7,3}
  int Lvert1 = faceL( fL, 2 ), Lvert2 = faceL( fL, 3 );
  // R face(2) node order: {0,2,6,4}
  int Rvert1 = faceR( fR, 1 ), Rvert2 = faceR( fR, 2 );

  nodeL( Lvert1, 0 ) += eps;  nodeR( Rvert1, 0 ) += eps;
  nodeL( Lvert2, 0 ) += eps;  nodeR( Rvert2, 0 ) += eps;
}

// create nonplanar case for distortion test
static inline void makeDistortedNonplanar( array2d< real64, nodes::REFERENCE_POSITION_PERM > & nodeL,
                                           array2d< real64, nodes::REFERENCE_POSITION_PERM > & nodeR,
                                           FaceManager::NodeMapType const & faceL,
                                           FaceManager::NodeMapType const & faceR,
                                           real64 eps )
{
  // indicate shared face
  localIndex const fL = 3, fR = 2;

  int Lvert = faceL( fL, 2 );
  int Rvert = faceR( fR, 2 );

  nodeL( Lvert, 0 ) += eps;
  nodeR( Rvert, 0 ) -= eps / 2;
}

// return the sum of all entries in row r of matrix T
static inline real64 rowSum( arraySlice2d< real64 const > T,
                             localIndex row )
{
  real64 sum = 0;
  // iterate over columns
  for( localIndex j = 0; j < T.size( 1 ); ++j )
  {
    sum += T( row, j );
  }
  return sum;
}

// recompute volume and center of distorted cell
static inline void computeDistortedVolumeAndCenter( array2d< real64, nodes::REFERENCE_POSITION_PERM > const & node,
                                                    real64 (& center)[3], real64 & vol )
{
  array1d< localIndex > order;
  order.resize( 8 );
  order( 0 ) = 0;
  order( 1 ) = 4;
  order( 2 ) = 2;
  order( 3 ) = 6;
  order( 4 ) = 1;
  order( 5 ) = 5;
  order( 6 ) = 3;
  order( 7 ) = 7;

  real64 X[8][3];
  center[0] = center[1] = center[2] = 0.0;
  for( int a = 0; a < 8; ++a )
  {
    X[a][0] = node( order( a ), 0 );
    X[a][1] = node( order( a ), 1 );
    X[a][2] = node( order( a ), 2 );
    center[0] += X[a][0];
    center[1] += X[a][1];
    center[2] += X[a][2];
  }
  vol = hexahedronVolume( X );
  center[0] /= 8.0; center[1] /= 8.0; center[2] /= 8.0;
}

static inline bool isBoundaryFace( FaceManager::NodeMapType const & faceTonode,
                                   array2d< real64, nodes::REFERENCE_POSITION_PERM > const & node,
                                   localIndex f )
{
  real64 fc[3], fn[3];
  centroid_3DPolygon( faceTonode[f], node.toViewConst(), fc, fn );
  constexpr real64 tol = 1e-12;

  return (std::abs( fc[0] - 0.0 ) < tol) || (std::abs( fc[0] - 2.0 ) < tol) ||
         (std::abs( fc[1] - 0.0 ) < tol) || (std::abs( fc[1] - 1.0 ) < tol) ||
         (std::abs( fc[2] - 0.0 ) < tol) || (std::abs( fc[2] - 1.0 ) < tol);
}

// calculate pressures on cell boundary (for Dirichlet boundary condition)
static inline real64 computeBoundaryPressure( FaceManager::NodeMapType const & faceTonode,
                                              array2d< real64, nodes::REFERENCE_POSITION_PERM > const & node,
                                              localIndex f,
                                              real64 alpha )
{
  real64 fc[3], fn[3];
  centroid_3DPolygon( faceTonode[f], node.toViewConst(), fc, fn );
  return fc[0] + alpha * fc[1];
}

// calculate pressure values based on the transmissbility matrix T and compare against the analytical values
template< int NF >
static double computeLinearPressure_error( int ipKind,
                                           DistortionMode mode = DistortionMode::None,
                                           real64 eps = 0.0 )
{
  array2d< real64, nodes::REFERENCE_POSITION_PERM > node_L;
  array2d< real64, nodes::REFERENCE_POSITION_PERM > node_R;
  FaceManager::NodeMapType faceTonode_L;
  FaceManager::NodeMapType faceTonode_R;
  array1d< localIndex > elemToface_L;
  array1d< localIndex > elemToface_R;
  real64 center_L[3], center_R[3], vol_L = 0, vol_R = 0;

  makeUnitCube( 0.0, node_L, faceTonode_L, elemToface_L, center_L, vol_L );
  makeUnitCube( 1.0, node_R, faceTonode_R, elemToface_R, center_R, vol_R );

  if( mode == DistortionMode::Planar )
  {
    makeDistortedPlanar( node_L, node_R, faceTonode_L, faceTonode_R, eps );
  }
  else if( mode == DistortionMode::NonPlanar )
  {
    makeDistortedNonplanar( node_L, node_R, faceTonode_L, faceTonode_R, eps );
  }

  computeDistortedVolumeAndCenter( node_L, center_L, vol_L );
  computeDistortedVolumeAndCenter( node_R, center_R, vol_R );

  // indicate the shared face
  localIndex const fL_int = 3;
  localIndex const fR_int = 2;

  real64 alpha_lin = 0.3; // pressure = x * alpha_lin * y
  auto isDirL = [&]( localIndex f ){ return isBoundaryFace( faceTonode_L, node_L, f ); };
  auto isDirR = [&]( localIndex f ){ return isBoundaryFace( faceTonode_R, node_R, f ); };

  // compute the transmissibility matrix T
  constexpr real64 ltol = 1e-12;
  real64 Kvec[3] = {4.0, 1.0, 0.5};
  stackArray2d< real64, NF * NF > TL( NF, NF ), TR( NF, NF );

  stackArray1d< real64, 3 > cLc( 3 ), cRc( 3 );
  for( int d = 0; d < 3; ++d )
  {
    cLc[d] = center_L[d];
    cRc[d] = center_R[d];
  }

  array1d< real64 > mult( NF );
  mult.setValues< parallelHostPolicy >( 1.0 );

  if( ipKind == InnerProductType::TPFA )
  {
    TPFAInnerProduct::compute< NF >( node_L.toViewConst(),
                                     mult.toViewConst(),
                                     faceTonode_L.toViewConst(),
                                     elemToface_L.toSliceConst(),
                                     cLc, vol_L, Kvec, ltol, TL.toSlice() );
    TPFAInnerProduct::compute< NF >( node_R.toViewConst(),
                                     mult.toViewConst(),
                                     faceTonode_R.toViewConst(),
                                     elemToface_R.toSliceConst(),
                                     cRc, vol_R, Kvec, ltol, TR.toSlice() );
  }
  else if( ipKind == InnerProductType::QUASI_TPFA )
  {
    QuasiTPFAInnerProduct::compute< NF >( node_L.toViewConst(),
                                          mult.toViewConst(),
                                          faceTonode_L.toViewConst(),
                                          elemToface_L.toSliceConst(),
                                          cLc, vol_L, Kvec, ltol, TL.toSlice() );
    QuasiTPFAInnerProduct::compute< NF >( node_R.toViewConst(),
                                          mult.toViewConst(),
                                          faceTonode_R.toViewConst(),
                                          elemToface_R.toSliceConst(),
                                          cRc, vol_R, Kvec, ltol, TR.toSlice() );
  }
  else if( ipKind == InnerProductType::SIMPLE )
  {
    SimpleInnerProduct::compute< NF >( node_L.toViewConst(),
                                       mult.toViewConst(),
                                       faceTonode_L.toViewConst(),
                                       elemToface_L.toSliceConst(),
                                       cLc, vol_L, Kvec, ltol, TL.toSlice() );
    SimpleInnerProduct::compute< NF >( node_R.toViewConst(),
                                       mult.toViewConst(),
                                       faceTonode_R.toViewConst(),
                                       elemToface_R.toSliceConst(),
                                       cRc, vol_R, Kvec, ltol, TR.toSlice() );
  }
  else if( ipKind == InnerProductType::BDVLM )
  {
    BdVLMInnerProduct::compute< NF >( node_L.toViewConst(),
                                      mult.toViewConst(),
                                      faceTonode_L.toViewConst(),
                                      elemToface_L.toSliceConst(),
                                      cLc, vol_L, Kvec, ltol, TL.toSlice() );
    BdVLMInnerProduct::compute< NF >( node_R.toViewConst(),
                                      mult.toViewConst(),
                                      faceTonode_R.toViewConst(),
                                      elemToface_R.toSliceConst(),
                                      cRc, vol_R, Kvec, ltol, TR.toSlice() );
  }

  // assemble linear system to compute pressure values
  auto assembleCell = [&]( arraySlice2d< real64 const > T,
                           FaceManager::NodeMapType const & faceTonode,
                           array2d< real64, nodes::REFERENCE_POSITION_PERM > const & node,
                           localIndex fInt,
                           bool leftCell,
                           real64 & a, real64 & alpha, real64 & r )
  {
    a = alpha = r = 0.0;
    for( localIndex i = 0; i < NF; ++i )
    {
      bool isDir = leftCell ? isDirL( i ) : isDirR( i );
      if( !isDir && i != fInt )
        continue;
      a += rowSum( T, i );
      alpha += T( i, fInt );

      for( localIndex j = 0; j < NF; ++j )
      {
        if( j == fInt )
          continue;
        bool jDir = leftCell ? isDirL( j ) : isDirR( j );
        if( jDir )
          r += T( i, j ) * computeBoundaryPressure( faceTonode, node, j, alpha_lin );
      }
    }
  };

  auto sumFaceRow = [&]( arraySlice2d< real64 const > T,
                         FaceManager::NodeMapType const & faceToNode,
                         array2d< real64, nodes::REFERENCE_POSITION_PERM > const & node,
                         localIndex fInt, bool leftCell ) -> real64
  {
    real64 sum = 0;
    for( localIndex j = 0; j < NF; ++j )
    {
      bool jDir = leftCell ? isDirL( j ) : isDirR( j );
      if( jDir )
        sum += T( fInt, j ) * computeBoundaryPressure( faceToNode, node, j, alpha_lin );
    }
    return sum;
  };

  real64 aL, alphaL, rL, aR, alphaR, rR;
  assembleCell( TL.toSliceConst(), faceTonode_L, node_L, fL_int, true, aL, alphaL, rL );
  assembleCell( TR.toSliceConst(), faceTonode_R, node_R, fR_int, false, aR, alphaR, rR );

  real64 betaL = rowSum( TL.toSliceConst(), fL_int );
  real64 betaR = rowSum( TR.toSliceConst(), fR_int );
  real64 gamma = TL( fL_int, fL_int ) + TR( fR_int, fR_int );
  real64 boundary_sum = sumFaceRow( TL.toSliceConst(), faceTonode_L, node_L, fL_int, true )
                        + sumFaceRow( TR.toSliceConst(), faceTonode_R, node_R, fR_int, false );

  real64 A11 = aL, A22 = aR, A13 = -alphaL, A23 = -alphaR, A31 = betaL, A32 = betaR, A33 = -gamma;
  real64 b1 = rL, b2 = rR, b3 = boundary_sum;

  real64 denom = A33 - A31 * A13 / A11 - A32 * A23 / A22;
  real64 lambdaI = ( b3 - A31 * (b1 / A11) - A32 * (b2 / A22) ) / denom;
  real64 pL = ( b1 - A13 * lambdaI ) / A11;
  real64 pR = ( b2 - A23 * lambdaI ) / A22;

  // compute analytical pressure values
  double pL_exact = center_L[0] + alpha_lin * center_L[1];
  double pR_exact = center_R[0] + alpha_lin * center_R[1];
  double err = std::max( std::abs( pL - pL_exact ) / std::abs( pL_exact ), std::abs( pR - pR_exact ) / std::abs( pR_exact ) );

  return err;
}

// ================= case 0: without distortion ===========================
TEST( MimeticIP_Linear, UnitCube_LinearPressure_TPFA )
{
  double err = computeLinearPressure_error< 6 >( InnerProductType::TPFA );
  EXPECT_LT( err, consistency_tol );
}

TEST( MimeticIP_Linear, UnitCube_LinearPressure_QuasiTPFA )
{
  double err = computeLinearPressure_error< 6 >( InnerProductType::QUASI_TPFA );
  EXPECT_LT( err, consistency_tol );
}

TEST( MimeticIP_Linear, UnitCube_LinearPressure_Simple )
{
  double err = computeLinearPressure_error< 6 >( InnerProductType::SIMPLE );
  EXPECT_LT( err, consistency_tol );
}

TEST( MimeticIP_Linear, UnitCube_LinearPressure_BdVLM )
{
  double err = computeLinearPressure_error< 6 >( InnerProductType::BDVLM );
  EXPECT_LT( err, consistency_tol );
}

// =================== case 1: with distortion (planar) ===========================
TEST( MimeticIP_Linear, Distortion_Planar_LinearPressure )
{
  int neps = 3;
  std::vector< double > eps_values( 3 );
  eps_values[0] = 0.0;   // no distortion
  eps_values[1] = 0.2;   // moderate distortion
  eps_values[2] = 0.9;   // severe distortion

  // all mimetic schemes are consistent
  for( int i = 0; i < neps; ++i )
  {
    double eps = eps_values[i];

    double errQTPFA = computeLinearPressure_error< 6 >( InnerProductType::QUASI_TPFA, DistortionMode::Planar, eps );
    EXPECT_LT( errQTPFA, consistency_tol );

    double errSIMPLE= computeLinearPressure_error< 6 >( InnerProductType::SIMPLE, DistortionMode::Planar, eps );
    EXPECT_LT( errSIMPLE, consistency_tol );

    double errBDVLM = computeLinearPressure_error< 6 >( InnerProductType::BDVLM, DistortionMode::Planar, eps );
    EXPECT_LT( errBDVLM, consistency_tol );
  }

  for( int i = 0; i < neps; ++i )
  {
    double eps = eps_values[i];

    double errTPFA  = computeLinearPressure_error< 6 >( InnerProductType::TPFA, DistortionMode::Planar, eps );

    if( i == 0 )
    {
      EXPECT_LT( errTPFA, consistency_tol );       // test TPFA is consistent with K-orthogonal grids
    }
    else
    {
      EXPECT_GT( errTPFA, consistency_tol );       // test TPFA is inconsistent with non K-orthogonal grids
    }
  }
}

// =================== case 2: with distortion (nonplanar) ===========================
TEST( MimeticIP_Linear, Distortion_NonPlanar_LinearPressure )
{
  int neps = 2;
  std::vector< double > eps_values( 2 );
  eps_values[0] = 0.2;   // moderate distortion
  eps_values[1] = 0.9;   // severe distortion

  // all schemes are inconsistent with nonplanar case
  for( int i = 0; i < neps; ++i )
  {
    double eps = eps_values[i];

    double errTPFA  = computeLinearPressure_error< 6 >( InnerProductType::TPFA, DistortionMode::NonPlanar, eps );
    EXPECT_GT( errTPFA, consistency_tol );

    double errQTPFA = computeLinearPressure_error< 6 >( InnerProductType::QUASI_TPFA, DistortionMode::NonPlanar, eps );
    EXPECT_GT( errQTPFA, consistency_tol );

    double errSIMPLE= computeLinearPressure_error< 6 >( InnerProductType::SIMPLE, DistortionMode::NonPlanar, eps );
    EXPECT_GT( errSIMPLE, consistency_tol );

    double errBDVLM = computeLinearPressure_error< 6 >( InnerProductType::BDVLM, DistortionMode::NonPlanar, eps );
    EXPECT_GT( errBDVLM, consistency_tol );
  }
}

//======================== Hydrostatic Equilibrium Consistency Test =============================
static inline void makeDistortedPlanar_gravity( array2d< real64, nodes::REFERENCE_POSITION_PERM > & nodeL,
                                                array2d< real64, nodes::REFERENCE_POSITION_PERM > & nodeU,
                                                FaceManager::NodeMapType const & faceL,
                                                FaceManager::NodeMapType const & faceU,
                                                localIndex fL, localIndex fU,
                                                real64 eps )
{
  int Lv1 = faceL( fL, 1 ), Lv2 = faceL( fL, 2 );
  int Rv1 = faceU( fU, 3 ), Rv2 = faceU( fU, 2 );

  nodeL( Lv1, 2 ) += eps;  nodeU( Rv1, 2 ) += eps;
  nodeL( Lv2, 2 ) += eps;  nodeU( Rv2, 2 ) += eps;
}

static inline void makeDistortedNonplanar_gravity( array2d< real64, nodes::REFERENCE_POSITION_PERM > & nodeL,
                                                   array2d< real64, nodes::REFERENCE_POSITION_PERM > & nodeU,
                                                   FaceManager::NodeMapType const & faceL,
                                                   FaceManager::NodeMapType const & faceU,
                                                   localIndex fL, localIndex fU,
                                                   real64 eps )
{
  int Lv = faceL( fL, 2 );
  int Rv = faceU( fU, 2 );

  nodeL( Lv, 2 ) += eps;
  nodeU( Rv, 2 ) += eps;
}

template< int NF >
static double computeGravityConsistency_error( int ipKind,
                                               DistortionMode mode = DistortionMode::None,
                                               real64 eps = 0.0,
                                               real64 rho = 1.0,
                                               real64 g = 1.0,
                                               real64 pD = 1.0 )
{
  // build two cells (z-stacked)
  array2d< real64, nodes::REFERENCE_POSITION_PERM > node_L, node_U;
  FaceManager::NodeMapType faceTonode_L, faceTonode_U;
  array1d< localIndex > elemToface_L, elemToface_U;
  real64 center_L[3], center_U[3], vol_L = 0, vol_U = 0;

  makeUnitCube( 0.0, node_L, faceTonode_L, elemToface_L, center_L, vol_L );
  makeUnitCube( 0.0, node_U, faceTonode_U, elemToface_U, center_U, vol_U );
  // stack upward of second cube (z + 1)
  for( int i = 0; i < 8; ++i )
  {
    node_U( i, 2 ) += 1.0;
  }

  // shared interface: lower top face fL = 5, upper bottom face fU = 0
  if( mode == DistortionMode::Planar )
  {
    makeDistortedPlanar_gravity( node_L, node_U, faceTonode_L, faceTonode_U, 5, 0, eps );
  }
  else if( mode == DistortionMode::NonPlanar )
  {
    makeDistortedNonplanar_gravity( node_L, node_U, faceTonode_L, faceTonode_U, 5, 0, eps );
  }

  computeDistortedVolumeAndCenter( node_L, center_L, vol_L );
  computeDistortedVolumeAndCenter( node_U, center_U, vol_U );

  // indicate the shared face for z-stacked cubes
  localIndex const fL_int = 5;   // top of lower cell
  localIndex const fU_int = 0;   // bottom of upper cell

  // build the matrix T per cell (reuse inner product computes)
  constexpr real64 ltol = 1e-12;
  real64 Kvec[3] = { 1.0, 1.0, 1.0 };
  stackArray2d< real64, NF * NF > TL( NF, NF ), TU( NF, NF );

  stackArray1d< real64, 3 > cLc( 3 ), cUc( 3 );
  for( int d = 0; d < 3; ++d )
  {
    cLc[d] = center_L[d];
    cUc[d] = center_U[d];
  }

  array1d< real64 > mult( NF );
  mult.setValues< parallelHostPolicy >( 1.0 );

  if( ipKind == InnerProductType::TPFA )
  {
    TPFAInnerProduct::compute< NF >( node_L.toViewConst(),
                                     mult.toViewConst(),
                                     faceTonode_L.toViewConst(),
                                     elemToface_L.toSliceConst(),
                                     cLc, vol_L, Kvec, ltol, TL.toSlice() );
    TPFAInnerProduct::compute< NF >( node_U.toViewConst(),
                                     mult.toViewConst(),
                                     faceTonode_U.toViewConst(),
                                     elemToface_U.toSliceConst(),
                                     cUc, vol_U, Kvec, ltol, TU.toSlice() );
  }
  else if( ipKind == InnerProductType::QUASI_TPFA )
  {
    QuasiTPFAInnerProduct::compute< NF >( node_L.toViewConst(),
                                          mult.toViewConst(),
                                          faceTonode_L.toViewConst(),
                                          elemToface_L.toSliceConst(),
                                          cLc, vol_L, Kvec, ltol, TL.toSlice() );
    QuasiTPFAInnerProduct::compute< NF >( node_U.toViewConst(),
                                          mult.toViewConst(),
                                          faceTonode_U.toViewConst(),
                                          elemToface_U.toSliceConst(),
                                          cUc, vol_U, Kvec, ltol, TU.toSlice() );
  }
  else if( ipKind == InnerProductType::SIMPLE )
  {
    SimpleInnerProduct::compute< NF >( node_L.toViewConst(),
                                       mult.toViewConst(),
                                       faceTonode_L.toViewConst(),
                                       elemToface_L.toSliceConst(),
                                       cLc, vol_L, Kvec, ltol, TL.toSlice() );
    SimpleInnerProduct::compute< NF >( node_U.toViewConst(),
                                       mult.toViewConst(),
                                       faceTonode_U.toViewConst(),
                                       elemToface_U.toSliceConst(),
                                       cUc, vol_U, Kvec, ltol, TU.toSlice() );
  }
  else if( ipKind == InnerProductType::BDVLM )
  {
    BdVLMInnerProduct::compute< NF >( node_L.toViewConst(),
                                      mult.toViewConst(),
                                      faceTonode_L.toViewConst(),
                                      elemToface_L.toSliceConst(),
                                      cLc, vol_L, Kvec, ltol, TL.toSlice() );
    BdVLMInnerProduct::compute< NF >( node_U.toViewConst(),
                                      mult.toViewConst(),
                                      faceTonode_U.toViewConst(),
                                      elemToface_U.toSliceConst(),
                                      cUc, vol_U, Kvec, ltol, TU.toSlice() );
  }

  auto faceCentroid = [&]( auto const & node,
                           auto const & f2n,
                           localIndex f,
                           real64 c[3] )
  {
    // fc: centroid, fn: area-weighted normal
    stackArray1d< real64, 3 > fc( 3 ), fn( 3 );
    computationalGeometry::centroid_3DPolygon( f2n[f], node.toViewConst(), fc.toSlice(), fn.toSlice() );

    c[0] = fc[0];
    c[1] = fc[1];
    c[2] = fc[2];
  };

  real64 xTop[3] = { 0.0, 0.0, 0.0 };
  {
    real64 zbest = std::numeric_limits< real64 >::lowest();
    for( localIndex f = 0; f < faceTonode_U.size(); ++f )
    {
      real64 c[3];
      faceCentroid( node_U, faceTonode_U, f, c );
      if( c[2] > zbest )
      {
        zbest = c[2];
        xTop[0] = c[0];
        xTop[1] = c[1];
        xTop[2] = c[2];
      }
    }
  }

  const real64 gvec[3] = { 0.0, 0.0, -g };
  const real64 p0 = pD + rho * ( gvec[0] * xTop[0] + gvec[1] * xTop[1] + gvec[2] * xTop[2] );

  auto buildCloc = [&]( auto const & node, auto const & f2n, real64 const Cc[3] )
  {
    stackArray2d< real64, NF * 3 > C( NF, 3 );
    for( localIndex k = 0; k < NF; ++k )
    {
      real64 cf[3];
      faceCentroid( node, f2n, k, cf );
      C( k, 0 ) = cf[0] - Cc[0];
      C( k, 1 ) = cf[1] - Cc[1];
      C( k, 2 ) = cf[2] - Cc[2];
    }
    return C;
  };
  auto CL = buildCloc( node_L, faceTonode_L, center_L );
  auto CU = buildCloc( node_U, faceTonode_U, center_U );

  // hydro flux per cell: u = T (pi - ep + C (rho g))
  auto hydroFlux = [&]( arraySlice2d< real64 const > Tloc,
                        auto const & C_loc,
                        auto const & node,
                        auto const & f2n,
                        real64 const Cc[3] )
  {
    // pE (cell center)
    real64 pE = p0 - rho * ( gvec[0] * Cc[0] + gvec[1] * Cc[1] + gvec[2] * Cc[2] );

    //piE (face center)
    stackArray1d< real64, NF > piE( NF );
    for( localIndex k = 0; k < NF; ++k )
    {
      real64 cf[3];
      faceCentroid( node, f2n, k, cf );
      piE[k] = p0 - rho * (gvec[0] * cf[0] + gvec[1] * cf[1] + gvec[2] * cf[2] );
    }

    // bE = C (rho g)
    real64 rg[3] = { rho * gvec[0], rho * gvec[1], rho * gvec[2] };
    stackArray1d< real64, NF > bE( NF );
    for( localIndex k = 0; k < NF; ++k )
    {
      bE[k] = C_loc( k, 0 ) * rg[0] + C_loc( k, 1 ) * rg[1] + C_loc( k, 2 ) * rg[2];
    }

    // rhs = piE - pE*1 + bE
    stackArray1d< real64, NF > rhs( NF ), uE( NF );
    for( localIndex k = 0; k < NF; ++k )
    {
      rhs[k] = piE[k] - pE + bE[k];
    }

    // uE = Tloc * rhs
    for( localIndex i = 0; i < NF; ++i )
    {
      real64 s = 0;
      for( localIndex j = 0; j < NF; ++j )
      {
        s += Tloc( i, j ) * rhs[j];
      }
      uE[i] = s;
    }

    return uE;
  };

  auto uL = hydroFlux( TL.toSliceConst(), CL.toSliceConst(), node_L, faceTonode_L, center_L );
  auto uU = hydroFlux( TU.toSliceConst(), CU.toSliceConst(), node_U, faceTonode_U, center_U );

  real64 maxErr = 0;
  for( localIndex i = 0; i < NF; ++i )
  {
    maxErr = std::max( maxErr, std::abs( uL[i] ) );
  }
  for( localIndex i = 0; i < NF; ++i )
  {
    maxErr = std::max( maxErr, std::abs( uU[i] ) );
  }

  // interface flux consistency (should be ~0)
  real64 faceErr = std::abs( uL[fL_int] + uU[fU_int] );

  return std::max( maxErr, faceErr );
}

// ================= case 0: without distortion ===========================
TEST( Hydrostatic, GravityConsistency_NoDistortion_TPFA )
{
  double err = computeGravityConsistency_error< 6 >( InnerProductType::TPFA, DistortionMode::None, 0.0, 1.0, 1.0, 1.0 );
  EXPECT_LE( err, consistency_tol );
}

TEST( Hydrostatic, GravityConsistency_NoDistortion_QuasiTPFA )
{
  double err = computeGravityConsistency_error< 6 >( InnerProductType::QUASI_TPFA, DistortionMode::None, 0.0, 1.0, 1.0, 1.0 );
  EXPECT_LE( err, consistency_tol );
}

TEST( Hydrostatic, GravityConsistency_NoDistortion_Simple )
{
  double err = computeGravityConsistency_error< 6 >( InnerProductType::SIMPLE, DistortionMode::None, 0.0, 1.0, 1.0, 1.0 );
  EXPECT_LE( err, consistency_tol );
}

TEST( Hydrostatic, GravityConsistency_NoDistortion_BdVLM )
{
  double err = computeGravityConsistency_error< 6 >( InnerProductType::BDVLM, DistortionMode::None, 0.0, 1.0, 1.0, 1.0 );
  EXPECT_LE( err, consistency_tol );
}

// =================== case 1: with distortion (planar) ===========================
TEST( Hydrostatic, GravityConsistency_Distortion_Planar )
{
  int neps = 2;
  std::vector< double > eps_values( 2 );
  eps_values[0] = 0.2;   // moderate distortion
  eps_values[1] = 0.9;   // severe distortion

  // all mimetic schemes are consistent
  for( int i = 0; i < neps; ++i )
  {
    double eps = eps_values[i];

    double errTPFA = computeGravityConsistency_error< 6 >( InnerProductType::TPFA, DistortionMode::Planar, eps );
    EXPECT_LE( errTPFA, consistency_tol );

    double errQTPFA = computeGravityConsistency_error< 6 >( InnerProductType::QUASI_TPFA, DistortionMode::Planar, eps );
    EXPECT_LE( errQTPFA, consistency_tol );

    double errSIMPLE = computeGravityConsistency_error< 6 >( InnerProductType::SIMPLE, DistortionMode::Planar, eps );
    EXPECT_LE( errSIMPLE, consistency_tol );

    double errBDVLM = computeGravityConsistency_error< 6 >( InnerProductType::BDVLM, DistortionMode::Planar, eps );
    EXPECT_LE( errBDVLM, consistency_tol );
  }
}

// =================== case 2: with distortion (nonplanar) ===========================
TEST( Hydrostatic, GravityConsistency_Distortion_NonPlanar )
{
  int neps = 2;
  std::vector< double > eps_values( 2 );
  eps_values[0] = 0.2;   // moderate distortion
  eps_values[1] = 0.9;   // severe distortion

  // all mimetic schemes are consistent
  for( int i = 0; i < neps; ++i )
  {
    double eps = eps_values[i];

    double errTPFA = computeGravityConsistency_error< 6 >( InnerProductType::TPFA, DistortionMode::NonPlanar, eps );
    EXPECT_LE( errTPFA, consistency_tol );

    double errQTPFA = computeGravityConsistency_error< 6 >( InnerProductType::QUASI_TPFA, DistortionMode::NonPlanar, eps );
    EXPECT_LE( errQTPFA, consistency_tol );

    double errSIMPLE = computeGravityConsistency_error< 6 >( InnerProductType::SIMPLE, DistortionMode::NonPlanar, eps );
    EXPECT_LE( errSIMPLE, consistency_tol );

    double errBDVLM = computeGravityConsistency_error< 6 >( InnerProductType::BDVLM, DistortionMode::NonPlanar, eps );
    EXPECT_LE( errBDVLM, consistency_tol );
  }
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  geos::basicSetup( argc, argv );

  int const result = RUN_ALL_TESTS();

  geos::basicCleanup();

  return result;
}
