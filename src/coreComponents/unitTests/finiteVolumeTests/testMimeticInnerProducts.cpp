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

template< localIndex NF, typename ARRAY_VIEW_T >
static void runConsistencyTest( array2d< real64, nodes::REFERENCE_POSITION_PERM > const & nodePosition,
                                FaceManager::NodeMapType const & faceToNodes,
                                array1d< localIndex > const & elemToFaces,
                                real64 const elemCenter[3],
                                real64 const elemPerm[3],
                                real64 elemVolume,
                                ARRAY_VIEW_T const & transMatrix,
                                std::string const & testName )
{
    real64 N[NF][3], C[NF][3], TC[NF][3], K[3][3];
    real64 faceCenter[3], faceNormal[3];

    // full tensor K
    real64 perm[3] = { elemPerm[0], elemPerm[1], elemPerm[2] };
    MimeticInnerProductHelpers::makeFullTensor( perm, K );

    // compute face normals
    for ( localIndex iface = 0; iface < NF; ++iface )
    {
        computationalGeometry::centroid_3DPolygon(
            faceToNodes[iface], nodePosition.toViewConst(),
            faceCenter, faceNormal );

        for ( int d = 0; d < 3; ++d )
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
    
    for ( localIndex i = 0; i < NF; ++i )
    {
      for ( localIndex j = 0; j < 3; ++j )
      {
        diffMat[i][j] = C[i][j] - TC[i][j];
      }
    }
    
    for( std::ptrdiff_t i = 0; i < NF; ++i )
    {
      diffNorm += LvArray::tensorOps::l2NormSquared<3>( diffMat[i] );
    }
    diffNorm = std::sqrt( diffNorm );
    
    EXPECT_LT( diffNorm, 1e-10 ) << testName << ": norm(NK - TC) = " << diffNorm;
    
    if ( diffNorm < 1e-10 )
    {
        std::cout << "[CONSISTENCY TEST PASSED] " << testName << " consistency test passed" <<  std::endl;
    }
    else
    {
        std::cout << "[CONSISTENCY TEST FAILED] " << testName << " consistency test failed: norm(NK - TC) = " << diffNorm << std::endl;
    }
}

static inline
void distortTopFaceNonPlanar( array2d< real64, nodes::REFERENCE_POSITION_PERM > & nodePosition,
                              real64 eps )
{
  nodePosition( 1, 2 ) += eps;   // vertex 1:  z -> z + eps
  nodePosition( 7, 2 ) += eps;   // vertex 7:  z -> z + eps

  // nodePosition( 5, 2 ) -= eps;
  // nodePosition( 3, 2 ) -= eps;
}

static inline
void computeDistortedVolumeAndCenter( array2d < real64, nodes::REFERENCE_POSITION_PERM > const & nodePosition,
                                      real64 ( & elemCenter )[3],
                                      real64 & elemVolume )
{
    array1d< localIndex > toNodes;
    toNodes.resize( 8 );
    toNodes( 0 ) = 0;
    toNodes( 1 ) = 4;
    toNodes( 2 ) = 2;
    toNodes( 3 ) = 6;
    toNodes( 4 ) = 1;
    toNodes( 5 ) = 5;
    toNodes( 6 ) = 3;
    toNodes( 7 ) = 7;
    
    computeVolumeAndCenter( nodePosition, toNodes, elemCenter, elemVolume);
}

// check if matrix A is SPD via Cholesky factorization.
// also tracks min/max pivot values (minP: smallest pivot, maxP = largest pivot)
// their ratio is approximately sqrt( condition number )
static inline
bool cholCheck( arraySlice2d< real64 const > const & A,
                real64 & minP, real64 & maxP )
{
  localIndex const n = A.size(0);
  stackArray2d< real64, 16 * 16 > R( n, n );
    
  for( localIndex i = 0; i < n; ++i )
    for( localIndex j = 0; j < n; ++j )
      R(i,j) = A(i,j);

  minP = std::numeric_limits< real64 >::infinity();
  maxP = 0.0;

  for( localIndex k = 0; k < n; ++k )
  {
    // symmetric
    for( localIndex j = k; j < n; ++j )
      R(j,k) = R(k,j) = 0.5 * ( R(k,j) + R(j,k) );

    real64 s = R(k,k);
    for( localIndex m = 0; m < k; ++m )
      s -= R(m,k) * R(m,k);

    if( s <= 0.0 || !std::isfinite(s) )
      return false;

    real64 const rkk = std::sqrt(s);
    R(k,k) = rkk;
    minP = std::min( minP, rkk );
    maxP = std::max( maxP, rkk );

    for( localIndex j = k + 1; j < n; ++j )
    {
      real64 t = R(k,j);
      for( localIndex m = 0; m < k; ++m )
          t -= R(m,k) * R(m,j);
        
      R(k,j) = t / rkk;
      R(j,k) = R(k,j);
    }
  }
  return true;
}

// returns standard GEOSX node ordering for hexahedral cell
static inline
void getHexaNodeOrder( array1d< localIndex > & toNodes )
{
  toNodes.resize(8);
  toNodes(0) = 0; toNodes(1) = 4; toNodes(2) = 2; toNodes(3) = 6;
  toNodes(4) = 1; toNodes(5) = 5; toNodes(6) = 3; toNodes(7) = 7;
}

// B = A^T A
static inline void mat3_ATxA( const real64 A[3][3], real64 B[3][3] )
{
  for( int i = 0; i < 3; ++i )
    for( int j = 0; j < 3; ++j )
    {
      real64 s = 0.0;
      for( int k = 0; k < 3; ++k )
          s += A[k][i] * A[k][j];
      B[i][j] = s;
    }
}

static inline real64 mat3_det( const real64 A[3][3] )
{
  return A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1])
       - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0])
       + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
}

static inline bool mat3_inv( const real64 A[3][3], real64 invA[3][3] )
{
  real64 d = mat3_det(A);
  if( std::abs(d) < 1e-30 ) return false;
  real64 id = 1.0/d;
  invA[0][0] =  (A[1][1] * A[2][2] - A[1][2] * A[2][1]) * id;
  invA[0][1] = -(A[0][1] * A[2][2] - A[0][2] * A[2][1]) * id;
  invA[0][2] =  (A[0][1] * A[1][2] - A[0][2] * A[1][1]) * id;
  invA[1][0] = -(A[1][0] * A[2][2] - A[1][2] * A[2][0]) * id;
  invA[1][1] =  (A[0][0] * A[2][2] - A[0][2] * A[2][0]) * id;
  invA[1][2] = -(A[0][0] * A[1][2] - A[0][2] * A[1][0]) * id;
  invA[2][0] =  (A[1][0] * A[2][1] - A[1][1] * A[2][0]) * id;
  invA[2][1] = -(A[0][0] * A[2][1] - A[0][1] * A[2][0]) * id;
  invA[2][2] =  (A[0][0] * A[1][1] - A[0][1] * A[1][0]) * id;
    
  return true;
}

// compute eigenvalues of 3 x 3 symmetric matrix by Jacobi eigenvalue algorithm
static inline void sym3x3_eigs( real64 S[3][3], real64 evals[3] )
{
  for( int sweep = 0; sweep < 20; ++sweep )
  {
    // find the largest off-diagonal element
    int p = 0, q = 1;
    real64 maxa = std::abs( S[0][1] );
    if( std::abs(S[0][2]) > maxa ){ p = 0; q = 2; maxa = std::abs(S[0][2]); }
    if( std::abs(S[1][2]) > maxa ){ p = 1; q = 2; maxa = std::abs(S[1][2]); }
    if( maxa < 1e-30 ) break;
    real64 app = S[p][p], aqq = S[q][q], apq = S[p][q];
    real64 phi = 0.5 * std::atan2( 2 * apq, (aqq - app) );
    real64 c = std::cos(phi), s = std::sin(phi);

    // apply the Jocobi rotation
    for( int k = 0; k < 3; ++k )
    {
      real64 skp = S[k][p], skq = S[k][q];
      S[k][p] = c * skp - s * skq;
      S[k][q] = s * skp + c * skq;
    }
    for( int k = 0; k < 3; ++k )
    {
      real64 spk = S[p][k], sqk = S[q][k];
      S[p][k] = c * spk - s * sqk;
      S[q][k] = s * spk + c * sqk;
    }
  }
  evals[0] = S[0][0];
  evals[1] = S[1][1];
  evals[2] = S[2][2];
}

struct AffineMetrics {
  real64 A[3][3];   // best-fit linear
  real64 b[3];      // translation
  real64 detA;      // det(A)
  real64 rms;       // RMS residual
};

static inline
AffineMetrics fitAffineAndMeasure(
  array2d< real64, nodes::REFERENCE_POSITION_PERM > const & X0, // base
  array2d< real64, nodes::REFERENCE_POSITION_PERM > const & X   // distorted
)
{
  array1d< localIndex > toNodes;
  getHexaNodeOrder( toNodes );

  // compute centroid
  real64 xbar[3] = { 0, 0, 0 }, Xbar[3] = { 0, 0, 0 };
  for( localIndex i = 0; i < toNodes.size(); ++i){
    localIndex a = toNodes(i);
    for( int d = 0; d < 3; ++d ){
      Xbar[d] += X0( a, d );
      xbar[d] += X( a, d );
    }
  }
  for( int d = 0; d < 3; ++d ){
    Xbar[d] /= toNodes.size();
    xbar[d] /= toNodes.size();
  }

  real64 SXX[3][3]={ {0,0,0}, {0,0,0}, {0,0,0} };
  real64 SXx[3][3]={ {0,0,0}, {0,0,0}, {0,0,0} };
  for( localIndex i = 0; i < toNodes.size(); ++i ){
    localIndex a = toNodes(i);
    real64 dX[3] = { X0(a,0) - Xbar[0], X0(a,1) - Xbar[1], X0(a,2) - Xbar[2] };
    real64 dx[3] = {  X(a,0) - xbar[0],  X(a,1) - xbar[1],  X(a,2) - xbar[2] };
    for( int r = 0; r < 3; ++r ){
      for( int c = 0; c < 3; ++c ){
        SXX[r][c] += dX[r] * dX[c];
        SXx[r][c] += dx[r] * dX[c];
      }
    }
  }

  real64 invSXX[3][3];
  AffineMetrics M;
  if( !mat3_inv(SXX, invSXX) ){
    for( int i = 0; i < 3; ++i ){ for( int j = 0; j < 3; ++j ) M.A[i][j] = ( i==j) ; M.b[i]=0; }
    M.detA = 1;
    M.rms = 0;
    return M;
  }
  
  for( int i = 0; i < 3; ++i )
    for( int j = 0; j < 3; ++j ){
      real64 s = 0.0;
      for( int k = 0; k < 3; ++k) s += SXx[i][k] * invSXX[k][j];
      M.A[i][j] = s;
    }
    
  for( int i = 0; i < 3; ++i ){
    real64 s = xbar[i] - ( M.A[i][0] * Xbar[0] + M.A[i][1] * Xbar[1] + M.A[i][2] * Xbar[2] );
    M.b[i] = s;
  }

  // sigma from eig(A^T A)
  real64 ATA[3][3]; mat3_ATxA(M.A, ATA);
 
  real64 S[3][3] = { {ATA[0][0], ATA[0][1], ATA[0][2]},
                   {ATA[1][0], ATA[1][1], ATA[1][2]},
                   {ATA[2][0], ATA[2][1], ATA[2][2]} };
  real64 evals[3]; sym3x3_eigs(S, evals);
 
  for( int i = 0; i < 3; ++i ) evals[i] = std::max<real64>( evals[i], 0.0 );
  real64 s1 = std::sqrt( std::max(evals[0], std::max(evals[1], evals[2])) );
  real64 s3 = std::sqrt( std::min(evals[0], std::min(evals[1], evals[2])) );

  M.detA = mat3_det(M.A);

  // RMS residual
  real64 se = 0.0;
  for( localIndex i = 0; i < toNodes.size(); ++i ){
    localIndex a = toNodes(i);
    real64 px[3] = {
      M.A[0][0] * X0(a,0) + M.A[0][1] * X0(a,1) + M.A[0][2] * X0(a,2) + M.b[0],
      M.A[1][0] * X0(a,0) + M.A[1][1] * X0(a,1) + M.A[1][2] * X0(a,2) + M.b[1],
      M.A[2][0] * X0(a,0) + M.A[2][1] * X0(a,1) + M.A[2][2] * X0(a,2) + M.b[2]
    };
    real64 rx = X(a,0) - px[0], ry = X(a,1) - px[1], rz = X(a,2) - px[2];
    se += rx * rx + ry * ry + rz * rz;
  }
  M.rms = std::sqrt( se / toNodes.size() );
  return M;
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
                            elemToFaces,
                            elemCenter,
                            elemPerm,
                            elemVolume,
                            transMatrix.toViewConst(),
                            "TPFA_hexa");
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
                            elemToFaces,
                            elemCenter,
                            elemPerm,
                            elemVolume,
                            transMatrix.toViewConst(),
                            "QTPFA_hexa");
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
                            elemToFaces,
                            elemCenter,
                            elemPerm,
                            elemVolume,
                            transMatrix.toViewConst(),
                            "Simple_hexa");
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
                            elemToFaces,
                            elemCenter,
                            elemPerm,
                            elemVolume,
                            transMatrix.toViewConst(),
                            "BdVLM_hexa");
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
                            elemToFaces,
                            elemCenter,
                            elemPerm,
                            elemVolume,
                            transMatrix.toViewConst(),
                            "TPFA_tetra");
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
                            elemToFaces,
                            elemCenter,
                            elemPerm,
                            elemVolume,
                            transMatrix.toViewConst(),
                            "QTPFA_tetra");
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
                            elemToFaces,
                            elemCenter,
                            elemPerm,
                            elemVolume,
                            transMatrix.toViewConst(),
                            "Simple_tetra");
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
                            elemToFaces,
                            elemCenter,
                            elemPerm,
                            elemVolume,
                            transMatrix.toViewConst(),
                            "BdVLM_tetra");
}

TEST( testMimeticInnerProducts, Hexa_DistortionTest )
{
  localIndex constexpr NF = 6;

  array2d< real64, nodes::REFERENCE_POSITION_PERM > nodePos_base;
  FaceManager::NodeMapType faceToNodes;
  array1d< localIndex > elemToFaces;
  real64 elemCenter_base[3] = { 0.0, 0.0, 0.0 };
  real64 elemPerm[3] = { 1.0, 1.0, 1.0 };
  real64 elemVolume_base = 0.0;
  real64 lengthTol = 0.0;

  // base hexa (undistorted)
  stackArray2d< real64, NF * NF > dummyRef( NF, NF );
  makeHexa(nodePos_base, faceToNodes, elemToFaces,
           elemCenter_base, elemVolume_base, elemPerm,
           lengthTol, InnerProductType::TPFA, dummyRef);

  // characteristic length for nondimensional RMS
  real64 Lchar = std::cbrt( std::max(elemVolume_base, 1e-30) );

  // eps list
  std::vector< real64 > epsilons = { 1e-6, 1e-3, 1e-1, 0.5, 1.0, 1.5, 2.0 };

//  std::string csvPath = "cond_number.csv";
//  std::ofstream fout(csvPath);
//  fout << "eps,volRatio,rms_affine,TPFA,QTPFA,Simple,BdVLM\n";

  for( real64 eps : epsilons )
  {
    std::cout << "\n=== eps = " << eps << " ===\n";

    auto nodePos = nodePos_base;

    // distort
    distortTopFaceNonPlanar(nodePos, eps);

    // update center & volume
    real64 elemCenter[3] = { 0, 0, 0 };
    real64 elemVolume = 0.0;
    computeDistortedVolumeAndCenter(nodePos, elemCenter, elemVolume);
    EXPECT_GT(elemVolume, 0.0);

    AffineMetrics metrics = fitAffineAndMeasure(nodePos_base, nodePos);
    real64 volRatio  = elemVolume / std::max(elemVolume_base, 1e-30);
    real64 rms_affine = metrics.rms / std::max(Lchar, 1e-30);

    std::cout << "volRatio = " << volRatio
              << "  rms_affine = " << rms_affine << "\n";

    array1d< real64 > transMult(NF);
    transMult.setValues<parallelHostPolicy>(1.0);

    stackArray1d< real64,3 > center(3);
    for( int i = 0; i < 3; ++i ) center[i] = elemCenter[i];

    stackArray2d< real64, NF * NF > T_tpfa(NF,NF), T_q(NF,NF), T_s(NF,NF), T_b(NF,NF);

    TPFAInnerProduct::compute<NF>( nodePos.toViewConst(), transMult.toViewConst(),
      faceToNodes.toViewConst(), elemToFaces.toSliceConst(),
      center, elemVolume, elemPerm, lengthTol, T_tpfa.toSlice() );

    QuasiTPFAInnerProduct::compute<NF>( nodePos.toViewConst(), transMult.toViewConst(),
      faceToNodes.toViewConst(), elemToFaces.toSliceConst(),
      center, elemVolume, elemPerm, lengthTol, T_q.toSlice() );

    SimpleInnerProduct::compute<NF>( nodePos.toViewConst(), transMult.toViewConst(),
      faceToNodes.toViewConst(), elemToFaces.toSliceConst(),
      center, elemVolume, elemPerm, lengthTol, T_s.toSlice() );

    BdVLMInnerProduct::compute<NF>( nodePos.toViewConst(), transMult.toViewConst(),
      faceToNodes.toViewConst(), elemToFaces.toSliceConst(),
      center, elemVolume, elemPerm, lengthTol, T_b.toSlice() );

    double cond_tpfa = 0, cond_qtpfa = 0, cond_simple = 0, cond_bdvlm = 0;
    {
      real64 minP, maxP;
      bool spd;

      spd = cholCheck(T_tpfa.toSliceConst(), minP, maxP);
      cond_tpfa = maxP / minP;
      std::cout << "TPFA: SPD="<<spd<<" cond~="<<cond_tpfa<<"\n";

      spd = cholCheck(T_q.toSliceConst(), minP, maxP);
      cond_qtpfa = maxP/minP;
      std::cout << "QTPFA: SPD="<<spd<<" cond~="<<cond_qtpfa<<"\n";

      spd = cholCheck(T_s.toSliceConst(), minP, maxP);
      cond_simple = maxP/minP;
      std::cout << "Simple: SPD="<<spd<<" cond~="<<cond_simple<<"\n";

      spd = cholCheck(T_b.toSliceConst(), minP, maxP);
      cond_bdvlm = maxP/minP;
      std::cout << "BdVLM: SPD="<<spd<<" cond~="<<cond_bdvlm<<"\n";
    }

//    fout << eps << "," << volRatio << ","  << rms_affine << ","
//         << cond_tpfa << "," << cond_qtpfa << ","
//         << cond_simple << "," << cond_bdvlm << "\n";
  }

//  fout.close();
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  geos::basicSetup( argc, argv );

  int const result = RUN_ALL_TESTS();

  geos::basicCleanup();

  return result;
}
