/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file ConformingVirtualElement_1.hpp
 */

#ifndef GEOSX_VIRTUALELEMENT_CONFORMINGVIRTUALELEMENTORDER1_HPP_
#define GEOSX_VIRTUALELEMENT_CONFORMINGVIRTUALELEMENTORDER1_HPP_

#include "VirtualElementBase.hpp"

namespace geosx
{
namespace virtualElement
{
template< localIndex MAXCELLNODES, localIndex MAXFACENODES >
class ConformingVirtualElementOrder1 final : public VirtualElementBase
{
public:
  static constexpr localIndex maxSupportPoints = MAXCELLNODES;
  static constexpr localIndex numQuadraturePoints = 1;

private:
  GEOSX_HOST_DEVICE
  static void
  ComputeFaceIntegrals( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodesCoords,
                        arraySlice1d< localIndex const > const & faceToNodes,
                        arraySlice1d< localIndex const > const & faceToEdges,
                        real64 const & faceArea,
                        arraySlice1d< real64 const > const & faceCenter,
                        arraySlice1d< real64 const > const & faceNormal,
                        EdgeManager::NodeMapType const & edgeToNodes,
                        real64 const & invCellDiameter,
                        arraySlice1d< real64 const > const & cellCenter,
                        real64 basisIntegrals[MAXFACENODES],
                        real64 threeDMonomialIntegrals[3] );

  static localIndex m_numSupportPoints;
  static real64 m_quadratureWeight;
  static real64 m_basisFunctionsIntegralMean[maxSupportPoints];
  static real64 m_stabilizationMatrix[maxSupportPoints][maxSupportPoints];
  static real64 m_basisDerivativesIntegralMean[maxSupportPoints][3];

public:
  virtual ~ConformingVirtualElementOrder1() override
  {}

  GEOSX_HOST_DEVICE
  static void
  ComputeProjectors( localIndex const & cellIndex,
                     arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodesCoords,
                     CellBlock::NodeMapType const & cellToNodeMap,
                     CellBlock::FaceMapType const & elementToFaceMap,
                     FaceManager::NodeMapType const & faceToNodeMap,
                     FaceManager::EdgeMapType const & faceToEdgeMap,
                     EdgeManager::NodeMapType const & edgeToNodeMap,
                     arrayView2d< real64 const > const faceCenters,
                     arrayView2d< real64 const > const faceNormals,
                     arrayView1d< real64 const > const faceAreas,
                     arraySlice1d< real64 const > const & cellCenter,
                     real64 const & cellVolume
                     );

  GEOSX_HOST_DEVICE
  static localIndex getNumQuadraturePoints()
  {
    return numQuadraturePoints;
  }

  GEOSX_HOST_DEVICE
  static localIndex getNumSupportPoints()
  {
    return m_numSupportPoints;
  }

  GEOSX_HOST_DEVICE
  static real64 calcStabilizationValue( localIndex const iBasisFunction,
                                        localIndex const jBasisFunction
                                        )
  { return m_stabilizationMatrix[iBasisFunction][jBasisFunction]; }

  GEOSX_HOST_DEVICE
  static real64 calcGradN( localIndex const q,
                           real64 ( & gradN )[maxSupportPoints][3] )
  {
    for( localIndex i = 0; i < maxSupportPoints; ++i )
      for( localIndex j = 0; j < 3; ++j )
        gradN[i][j] = m_basisDerivativesIntegralMean[i][j];
    return transformedQuadratureWeight( q );
  }

  GEOSX_HOST_DEVICE
  static void calcN( localIndex const q,
                     real64 ( & N )[maxSupportPoints] )
  {
    GEOSX_UNUSED_VAR( q );
    for( localIndex i = 0; i < maxSupportPoints; ++i )
      N[i] = m_basisFunctionsIntegralMean[i];
  }

  GEOSX_HOST_DEVICE
  static real64 transformedQuadratureWeight( localIndex const GEOSX_UNUSED_PARAM( q ) )
  {
    return m_quadratureWeight;
  }


public:

  template< localIndex DIMENSION, typename POINT_COORDS_TYPE >
  GEOSX_HOST_DEVICE
  static real64 ComputeDiameter( POINT_COORDS_TYPE points,
                                 localIndex const & numPoints )
  {
    array1d< localIndex > selectAllPoints( numPoints );
    for( localIndex i = 0; i < numPoints; ++i )
      selectAllPoints[i] = i;
    return ComputeDiameter< DIMENSION, POINT_COORDS_TYPE,
                            array1d< localIndex > const & >( points,
                                                             selectAllPoints,
                                                             numPoints );
  }

  template< localIndex DIMENSION, typename POINT_COORDS_TYPE, typename POINT_SELECTION_TYPE >
  GEOSX_HOST_DEVICE
  static real64 ComputeDiameter( POINT_COORDS_TYPE points,
                                 POINT_SELECTION_TYPE selectedPoints,
                                 localIndex const & numSelectedPoints )
  {
    real64 diameter = 0;
    for( localIndex numPoint = 0; numPoint < numSelectedPoints; ++numPoint )
    {
      for( localIndex numOthPoint = 0; numOthPoint < numPoint; ++numOthPoint )
      {
        real64 candidateDiameter = 0.0;
        for( localIndex i = 0; i < DIMENSION; ++i )
        {
          real64 coordDiff = points[selectedPoints[numPoint]][i] -
                             points[selectedPoints[numOthPoint]][i];
          candidateDiameter += coordDiff * coordDiff;
        }
        if( diameter < candidateDiameter )
          diameter = candidateDiameter;
      }
      diameter = LvArray::math::sqrt< real64 >( diameter );
    }
    return diameter;
  }

};
}
}

#include "ConformingVirtualElementOrder1_impl.hpp"

#endif // GEOSX_VIRTUALELEMENT_CONFORMINGVIRTUALELEMENTORDER1_HPP_
