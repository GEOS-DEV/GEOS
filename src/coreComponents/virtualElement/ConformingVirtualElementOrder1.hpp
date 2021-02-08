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
template< localIndex MAXCELLNODES >
class ConformingVirtualElementOrder1 final : public VirtualElementBase
{
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
                        array1d< real64 > & basisIntegrals,
                        real64 * threeDMonomialIntegrals );

  static localIndex m_numQuadraturePoints;
  static localIndex m_numSupportPoints;
  static array1d< real64 > m_quadratureWeights;
  static array1d< real64 > m_basisFunctionsIntegralMean;
  static array2d< real64 > m_stabilizationMatrix;
  static array2d< real64 > m_basisDerivativesIntegralMean;

public:

  GEOSX_HOST_DEVICE
  template< localIndex DIMENSION, typename POINT_COORDS_TYPE >
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

  GEOSX_HOST_DEVICE
  template< localIndex DIMENSION, typename POINT_COORDS_TYPE, typename POINT_SELECTION_TYPE >
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
  virtual localIndex getNumQuadraturePoints() const override
  {
    return m_numQuadraturePoints;
  }

  GEOSX_HOST_DEVICE
  virtual localIndex getNumSupportPoints() const override
  {
    return m_numSupportPoints;
  }

  GEOSX_HOST_DEVICE
  real64 getStabilizationValue( localIndex const iBasisFunction,
                                localIndex const jBasisFunction
                                ) const override
  { return m_stabilizationMatrix[iBasisFunction][jBasisFunction]; }

  /**
   * @brief Get the shape function projected derivatives at a given quadrature point.
   * @details The output contains the integral mean of derivatives.
   * @param q The quadrature point index.
   * @param gradN Return array of the shape function projected derivatives. Size will be @ref
   * getNumSupportPoints() x 3.
   */
  GEOSX_HOST_DEVICE
  void getGradN( localIndex const q,
                 arrayView2d< real64 const > & gradN ) const override
  {
    GEOSX_UNUSED_VAR( q );
    gradN = m_basisDerivativesIntegralMean;
  }

  /**
   * @brief Get the shape function projections at a given quadrature point.
   * @details The output contains the integral mean of functions
   * @param q The quadrature point index.
   * @param gradN Return array of the shape function projections. Size will be @ref getNumSupportPoints().
   */
  GEOSX_HOST_DEVICE
  void getN( localIndex const q,
             arrayView1d< real64 const > & N ) const override
  {
    GEOSX_UNUSED_VAR( q );
    N = m_basisFunctionsIntegralMean;
  }

  GEOSX_HOST_DEVICE
  virtual real64 transformedQuadratureWeight( localIndex const q ) const
  {
    return m_quadratureWeights[q];
  }

};
}
}

#include "ConformingVirtualElementOrder1_impl.hpp"

#endif // GEOSX_VIRTUALELEMENT_CONFORMINGVIRTUALELEMENTORDER1_HPP_
