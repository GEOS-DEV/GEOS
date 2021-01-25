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
class ConformingVirtualElementOrder1 final : public VirtualElementBase
{
  private:
  void
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
  public:

  localIndex m_numQuadraturePoints;
  localIndex m_numSupportPoints;
  array1d< real64 > m_basisFunctionsIntegralMean;
  array2d< real64 > m_basisDerivativesIntegralMean;
  array2d< real64 > m_stabilizationMatrix;

  template< localIndex DIMENSION, typename POINT_COORDS_TYPE >
  real64 ComputeDiameter( POINT_COORDS_TYPE points,
                          localIndex const & numPoints ) const
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
  real64 ComputeDiameter( POINT_COORDS_TYPE points,
                          POINT_SELECTION_TYPE selectedPoints,
                          localIndex const & numSelectedPoints ) const
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

  void ComputeProjectors( localIndex const & cellIndex,
                          arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodesCoords,
                          CellBlock::NodeMapType const & cellToNodes,
                          CellBlock::FaceMapType const & elementToFaceMap,
                          FaceManager::NodeMapType const & faceToNodeMap,
                          FaceManager::EdgeMapType const & faceToEdgeMap,
                          EdgeManager::NodeMapType const & edgeToNodeMap,
                          arrayView2d< real64 const > const faceCenters,
                          arrayView2d< real64 const > const faceNormals,
                          arrayView1d< real64 const> const faceAreas,
                          arraySlice1d< real64 const > const & cellCenter,
                          real64 const & cellVolume
                          );

  virtual localIndex getNumQuadraturePoints() const override
  {
    return m_numQuadraturePoints;
  }
  virtual localIndex getNumSupportPoints() const override
  {
    return m_numSupportPoints;
  }
};
}
}

#endif // GEOSX_VIRTUALELEMENT_CONFORMINGVIRTUALELEMENTORDER1_HPP_
