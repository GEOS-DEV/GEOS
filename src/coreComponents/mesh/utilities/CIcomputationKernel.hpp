/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */


/**
 * @file CIcomputationKernel.hpp
 */
#ifndef GEOSX_MESH_UTILITIES_CICOMPUTATIONKERNEL_HPP_
#define GEOSX_MESH_UTILITIES_CICOMPUTATIONKERNEL_HPP_

#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "finiteElement/FiniteElementDispatch.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "mesh/EmbeddedSurfaceSubRegion.hpp"
#include "mesh/CellElementSubRegion.hpp"
#include "mesh/ElementType.hpp"


namespace geosx
{

template< typename FE_TYPE >
class CIcomputationKernel
{
public:
  CIcomputationKernel( FE_TYPE const & finiteElementSpace,
                       NodeManager const & nodeManager,
                       CellElementSubRegion const & elementSubRegion,
                       EmbeddedSurfaceSubRegion & embeddedSurfSubRegion ):
    m_finiteElementSpace( finiteElementSpace ),
    m_elementType( elementSubRegion.getElementType() ),
    m_X( nodeManager.referencePosition() ),
    m_elemsToNodes( elementSubRegion.nodeList().toViewConst() ),
    m_fracturedElems( elementSubRegion.fracturedElementsList().toViewConst()),
    m_cellsToEmbeddedSurfaces( elementSubRegion.embeddedSurfacesList().toViewConst() ),
    m_normalVector( embeddedSurfSubRegion.getNormalVector().toViewConst()),
    m_fracCenter( embeddedSurfSubRegion.getElementCenter().toViewConst()),
    m_fractureSurfaceArea( embeddedSurfSubRegion.getElementArea().toViewConst() ),
    m_connectivityIndex( embeddedSurfSubRegion.getConnectivityIndex() )
  {}

  static constexpr int numNodesPerElem = FE_TYPE::maxSupportPoints;

  static constexpr int numSamplingPoints = FE_TYPE::numSamplingPoints;

  struct StackVariables
  {
public:

    /// Constructor.
    GEOSX_HOST_DEVICE
    StackVariables():
      xLocal(),
      samplingPointCoord()
    {}

    /// C-array stack storage for element local the nodal positions.
    real64 xLocal[ numNodesPerElem ][ 3 ];

    real64 samplingPointCoord[3];
  };


  template< typename POLICY,
            typename KERNEL_TYPE >
  static
  void launchCICompuationKernel( KERNEL_TYPE & kernelComponent )
  {
    GEOSX_MARK_FUNCTION;
    forAll< POLICY >( kernelComponent.m_fracturedElems.size(),
                      [=] GEOSX_HOST_DEVICE ( localIndex const i )
    {

      localIndex k = kernelComponent.m_fracturedElems[i];

      StackVariables stack;
      kernelComponent.setup( k, stack );

      real64 averageDistance = 0.0;
      for( integer np=0; np<numSamplingPoints; ++np )
      {
        kernelComponent.samplingPointCoord( np, stack );
        averageDistance += kernelComponent.computeDistance( k, stack.samplingPointCoord );
      }
      averageDistance /= numSamplingPoints;
      kernelComponent.setConnectivityIndex( k, averageDistance );
    } );
  }

  /**
   * @brief Copy global values from primary field to a local stack array.
   * @copydoc
   *
   */
  GEOSX_HOST_DEVICE
  void setup( localIndex const k,
              StackVariables & stack ) const
  {
    for( localIndex a=0; a<numNodesPerElem; ++a )
    {
      localIndex const localNodeIndex = m_elemsToNodes( k, a );

      for( int i=0; i<3; ++i )
      {
        stack.xLocal[ a ][ i ] = m_X[ localNodeIndex ][ i ];
      }
    }
  }

  GEOSX_HOST_DEVICE
  real64 computeDistance( localIndex const k,
                          real64 const (&point)[3] ) const
  {
    localIndex const embSurfIndex = m_cellsToEmbeddedSurfaces[k][0];
    real64 pointToFracCenter[3];
    LvArray::tensorOps::copy< 3 >( pointToFracCenter, point );
    LvArray::tensorOps::subtract< 3 >( pointToFracCenter, m_fracCenter[embSurfIndex] );
    return LvArray::math::abs( LvArray::tensorOps::AiBi< 3 >( pointToFracCenter, m_normalVector[embSurfIndex] ));
  }

  GEOSX_HOST_DEVICE
  void samplingPointCoord( integer const np,
                           StackVariables & stack ) const
  {
    // Get sampling point coord in parent space.
    real64 parentSamplingPointCoord[3] = {0.0, 0.0, 0.0};
    FE_TYPE::getSamplingPointCoordInParentSpace( np, parentSamplingPointCoord );

    // Compute shape function values at sampling point
    real64 N[numNodesPerElem];
    FE_TYPE::calcN( parentSamplingPointCoord, N );

    LvArray::tensorOps::fill< 3 >( stack.samplingPointCoord, 0.0 );

    // Compute sampling point coord in the physical space
    for( localIndex a=0; a<numNodesPerElem; a++ )
    {
      for( int i =0; i < 3; i++ )
      {
        stack.samplingPointCoord[i] += stack.xLocal[a][i] * N[a];
      }
    }
  }

  GEOSX_HOST_DEVICE
  void setConnectivityIndex( localIndex const k,
                             real64 const averageDistance ) const
  {
    localIndex const embSurfIndex = m_cellsToEmbeddedSurfaces[k][0];
    m_connectivityIndex[embSurfIndex] = m_fractureSurfaceArea[embSurfIndex] / averageDistance;
  }

private:

  FE_TYPE const & m_finiteElementSpace;

  ElementType const m_elementType;

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const m_X;

  /// The element to nodes map.
  traits::ViewTypeConst< typename CellElementSubRegion::NodeMapType::base_type > const m_elemsToNodes;

  SortedArrayView< localIndex const > const m_fracturedElems;

  ArrayOfArraysView< localIndex const > const m_cellsToEmbeddedSurfaces;

  arrayView2d< real64 const > const m_normalVector;

  arrayView2d< real64 const > const m_fracCenter;

  arrayView1d< real64 const > const m_fractureSurfaceArea;

  arrayView1d< real64 > const m_connectivityIndex;
};

}

#endif /* GEOSX_MESH_UTILITIES_CICOMPUTATIONKERNEL_HPP_ */
