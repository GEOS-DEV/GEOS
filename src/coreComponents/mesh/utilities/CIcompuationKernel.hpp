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
 * @file KernelBase.hpp
 */

#ifndef GEOSX_MESH_UTILITIES_CICOMPUTATIONKERNEL_HPP_
#define GEOSX_MESH_UTILITIES_CICOMPUTATIONKERNEL_HPP_

#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "finiteElement/FiniteElementDispatch.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

namespace geosx
{

template< typename FE_TYPE >
class CIcomputationKernel
{
 public:  
  CIcomputationKernel( FE_TYPE const & finiteElementSpace,
                       NodeManager const & nodeManager, 
                       CellElementSubRegion const & elementSubRegion,
                       EmbeddedSurfaceSubregion const & embeddedSurfSubRegion ):
    m_finiteElementSpace( finiteElementSpace ),
    m_X( nodeManger.referencePosition() ),
    m_elemsToNodes( elementSubRegion.nodeList().toViewConst() ),
    m_fracturedElems( elementSubRegion.fracturedElementsList().toViewConst()),
    m_cellsToEmbeddedSurfaces( elementSubRegion.embeddedSurfacesList().toViewConst() ),
    m_normalVector( embeddedSurfSubRegion.getNormalVector().toViewConst()),
    m_elemCenter( embeddedSurfSubRegion.getElementCenter().toViewConst())
  {}
   
  static constexpr int numNodesPerElem = FE_TYPE::maxSupportPoints;

  static constexpr int numSamplingPoints = FE_TYPE::numSampligPoints;

  struct StackVariables
  {
    public:

    /// Constructor.
    GEOSX_HOST_DEVICE
    StackVariables():
    xLocal(),
    samplingPointCoord(0.0)
    {}

    /// C-array stack storage for element local the nodal positions.
    real64 xLocal[ numNodesPerElem ][ 3 ];

    real64 samplingPointCoord[3];
  };


  template< typename POLICY >
  void computeCIs()
  {
    GEOSX_MARK_FUNCTION;

    forAll< POLICY >(  m_fracturedElems.size(),
                       [=] GEOSX_HOST_DEVICE ( localIndex const i )
    {
      real64 averageDistance = 0.0;
      
      StackVariables stack;
      setup( k, stack );
      for( integer np=0; q<numOfSamplingPoints; ++np )
      {
        samplingPointCoord( k, np, stack );
        averageDistance += computeDistance( point )
      }
      
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
    real64 pointToFracCenter[3];
    LvArray::tensorOps::copy< 3 >( pointToFracCenter, point );
    LvArray::tensorOps::subtract< 3 >( pointToFracCenter, m_fracCenter[k] );
    return LvArray::tensorOps::AiBi< 3 >( pointToFracCenter, m_normalVector[k] );   
  }

  GEOSX_HOST_DEVICE
  samplingPointCoord( localIndex const k, 
                      integer const np, 
                      StackVariable & stack ) const
  {
    // Get sampling point coord in parent space.

    // Compute shape function values at sampling point
    for (localIndex a=0; a<numNodesPerElem; a++)
    {
      LvArray::tensorOps::scaledAdd(stack.samplingPointCoord, stack.xLocal[a], N[a] );
    }
  }

  private:

    FE_TYPE const & m_finiteElementSpace;

    arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const m_X;

    /// The element to nodes map.
    traits::ViewTypeConst< typename CellElementSubRegion::NodeMapType::base_type > const m_elemsToNodes;

    SortedArrayView< localIndex const > const m_fracturedElems;

    ArrayOfArraysView< localIndex const > const m_cellsToEmbeddedSurfaces;

    arrayView2d< real64 const > const m_normalVector;

    arrayView2d< real64 const > const m_fracCenter;
}

}