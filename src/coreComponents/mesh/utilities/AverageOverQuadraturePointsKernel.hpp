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
 * @file AverageOverQuadraturePointsKernel.hpp
 */

#ifndef GEOS_MESH_UTILITIES_AVERAGEOVERQUADRATUREPOINTSKERNEL_HPP_
#define GEOS_MESH_UTILITIES_AVERAGEOVERQUADRATUREPOINTSKERNEL_HPP_

#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "finiteElement/FiniteElementDispatch.hpp"
#include "mesh/CellElementSubRegion.hpp"

namespace geos
{

template< typename SUBREGION_TYPE,
          typename FE_TYPE >
class AverageOverQuadraturePointsBase
{
public:

  AverageOverQuadraturePointsBase( NodeManager & nodeManager,
                                   EdgeManager const & edgeManager,
                                   FaceManager const & faceManager,
                                   SUBREGION_TYPE const & elementSubRegion,
                                   FE_TYPE const & finiteElementSpace ):
    m_finiteElementSpace( finiteElementSpace ),
    m_elemsToNodes( elementSubRegion.nodeList().toViewConst() ),
    m_X( nodeManager.referencePosition() ),
    m_elementVolume( elementSubRegion.getElementVolume() )
  {
    finiteElement::FiniteElementBase::
      initialize< FE_TYPE >( nodeManager,
                             edgeManager,
                             faceManager,
                             elementSubRegion,
                             m_meshData );
  }

  //*****************************************************************************
  struct StackVariables
  {
public:
    GEOS_HOST_DEVICE
    StackVariables():
      xLocal()
    {}

    /// C-array stack storage for element local the nodal positions.
    real64 xLocal[ FE_TYPE::maxSupportPoints ][ 3 ];

    /// Stack variables needed for the underlying FEM type
    typename FE_TYPE::StackVariables feStack;
  };
  //***************************************************************************

  GEOS_HOST_DEVICE
  void setup( localIndex const k,
              StackVariables & stack ) const
  {
    m_finiteElementSpace.template setup< FE_TYPE >( k, m_meshData, stack.feStack );

    for( localIndex a = 0; a < FE_TYPE::maxSupportPoints; ++a )
    {
      localIndex const localNodeIndex = m_elemsToNodes( k, a );

      for( integer i = 0; i < 3; ++i )
      {
        stack.xLocal[a][i] = m_X[localNodeIndex][i];
      }
    }
  }

protected:

  /// The finite element space/discretization object for the element type in
  /// the SUBREGION_TYPE.
  FE_TYPE const & m_finiteElementSpace;

  /// The element to nodes map.
  traits::ViewTypeConst< typename SUBREGION_TYPE::NodeMapType::base_type > const m_elemsToNodes;

  /// The reference position of the nodes
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const m_X;

  /// The volume of the elements
  arrayView1d< real64 const > const m_elementVolume;

  /// Data structure containing mesh data used to setup the finite element
  typename FE_TYPE::template MeshData< SUBREGION_TYPE > m_meshData;
};

template< typename SUBREGION_TYPE,
          typename FE_TYPE >
class AverageOverQuadraturePoints1D :
  public AverageOverQuadraturePointsBase< SUBREGION_TYPE,
                                          FE_TYPE >
{
public:

  /// Alias for the base class;
  using Base = AverageOverQuadraturePointsBase< SUBREGION_TYPE,
                                                FE_TYPE >;

  using Base::m_elementVolume;

  AverageOverQuadraturePoints1D( NodeManager & nodeManager,
                                 EdgeManager const & edgeManager,
                                 FaceManager const & faceManager,
                                 SUBREGION_TYPE const & elementSubRegion,
                                 FE_TYPE const & finiteElementSpace,
                                 arrayView2d< real64 const > const property,
                                 arrayView1d< real64 > const averageProperty ):
    Base( nodeManager,
          edgeManager,
          faceManager,
          elementSubRegion,
          finiteElementSpace ),
    m_property( property ),
    m_averageProperty( averageProperty )
  {}

  struct StackVariables : Base::StackVariables
  {};

  GEOS_HOST_DEVICE
  void setup( localIndex const k,
              StackVariables & stack ) const
  {
    Base::setup( k, stack );
    m_averageProperty[k] = 0.0;
  }

  GEOS_HOST_DEVICE
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack ) const
  {
    real64 const weight = FE_TYPE::transformedQuadratureWeight( q, stack.xLocal, stack.feStack ) / m_elementVolume[k];
    m_averageProperty[k] += weight * m_property[k][q];
  }

  template< typename POLICY,
            typename KERNEL_TYPE >
  static real64
  kernelLaunch( localIndex const numElems,
                KERNEL_TYPE const & kernelComponent )
  {
    forAll< POLICY >( numElems,
                      [=] GEOS_HOST_DEVICE ( localIndex const k )
    {
      typename KERNEL_TYPE::StackVariables stack;

      kernelComponent.setup( k, stack );
      for( integer q = 0; q < FE_TYPE::numQuadraturePoints; ++q )
      {
        kernelComponent.quadraturePointKernel( k, q, stack );
      }
    } );
    return 0.0;
  }

protected:

  /// The property living on quadrature points
  arrayView2d< real64 const > const m_property;

  /// The average property
  arrayView1d< real64 > const m_averageProperty;

};


/**
 * @class AverageOverQuadraturePointsKernelFactory
 */
class AverageOverQuadraturePoints1DKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   */
  template< typename SUBREGION_TYPE,
            typename FE_TYPE,
            typename POLICY >
  static void
  createAndLaunch( NodeManager & nodeManager,
                   EdgeManager const & edgeManager,
                   FaceManager const & faceManager,
                   SUBREGION_TYPE const & elementSubRegion,
                   FE_TYPE const & finiteElementSpace,
                   arrayView2d< real64 const > const property,
                   arrayView1d< real64 > const averageProperty )
  {
    AverageOverQuadraturePoints1D< SUBREGION_TYPE, FE_TYPE >
    kernel( nodeManager, edgeManager, faceManager, elementSubRegion, finiteElementSpace,
            property, averageProperty );

    AverageOverQuadraturePoints1D< SUBREGION_TYPE, FE_TYPE >::template
    kernelLaunch< POLICY >( elementSubRegion.size(), kernel );
  }
};

}

#endif /* GEOS_MESH_UTILITIES_AVERAGEOVERQUADRATUREPOINTSKERNEL_HPP_ */
