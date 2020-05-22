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
 * @file RegionLoop.hpp
 */

#ifndef GEOSX_FINITEELEMENT_REGIONLOOP_HPP_
#define GEOSX_FINITEELEMENT_REGIONLOOP_HPP_

#include "common/DataTypes.hpp"
#include "constitutive/solid/LinearElasticIsotropic.hpp"
#include "constitutive/ConstitutivePassThru.hpp"
#include "finiteElement/FiniteElementDiscretization.hpp"
#include "mesh/ElementRegionManager.hpp"


namespace geosx
{

/**
 * @namespace Contains implementations of physics loops.
 */
namespace finiteElement
{

/**
 * @brief Function to take runtime integers and call generic lambda with
 *        compile time integers.
 * @tparam LAMBDA the type of the generic lambda/function
 * @param[in] NUM_NODES_PER_ELEM The number of nodes per element
 * @param[in] NUM_QUADRATURE_POINTS The number of quadrature points per element.
 * @param[in] lambda The generic lambda that will be passed the integral_constant
 *            representation of
 */
template< typename LAMBDA >
void
discretizationLaunchSelector( localIndex NUM_NODES_PER_ELEM,
                              localIndex NUM_QUADRATURE_POINTS,
                              LAMBDA && lambda )
{
  if( NUM_NODES_PER_ELEM==8 && NUM_QUADRATURE_POINTS==8 )
  {
    lambda( std::integral_constant< int, 8 >(), std::integral_constant< int, 8 >() );
  }
  else if( NUM_NODES_PER_ELEM==8 && NUM_QUADRATURE_POINTS==1 )
  {
    lambda( std::integral_constant< int, 8 >(), std::integral_constant< int, 1 >() );
  }
  else if( NUM_NODES_PER_ELEM==4 && NUM_QUADRATURE_POINTS==1 )
  {
    lambda( std::integral_constant< int, 4 >(), std::integral_constant< int, 1 >() );
  }
  else
  {
    GEOSX_ERROR( "Valid Branch not found." );
  }
}


/**
 * @class FiniteElementRegionLoop
 *
 * This class encapsulates base components and interface for applying a finite
 * element method over a loop of element regions.
 *
 */
class KernelBase
{
public:

  //***************************************************************************
  /**
   * @class Parameters
   *
   * Contains non-mesh parameters passed in from the calling scope.
   */
  struct Parameters
  {};


  //***************************************************************************
  /**
   * @class StackVariables
   * @tparam NUM_ROWS The number rows to allocate for the residual/jacobian.
   * @tparam NUM_COLS The number or columns to allocate for the jacobian.
   * Contains variables that will be allocated on the stack of the main kernel.
   * This will typically consist of local arrays to hold data mapped from the
   * global data arrays, and local storage for the residual and jacobian
   * contributions.
   */
  template< int NUM_ROWS,
            int NUM_COLS >
  struct StackVariables
  {};

  //***************************************************************************
  /**
   * @class Kernels
   */
  template< typename SUBREGION_TYPE,
            typename CONSTITUTIVE_TYPE,
            int NUM_TEST_SUPPORT_POINTS_PER_ELEM,
            int NUM_TRIAL_SUPPORT_POINTS_PER_ELEM,
            int NUM_DOF_PER_TEST_SP,
            int NUM_DOF_PER_TRIAL_SP >
  class Components
  {
public:
    Components( SUBREGION_TYPE const & elementSubRegion,
                FiniteElementBase const * const finiteElementSpace,
                CONSTITUTIVE_TYPE * const inputConstitutiveType,
                Parameters const & GEOSX_UNUSED_PARAM( parameters ) ):
      Components( elementSubRegion,
                  finiteElementSpace,
                  inputConstitutiveType,
                  typename CONSTITUTIVE_TYPE::KernelWrapper() )
    {}


    Components( SUBREGION_TYPE const & elementSubRegion,
                FiniteElementBase const * const finiteElementSpace,
                CONSTITUTIVE_TYPE * const GEOSX_UNUSED_PARAM( inputConstitutiveType ),
                typename CONSTITUTIVE_TYPE::KernelWrapper const & inputConstitutiveUpdate ):
      elemsToNodes( elementSubRegion.nodeList().toViewConst() ),
      elemGhostRank( elementSubRegion.ghostRank() ),
      constitutiveUpdate( inputConstitutiveUpdate ),
      m_finiteElementSpace( finiteElementSpace )
    {}

    template< typename STACK_VARIABLE_TYPE >
    GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    void setup( localIndex const GEOSX_UNUSED_PARAM( k ),
                STACK_VARIABLE_TYPE & GEOSX_UNUSED_PARAM( stack ) ) const
    {}

    template< typename STACK_VARIABLE_TYPE >
    GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    void quadraturePointStateUpdate( localIndex const GEOSX_UNUSED_PARAM( k ),
                                     localIndex const GEOSX_UNUSED_PARAM( q ),
                                     STACK_VARIABLE_TYPE & GEOSX_UNUSED_PARAM( stack ) ) const
    {}

    template< typename STACK_VARIABLE_TYPE >
    GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    void quadraturePointJacobianContribution( localIndex const GEOSX_UNUSED_PARAM( k ),
                                              localIndex const GEOSX_UNUSED_PARAM( q ),
                                              STACK_VARIABLE_TYPE & GEOSX_UNUSED_PARAM( stack ) ) const
    {}

    template< typename STACK_VARIABLE_TYPE >
    GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    void quadraturePointResidualContribution( localIndex const GEOSX_UNUSED_PARAM( k ),
                                              localIndex const GEOSX_UNUSED_PARAM( q ),
                                              STACK_VARIABLE_TYPE & GEOSX_UNUSED_PARAM( stack ) ) const
    {}

    template< typename STACK_VARIABLE_TYPE >
    GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    real64 complete( localIndex const GEOSX_UNUSED_PARAM( k ),
                     STACK_VARIABLE_TYPE & GEOSX_UNUSED_PARAM( stack ) ) const
    {
      return 0;
    }



    template< typename POLICY,
              int NUM_QUADRATURE_POINTS,
              typename STACK_VARIABLES,
              typename COMPONENT_TYPE >
    static
    typename std::enable_if< std::is_same< POLICY, serialPolicy >::value ||
                             std::is_same< POLICY, parallelHostPolicy >::value, real64 >::type
    Launch( localIndex const numElems,
            COMPONENT_TYPE const & kernelComponent )
    {
      GEOSX_MARK_FUNCTION;
      RAJA::ReduceMax< typename ReducePolicy< POLICY >::type, real64 > maxResidual( 0 );

      forAll< POLICY >( numElems,
                        [=] ( localIndex const k )
      {
        STACK_VARIABLES stack;

        kernelComponent.setup( k, stack );
        for( integer q=0; q<NUM_QUADRATURE_POINTS; ++q )
        {
          kernelComponent.quadraturePointStateUpdate( k, q, stack );

          kernelComponent.quadraturePointJacobianContribution( k, q, stack );

          kernelComponent.quadraturePointResidualContribution( k, q, stack );
        }
        if( kernelComponent.elemGhostRank[k] < 0 )
        {
          maxResidual.max( kernelComponent.complete( k, stack ) );
        }
      } );
      return maxResidual.get();
    }

    template< typename POLICY,
              int NUM_QUADRATURE_POINTS,
              typename STACK_VARIABLES,
              typename COMPONENT_TYPE >
    static
    typename std::enable_if< !( std::is_same< POLICY, serialPolicy >::value ||
                                std::is_same< POLICY, parallelHostPolicy >::value ), real64 >::type
    Launch( localIndex const numElems,
            COMPONENT_TYPE const & kernelComponent )
    {
      GEOSX_MARK_FUNCTION;
      RAJA::ReduceMax< typename ReducePolicy< POLICY >::type, real64 > maxResidual( 0 );

      forAll< POLICY >( numElems,
                        [=] GEOSX_DEVICE ( localIndex const k )
      {
        STACK_VARIABLES stack;

        kernelComponent.setup( k, stack );
        for( integer q=0; q<NUM_QUADRATURE_POINTS; ++q )
        {
          kernelComponent.quadraturePointStateUpdate( k, q, stack );

          kernelComponent.quadraturePointJacobianContribution( k, q, stack );

          kernelComponent.quadraturePointResidualContribution( k, q, stack );
        }
        if( kernelComponent.elemGhostRank[k] < 0 )
        {
          maxResidual.max( kernelComponent.complete( k, stack ) );
        }
      } );
      return maxResidual.get();
    }



//protected:
    typename SUBREGION_TYPE::NodeMapType::base_type::ViewTypeConst const elemsToNodes;
    arrayView1d< integer const > const elemGhostRank;
    typename CONSTITUTIVE_TYPE::KernelWrapper const constitutiveUpdate;
    FiniteElementBase const * m_finiteElementSpace;

  };
};


//***************************************************************************
template< typename POLICY,
          typename UPDATE_CLASS,
          typename CONSTITUTIVE_BASE,
          typename REGION_TYPE,
          template< typename SUBREGION_TYPE,
                    typename CONSTITUTIVE_TYPE,
                    int NUM_TEST_SUPPORT_POINTS_PER_ELEM,
                    int NUM_TRIAL_SUPPORT_POINTS_PER_ELEM > class COMPONENTS_TYPE>
static
real64 RegionBasedKernelApplication( MeshLevel & mesh,
                                     arrayView1d< string const > const & targetRegions,
                                     arrayView1d< string const > const & constitutiveNames,
                                     FiniteElementDiscretization const * const feDiscretization,
                                     arrayView1d< globalIndex const > const & inputDofNumber,
                                     ParallelMatrix & inputMatrix,
                                     ParallelVector & inputRhs,
                                     typename UPDATE_CLASS::Parameters const & parameters )
{

  real64 maxResidual = 0;

  NodeManager & nodeManager = *(mesh.getNodeManager());
  ElementRegionManager & elementRegionManager = *(mesh.getElemManager());


  elementRegionManager.forElementSubRegions< REGION_TYPE >( targetRegions,
                                                            [&] ( localIndex const targetRegionIndex,
                                                                  auto & elementSubRegion )
  {
    localIndex const numElems = elementSubRegion.size();
    typedef TYPEOFREF( elementSubRegion ) SUBREGIONTYPE;


    FiniteElementBase const * finiteElementSpace = nullptr;
    if( feDiscretization != nullptr )
    {
      finiteElementSpace = ( feDiscretization->getFiniteElement( elementSubRegion.GetElementTypeString() ) ).get();
    }

    localIndex const
    numQuadraturePointsPerElem = finiteElementSpace == nullptr ?
                                 1 :
                                 finiteElementSpace->n_quadrature_points();

    discretizationLaunchSelector( elementSubRegion.numNodesPerElement(),
                                  numQuadraturePointsPerElem,
                                  [&]( auto constNNPE,
                                       auto constNQPPE )
    {
      constexpr int NUM_NODES_PER_ELEM = decltype( constNNPE )::value;
      constexpr int NUM_QUADRATURE_POINTS = decltype( constNQPPE )::value;

      constitutive::ConstitutiveBase * const
      constitutiveRelation = ( targetRegionIndex <= constitutiveNames.size()-1 ) ?
                             elementSubRegion.template getConstitutiveModel( constitutiveNames[targetRegionIndex] ) : nullptr;

      constitutive::ConstitutivePassThru< CONSTITUTIVE_BASE >::Execute( constitutiveRelation,
                                                                        [&]( auto * const castedConstitutiveRelation )
      {
        using CONSTITUTIVE_TYPE = TYPEOFPTR( castedConstitutiveRelation );

        using KERNEL_TYPE = COMPONENTS_TYPE< SUBREGIONTYPE,
                                             CONSTITUTIVE_TYPE,
                                             NUM_NODES_PER_ELEM,
                                             NUM_NODES_PER_ELEM >;
        KERNEL_TYPE kernelComponent( inputDofNumber,
                                     inputMatrix,
                                     inputRhs,
                                     nodeManager,
                                     elementSubRegion,
                                     finiteElementSpace,
                                     castedConstitutiveRelation,
                                     parameters );

        maxResidual = std::max( maxResidual,
                                KERNEL_TYPE::template Launch< POLICY,
                                                              NUM_QUADRATURE_POINTS,
                                                              typename decltype(kernelComponent)::StackVars >( numElems,
                                                                                                               kernelComponent ) );
      } );
    } );
  } );

  return maxResidual;
}
}
}



#endif /* GEOSX_FINITEELEMENT_REGIONLOOP_HPP_ */
