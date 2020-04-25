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
 * @file PhysicsLoopInterface.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_PHYSICSLOOPINTERFACE_HPP_
#define GEOSX_PHYSICSSOLVERS_PHYSICSLOOPINTERFACE_HPP_

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
namespace physicsLoopInterface
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
class FiniteElementRegionLoop
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
            int NUM_COLS,
            bool ROWS_EQ_COL = (NUM_ROWS)==(NUM_COLS) >
  struct StackVariables
  {
public:
    static constexpr int numRows = NUM_ROWS;
    static constexpr int numCols = NUM_COLS;

    //    GEOSX_HOST_DEVICE
    StackVariables():
      localRowDofIndex{ 0 },
      localColDofIndex{ 0 },
      localResidual{ 0.0 },
      localJacobian{ {0.0} }
    {}

    globalIndex localRowDofIndex[numRows];
    globalIndex localColDofIndex[numRows];
    real64 localResidual[numRows];
    real64 localJacobian[numRows][numCols];
  };


  template< int NUM_ROWS,
            int NUM_COLS >
  struct StackVariables< NUM_ROWS,
                         NUM_COLS,
                         true >
  {
public:
    static constexpr int numRows = NUM_ROWS;
    static constexpr int numCols = NUM_COLS;

    //    GEOSX_HOST_DEVICE
    StackVariables():
      localRowDofIndex{ 0 },
      localColDofIndex( localRowDofIndex ),
      localResidual{ 0.0 },
      localJacobian{ {0.0} }
    {}

    globalIndex localRowDofIndex[numRows];
    globalIndex * const localColDofIndex;   // non-memcopyable...problem if we are capturing.
    real64 localResidual[numRows];
    real64 localJacobian[numRows][numCols];
  };


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
  class Kernels
  {
public:

    static constexpr int numTestSupportPointsPerElem  = NUM_TEST_SUPPORT_POINTS_PER_ELEM;
    static constexpr int numTrialSupportPointsPerElem = NUM_TRIAL_SUPPORT_POINTS_PER_ELEM;
    static constexpr int numDofPerTestSupportPoint    = NUM_DOF_PER_TEST_SP;
    static constexpr int numDofPerTrialSupportPoint   = NUM_DOF_PER_TRIAL_SP;

    using StackVars = StackVariables< numTestSupportPointsPerElem*numDofPerTestSupportPoint,
                                      numTrialSupportPointsPerElem*numDofPerTrialSupportPoint >;

    Kernels( arrayView1d< globalIndex const > const & inputDofNumber,
             ParallelMatrix & inputMatrix,
             ParallelVector & inputRhs,
             NodeManager const & nodeManager,
             SUBREGION_TYPE const & elementSubRegion,
             FiniteElementBase const * const finiteElementSpace,
             CONSTITUTIVE_TYPE * const inputConstitutiveType,
             Parameters const & GEOSX_UNUSED_PARAM( parameters ) ):
      Kernels( inputDofNumber,
               inputMatrix,
               inputRhs,
               nodeManager,
               elementSubRegion,
               finiteElementSpace,
               inputConstitutiveType,
               typename CONSTITUTIVE_TYPE::KernelWrapper() )
    {}

//    Kernels() = delete;
////    Kernels( Kernels const & ) = delete;
//    Kernels & operator=( Kernels const & ) = delete;
//    Kernels( Kernels && ) = delete;
//    Kernels & operator=( Kernels && ) = delete;


    template< typename STACK_VARIABLE_TYPE >
    //  GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    void preKernel( localIndex const k,
                    STACK_VARIABLE_TYPE & stack ) const
    {
      for( localIndex a=0; a<numTestSupportPointsPerElem; ++a )
      {
        localIndex const localNodeIndex = elemsToNodes[ k][ a ];
        for( int i=0; i<numDofPerTestSupportPoint; ++i )
        {
          stack.localRowDofIndex[a*numDofPerTestSupportPoint+i] = m_dofNumber[localNodeIndex]+i;
        }
      }

      // TODO This is incorrect. The support points of the trial space is not necessarily the nodes.
      for( localIndex a=0; a<numTrialSupportPointsPerElem; ++a )
      {
        localIndex const localNodeIndex = elemsToNodes[ k][ a];
        for( int i=0; i<numDofPerTrialSupportPoint; ++i )
        {
          stack.localColDofIndex[a*numDofPerTrialSupportPoint+i] = m_dofNumber[localNodeIndex]+i;
        }
      }

    }

    template< typename STACK_VARIABLE_TYPE >
    //  GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    void updateKernel( localIndex const GEOSX_UNUSED_PARAM( k ),
                       localIndex const GEOSX_UNUSED_PARAM( q ),
                       STACK_VARIABLE_TYPE & GEOSX_UNUSED_PARAM( stack ) ) const
    {}

    template< typename PARAMETERS_TYPE, typename STACK_VARIABLE_TYPE >
    //  GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    void stiffnessKernel( localIndex const GEOSX_UNUSED_PARAM( k ),
                          localIndex const GEOSX_UNUSED_PARAM( q ),
                          PARAMETERS_TYPE const & GEOSX_UNUSED_PARAM( parameters ),
                          STACK_VARIABLE_TYPE & GEOSX_UNUSED_PARAM( stack ) ) const
    {}

    template< typename PARAMETERS_TYPE, typename STACK_VARIABLE_TYPE >
    //  GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    void integrationKernel( localIndex const GEOSX_UNUSED_PARAM( k ),
                            localIndex const GEOSX_UNUSED_PARAM( q ),
                            PARAMETERS_TYPE const & GEOSX_UNUSED_PARAM( parameters ),
                            STACK_VARIABLE_TYPE & GEOSX_UNUSED_PARAM( stack ) ) const
    {}

    template< typename PARAMETERS_TYPE, typename STACK_VARIABLE_TYPE >
    //  GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    real64 postKernel( PARAMETERS_TYPE const & GEOSX_UNUSED_PARAM( parameters ),
                       STACK_VARIABLE_TYPE & stack ) const
    {
      m_matrix.insert( stack.localRowDofIndex,
                       stack.localColDofIndex,
                       &(stack.localJacobian[0][0]),
                       stack.numRows,
                       stack.numCols );
      return 0;
    }

    arrayView1d< globalIndex const > const m_dofNumber;
    ParallelMatrix & m_matrix;
    ParallelVector & m_rhs;
    typename SUBREGION_TYPE::NodeMapType::base_type::ViewTypeConst const elemsToNodes;
    arrayView1d< integer const > const elemGhostRank;
    typename CONSTITUTIVE_TYPE::KernelWrapper const constitutiveUpdate;
    FiniteElementBase const * m_finiteElementSpace;

protected:
    Kernels( arrayView1d< globalIndex const > const & inputDofNumber,
             ParallelMatrix & inputMatrix,
             ParallelVector & inputRhs,
             NodeManager const & GEOSX_UNUSED_PARAM( nodeManager ),
             SUBREGION_TYPE const & elementSubRegion,
             FiniteElementBase const * const finiteElementSpace,
             CONSTITUTIVE_TYPE * const GEOSX_UNUSED_PARAM( inputConstitutiveType ),
             typename CONSTITUTIVE_TYPE::KernelWrapper const & inputConstitutiveUpdate ):
      m_dofNumber( inputDofNumber ),
      m_matrix( inputMatrix ),
      m_rhs( inputRhs ),
      elemsToNodes( elementSubRegion.nodeList().toViewConst() ),
      elemGhostRank( elementSubRegion.ghostRank() ),
      constitutiveUpdate( inputConstitutiveUpdate ),
      m_finiteElementSpace( finiteElementSpace )
    {}

  };

  //***************************************************************************
  template< typename POLICY,
            int NUM_QUADRATURE_POINTS,
            typename STACK_VARIABLES,
            typename PARAMETERS_TYPE,
            typename KERNEL_CLASS >
  static
  real64 Launch( localIndex const numElems,
                 PARAMETERS_TYPE const & parameters,
                 KERNEL_CLASS const & kernelClass )
  {
    RAJA::ReduceMax< serialReduce, real64 > maxResidual( 0 );

    forAll< POLICY >( numElems,
                      [=] ( localIndex const k )
    {
      STACK_VARIABLES stack;

      kernelClass.preKernel( k, stack );
      for( integer q=0; q<NUM_QUADRATURE_POINTS; ++q )
      {
        kernelClass.updateKernel( k, q, stack );

        kernelClass.stiffnessKernel( k, q, parameters, stack );

        kernelClass.integrationKernel( k, q, parameters, stack );
      }
      if( kernelClass.elemGhostRank[k] < 0 )
      {
        maxResidual.max( kernelClass.postKernel( parameters, stack ) );
      }
    } );
    return maxResidual.get();
  }


  //***************************************************************************
  template< typename POLICY,
            typename UPDATE_CLASS,
            typename CONSTITUTIVE_BASE,
            typename REGION_TYPE,
            typename PARAMETER_CLASS,
            template< typename SUBREGION_TYPE,
                      typename CONSTITUTIVE_TYPE,
                      int NUM_TEST_SUPPORT_POINTS_PER_ELEM,
                      int NUM_TRIAL_SUPPORT_POINTS_PER_ELEM > class KERNEL_CLASS = UPDATE_CLASS::template Kernels >
  static
  real64 Execute( MeshLevel & mesh,
                  arrayView1d< string const > const & targetRegions,
                  arrayView1d< string const > const & constitutiveNames,
                  FiniteElementDiscretization const * const feDiscretization,
                  arrayView1d< globalIndex const > const & inputDofNumber,
                  ParallelMatrix & inputMatrix,
                  ParallelVector & inputRhs,
                  PARAMETER_CLASS const & parameters )
  {

    real64 maxResidual = 0;

    NodeManager const & nodeManager = *(mesh.getNodeManager());
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

          KERNEL_CLASS< SUBREGIONTYPE,
                        CONSTITUTIVE_TYPE,
                        NUM_NODES_PER_ELEM,
                        NUM_NODES_PER_ELEM > kernelClass( inputDofNumber,
                                                          inputMatrix,
                                                          inputRhs,
                                                          nodeManager,
                                                          elementSubRegion,
                                                          finiteElementSpace,
                                                          castedConstitutiveRelation,
                                                          parameters );

          maxResidual = std::max( maxResidual,
                                  Launch< POLICY,
                                          NUM_QUADRATURE_POINTS,
                                          typename decltype(kernelClass)::StackVars >( numElems,
                                                                                       parameters,
                                                                                       kernelClass ) );
        } );
      } );
    } );

    return maxResidual;
  }



  //***************************************************************************
  template< typename POLICY,
            typename UPDATE_CLASS,
            typename REGION_TYPE >
  static
  real64 FillSparsity( MeshLevel & mesh,
                       arrayView1d< string const > const & targetRegions,
                       FiniteElementDiscretization const * const feDiscretization,
                       arrayView1d< globalIndex const > const & inputDofNumber,
                       ParallelMatrix & inputMatrix,
                       ParallelVector & inputRhs )
  {
    return Execute< POLICY,
                    UPDATE_CLASS,
                    constitutive::Dummy,
                    REGION_TYPE,
                    FiniteElementRegionLoop::Parameters,
                    UPDATE_CLASS::template SparsityKernels >( mesh,
                                                              targetRegions,
                                                              array1d<string>(),
                                                              feDiscretization,
                                                              inputDofNumber,
                                                              inputMatrix,
                                                              inputRhs,
                                                              FiniteElementRegionLoop::Parameters() );
  }


};

}
}



#endif /* GEOSX_PHYSICSSOLVERS_PHYSICSLOOPINTERFACE_HPP_ */
