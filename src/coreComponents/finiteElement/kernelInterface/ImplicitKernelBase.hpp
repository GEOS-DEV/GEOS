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

#include "KernelBase.hpp"

/**
 * @file RegionLoopSparsity.hpp
 */

#ifndef GEOSX_FINITEELEMENT_REGIONLOOPSPARSITY_HPP_
#define GEOSX_FINITEELEMENT_REGIONLOOPSPARSITY_HPP_




namespace geosx
{

/**
 * @namespace Contains implementations of physics loops.
 */
namespace finiteElement
{

/**
 * @class FiniteElementRegionLoop
 *
 * This class encapsulates base components and interface for applying a finite
 * element method over a loop of element regions.
 *
 */
class ImplicitKernelBase : public KernelBase
{
public:

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
  struct StackVariables : KernelBase::StackVariables<NUM_ROWS,NUM_COLS>
  {
public:
    static constexpr int numRows = NUM_ROWS;
    static constexpr int numCols = NUM_COLS;

    GEOSX_HOST_DEVICE
    StackVariables():
      localRowDofIndex{ 0 },
      localColDofIndex{ 0 },
      localResidual{ 0.0 },
      localJacobian{ {0.0} }
    {}

    globalIndex localRowDofIndex[numRows];
    globalIndex localColDofIndex[numCols];
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
  class Components : KernelBase::Components< SUBREGION_TYPE,
                                      CONSTITUTIVE_TYPE,
                                      NUM_TEST_SUPPORT_POINTS_PER_ELEM,
                                      NUM_TRIAL_SUPPORT_POINTS_PER_ELEM,
                                      NUM_DOF_PER_TEST_SP,
                                      NUM_DOF_PER_TRIAL_SP >
  {
public:

    using ComponentsBase = KernelBase::Components< SUBREGION_TYPE,
                                             CONSTITUTIVE_TYPE,
                                             NUM_TEST_SUPPORT_POINTS_PER_ELEM,
                                             NUM_TRIAL_SUPPORT_POINTS_PER_ELEM,
                                             NUM_DOF_PER_TEST_SP,
                                             NUM_DOF_PER_TRIAL_SP >;

    static constexpr int numTestSupportPointsPerElem  = NUM_TEST_SUPPORT_POINTS_PER_ELEM;
    static constexpr int numTrialSupportPointsPerElem = NUM_TRIAL_SUPPORT_POINTS_PER_ELEM;
    static constexpr int numDofPerTestSupportPoint    = NUM_DOF_PER_TEST_SP;
    static constexpr int numDofPerTrialSupportPoint   = NUM_DOF_PER_TRIAL_SP;

    using StackVars = StackVariables< numTestSupportPointsPerElem*numDofPerTestSupportPoint,
                                      numTrialSupportPointsPerElem*numDofPerTrialSupportPoint >;

    using ComponentsBase::elemsToNodes;
    using ComponentsBase::elemGhostRank;
    using ComponentsBase::constitutiveUpdate;
    using ComponentsBase::m_finiteElementSpace;
    using ComponentsBase::Launch;


    Components( arrayView1d< globalIndex const > const & inputDofNumber,
             ParallelMatrix & inputMatrix,
             ParallelVector & inputRhs,
             NodeManager const & GEOSX_UNUSED_PARAM( nodeManager ),
             SUBREGION_TYPE const & elementSubRegion,
             FiniteElementBase const * const finiteElementSpace,
             CONSTITUTIVE_TYPE * const inputConstitutiveType,
             Parameters const & GEOSX_UNUSED_PARAM( parameters ) ):
      ComponentsBase( elementSubRegion,
                   finiteElementSpace,
                   inputConstitutiveType,
                   inputConstitutiveType->createKernelWrapper() ),
      m_dofNumber( inputDofNumber ),
      m_matrix( inputMatrix ),
      m_rhs( inputRhs )
    {}

    template< typename STACK_VARIABLE_TYPE >
    GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    void setup( localIndex const k,
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
    GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    void quadraturePointStateUpdate( localIndex const GEOSX_UNUSED_PARAM( k ),
                       localIndex const GEOSX_UNUSED_PARAM( q ),
                       STACK_VARIABLE_TYPE & GEOSX_UNUSED_PARAM( stack ) ) const
    {}

    template< typename PARAMETERS_TYPE, typename STACK_VARIABLE_TYPE >
    GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    void quadraturePointJacobianContribution( localIndex const GEOSX_UNUSED_PARAM( k ),
                          localIndex const GEOSX_UNUSED_PARAM( q ),
                          PARAMETERS_TYPE const & GEOSX_UNUSED_PARAM( parameters ),
                          STACK_VARIABLE_TYPE & GEOSX_UNUSED_PARAM( stack ) ) const
    {}

    template< typename PARAMETERS_TYPE, typename STACK_VARIABLE_TYPE >
    GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    void quadraturePointResidualContribution( localIndex const GEOSX_UNUSED_PARAM( k ),
                            localIndex const GEOSX_UNUSED_PARAM( q ),
                            PARAMETERS_TYPE const & GEOSX_UNUSED_PARAM( parameters ),
                            STACK_VARIABLE_TYPE & GEOSX_UNUSED_PARAM( stack ) ) const
    {}

    template< typename PARAMETERS_TYPE, typename STACK_VARIABLE_TYPE >
//    GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    real64 complete( localIndex const GEOSX_UNUSED_PARAM( k ),
                       PARAMETERS_TYPE const & GEOSX_UNUSED_PARAM( parameters ),
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
  };
};

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
  return RegionBasedKernelApplication< POLICY,
                  ImplicitKernelBase,
                  constitutive::Dummy,
                  REGION_TYPE,
                  UPDATE_CLASS::template SparsityComponents >( mesh,
                                                            targetRegions,
                                                            array1d<string>(),
                                                            feDiscretization,
                                                            inputDofNumber,
                                                            inputMatrix,
                                                            inputRhs,
                                                            ImplicitKernelBase::Parameters() );
}

}
}



#endif /* GEOSX_FINITEELEMENT_REGIONLOOPSPARSITY_HPP_ */
