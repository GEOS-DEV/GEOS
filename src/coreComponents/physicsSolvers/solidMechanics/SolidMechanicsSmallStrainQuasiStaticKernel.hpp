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
 * @file SolidMechanicsSmallStrainQuasiStaticKernels.hpp
 */

#pragma once

#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveBase.hpp"
#include "finiteElement/ElementLibrary/FiniteElementBase.h"
#include "finiteElement/FiniteElementShapeFunctionKernel.hpp"
#include "finiteElement/Kinematics.h"
#include "finiteElement/kernelInterface/ImplicitKernelBase.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"
#include "TimeIntegrationOption.hpp"

namespace geosx
{

namespace SolidMechanicsLagrangianFEMKernels
{

class QuasiStatic : public finiteElement::ImplicitKernelBase
{
public:
  using BaseKernel = finiteElement::ImplicitKernelBase;
  static constexpr int numTestDofPerSP = 3;
  static constexpr int numTrialDofPerSP = 3;

  struct Parameters : public BaseKernel::Parameters
  {
    Parameters( real64 const inputGravityVector[3] ):
      m_gravityVector{ inputGravityVector[0], inputGravityVector[1], inputGravityVector[2] }
    {}

    real64 const m_gravityVector[3];
  };


  template< int NUM_TEST_SUPPORT_POINTS_PER_ELEM,
            int NUM_TRIAL_SUPPORT_POINTS_PER_ELEM >
  struct StackVariables : BaseKernel::StackVariables< NUM_TEST_SUPPORT_POINTS_PER_ELEM *numTestDofPerSP,
                                                      NUM_TRIAL_SUPPORT_POINTS_PER_ELEM *numTrialDofPerSP >
  {
public:
    using StackVariablesBase = BaseKernel::StackVariables< NUM_TEST_SUPPORT_POINTS_PER_ELEM *numTestDofPerSP,
                                                           NUM_TRIAL_SUPPORT_POINTS_PER_ELEM *numTrialDofPerSP >;
    using StackVariablesBase::numRows;
    using StackVariablesBase::numCols;
    static constexpr int numNodes = NUM_TEST_SUPPORT_POINTS_PER_ELEM;


    GEOSX_HOST_DEVICE
    StackVariables():
      StackVariablesBase(),
      u_local(),
      uhat_local(),
      constitutiveStiffness{ {0.0} }
    {}

    R1Tensor u_local[numNodes];
    R1Tensor uhat_local[numNodes];
    real64 constitutiveStiffness[6][6];
  };

  template< typename SUBREGION_TYPE,
            typename CONSTITUTIVE_TYPE,
            int NUM_NODES_PER_ELEM,
            int >
  using SparsityComponents = BaseKernel::Components< SUBREGION_TYPE,
                                                     CONSTITUTIVE_TYPE,
                                                     NUM_NODES_PER_ELEM,
                                                     NUM_NODES_PER_ELEM,
                                                     numTestDofPerSP,
                                                     numTrialDofPerSP >;

  template< typename SUBREGION_TYPE,
            typename CONSTITUTIVE_TYPE,
            int NUM_NODES_PER_ELEM,
            int >
  class Components : public BaseKernel::Components< SUBREGION_TYPE,
                                                    CONSTITUTIVE_TYPE,
                                                    NUM_NODES_PER_ELEM,
                                                    NUM_NODES_PER_ELEM,
                                                    numTestDofPerSP,
                                                    numTrialDofPerSP >
  {
public:
    using ComponentsBase = BaseKernel::Components< SUBREGION_TYPE,
                                                   CONSTITUTIVE_TYPE,
                                                   NUM_NODES_PER_ELEM,
                                                   NUM_NODES_PER_ELEM,
                                                   numTestDofPerSP,
                                                   numTrialDofPerSP >;

    static constexpr int numNodesPerElem = NUM_NODES_PER_ELEM;

    using ComponentsBase::m_dofNumber;
    using ComponentsBase::m_matrix;
    using ComponentsBase::m_rhs;
    using ComponentsBase::elemsToNodes;
    using ComponentsBase::constitutiveUpdate;
    using ComponentsBase::m_finiteElementSpace;

    using StackVars = StackVariables< numNodesPerElem,
                                      numNodesPerElem >;



    Components( arrayView1d< globalIndex const > const & inputDofNumber,
                ParallelMatrix & inputMatrix,
                ParallelVector & inputRhs,
                NodeManager const & nodeManager,
                SUBREGION_TYPE const & elementSubRegion,
                FiniteElementBase const * const finiteElementSpace,
                CONSTITUTIVE_TYPE * const inputConstitutiveType,
                Parameters const & parameters ):
      ComponentsBase( inputDofNumber,
                      inputMatrix,
                      inputRhs,
                      nodeManager,
                      elementSubRegion,
                      finiteElementSpace,
                      inputConstitutiveType,
                      parameters ),
      m_disp( nodeManager.totalDisplacement()),
      m_uhat( nodeManager.incrementalDisplacement()),
      dNdX( elementSubRegion.template getReference< array3d< R1Tensor > >( dataRepository::keys::dNdX )),
      detJ( elementSubRegion.template getReference< array2d< real64 > >( dataRepository::keys::detJ ) ),
      m_gravityVector{ parameters.m_gravityVector[0], parameters.m_gravityVector[1], parameters.m_gravityVector[2] }
    {}

    arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const m_disp;
    arrayView2d< real64 const, nodes::INCR_DISPLACEMENT_USD > const m_uhat;

    arrayView3d< R1Tensor const > const dNdX;
    arrayView2d< real64 const > const detJ;
    real64 const m_gravityVector[3];


    template< typename STACK_VARIABLE_TYPE >
    GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    void setup( localIndex const k,
                STACK_VARIABLE_TYPE & stack ) const
    {
      for( localIndex a=0; a<NUM_NODES_PER_ELEM; ++a )
      {
        localIndex const localNodeIndex = elemsToNodes( k, a );

        stack.u_local[ a ] = m_disp[ localNodeIndex ];
        stack.uhat_local[ a ] = m_uhat[ localNodeIndex ];

        for( int i=0; i<3; ++i )
        {
          stack.localRowDofIndex[a*3+i] = m_dofNumber[localNodeIndex]+i;
          stack.localColDofIndex[a*3+i] = m_dofNumber[localNodeIndex]+i;
        }
      }

    }

    template< typename STACK_VARIABLE_TYPE >
    GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    void quadraturePointStateUpdate( localIndex const k,
                                     localIndex const q,
                                     STACK_VARIABLE_TYPE & stack ) const
    {
      real64 strainInc[6] = {0};
      for( localIndex a = 0; a < NUM_NODES_PER_ELEM; ++a )
      {
        strainInc[0] = strainInc[0] + dNdX( k, q, a )[0] * stack.uhat_local[a][0];
        strainInc[1] = strainInc[1] + dNdX( k, q, a )[1] * stack.uhat_local[a][1];
        strainInc[2] = strainInc[2] + dNdX( k, q, a )[2] * stack.uhat_local[a][2];
        strainInc[3] = strainInc[3] + dNdX( k, q, a )[2] * stack.uhat_local[a][1] +
                       dNdX( k, q, a )[1] * stack.uhat_local[a][2];
        strainInc[4] = strainInc[4] + dNdX( k, q, a )[2] * stack.uhat_local[a][0] +
                       dNdX( k, q, a )[0] * stack.uhat_local[a][2];
        strainInc[5] = strainInc[5] + dNdX( k, q, a )[1] * stack.uhat_local[a][0] +
                       dNdX( k, q, a )[0] * stack.uhat_local[a][1];
      }

      constitutiveUpdate.SmallStrain( k, q, strainInc );

      GEOSX_UNUSED_VAR( q )
      constitutiveUpdate.GetStiffness( k, stack.constitutiveStiffness );
    }


    template< typename STACK_VARIABLE_TYPE /*,
                                              typename DYNAMICS_LAMBDA = nvstd::function< void( localIndex, localIndex)
                                                 >*/>
    GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    void quadraturePointJacobianContribution( localIndex const k,
                                              localIndex const q,
                                              STACK_VARIABLE_TYPE & stack /*,
                                                                             DYNAMICS_LAMBDA && dynamicsTerms = []
                                                                                GEOSX_DEVICE ( localIndex,
                                                                                localIndex){}*/) const
    {
      for( localIndex a=0; a<NUM_NODES_PER_ELEM; ++a )
      {
        for( localIndex b=0; b<NUM_NODES_PER_ELEM; ++b )
        {
          real64 const (&c)[6][6] = stack.constitutiveStiffness;
          stack.localJacobian[ a*3+0 ][ b*3+0 ] -= ( c[0][0]*dNdX( k, q, a )[0]*dNdX( k, q, b )[0] +
                                                     c[5][5]*dNdX( k, q, a )[1]*dNdX( k, q, b )[1] +
                                                     c[4][4]*dNdX( k, q, a )[2]*dNdX( k, q, b )[2] ) * detJ( k, q );

          stack.localJacobian[ a*3+0 ][ b*3+1 ] -= ( c[5][5]*dNdX( k, q, a )[1]*dNdX( k, q, b )[0] +
                                                     c[0][1]*dNdX( k, q, a )[0]*dNdX( k, q, b )[1] ) * detJ( k, q );

          stack.localJacobian[ a*3+0 ][ b*3+2 ] -= ( c[4][4]*dNdX( k, q, a )[2]*dNdX( k, q, b )[0] +
                                                     c[0][2]*dNdX( k, q, a )[0]*dNdX( k, q, b )[2] ) * detJ( k, q );

          stack.localJacobian[ a*3+1 ][ b*3+1 ] -= ( c[5][5]*dNdX( k, q, a )[0]*dNdX( k, q, b )[0] +
                                                     c[1][1]*dNdX( k, q, a )[1]*dNdX( k, q, b )[1] +
                                                     c[3][3]*dNdX( k, q, a )[2]*dNdX( k, q, b )[2] ) * detJ( k, q );

          stack.localJacobian[ a*3+1 ][ b*3+0 ] -= ( c[0][1]*dNdX( k, q, a )[1]*dNdX( k, q, b )[0] +
                                                     c[5][5]*dNdX( k, q, a )[0]*dNdX( k, q, b )[1] ) * detJ( k, q );

          stack.localJacobian[ a*3+1 ][ b*3+2 ] -= ( c[3][3]*dNdX( k, q, a )[2]*dNdX( k, q, b )[1] +
                                                     c[1][2]*dNdX( k, q, a )[1]*dNdX( k, q, b )[2] ) * detJ( k, q );

          stack.localJacobian[ a*3+2 ][ b*3+0 ] -= ( c[0][2]*dNdX( k, q, a )[2]*dNdX( k, q, b )[0] +
                                                     c[4][4]*dNdX( k, q, a )[0]*dNdX( k, q, b )[2] ) * detJ( k, q );

          stack.localJacobian[ a*3+2 ][ b*3+1 ] -= ( c[1][2]*dNdX( k, q, a )[2]*dNdX( k, q, b )[1] +
                                                     c[3][3]*dNdX( k, q, a )[1]*dNdX( k, q, b )[2] ) * detJ( k, q );

          stack.localJacobian[ a*3+2 ][ b*3+2 ] -= ( c[4][4]*dNdX( k, q, a )[0]*dNdX( k, q, b )[0] +
                                                     c[3][3]*dNdX( k, q, a )[1]*dNdX( k, q, b )[1] +
                                                     c[2][2]*dNdX( k, q, a )[2]*dNdX( k, q, b )[2] ) * detJ( k, q );

//          dynamicsTerms( a, b );
        }
      }
    }

    template< typename STACK_VARIABLE_TYPE /*,
                                              typename DYNAMICS_LAMBDA = STD_FUNCTION< void( real64 * ) >*/>
    GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    void quadraturePointResidualContribution( localIndex const k,
                                              localIndex const q,
                                              STACK_VARIABLE_TYPE & stack /*,
                                                                             DYNAMICS_LAMBDA && stressModifier = []
                                                                                GEOSX_DEVICE ( real64 * ) {}*/) const
    {
      real64 stress[6] = { constitutiveUpdate.m_stress( k, q, 0 ),
                           constitutiveUpdate.m_stress( k, q, 1 ),
                           constitutiveUpdate.m_stress( k, q, 2 ),
                           constitutiveUpdate.m_stress( k, q, 3 ),
                           constitutiveUpdate.m_stress( k, q, 4 ),
                           constitutiveUpdate.m_stress( k, q, 5 ) };

//      stressModifier( stress );

      for( localIndex a = 0; a < NUM_NODES_PER_ELEM; ++a )
      {
        stack.localResidual[ a * 3 + 0 ] -= ( stress[ 0 ] * dNdX( k, q, a )[ 0 ] +
                                              stress[ 5 ] * dNdX( k, q, a )[ 1 ] +
                                              stress[ 4 ] * dNdX( k, q, a )[ 2 ] -
                                              m_gravityVector[0] ) * detJ( k, q );
        stack.localResidual[ a * 3 + 1 ] -= ( stress[ 5 ] * dNdX( k, q, a )[ 0 ] +
                                              stress[ 1 ] * dNdX( k, q, a )[ 1 ] +
                                              stress[ 3 ] * dNdX( k, q, a )[ 2 ] -
                                              m_gravityVector[1] ) * detJ( k, q );
        stack.localResidual[ a * 3 + 2 ] -= ( stress[ 4 ] * dNdX( k, q, a )[ 0 ] +
                                              stress[ 3 ] * dNdX( k, q, a )[ 1 ] +
                                              stress[ 2 ] * dNdX( k, q, a )[ 2 ] -
                                              m_gravityVector[2] ) * detJ( k, q );
      }
    }

    template< typename STACK_VARIABLE_TYPE >
    //GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    real64 complete( localIndex const GEOSX_UNUSED_PARAM( k ),
                     STACK_VARIABLE_TYPE & stack ) const
    {
      real64 meanForce = 0;
      for( localIndex a=0; a<stack.numRows; ++a )
      {
//        RAJA::atomicMax< RAJA::auto_atomic >( &meanForce, stack.localResidual[a] );
        meanForce = std::max( meanForce, stack.localResidual[a] );
//                meanForce += fabs( stack.localResidual[a] );
      }
//            meanForce /= stack.ndof;

      m_matrix.add( stack.localRowDofIndex,
                    stack.localColDofIndex,
                    &(stack.localJacobian[0][0]),
                    stack.numRows,
                    stack.numCols );

      m_rhs.add( stack.localRowDofIndex,
                 stack.localResidual,
                 stack.numRows );

      return meanForce;
    }

  };
};

} // namespace SolidMechanicsLagrangianFEMKernels

} // namespace geosx
