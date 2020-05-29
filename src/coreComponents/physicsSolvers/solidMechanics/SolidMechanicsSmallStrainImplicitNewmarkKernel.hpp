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
 * @file SolidMechanicsSmallStrainImplicitNewmarkKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSSMALLSTRAINIMPLICITNEWMARK_HPP_
#define GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSSMALLSTRAINIMPLICITNEWMARK_HPP_

#include "SolidMechanicsSmallStrainQuasiStaticKernel.hpp"


namespace geosx
{

namespace SolidMechanicsLagrangianFEMKernels
{

struct ImplicitNewmarkConstructorParams : QuasiStaticConstructorParams
{
  ImplicitNewmarkConstructorParams( arrayView1d< globalIndex const > const & inputDofNumber,
                     ParallelMatrix & inputMatrix,
                     ParallelVector & inputRhs,
                     real64 const inputGravityVector[3],
                     real64 const inputNewmarkGamma,
                     real64 const inputNewmarkBeta,
                     real64 const inputMassDamping,
                     real64 const inputStiffnessDamping,
                     real64 const inputDt):
                       QuasiStaticConstructorParams( inputDofNumber,
                             inputMatrix,
                             inputRhs,
                             inputGravityVector ),
    m_newmarkGamma(inputNewmarkGamma),
    m_newmarkBeta(inputNewmarkBeta),
    m_massDamping(inputMassDamping),
    m_stiffnessDamping(inputStiffnessDamping),
    m_dt(inputDt)
  {}

  real64 const m_newmarkGamma;
  real64 const m_newmarkBeta;
  real64 const m_massDamping;
  real64 const m_stiffnessDamping;
  real64 const m_dt;
};

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          int NUM_NODES_PER_ELEM,
          int >
class ImplicitNewmark : public QuasiStatic< SUBREGION_TYPE,
                                            CONSTITUTIVE_TYPE,
                                            NUM_NODES_PER_ELEM,
                                            NUM_NODES_PER_ELEM >
{
public:
  using Base = QuasiStatic< SUBREGION_TYPE,
                            CONSTITUTIVE_TYPE,
                            NUM_NODES_PER_ELEM,
                            NUM_NODES_PER_ELEM >;

  using Base::numNodesPerElem;
  using Base::numTestSupportPointsPerElem;
  using Base::numTrialSupportPointsPerElem;
  using Base::numDofPerTestSupportPoint;
  using Base::numDofPerTrialSupportPoint;

  using Base::m_dofNumber;
  using Base::m_matrix;
  using Base::m_rhs;
  using Base::m_elemsToNodes;
  using Base::m_disp;
  using Base::m_uhat;
  using Base::m_detJ;

  /// Alias for the struct that holds the constructor parameters
  using ConstructorParams = ImplicitNewmarkConstructorParams;

  ImplicitNewmark( NodeManager const & nodeManager,
                   EdgeManager const & edgeManager,
                   FaceManager const & faceManager,
                   SUBREGION_TYPE const & elementSubRegion,
                   FiniteElementBase const * const finiteElementSpace,
                   CONSTITUTIVE_TYPE * const inputConstitutiveType,
                   arrayView1d< globalIndex const > const & inputDofNumber,
                   ParallelMatrix & inputMatrix,
                   ParallelVector & inputRhs,
                   real64 const inputGravityVector[3],
                   real64 const inputNewmarkGamma,
                   real64 const inputNewmarkBeta,
                   real64 const inputMassDamping,
                   real64 const inputStiffnessDamping,
                   real64 const inputDt ):
    Base( nodeManager,
          edgeManager,
          faceManager,
          elementSubRegion,
          finiteElementSpace,
          inputConstitutiveType,
          inputDofNumber,
          inputMatrix,
          inputRhs,
          inputGravityVector ),
    m_vtilde( nodeManager.totalDisplacement()),
    m_uhattilde( nodeManager.totalDisplacement()),
    m_density( inputConstitutiveType->getDensity()),
    m_newmarkGamma( inputNewmarkGamma ),
    m_newmarkBeta( inputNewmarkBeta ),
    m_massDamping( inputMassDamping ),
    m_stiffnessDamping( inputStiffnessDamping ),
    m_dt( inputDt )
  {}

  ImplicitNewmark( NodeManager const & nodeManager,
               EdgeManager const & edgeManager,
               FaceManager const & faceManager,
               SUBREGION_TYPE const & elementSubRegion,
               FiniteElementBase const * const finiteElementSpace,
               CONSTITUTIVE_TYPE * const inputConstitutiveType,
               ConstructorParams & params ):
    ImplicitNewmark( nodeManager,
                     edgeManager,
                     faceManager,
                     elementSubRegion,
                     finiteElementSpace,
                     inputConstitutiveType,
                     params.m_dofNumber,
                     params.m_matrix,
                     params.m_rhs,
                     params.m_gravityVector,
                     params.m_newmarkGamma,
                     params.m_newmarkBeta,
                     params.m_massDamping,
                     params.m_stiffnessDamping,
                     params.m_dt )
  {}

  template< typename STACK_VARIABLE_TYPE >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void setup( localIndex const k,
              STACK_VARIABLE_TYPE & stack ) const
  {
    for( localIndex a=0; a<numNodesPerElem; ++a )
    {
      localIndex const localNodeIndex = m_elemsToNodes( k, a );

      stack.u_local[ a ] = m_disp[ localNodeIndex ];
      stack.uhat_local[ a ] = m_uhat[ localNodeIndex ];
      stack.vtilde_local[ a ] = m_vtilde[ localNodeIndex ];
      stack.uhattilde_local[ a ] = m_uhattilde[ localNodeIndex ];
    }
    Base::setup( k, stack );
  }

  template< typename STACK_VARIABLE_TYPE >
  GEOSX_DEVICE
  GEOSX_FORCE_INLINE
  void quadraturePointJacobianContribution( localIndex const k,
                                            localIndex const q,
                                            STACK_VARIABLE_TYPE & stack ) const
  {

    real64 N[numNodesPerElem];
    FiniteElementShapeKernel::shapeFunctionValues( q, N );

    Base::quadraturePointJacobianContribution( k, q, stack, [&] GEOSX_DEVICE ( localIndex const a,
                                                                               localIndex const b ) mutable
    {
      real64 const integrationFactor = m_density( k, q ) * N[a] * N[b] * m_detJ( k, q );
      real64 const temp1 = ( m_massDamping * m_newmarkGamma/( m_newmarkBeta * m_dt )
                       + 1.0 / ( m_newmarkBeta * m_dt * m_dt ) )* integrationFactor;

      constexpr int nsdof = numDofPerTestSupportPoint;
      for( int i=0; i<nsdof; ++i )
      {
        realT const acc = 1.0 / ( m_newmarkBeta * m_dt * m_dt ) * ( stack.uhat_local[b][i] - stack.uhattilde_local[b][i] );
        realT const vel = stack.vtilde_local[b][i] +
                          m_newmarkGamma/( m_newmarkBeta * m_dt ) *( stack.uhat_local[b][i]
                                                                     - stack.uhattilde_local[b][i] );

        stack.dRdU_InertiaMassDamping[ a*nsdof+i][ b*nsdof+i ] -= temp1;
        stack.localResidual[ a*nsdof+i ] -= ( m_massDamping * vel + acc ) * integrationFactor;
      }
    } );
  }

  template< typename STACK_VARIABLE_TYPE >
  //    GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  real64 complete( localIndex const k,
                   STACK_VARIABLE_TYPE & stack ) const
  {

    for( int a=0; a<numNodesPerElem; ++a )
    {
      for( int b=0; b<numNodesPerElem; ++b )
      {
        for( int i=0; i<numDofPerTestSupportPoint; ++i )
        {
          for( int j=0; j<numDofPerTrialSupportPoint; ++j )
          {
            stack.localResidual[ a*numDofPerTestSupportPoint+i ] =
              stack.localResidual[ a*numDofPerTestSupportPoint+i ] +
              m_stiffnessDamping * stack.localJacobian[ a*numDofPerTestSupportPoint+i][ b*numDofPerTrialSupportPoint+j ] *
              ( stack.vtilde_local[b][j] + m_newmarkGamma/(m_newmarkBeta * m_dt)*(stack.uhat_local[b][j]-stack.uhattilde_local[b][j]) );

            stack.localJacobian[a*numDofPerTestSupportPoint+i][b*numDofPerTrialSupportPoint+j] =
              stack.localJacobian[a*numDofPerTestSupportPoint+i][b*numDofPerTrialSupportPoint+j] +
              stack.localJacobian[a][b] * (1.0 + m_stiffnessDamping * m_newmarkGamma / ( m_newmarkBeta * m_dt ) ) +
              stack.dRdU_InertiaMassDamping[ a*numDofPerTestSupportPoint+i ][ b*numDofPerTrialSupportPoint+j ];
          }
        }
      }
    }

    for( int a=0; a<stack.numRows; ++a )
    {
      for( int b=0; b<stack.numCols; ++b )
      {
        stack.localJacobian[a][b] += stack.localJacobian[a][b] * (1.0 + m_stiffnessDamping * m_newmarkGamma / ( m_newmarkBeta * m_dt ) )
                                     + stack.dRdU_InertiaMassDamping[ a ][ b ];
      }
    }

    return Base::complete( k, stack );
  }



  struct StackVariables : public Base::StackVariables
  {
public:
    using Base::StackVariables::numRows;
    using Base::StackVariables::numCols;

    GEOSX_HOST_DEVICE
    StackVariables():
      Base::StackVariables(),
            dRdU_InertiaMassDamping{ {0.0} },
      vtilde_local(),
      uhattilde_local()
    {}

    real64 dRdU_InertiaMassDamping[ numRows ][ numCols ];
    R1Tensor vtilde_local[numNodesPerElem];
    R1Tensor uhattilde_local[numNodesPerElem];
  };


protected:
  arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const m_vtilde;
  arrayView2d< real64 const, nodes::INCR_DISPLACEMENT_USD > const m_uhattilde;
  arrayView2d< real64 const > const m_density;
  real64 const m_newmarkGamma;
  real64 const m_newmarkBeta;
  real64 const m_massDamping;
  real64 const m_stiffnessDamping;
  real64 const m_dt;

};

} // namespace SolidMechanicsLagrangianFEMKernels

} // namespace geosx

#endif //GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSSMALLSTRAINIMPLICITNEWMARK_HPP_
