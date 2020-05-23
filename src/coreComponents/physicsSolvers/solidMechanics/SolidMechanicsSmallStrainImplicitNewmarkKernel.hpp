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

#pragma once

#include "SolidMechanicsSmallStrainQuasiStaticKernel.hpp"


namespace geosx
{

namespace SolidMechanicsLagrangianFEMKernels
{

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
  using Base::numDofPerNode;
  using Base::ndof;


  struct StackVariables : Base::StackVariables
  {
public:

    GEOSX_HOST_DEVICE
    StackVariables():
      Base::StackVariables(),
            dRdU_InertiaMassDamping{ {0.0} },
      vtilde_local{ { 0.0, 0.0, 0.0} },
      uhattilde_local{ { 0.0, 0.0, 0.0} }
    {}

    real64 dRdU_InertiaMassDamping[ ndof ][ ndof ];
    R1Tensor vtilde_local[NUM_NODES_PER_ELEM];
    R1Tensor uhattilde_local[NUM_NODES_PER_ELEM];
  };

  using Base::m_dofNumber;
  using Base::m_matrix;
  using Base::m_rhs;
  using Base::elemsToNodes;
  using Base::constitutiveUpdate;
  using Base::m_disp;
  using Base::m_uhat;
  using Base::dNdX;
  using Base::detJ;
  using Base::m_finiteElementSpace;



  ImplicitNewmark( arrayView1d< globalIndex const > const & inputDofNumber,
                   ParallelMatrix & inputMatrix,
                   ParallelVector & inputRhs,
                   NodeManager const & nodeManager,
                   SUBREGION_TYPE const & elementSubRegion,
                   FiniteElementBase const * const finiteElementSpace,
                   CONSTITUTIVE_TYPE & constitutiveModel,
                   real64 const inputGravityVector[3],
                   real64 const inputNewmarkGamma,
                   real64 const inputNewmarkBeta,
                   real64 const inputMassDamping,
                   real64 const inputStiffnessDamping,
                   real64 const inputDt ):
    Base( inputDofNumber,
          inputMatrix,
          inputRhs,
          nodeManager,
          elementSubRegion,
          finiteElementSpace,
          constitutiveModel,
          inputGravityVector ),
    m_vtilde( nodeManager.totalDisplacement()),
    m_uhattilde( nodeManager.totalDisplacement()),
    m_density( constitutiveModel.getDensity()),
    m_newmarkGamma( inputNewmarkGamma ),
    m_newmarkBeta( inputNewmarkBeta ),
    m_massDamping( inputMassDamping ),
    m_stiffnessDamping( inputStiffnessDamping ),
    m_dt( inputDt )
  {}

  arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const m_vtilde;
  arrayView2d< real64 const, nodes::INCR_DISPLACEMENT_USD > const m_uhattilde;
  arrayView2d< real64 const > const m_density;
  real64 const m_newmarkGamma;
  real64 const m_newmarkBeta;
  real64 const m_massDamping;
  real64 const m_stiffnessDamping;
  real64 const m_dt;

  template< typename STACK_VARIABLE_TYPE >
  //    GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void setup( localIndex const k,
              STACK_VARIABLE_TYPE & stack ) const
  {
    for( localIndex a=0; a<STACK_VARIABLE_TYPE::numNodesPerElem; ++a )
    {
      localIndex const localNodeIndex = elemsToNodes( k, a );

      stack.u_local[ a ] = m_disp[ localNodeIndex ];
      stack.uhat_local[ a ] = m_uhat[ localNodeIndex ];
      stack.vtilde_local[ a ] = m_vtilde[ localNodeIndex ];
      stack.uhattilde_local[ a ] = m_uhattilde[ localNodeIndex ];

      for( int i=0; i<3; ++i )
      {
        stack.elementLocalDofIndex[a*3+i] = m_dofNumber[localNodeIndex]+i;
      }
    }

  }

  template< typename STACK_VARIABLE_TYPE >
  //    GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void quadraturePointJacobianContribution( localIndex const k,
                                            localIndex const q,
                                            STACK_VARIABLE_TYPE & stack ) const
  {

    std::vector< double > const & N = m_finiteElementSpace->values( q );

//      real64 N[STACK_VARIABLE_TYPE::numNodesPerElem];

    Base::quadraturePointJacobianContribution( k, q, stack, [&]( localIndex const a, localIndex const b )
    {
      real64 integrationFactor = m_density( k, q ) * N[a] * N[b] * detJ( k, q );
      real64 temp1 = ( m_massDamping * m_newmarkGamma/( m_newmarkBeta * m_dt )
                       + 1.0 / ( m_newmarkBeta * m_dt * m_dt ) )* integrationFactor;

      constexpr int nsdof = STACK_VARIABLE_TYPE::numDofPerNode;
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
  real64 complete( STACK_VARIABLE_TYPE & stack ) const
  {
    constexpr int nsdof = STACK_VARIABLE_TYPE::numDofPerNode;

    for( localIndex a=0; a<STACK_VARIABLE_TYPE::numNodesPerElem; ++a )
    {
      for( localIndex b=0; b<STACK_VARIABLE_TYPE::numNodesPerElem; ++b )
      {
        for( int i=0; i<nsdof; ++i )
        {
          for( int j=0; j<nsdof; ++j )
          {
            stack.localResidual[ a*nsdof+i ] += m_stiffnessDamping * stack.localJacobian[ a*nsdof+i][ b*nsdof+j ] *
                                                ( stack.vtilde_local[b][j] + m_newmarkGamma/(m_newmarkBeta * m_dt)*(stack.uhat_local[b][j]-stack.uhattilde_local[b][j]) );

            stack.localJacobian[a*nsdof+i][b*nsdof+j] += stack.localJacobian[a][b] * (1.0 + m_stiffnessDamping * m_newmarkGamma / ( m_newmarkBeta * m_dt ) )
                                                         + stack.dRdU_InertiaMassDamping[ a ][ b ];
          }
        }
      }
    }

    for( localIndex a=0; a<STACK_VARIABLE_TYPE::ndof; ++a )
    {
      for( localIndex b=0; b<STACK_VARIABLE_TYPE::ndof; ++b )
      {
        stack.localJacobian[a][b] += stack.localJacobian[a][b] * (1.0 + m_stiffnessDamping * m_newmarkGamma / ( m_newmarkBeta * m_dt ) )
                                     + stack.dRdU_InertiaMassDamping[ a ][ b ];
      }
    }

    return Base::complete( stack );
  }
};

} // namespace SolidMechanicsLagrangianFEMKernels

} // namespace geosx
