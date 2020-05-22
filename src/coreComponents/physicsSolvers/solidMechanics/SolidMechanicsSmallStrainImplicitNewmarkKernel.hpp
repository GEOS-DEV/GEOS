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

class ImplicitNewmark
{
public:
  using Base = QuasiStatic;

  struct Parameters : public Base::Parameters
  {
    Parameters( real64 const inputGravityVector[3],
                real64 const inputNewmarkGamma,
                real64 const inputNewmarkBeta,
                real64 const inputMassDamping,
                real64 const inputStiffnessDamping,
                real64 const inputDt ):
      Base::Parameters( inputGravityVector ),
      newmarkGamma( inputNewmarkGamma ),
      newmarkBeta( inputNewmarkBeta ),
      massDamping( inputMassDamping ),
      stiffnessDamping( inputStiffnessDamping ),
      dt( inputDt )
    {}

    real64 const newmarkGamma;
    real64 const newmarkBeta;
    real64 const massDamping;
    real64 const stiffnessDamping;
    real64 const dt;
  };

  template< int NUM_NODES_PER_ELEM, int NUM_DOF_PER_NODE >
  struct StackVariables : Base::StackVariables< NUM_NODES_PER_ELEM, NUM_DOF_PER_NODE >
  {
public:
    using StackVariablesBase = Base::StackVariables< NUM_NODES_PER_ELEM, NUM_DOF_PER_NODE >;

    using StackVariablesBase::numNodesPerElem;
    using StackVariablesBase::numDofPerNode;
    using StackVariablesBase::ndof;

//      GEOSX_HOST_DEVICE
    StackVariables():
      StackVariablesBase(),
      dRdU_InertiaMassDamping{ {0.0} },
      vtilde_local{ { 0.0, 0.0, 0.0} },
      uhattilde_local{ { 0.0, 0.0, 0.0} }
    {}

    real64 dRdU_InertiaMassDamping[ ndof ][ ndof ];

    R1Tensor vtilde_local[NUM_NODES_PER_ELEM];
    R1Tensor uhattilde_local[NUM_NODES_PER_ELEM];
  };

  template< typename SUBREGION_TYPE,
            typename CONSTITUTIVE_TYPE,
            int NUM_NODES_PER_ELEM,
            int >
  class Components : public Base::Components< SUBREGION_TYPE, CONSTITUTIVE_TYPE, NUM_NODES_PER_ELEM, NUM_NODES_PER_ELEM >
  {
public:

    using ComponentsBase = Base::Components< SUBREGION_TYPE, CONSTITUTIVE_TYPE, NUM_NODES_PER_ELEM, NUM_NODES_PER_ELEM >;
    using ComponentsBase::m_dofNumber;
    using ComponentsBase::m_matrix;
    using ComponentsBase::m_rhs;
    using ComponentsBase::elemsToNodes;
    using ComponentsBase::constitutiveUpdate;
    using ComponentsBase::m_disp;
    using ComponentsBase::m_uhat;
    using ComponentsBase::dNdX;
    using ComponentsBase::detJ;
    using ComponentsBase::m_finiteElementSpace;



    Components( arrayView1d< globalIndex const > const & inputDofNumber,
                ParallelMatrix & inputMatrix,
                ParallelVector & inputRhs,
                NodeManager const & nodeManager,
                SUBREGION_TYPE const & elementSubRegion,
                FiniteElementBase const * const finiteElementSpace,
                CONSTITUTIVE_TYPE & constitutiveModel,
                Parameters const & parameters ):
      ComponentsBase( inputDofNumber,
                      inputMatrix,
                      inputRhs,
                      nodeManager,
                      elementSubRegion,
                      finiteElementSpace,
                      constitutiveModel,
                      parameters ),
      m_vtilde( nodeManager.totalDisplacement()),
      m_uhattilde( nodeManager.totalDisplacement()),
      m_density( constitutiveModel.getDensity()),
      m_newmarkGamma( parameters.newmarkGamma ),
      m_newmarkBeta( parameters.newmarkBeta ),
      m_massDamping( parameters.massDamping ),
      m_stiffnessDamping( parameters.stiffnessDamping ),
      m_dt( parameters.dt )
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

      ComponentsBase::quadraturePointJacobianContribution( k, q, stack, [&]( localIndex const a, localIndex const b )
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

      return ComponentsBase::complete( stack );
    }
  };

};

} // namespace SolidMechanicsLagrangianFEMKernels

} // namespace geosx
