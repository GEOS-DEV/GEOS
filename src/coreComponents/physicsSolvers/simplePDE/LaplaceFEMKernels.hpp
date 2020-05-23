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
 * @file SolidMechanicsLagrangianFEMKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_SIMPLEPDE_LAPLACE_KERNELS_FEM_HPP_
#define GEOSX_PHYSICSSOLVERS_SIMPLEPDE_LAPLACE_KERNELS_FEM_HPP_

#include "finiteElement/kernelInterface/ImplicitKernelBase.hpp"

namespace geosx
{
//***************************************************************************
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          int NUM_NODES_PER_ELEM,
          int >
class LaplaceFEMKernel :
  public finiteElement::ImplicitKernelBase< SUBREGION_TYPE,
                                            CONSTITUTIVE_TYPE,
                                            NUM_NODES_PER_ELEM,
                                            NUM_NODES_PER_ELEM,
                                            1,
                                            1 >
{
public:
  using Base = finiteElement::ImplicitKernelBase< SUBREGION_TYPE,
                                                  CONSTITUTIVE_TYPE,
                                                  NUM_NODES_PER_ELEM,
                                                  NUM_NODES_PER_ELEM,
                                                  1,
                                                  1 >;

  static constexpr int numNodesPerElem = NUM_NODES_PER_ELEM;
  static constexpr int numTestDofPerSP = 1;
  static constexpr int numTrialDofPerSP = 1;



  //***************************************************************************
  struct StackVariables : Base::StackVariables
  {
public:

    using Base::StackVariables::numRows;
    using Base::StackVariables::numCols;

    GEOSX_HOST_DEVICE
    StackVariables():
      Base::StackVariables(),
            primaryField_local{ 0.0 }
    {}

    real64 primaryField_local[numNodesPerElem];
  };



  using Base::m_dofNumber;
  using Base::m_matrix;
  using Base::m_rhs;
  using Base::elemsToNodes;
  using Base::elemGhostRank;
  using Base::constitutiveUpdate;
  using Base::m_finiteElementSpace;
  using Base::Launch;

  LaplaceFEMKernel( NodeManager const & nodeManager,
                    EdgeManager const & edgeManager,
                    FaceManager const & faceManager,
                    SUBREGION_TYPE const & elementSubRegion,
                    FiniteElementBase const * const finiteElementSpace,
                    CONSTITUTIVE_TYPE * const inputConstitutiveType,
                    arrayView1d< globalIndex const > const & inputDofNumber,
                    ParallelMatrix & inputMatrix,
                    ParallelVector & inputRhs,
                    string const & fieldName ):
    Base( nodeManager,
          edgeManager,
          faceManager,
          elementSubRegion,
          finiteElementSpace,
          inputConstitutiveType,
          inputDofNumber,
          inputMatrix,
          inputRhs ),
    m_primaryField( nodeManager.template getReference< array1d< real64 > >( fieldName )),
    dNdX( elementSubRegion.template getReference< array3d< R1Tensor > >( dataRepository::keys::dNdX )),
    detJ( elementSubRegion.template getReference< array2d< real64 > >( dataRepository::keys::detJ ) )  //,
  {}

  arrayView1d< real64 const > const m_primaryField;
  arrayView3d< R1Tensor const > const dNdX;
  arrayView2d< real64 const > const detJ;


  template< typename STACK_VARIABLE_TYPE >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void setup( localIndex const k,
              STACK_VARIABLE_TYPE & stack ) const
  {
    for( localIndex a=0; a<NUM_NODES_PER_ELEM; ++a )
    {
      localIndex const localNodeIndex = elemsToNodes( k, a );

      stack.primaryField_local[ a ] = m_primaryField[ localNodeIndex ];
      stack.localRowDofIndex[a] = m_dofNumber[localNodeIndex];
      stack.localColDofIndex[a] = m_dofNumber[localNodeIndex];
    }
  }

  template< typename STACK_VARIABLE_TYPE,
            typename DYNAMICS_LAMBDA = std::function< void( localIndex, localIndex) > >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void quadraturePointJacobianContribution( localIndex const k,
                                            localIndex const q,
                                            STACK_VARIABLE_TYPE & stack ) const
  {
    for( localIndex a=0; a<NUM_NODES_PER_ELEM; ++a )
    {
      for( localIndex b=0; b<NUM_NODES_PER_ELEM; ++b )
      {
        stack.localJacobian[ a ][ b ] += Dot( dNdX( k, q, a ), dNdX( k, q, b ) ) * detJ( k, q );
      }
    }
  }

  template< typename STACK_VARIABLE_TYPE >
  //GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  real64 complete( localIndex const GEOSX_UNUSED_PARAM( k ),
                   STACK_VARIABLE_TYPE & stack ) const
  {
    for( localIndex a = 0; a < NUM_NODES_PER_ELEM; ++a )
    {
      for( localIndex b = 0; b < NUM_NODES_PER_ELEM; ++b )
      {
        stack.localResidual[ a ] += stack.localJacobian[ a ][ b ] * stack.primaryField_local[ b ];
      }
    }

    m_matrix.add( stack.localRowDofIndex,
                  stack.localColDofIndex,
                  &(stack.localJacobian[0][0]),
                  stack.numRows,
                  stack.numCols );

    m_rhs.add( stack.localRowDofIndex,
               stack.localResidual,
               stack.numRows );

    return 1.0;
  }
};


template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          int NUM_NODES_PER_ELEM,
          int >
using LaplaceFEMSparsity = finiteElement::ImplicitKernelBase< SUBREGION_TYPE,
                                                              CONSTITUTIVE_TYPE,
                                                              NUM_NODES_PER_ELEM,
                                                              NUM_NODES_PER_ELEM,
                                                              LaplaceFEMKernel< SUBREGION_TYPE,
                                                                                CONSTITUTIVE_TYPE,
                                                                                NUM_NODES_PER_ELEM,
                                                                                NUM_NODES_PER_ELEM >::numTestDofPerSP,
                                                              LaplaceFEMKernel< SUBREGION_TYPE,
                                                                                CONSTITUTIVE_TYPE,
                                                                                NUM_NODES_PER_ELEM,
                                                                                NUM_NODES_PER_ELEM >::numTestDofPerSP >;


} // namespace geosx
#endif // GEOSX_PHYSICSSOLVERS_SIMPLEPDE_LAPLACE_KERNELS_FEM_HPP_
