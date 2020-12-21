/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#ifndef GEOSX_PHYSICSSOLVERS_SIMPLEPDE_DIFFUSIONFEMKERNELS_HPP_
#define GEOSX_PHYSICSSOLVERS_SIMPLEPDE_DIFFUSIONFEMKERNELS_HPP_

#include "finiteElement/kernelInterface/ImplicitKernelBase.hpp"
#include "finiteElement/kernelInterface/SparsityKernelBase.hpp"

namespace geosx
{

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class DiffusionFEMKernel :
  public finiteElement::ImplicitKernelBase< SUBREGION_TYPE,
                                            CONSTITUTIVE_TYPE,
                                            FE_TYPE,
                                            1,
                                            1 >
{
public:

  using Base = finiteElement::ImplicitKernelBase< SUBREGION_TYPE,
                                                  CONSTITUTIVE_TYPE,
                                                  FE_TYPE,
                                                  1,
                                                  1 >;

  using Base::numDofPerTestSupportPoint;
  using Base::numDofPerTrialSupportPoint;
  using Base::m_dofNumber;
  using Base::m_dofRankOffset;
  using Base::m_matrix;
  using Base::m_rhs;
  using Base::m_elemsToNodes;
  using Base::m_finiteElementSpace;

  static constexpr int numNodesPerElem = Base::numTestSupportPointsPerElem;

  DiffusionFEMKernel( NodeManager const & nodeManager,
                    EdgeManager const & edgeManager,
                    FaceManager const & faceManager,
                    SUBREGION_TYPE const & elementSubRegion,
                    FE_TYPE const & finiteElementSpace,
                    CONSTITUTIVE_TYPE * const inputConstitutiveType,
                    arrayView1d< globalIndex const > const & inputDofNumber,
                    globalIndex const rankOffset,
                    CRSMatrixView< real64, globalIndex const > const & inputMatrix,
                    arrayView1d< real64 > const & inputRhs,
                    real64 const & dt,
                    string const & fieldName,
                    string const & dFieldName,
                    real64 const & diffusion ):
    Base( nodeManager,
          edgeManager,
          faceManager,
          elementSubRegion,
          finiteElementSpace,
          inputConstitutiveType,
          inputDofNumber,
          rankOffset,
          inputMatrix,
          inputRhs ),
    m_X( nodeManager.referencePosition() ),
    m_primaryField( nodeManager.template getReference< array1d< real64 > >( fieldName ) ),
    m_dPrimaryField( nodeManager.template getReference< array1d< real64 > >( dFieldName ) ),
    m_dt( dt ),
    m_diffusion( diffusion )
  {}

  struct StackVariables : Base::StackVariables
  {
public:

    GEOSX_HOST_DEVICE
    StackVariables():
      Base::StackVariables(),
            xLocal(),
            localPrimaryField{ 0.0 },
            dLocalPrimaryField{ 0.0 }
    {}

#if !defined(CALC_FEM_SHAPE_IN_KERNEL)
    // Dummy
    real64 xLocal;
#else
    // Local nodal positions.
    real64 xLocal[ numNodesPerElem ][ 3 ];
#endif

    // Local primary field variable and it variation at the last time step.
    real64 localPrimaryField[ numNodesPerElem ];
    real64 dLocalPrimaryField[ numNodesPerElem ];
  };

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void setup( localIndex const k,
              StackVariables & stack ) const
  {
    for( localIndex a=0; a<numNodesPerElem; ++a )
    {
      localIndex const localNodeIndex = m_elemsToNodes( k, a );

#if defined(CALC_FEM_SHAPE_IN_KERNEL)
      for( int i=0; i<3; ++i )
      {
        stack.xLocal[ a ][ i ] = m_X[ localNodeIndex ][ i ];
      }
#endif

      stack.localPrimaryField[ a ] = m_primaryField[ localNodeIndex ];
      stack.dLocalPrimaryField[ a ] = m_dPrimaryField[ localNodeIndex ];
      stack.localRowDofIndex[ a ] = m_dofNumber[ localNodeIndex ];
      stack.localColDofIndex[ a ] = m_dofNumber[ localNodeIndex ];
    }
  }

  /**
   * @copydoc geosx::finiteElement::ImplicitKernelBase::quadraturePointKernel
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack ) const
  {
    real64 N[ numNodesPerElem ]; 
    real64 dNdX[ numNodesPerElem ][ 3 ];
    real64 const detJ = m_finiteElementSpace.template getGradN< FE_TYPE >( k, q, stack.xLocal, dNdX );
    FE_TYPE::calcN( q, N );

    for( localIndex a=0; a<numNodesPerElem; ++a )
    {
      for( localIndex b=0; b<numNodesPerElem; ++b )
      {
        real64 const val1 = - 1.0 / m_dt * N[ a ] * N[ b ] * detJ;
        real64 const val2 = - m_diffusion * LvArray::tensorOps::AiBi< 3 >( dNdX[a], dNdX[b] ) * detJ;

        stack.localJacobian[ a ][ b ] += val1 + val2;  
        stack.localResidual[ a ] += val1 * stack.dLocalPrimaryField[ b ] + val2 * stack.localPrimaryField[ b ];
      }
    }
  }

  /**
   * @copydoc geosx::finiteElement::ImplicitKernelBase::complete
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  real64 complete( localIndex const k,
                   StackVariables & stack ) const
  {
    GEOSX_UNUSED_VAR( k );
    real64 maxForce = 0;

    for( int a = 0; a < numNodesPerElem; ++a )
    {
      localIndex const dof = LvArray::integerConversion< localIndex >( stack.localRowDofIndex[ a ] - m_dofRankOffset );
      if( dof < 0 || dof >= m_matrix.numRows() ) continue;
      m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                              stack.localColDofIndex,
                                                                              stack.localJacobian[ a ],
                                                                              numNodesPerElem );

      RAJA::atomicAdd< parallelDeviceAtomic >( &m_rhs[ dof ], stack.localResidual[ a ] );
      maxForce = fmax( maxForce, fabs( stack.localResidual[ a ] ) );
    }

    return maxForce;
  }



protected:
  // Nodal positions.
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const m_X;

  // Primary field (pore pressure or temperature).
  arrayView1d< real64 const > const m_primaryField;

  // Primary field variation at the previous time step.
  arrayView1d< real64 const > const m_dPrimaryField;

  // Time step.
  real64 const m_dt;

  // Diffusion coefficient.
  real64 const m_diffusion;
};

} // namespace geosx

#endif // GEOSX_PHYSICSSOLVERS_SIMPLEPDE_DIFFUSIONFEMKERNELS_HPP_
