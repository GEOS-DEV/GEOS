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

#ifndef GEOSX_PHYSICSSOLVERS_SIMPLEPDE_ELASTICSTATICFEMKERNELS_HPP_
#define GEOSX_PHYSICSSOLVERS_SIMPLEPDE_ELASTICSTATICFEMKERNELS_HPP_

#include "finiteElement/kernelInterface/ImplicitKernelBase.hpp"
#include "finiteElement/kernelInterface/SparsityKernelBase.hpp"

namespace geosx
{

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class ElasticStaticFEMKernels :
  public finiteElement::ImplicitKernelBase< SUBREGION_TYPE,
                                            CONSTITUTIVE_TYPE,
                                            FE_TYPE,
                                            3,
                                            3 >
{
public:
  using Base = finiteElement::ImplicitKernelBase< SUBREGION_TYPE,
                                                  CONSTITUTIVE_TYPE,
                                                  FE_TYPE,
                                                  3,
                                                  3 >;

  using Base::m_dofNumber;
  using Base::m_dofRankOffset;
  using Base::m_matrix;
  using Base::m_rhs;
  using Base::m_elemsToNodes;
  using Base::m_finiteElementSpace;
  using Base::m_constitutiveUpdate;

  static constexpr int numNodesPerElem = Base::numTestSupportPointsPerElem;

  ElasticStaticFEMKernels( NodeManager const & nodeManager,
                           EdgeManager const & edgeManager,
                           FaceManager const & faceManager,
                           SUBREGION_TYPE const & elementSubRegion,
                           FE_TYPE const & finiteElementSpace,
                           CONSTITUTIVE_TYPE * const inputConstitutiveType,
                           arrayView1d< globalIndex const > const & inputDofNumber,
                           globalIndex const rankOffset,
                           CRSMatrixView< real64, globalIndex const > const & inputMatrix,
                           arrayView1d< real64 > const & inputRhs):
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
    m_X( nodeManager.referencePosition()),
    m_deltaDispl( nodeManager.incrementalDisplacement())
  {}

  struct StackVariables : public Base::StackVariables
  {
public:
    GEOSX_HOST_DEVICE
    StackVariables():
      Base::StackVariables(),
      xLocal(),
      deltaDispl_local()
    {}

#if !defined(CALC_FEM_SHAPE_IN_KERNEL)
    // Dummy
    real64 xLocal;
#else
    // Stack storage for element local the nodal positions.
    real64 xLocal[ numNodesPerElem ][ 3 ];
#endif

    // Stack storage for the element local nodal incremental displacement
    real64 deltaDispl_local[ numNodesPerElem ][ 3 ];
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

      for( int i=0; i<3; ++i )
      {
        stack.deltaDispl_local[ a ][ i ] = m_deltaDispl[ localNodeIndex ][ i ];
        stack.localRowDofIndex[ a * 3 + i ] = m_dofNumber[ localNodeIndex ] + i;
        stack.localColDofIndex[ a * 3 + i ] = m_dofNumber[ localNodeIndex ] + i;
      }
    }
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack) const
  {
    real64 dNdX[ numNodesPerElem ][ 3 ];
    real64 const detJ = m_finiteElementSpace.template getGradN< FE_TYPE >( k, q, stack.xLocal, dNdX );

    // Upper part of the local symmetric Jacobian matrix.
    typename CONSTITUTIVE_TYPE::KernelWrapper::DiscretizationOps stiffnessHelper;
    m_constitutiveUpdate.setDiscretizationOps( k, q, stiffnessHelper );
    stiffnessHelper.template upperBTDB< numNodesPerElem >( dNdX, -detJ, stack.localJacobian );

    // Contribution of stress to residual.
    real64 strainInc[ 6 ] = { 0 };
    real64 stress[ 6 ];

    FE_TYPE::symmetricGradient( dNdX, stack.deltaDispl_local, strainInc );  // Compute local strain increment
    m_constitutiveUpdate.SmallStrain( k, q, strainInc );                // Compute local stress increment
    m_constitutiveUpdate.getStress( k, q, stress );                     // Get updated local stress                              

    for( localIndex a=0; a<numNodesPerElem; ++a )
    {
      stack.localResidual[ 3 * a ]     -= ( stress[ 0 ] * dNdX[ a ][ 0 ] + stress[ 5 ] * dNdX[ a ][ 1 ] + stress[ 4 ] * dNdX[ a ][ 2 ] ) * detJ;
      stack.localResidual[ 3 * a + 1 ] -= ( stress[ 5 ] * dNdX[ a ][ 0 ] + stress[ 1 ] * dNdX[ a ][ 1 ] + stress[ 3 ] * dNdX[ a ][ 2 ] ) * detJ;
      stack.localResidual[ 3 * a + 2 ] -= ( stress[ 4 ] * dNdX[ a ][ 0 ] + stress[ 3 ] * dNdX[ a ][ 1 ] + stress[ 2 ] * dNdX[ a ][ 2 ] ) * detJ;
    }
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  real64 complete( localIndex const k,
                   StackVariables & stack ) const
  {
    GEOSX_UNUSED_VAR( k );
    real64 maxForce = 0;

    CONSTITUTIVE_TYPE::KernelWrapper::DiscretizationOps::template fillLowerBTDB< numNodesPerElem >( stack.localJacobian );

    for( localIndex a=0; a<numNodesPerElem; ++a )
    {
      for( int dim = 0; dim < 3; ++dim )
      {
        localIndex const dof = LvArray::integerConversion< localIndex >( stack.localRowDofIndex[ 3 * a + dim ] - m_dofRankOffset );
        if( dof < 0 || dof >= m_matrix.numRows() ) continue;

        m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                                stack.localRowDofIndex,
                                                                                stack.localJacobian[ 3 * a + dim ],
                                                                                numNodesPerElem * 3 );

        RAJA::atomicAdd< parallelDeviceAtomic >( &m_rhs[ dof ], stack.localResidual[ 3 * a + dim ] );

        maxForce = fmax( maxForce, fabs( stack.localResidual[ 3 * a + dim ] ) );
      }
    }

    return maxForce;
  }

protected:

  // The array containing the nodal position array.
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const m_X;

  // The rank-global incremental displacement array.
  arrayView2d< real64 const, nodes::INCR_DISPLACEMENT_USD > const m_deltaDispl;
};

} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_SIMPLEPDE_ELASTICSTATICFEMKERNELS_HPP_ */
