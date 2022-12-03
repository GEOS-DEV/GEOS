/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */


/**
 * @file TeamSolidMechanicsFEMKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_TEAMSOLIDMECHANICSFEMKERNELS_HPP_
#define GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_TEAMSOLIDMECHANICSFEMKERNELS_HPP_

#define GEOSX_DISPATCH_VEM /// enables VEM in FiniteElementDispatch

#include "finiteElement/kernelInterface/KernelBase.hpp"
#include "finiteElement/TeamKernelInterface/TeamKernelBase.hpp"
#include "finiteElement/TeamKernelInterface/TeamKernelFunctions/TeamKernelFunctions.hpp"
#include "finiteElement/TeamKernelInterface/QuadraturePointKernelFunctions/QuadratureFunctionsHelper.hpp"
#include "finiteElement/TeamKernelInterface/QuadraturePointKernelFunctions/qLocalLoad.hpp"
#include "finiteElement/TeamKernelInterface/QuadraturePointKernelFunctions/qLocalWrite.hpp"
#include "finiteElement/TeamKernelInterface/StackVariables/MeshStackVariables.hpp"
#include "finiteElement/TeamKernelInterface/StackVariables/VectorElementStackVariables.hpp"
#include "finiteElement/TeamKernelInterface/StackVariables/BasisStackVariables.hpp"
#include "finiteElement/TeamKernelInterface/StackVariables/QuadratureWeightsStackVariables.hpp"
#include "finiteElement/TeamKernelInterface/StackVariables/SharedStackVariables.hpp"
#include "finiteElement/TeamKernelInterface/StackVariables/SharedMemStackVariables.hpp"

#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutivePassThru.hpp"
#include "finiteElement/FiniteElementDispatch.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "tensor/tensor_types.hpp"

namespace geosx
{

//*****************************************************************************
/**
 * @brief Implements kernels for solving Laplace's equation.
 * @copydoc geosx::finiteElement::TeamKernelBase
 * @tparam NUM_NODES_PER_ELEM The number of nodes per element for the
 *                            @p SUBREGION_TYPE.
 * @tparam UNUSED An unused parameter since we are assuming that the test and
 *                trial space have the same number of support points.
 *
 * ### TeamSolidMechanicsFEMKernel Description
 * Implements the KernelBase interface functions required for solving the
 * quasi-static equilibrium equations using one of the
 * "finite element kernel application" functions such as
 * geosx::finiteElement::RegionBasedKernelApplication.
 *
 * In this implementation, the template parameter @p NUM_NODES_PER_ELEM is used
 * in place of both @p NUM_TEST_SUPPORT_POINTS_PER_ELEM and
 * @p NUM_TRIAL_SUPPORT_POINTS_PER_ELEM, which are assumed to be equal. This
 * results in the @p UNUSED template parameter as only the NUM_NODES_PER_ELEM
 * is passed to the ImplicitKernelBase template to form the base class.
 *
 * Additionally, the number of degrees of freedom per support point for both
 * the test and trial spaces are specified as `3` when specifying the base
 * class.
 */
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class TeamSolidMechanicsFEMKernel :
  public finiteElement::TeamKernelBase< SUBREGION_TYPE,
                                        CONSTITUTIVE_TYPE,
                                        FE_TYPE,
                                        3,
                                        3 >
{
public:
  /// An alias for the base class.
  using Base = finiteElement::TeamKernelBase< SUBREGION_TYPE,
                                              CONSTITUTIVE_TYPE,
                                              FE_TYPE,
                                              3,
                                              3 >;

  static constexpr int numNodesPerElem = Base::maxNumTestSupportPointsPerElem;

  using Base::numDofPerTestSupportPoint;
  using Base::numDofPerTrialSupportPoint;
  using Base::m_elemsToNodes;
  using Base::m_elemGhostRank;
  using Base::m_constitutiveUpdate;
  using Base::m_finiteElementSpace;

  /**
   * @brief Constructor
   * @copydoc geosx::finiteElement::TeamKernelBase::TeamKernelBase
   * @param fieldName The name of the primary field
   *                  (i.e. Temperature, Pressure, etc.)
   */
  TeamSolidMechanicsFEMKernel( NodeManager const & nodeManager,
                               EdgeManager const & edgeManager,
                               FaceManager const & faceManager,
                               localIndex const targetRegionIndex,
                               SUBREGION_TYPE const & elementSubRegion,
                               FE_TYPE const & finiteElementSpace,
                               CONSTITUTIVE_TYPE & inputConstitutiveType, // end of default args
                               arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const inputSrc,
                               arrayView2d< real64, nodes::TOTAL_DISPLACEMENT_USD > const inputDst ):
    Base( elementSubRegion,
          finiteElementSpace,
          inputConstitutiveType ),
    m_X( nodeManager.referencePosition() ),
    m_src( inputSrc ),
    m_dst( inputDst )
  {
    GEOSX_UNUSED_VAR( edgeManager );
    GEOSX_UNUSED_VAR( faceManager );
    GEOSX_UNUSED_VAR( targetRegionIndex );
  }

  //***************************************************************************
  /**
   * @class StackVariables
   * @copydoc geosx::finiteElement::TeamKernelBase::StackVariables
   *
   * Adds a stack array for the primary field.
   */
  template < typename KernelConfig >
  struct StackVariables : public Base::template StackVariables< KernelConfig >
  {
    static constexpr localIndex dim = 3;
    static constexpr localIndex num_dofs_mesh_1d = 2; // TODO
    static constexpr localIndex num_dofs_1d = 2; // TODO
    static constexpr localIndex num_quads_1d = 2; // TODO

    using Basis = stackVariables::StackBasis<num_dofs_1d,num_quads_1d>;

    /**
     * @brief Constructor
     */
    GEOSX_HOST_DEVICE
    StackVariables( LaunchContext & ctx ):
      Base::template StackVariables<KernelConfig>( ctx )
    {}

    /// Shared memory buffers, using buffers allows to avoid using too much shared memory.
    // static constexpr localIndex buffer_size = num_quads_1d * num_quads_1d * num_quads_1d;
    // static constexpr localIndex num_buffers = 2 * dim;
    // static constexpr localIndex buffer_size = num_quads_1d * num_quads_1d * num_quads_1d * dim * dim;
    // static constexpr localIndex num_buffers = 1;
    // stackVariables::SharedMemBuffers< buffer_size, num_buffers, batch_size > shared_mem;
  };

  template< typename POLICY, // ignored
            typename KERNEL_TYPE >
  static
  real64
  kernelLaunch( localIndex const numElems,
                KERNEL_TYPE const & fields )
  {
    GEOSX_MARK_FUNCTION;

    fields.m_X.move( LvArray::MemorySpace::cuda, false );
    fields.m_src.move( LvArray::MemorySpace::cuda, false );
    fields.m_dst.move( LvArray::MemorySpace::cuda, true );
    fields.m_elemsToNodes.move( LvArray::MemorySpace::cuda, false );

    // FE Config
    constexpr localIndex num_dofs_1d = 2;
    constexpr localIndex num_dofs_mesh_1d = 2;
    constexpr localIndex num_quads_1d = 2;

    // Kernel Config
    constexpr ThreadingModel threading_model = ThreadingModel::Serial;
    constexpr localIndex num_threads_1d = num_quads_1d;
    constexpr localIndex batch_size = 32;

    using KernelConf = finiteElement::KernelConfiguration< threading_model, num_threads_1d, batch_size >;
    using Stack = StackVariables< KernelConf >;

    finiteElement::forallElements< KernelConf >( numElems, fields, [=] GEOSX_HOST_DEVICE ( Stack & stack )
    {
      typename Stack::template Tensor< real64, num_dofs_1d, num_dofs_1d, num_dofs_1d, 3 > dofs_in;
      typename Stack::template Tensor< real64, num_dofs_mesh_1d, num_dofs_mesh_1d, num_dofs_mesh_1d, 3 > nodes;
      readField( stack, fields.m_elemsToNodes, fields.m_X, nodes );
      readField( stack, fields.m_elemsToNodes, fields.m_src, dofs_in );

      typename Stack::Basis basis( stack );
      /// Computation of the Jacobians
      typename Stack::template Tensor< real64, num_dofs_mesh_1d, num_dofs_mesh_1d, num_dofs_mesh_1d, 3, 3 > Jac;
      interpolateGradientAtQuadraturePoints( stack, basis, nodes, Jac );

      /// Computation of the Gradient of the solution field
      typename Stack::template Tensor< real64, num_dofs_1d, num_dofs_1d, num_dofs_1d, 3, 3 > Gu;
      interpolateGradientAtQuadraturePoints( stack, basis, dofs_in, Gu );

      /// QFunction
      forallQuadratureIndices( stack, [&]( auto const & quad_index )
      {
        constexpr localIndex dim = 3;
        // Load q-local gradients
        real64 grad_ref_u[ dim ][ dim ];
        qLocalLoad( quad_index, Gu, grad_ref_u );

        // load q-local jacobian
        real64 J[ dim ][ dim ];
        qLocalLoad( quad_index, Jac, J );

        // Compute D_q u_q = - w_q * det(J_q) * J_q^-1 * sigma ( J_q^-1, grad( u_q ) )
        real64 Jinv[ dim ][ dim ];
        real64 const detJ = computeInverseAndDeterminant( J, Jinv );

        real64 stress[ dim ][ dim ];
        computeStress( stack, fields.m_constitutiveUpdate, 0, Jinv, grad_ref_u, stress );

        // Apply - w_q * det(J_q) * J_q^-1
        real64 D[ dim ][ dim ];
        computeReferenceGradient( Jinv, stress, D );

        real64 const weight = 1.0; //stack.weights( quad_index );
        applyQuadratureWeights( weight, detJ, D );

        qLocalWrite( quad_index, D, Gu );
      } );

      /// Application of the test functions
      typename Stack::template Tensor< real64, num_dofs_1d, num_dofs_1d, num_dofs_1d, 3 > dofs_out;
      applyGradientTestFunctions( stack, basis, Gu, dofs_out );

      writeAddField( stack, fields.m_elemsToNodes, dofs_out, fields.m_dst );

    } );
    return 0.0;
  }


  #if 0
  template< typename POLICY,
            typename KERNEL_TYPE >
  static real64
  kernelLaunch( localIndex const numElems,
                KERNEL_TYPE const & kernelComponent )
  {
    GEOSX_MARK_FUNCTION;
    constexpr int dims=3;

  traits::ViewTypeConst< typename SUBREGION_TYPE::NodeMapType::base_type > const elemsToNodes = kernelComponent.m_elemsToNodes;
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X = kernelComponent.m_X;
  arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const src = kernelComponent.m_src;
  arrayView2d< real64, nodes::TOTAL_DISPLACEMENT_USD > const dst = kernelComponent.m_dst;
  typename CONSTITUTIVE_TYPE::KernelWrapper const & constitutiveUpdate = kernelComponent.m_constitutiveUpdate;


    forAll< parallelDevicePolicy<> >( numElems,
                      [=] GEOSX_DEVICE ( localIndex const k )
    {

      real64 xLocal[numNodesPerElem][dims];
      real64 srcLocal[numNodesPerElem][dims];
      real64 dstLocal[numNodesPerElem][dims] = {};
      localIndex nodeIndices[numNodesPerElem];

      for( localIndex a=0; a<numNodesPerElem; ++a )
      {
        nodeIndices[a] = elemsToNodes( k, a );
      }

      for( localIndex a=0; a<numNodesPerElem; ++a )
      {
        for( int i=0; i<dims; ++i )
        {
          xLocal[ a ][ i ] = X( nodeIndices[a],i);
          srcLocal[ a ][ i ] = src( nodeIndices[a],i);
        }
      }


      for( integer q=0; q<KERNEL_TYPE::numQuadraturePointsPerElem; ++q )
      {
        // real64 dNdX[ numNodesPerElem ][ dims ];
        // real64 const detJ = FE_TYPE::calcGradN( q, xLocal, dNdX );
        // /// Macro to substitute in the shape function derivatives.
        // real64 strain[6] = {0};
        // FE_TYPE::symmetricGradient( dNdX, srcLocal, strain );

        // real64 stressLocal[ 6 ] = {0};
        // constitutiveUpdate.smallStrainNoStateUpdate_StressOnly( k, q, strain, stressLocal );

        // for( localIndex c = 0; c < 6; ++c )
        // {
        //   stressLocal[ c ] *= -detJ;
        // }
        // FE_TYPE::plusGradNajAij( dNdX, stressLocal, dstLocal );
      }
      for( localIndex a = 0; a < numNodesPerElem; ++a )
      {
        for( int i = 0; i < numDofPerTestSupportPoint; ++i )
        {
          dstLocal[a][i] += xLocal[ a ][ i ] ;
        }
      }

      
      for( localIndex a = 0; a < numNodesPerElem; ++a )
      {
//        localIndex const nodeIndex = elemsToNodes( k, a );
        for( int i = 0; i < numDofPerTestSupportPoint; ++i )
        {
          RAJA::atomicAdd< parallelDeviceAtomic >( &(dst[nodeIndices[a]][i]), dstLocal[ a ][ i ] );
        }
      }

    } );
    return 0;
  }
  #endif

  /// The array containing the nodal position array.
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const m_X;
  /// The arrays containing input/output displacement.
  arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const m_src;
  arrayView2d< real64, nodes::TOTAL_DISPLACEMENT_USD > const m_dst;
};

/// The factory used to construct a TeamSolidMechanicsFEMKernel.
using TeamSolidMechanicsFEMKernelFactory = finiteElement::KernelFactory< TeamSolidMechanicsFEMKernel,
                                                                         arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const,
                                                                         arrayView2d< real64, nodes::TOTAL_DISPLACEMENT_USD > const >;

} // namesapce geosx

#endif // GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_TEAMSOLIDMECHANICSFEMKERNELS_HPP_
