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
 * @file TeamLaplaceFEMKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_SIMPLEPDE_TEAMLAPLACEFEMKERNELS_HPP_
#define GEOSX_PHYSICSSOLVERS_SIMPLEPDE_TEAMLAPLACEFEMKERNELS_HPP_

#define GEOSX_DISPATCH_VEM /// enables VEM in FiniteElementDispatch

#include "finiteElement/TeamKernelInterface/TeamKernelInterface.hpp"

#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutivePassThru.hpp"
#include "finiteElement/FiniteElementDispatch.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

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
 * ### TeamLaplaceFEMKernel Description
 * Implements the TeamKernelBase interface functions required for solving Laplace's
 * equation using the finite element kernel application functions such as
 * geosx::finiteElement::RegionBasedKernelApplication.
 *
 * In this implementation, the template parameter @p NUM_NODES_PER_ELEM is used
 * in place of both @p NUM_TEST_SUPPORT_POINTS_PER_ELEM and
 * @p NUM_TRIAL_SUPPORT_POINTS_PER_ELEM, which are assumed to be equal. This
 * results in the @p UNUSED template parameter as only the NUM_NODES_PER_ELEM
 * is passed to the TeamKernelBase template to form the base class.
 *
 * Additionally, the number of degrees of freedom per support point for both
 * the test and trial spaces are specified as `1` when specifying the base
 * class.
 */
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class TeamLaplaceFEMKernel :
  public finiteElement::TeamKernelBase< SUBREGION_TYPE,
                                        CONSTITUTIVE_TYPE,
                                        FE_TYPE,
                                        1,
                                        1 >
{
public:
  /// An alias for the base class.
  using Base = finiteElement::TeamKernelBase< SUBREGION_TYPE,
                                              CONSTITUTIVE_TYPE,
                                              FE_TYPE,
                                              1,
                                              1 >;

  /**
   * @brief Constructor
   * @copydoc geosx::finiteElement::TeamKernelBase::TeamKernelBase
   * @param fieldName The name of the primary field
   *                  (i.e. Temperature, Pressure, etc.)
   */
  TeamLaplaceFEMKernel( NodeManager const & nodeManager,
                        EdgeManager const & edgeManager,
                        FaceManager const & faceManager,
                        localIndex const targetRegionIndex,
                        SUBREGION_TYPE const & elementSubRegion,
                        FE_TYPE const & finiteElementSpace,
                        CONSTITUTIVE_TYPE & inputConstitutiveType, // end of default args
                        arrayView1d< real64 const > const inputSrc,
                        arrayView1d< real64 > const inputDst ):
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

  TeamLaplaceFEMKernel( NodeManager const & nodeManager,
                        EdgeManager const & edgeManager,
                        FaceManager const & faceManager,
                        localIndex const targetRegionIndex,
                        SUBREGION_TYPE const & elementSubRegion,
                        FE_TYPE const & finiteElementSpace,
                        CONSTITUTIVE_TYPE & inputConstitutiveType, // end of default args
                        arrayView1d< real64 > const inputRhs,
                        string const fieldName ):
    Base( elementSubRegion,
          finiteElementSpace,
          inputConstitutiveType ),
    m_X( nodeManager.referencePosition() ),
    m_src( nodeManager.template getReference< array1d< real64 > >( fieldName )),
    m_dst( inputRhs )
  {
    GEOSX_UNUSED_VAR( edgeManager );
    GEOSX_UNUSED_VAR( faceManager );
    GEOSX_UNUSED_VAR( targetRegionIndex );
  }

  template< typename POLICY, // ignored
            typename KERNEL_TYPE >
  static
  real64
  kernelLaunch( localIndex const numElems,
                KERNEL_TYPE const & fields )
  {
    GEOSX_MARK_FUNCTION;

    fields.m_X.move( parallelDeviceMemorySpace, false );
    fields.m_src.move( parallelDeviceMemorySpace, false );
    fields.m_dst.move( parallelDeviceMemorySpace, true );
    fields.m_elemsToNodes.move( parallelDeviceMemorySpace, false );

    // FE Config
    constexpr localIndex num_dofs_1d = 2;
    constexpr localIndex num_dofs_mesh_1d = 2;
    constexpr localIndex num_quads_1d = 2;

    // Kernel Config
    constexpr ThreadingModel threading_model = ThreadingModel::Serial;
    constexpr localIndex num_threads_1d = num_quads_1d;
    constexpr localIndex batch_size = 32;

    using KernelConf = finiteElement::KernelConfiguration< threading_model, num_threads_1d, batch_size >;
    using Stack = finiteElement::KernelContext< KernelConf >;

    finiteElement::forallElements<KernelConf>( numElems, [=] GEOSX_HOST_DEVICE ( Stack & stack )
    {
      typename Stack::template Tensor< real64, num_dofs_1d, num_dofs_1d, num_dofs_1d, 3 > nodes;
      typename Stack::template Tensor< real64, num_dofs_1d, num_dofs_1d, num_dofs_1d > dofs_in;
      readField( stack, fields.m_elemsToNodes, fields.m_X, nodes );
      readField( stack, fields.m_elemsToNodes, fields.m_src, dofs_in );

      typename Stack::template Basis< num_dofs_1d, num_quads_1d > basis( stack );
      /// Computation of the Jacobians
      typename Stack::template Tensor< real64, num_dofs_mesh_1d, num_dofs_mesh_1d, num_dofs_mesh_1d, 3, 3 > Jac;
      interpolateGradientAtQuadraturePoints( stack, basis, nodes, Jac );

      /// Computation of the Gradient of the solution field
      typename Stack::template Tensor< real64, num_dofs_1d, num_dofs_1d, num_dofs_1d, 3 > Gu;
      interpolateGradientAtQuadraturePoints( stack, basis, dofs_in, Gu );

      /// QFunction
      // forallQuadratureIndices( stack, [&]( auto const & quad_index )
      // {
      //   /// QFunction for Laplace operator
      //   constexpr int dim = Stack::dim;

      //   // Load q-local gradients
      //   real64 grad_ref_u[ dim ];
      //   qLocalLoad( quad_index, Gu, grad_ref_u );

      //   // load q-local jacobian
      //   real64 J[ dim ][ dim ];
      //   qLocalLoad( quad_index, Jac, J );

      //   // Compute D_q = w_q * det(J_q) * J_q^-1 * J_q^-T = w_q / det(J_q) * adj(J_q) * adj(J_q)^T
      //   real64 const weight = 1.0; // stack.weights( quad_index ); // FIXME

      //   real64 Jinv[ dim ][ dim ];
      //   real64 const detJ = computeInverseAndDeterminant( J, Jinv );

      //   real64 grad_phys_u[ dim ];
      //   computePhysicalGradient( Jinv, grad_ref_u, grad_phys_u );

      //   real64 Du[ dim ];
      //   computeReferenceGradient( Jinv, grad_phys_u, Du );

      //   applyQuadratureWeights( weight, detJ, Du );

      //   qLocalWrite( quad_index, Du, Gu );
      // } );

      /// Application of the test functions
      typename Stack::template Tensor< real64, num_dofs_1d, num_dofs_1d, num_dofs_1d > dofs_out;
      applyGradientTestFunctions( stack, basis, Gu, dofs_out );

      writeAddField( stack, fields.m_elemsToNodes, dofs_out, fields.m_dst );

    } );
    return 0.0;
  }

  /// The array containing the nodal position array.
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const m_X;

  arrayView1d< real64 const > const m_src;
  arrayView1d< real64 > const m_dst;
};

/// The factory used to construct a TeamLaplaceFEMKernel.
using TeamLaplaceFEMKernelFactory = finiteElement::KernelFactory< TeamLaplaceFEMKernel,
                                                                  arrayView1d< real64 > const,
                                                                  string const >;
using TeamLaplaceFEMKernelFactory2 = finiteElement::KernelFactory< TeamLaplaceFEMKernel,
                                                                  arrayView1d< real64 const > const,
                                                                  arrayView1d< real64 > const >;

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class TeamLaplaceFEMDiagonalKernel :
  public finiteElement::TeamKernelBase< SUBREGION_TYPE,
                                        CONSTITUTIVE_TYPE,
                                        FE_TYPE,
                                        1,
                                        1 >
{
public:
  /// An alias for the base class.
  using Base = finiteElement::TeamKernelBase< SUBREGION_TYPE,
                                              CONSTITUTIVE_TYPE,
                                              FE_TYPE,
                                              1,
                                              1 >;

  /**
   * @brief Constructor
   * @copydoc geosx::finiteElement::TeamKernelBase::TeamKernelBase
   * @param fieldName The name of the primary field
   *                  (i.e. Temperature, Pressure, etc.)
   */
  TeamLaplaceFEMDiagonalKernel( NodeManager const & nodeManager,
                                EdgeManager const & edgeManager,
                                FaceManager const & faceManager,
                                localIndex const targetRegionIndex,
                                SUBREGION_TYPE const & elementSubRegion,
                                FE_TYPE const & finiteElementSpace,
                                CONSTITUTIVE_TYPE & inputConstitutiveType, // end of default args
                                arrayView1d< real64 > const inputDiag ):
    Base( elementSubRegion,
          finiteElementSpace,
          inputConstitutiveType ),
    m_X( nodeManager.referencePosition() ),
    m_diag( inputDiag )
  {
    GEOSX_UNUSED_VAR( edgeManager );
    GEOSX_UNUSED_VAR( faceManager );
    GEOSX_UNUSED_VAR( targetRegionIndex );
  }

  template< typename POLICY, // ignored
            typename KERNEL_TYPE >
  static
  real64
  kernelLaunch( localIndex const numElems,
                KERNEL_TYPE const & fields )
  {
    GEOSX_MARK_FUNCTION;

    fields.m_X.move( parallelDeviceMemorySpace, false );
    fields.m_diag.move( parallelDeviceMemorySpace, true );
    fields.m_elemsToNodes.move( parallelDeviceMemorySpace, false );

    // FE Config
    constexpr localIndex num_dofs_1d = 2;
    constexpr localIndex num_dofs_mesh_1d = 2;
    constexpr localIndex num_quads_1d = 2;

    // Kernel Config
    constexpr ThreadingModel threading_model = ThreadingModel::Serial;
    constexpr localIndex num_threads_1d = num_quads_1d;
    constexpr localIndex batch_size = 32;

    using KernelConf = finiteElement::KernelConfiguration< threading_model, num_threads_1d, batch_size >;
    using Stack = finiteElement::KernelContext< KernelConf >;

    finiteElement::forallElements<KernelConf>( numElems, [=] GEOSX_HOST_DEVICE ( Stack & stack )
    {
      // typename Stack::template Tensor< real64, num_dofs_mesh_1d, num_dofs_mesh_1d, num_dofs_mesh_1d, 3 > nodes;
      // readField( stack, fields.m_elemsToNodes, fields.m_X, nodes );

      // typename Stack::template Basis< num_dofs_mesh_1d, num_quads_1d > basis( stack );
      // /// Computation of the Jacobians
      // typename Stack::template Tensor< real64, num_quads_1d, num_quads_1d, num_quads_1d, 3, 3 > Jac;
      // interpolateGradientAtQuadraturePoints( stack, basis, nodes, Jac );

      // /// QFunction
      // forallQuadratureIndices( stack, [&]( auto const & quad_index )
      // {
      //   /// QFunction for Laplace operator
      //   constexpr int dim = Stack::dim;

      //   // load q-local jacobian
      //   real64 J[ dim ][ dim ];
      //   qLocalLoad( quad_index, Jac, J );

      //   // Compute D_q = w_q * det(J_q) * J_q^-1 * J_q^-T = w_q / det(J_q) * adj(J_q) * adj(J_q)^T
      //   real64 const weight = 1.0; // stack.weights( quad_index ); // FIXME

      //   real64 const detJ = determinant( J );
      //   real64 const detJinv = 1.0 / detJ;

      //   real64 AdjJ[ dim ][ dim ];
      //   adjugate( J, AdjJ );

      //   real64 D[ dim ][ dim ];
      //   D[0][0] = weight * detJinv * (AdjJ[0][0]*AdjJ[0][0] + AdjJ[0][1]*AdjJ[0][1] + AdjJ[0][2]*AdjJ[0][2]);
      //   D[1][0] = weight * detJinv * (AdjJ[0][0]*AdjJ[1][0] + AdjJ[0][1]*AdjJ[1][1] + AdjJ[0][2]*AdjJ[1][2]);
      //   D[2][0] = weight * detJinv * (AdjJ[0][0]*AdjJ[2][0] + AdjJ[0][1]*AdjJ[2][1] + AdjJ[0][2]*AdjJ[2][2]);
      //   D[0][1] = D[1][0];
      //   D[1][1] = weight * detJinv * (AdjJ[1][0]*AdjJ[1][0] + AdjJ[1][1]*AdjJ[1][1] + AdjJ[1][2]*AdjJ[1][2]);
      //   D[2][1] = weight * detJinv * (AdjJ[1][0]*AdjJ[2][0] + AdjJ[1][1]*AdjJ[2][1] + AdjJ[1][2]*AdjJ[2][2]);
      //   D[0][2] = D[2][0];
      //   D[1][2] = D[2][1];
      //   D[2][2] = weight * detJinv * (AdjJ[2][0]*AdjJ[2][0] + AdjJ[2][1]*AdjJ[2][1] + AdjJ[2][2]*AdjJ[2][2]);

      //   qLocalWrite( quad_index, D, Jac );
      // } );

      // // Computation of the local diagonal
      // typename Stack::template Tensor< real64, num_dofs_1d, num_dofs_1d, num_dofs_1d > diag;
      // computeGradGradLocalDiagonal( stack,
      //                               basis.getValuesAtQuadPts(),
      //                               basis.getGradientValuesAtQuadPts(),
      //                               Jac,
      //                               diag );

      // writeAddField( stack, fields.m_elemsToNodes, diag, fields.m_diag );

    } );
    return 0.0;
  }

  /// The array containing the nodal position array.
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const m_X;

  arrayView1d< real64 > const m_diag;
};

using TeamLaplaceFEMDiagonalKernelFactory = finiteElement::KernelFactory< TeamLaplaceFEMDiagonalKernel,
                                                                          arrayView1d< real64 > const >;

} // namesapce geosx

#endif // GEOSX_PHYSICSSOLVERS_SIMPLEPDE_TEAMLAPLACEFEMKERNELS_HPP_
