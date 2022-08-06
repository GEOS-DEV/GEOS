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
 * @file TeamKernelBase.hpp
 */

#ifndef GEOSX_FINITEELEMENT_TEAMKERNELBASE_HPP_
#define GEOSX_FINITEELEMENT_TEAMKERNELBASE_HPP_

#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutivePassThru.hpp"
#include "finiteElement/FiniteElementDispatch.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "finiteElement/TeamKernelInterface/common.hpp"

namespace geosx
{

/**
 * @namespace finiteElement Contains the finite element implementation.
 */
namespace finiteElement
{

/**
 * @class TeamKernelBase
 * @brief Define the base interface for finite element kernels using block
 *        of threads.
 * @tparam SUBREGION_TYPE The type of subregion that the kernel will act on.
 * @tparam CONSTITUTIVE_TYPE The type of constitutive model present in the
 *                           subregion.
 * @tparam NUM_TEST_SUPPORT_POINTS_PER_ELEM The number of test space support
 *                                          points per element.
 * @tparam NUM_TRIAL_SUPPORT_POINTS_PER_ELEM The number of trial space support
 *                                           points per element.
 * @tparam NUM_DOF_PER_TEST_SP The number of DOF per test support point.
 * @tparam NUM_DOF_PER_TRIAL_SP The number of DOF per trial support point.
 *
 * ### General TeamKernelBase Description
 *
 * TeamKernelBase defines an interface for implementing finite element kernels
 * that will be callable by the family of kernel launch functions. Specific
 * physics kernels may or may not derive from TeamKernelBase, but must follow
 * the same interface in order to be callable from the generic launching
 * functions.
 *
 * The template parameters of TeamKernelBase should be duplicated as part of the
 * interface, EXCEPT for @p NUM_DOF_PER_TEST_SP and @p NUM_DOF_PER_TRIAL_SP.
 * These values should be set internally by the physics solver since each
 * physics discretization will have a constant intrinsic value for these
 * quantities. For example, when solving or the heat equation with scalar
 * temperature as the primary variable at the support point, these will have
 * a value of 1. In contrast, when solving a solid mechanics problem, with
 * vector displacement as the primary variable at the support point, these
 * will have a value of 3. Note that the interface provided by
 * geosx::finiteElement::RegionBasedKernelApplication will construct a
 * kernel assuming only the first 4 template arguments.
 */
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE,
          int NUM_DOF_PER_TEST_SP,
          int NUM_DOF_PER_TRIAL_SP >
class TeamKernelBase
{
public:
  /// Compile time value for the number of test function support points per
  /// element.
  static constexpr int maxNumTestSupportPointsPerElem  = FE_TYPE::maxSupportPoints;

  /// Compile time value for the number of trial function support points per
  /// element.
  static constexpr int maxNumTrialSupportPointsPerElem = FE_TYPE::maxSupportPoints;

  /// Compile time value for the number of degrees of freedom per test function
  /// support point.
  static constexpr int numDofPerTestSupportPoint    = NUM_DOF_PER_TEST_SP;

  /// Compile time value for the number of degrees of freedom per trial
  /// function support point.
  static constexpr int numDofPerTrialSupportPoint   = NUM_DOF_PER_TRIAL_SP;

  /// Compile time value for the number of quadrature points per element.
  static constexpr int numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;

  /**
   * @brief Constructor
   * @param elementSubRegion Reference to the SUBREGION_TYPE(class template
   *                         parameter) object.
   * @param inputConstitutiveType The constitutive object.
   * @param finiteElementSpace Placeholder for the finite element space object,
   *                           which currently doesn't do much.
   */
  TeamKernelBase( SUBREGION_TYPE const & elementSubRegion,
                  FE_TYPE const & finiteElementSpace,
                  CONSTITUTIVE_TYPE & inputConstitutiveType ):
    m_elemsToNodes( elementSubRegion.nodeList().toViewConst() ),
    m_elemGhostRank( elementSubRegion.ghostRank() ),
    m_constitutiveUpdate( inputConstitutiveType.createKernelUpdates() ),
    m_finiteElementSpace( finiteElementSpace )
  {}

  /**
   * @struct StackVariables
   * @brief Kernel variables allocated on the stack.
   *
   * ### ImplicitKernelBase::StackVariables Description
   *
   * Contains variables that will be allocated on the stack of the main kernel.
   * This will typically consist of local arrays to hold data mapped from the
   * global data arrays, and/or local storage for the residual and jacobian
   * contributions.
   */
  struct StackVariables
  {
public:

    /**
     * @brief Constructor
     */
    GEOSX_HOST_DEVICE
    StackVariables( LaunchContext & ctx )
      : ctx(ctx)
    {}

    /// RAJA launch context
    LaunchContext & ctx;

    /// Index of the finite element
    localIndex element_index;

    static constexpr localIndex batch_size = 1;
    static constexpr localIndex num_quads_1d = 2; // TODO
  };

  /**
   * @brief Performs the setup phase for the kernel (common to all elements).
   * @tparam STACK_VARIABLE_TYPE The type of StackVariable that holds the stack
   *                             variables. This is most likely a defined in a
   *                             type that derives from TeamKernelBase.
   * @param stack The StackVariable object that hold the stack variables.
   *
   * ### TeamKernelBase::setup() Description
   *
   * The operations typically found in setup are thing such as the collection
   * of global data into local stack storage.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void kernelSetup( StackVariables & stack ) const
  {
    GEOSX_UNUSED_VAR( stack );
  }

  /**
   * @brief Performs the setup phase for the kernel.
   * @tparam STACK_VARIABLE_TYPE The type of StackVariable that holds the stack
   *                             variables. This is most likely a defined in a
   *                             type that derives from TeamKernelBase.
   * @param k The element index.
   * @param stack The StackVariable object that hold the stack variables.
   *
   * ### TeamKernelBase::setup() Description
   *
   * The operations typically found in setup are thing such as the collection
   * of global data into local stack storage.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void setup( StackVariables & stack,
              localIndex const element_index ) const
  {
    stack.element_index = element_index;
  }

  /**
   * @brief Performs a state update at a quadrature point.
   * @tparam STACK_VARIABLE_TYPE The type of StackVariable that holds the stack
   *                             variables. This is most likely a defined in a
   *                             type that derives from TeamKernelBase.
   * @param k The element index.
   * @param q The quadrature point index.
   * @param stack The StackVariable object that hold the stack variables.
   *
   * ### TeamKernelBase::quadraturePointKernel() Description
   *
   * The operations found here are the mapping from the support points to the
   * quadrature point, calculation of gradients, etc. From this data the
   * state of the constitutive model is updated if required by the physics
   * package.
   */
  template < typename QuadraturePointIndex >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void quadraturePointKernel( StackVariables & stack,
                              QuadraturePointIndex const & quad_index ) const
  {
    GEOSX_UNUSED_VAR( stack );
    GEOSX_UNUSED_VAR( quad_index );
  }

  /**
   * @brief Performs the complete phase for the kernel.
   * @tparam STACK_VARIABLE_TYPE The type of StackVariable that holds the stack
   *                             variables. This is most likely a defined in a
   *                             type that derives from TeamKernelBase.
   * @param k The element index.
   * @param stack The StackVariable object that hold the stack variables.
   * @return The maximum contribution to the residual.
   *
   * ### TeamKernelBase::complete() Description
   *
   * The operations typically found in complete are the mapping of the local
   * Jacobian and Residual into the global Jacobian and Residual.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  real64 complete( StackVariables & stack ) const
  {
    GEOSX_UNUSED_VAR( stack );
    return 0;
  }


  /**
   * @brief Kernel Launcher.
   * @tparam POLICY The RAJA policy to use for the launch.
   * @tparam NUM_QUADRATURE_POINTS The number of quadrature points per element.
   * @tparam KERNEL_TYPE The type of Kernel to execute.
   * @param numElems The number of elements to process in this launch.
   * @param kernelComponent The instantiation of KERNEL_TYPE to execute.
   * @return The maximum residual contribution.
   *
   * This is a generic launching function for all of the finite element kernels
   * that follow the interface set by TeamKernelBase.
   */
  //START_kernelLauncher
  template< typename POLICY,
            typename KERNEL_TYPE >
  static
  real64
  kernelLaunch( localIndex const numElems,
                KERNEL_TYPE const & kernelComponent )
  {
    GEOSX_MARK_FUNCTION;
    using RAJA::RangeSegment;

    // Define a RAJA reduction variable to get the maximum residual contribution.
    // RAJA::ReduceMax< ReducePolicy< POLICY >, real64 > maxResidual( 0 );
    // RAJA::ReduceMax< RAJA::seq_reduce, real64 > maxResidual( 0 );
#if defined(RAJA_ENABLE_CUDA)
    RAJA::ReduceMax< RAJA::cuda_reduce, real64 > maxResidual( 0 );
#else
    RAJA::ReduceMax< RAJA::seq_reduce, real64 > maxResidual( 0 );
#endif

    constexpr localIndex batch_size = KERNEL_TYPE::StackVariables::batch_size;

    const localIndex num_blocks = ( numElems + batch_size - 1 ) / batch_size;
    // const localIndex num_SM = 80; // For V100
    // const localIndex num_blocks = 2 * num_SM;

    constexpr localIndex num_quads_1d = KERNEL_TYPE::StackVariables::num_quads_1d;

    launch< POLICY >
    ( GEOSX_RAJA_DEVICE, Resources( Teams( num_blocks ), Threads( num_quads_1d, num_quads_1d, batch_size ) ),
    [=] GEOSX_HOST_DEVICE ( LaunchContext ctx )
    {
      typename KERNEL_TYPE::StackVariables stack( kernelComponent, ctx );

      kernelComponent.kernelSetup( stack );

      // Each block of threads treats "batch_size" elements.
      loop<team_x>( ctx, RangeSegment( 0, num_blocks ), [&]( const int block_index )
      {
        // We batch elements over the z-thread dimension
        loop<thread_z>( ctx, RangeSegment( 0, batch_size ), [&]( const int thread_index_z )
        {
          const localIndex element_index = block_index * batch_size + thread_index_z;
          if ( element_index >= numElems ) { return; }

          kernelComponent.setup( stack, element_index );
          loop<thread_y>( ctx, RangeSegment( 0, num_quads_1d ), [&] ( localIndex quad_y )
          {
            loop<thread_x>( ctx, RangeSegment( 0, num_quads_1d ), [&] ( localIndex quad_x )
            {
              for( localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++ )
              {
                TensorIndex quad_index { quad_x, quad_y, quad_z };
                kernelComponent.quadraturePointKernel( stack, quad_index );
              }
            } );
          } );
          
          maxResidual.max( kernelComponent.complete( stack ) );
        } );
      } );
    } );
    return maxResidual.get();
  }
  //END_kernelLauncher

  /// The element to nodes map.
  traits::ViewTypeConst< typename SUBREGION_TYPE::NodeMapType::base_type > const m_elemsToNodes;

  /// The element ghost rank array.
  arrayView1d< integer const > const m_elemGhostRank;

  /// The constitutive update object used to update the constitutive state,
  /// and extract constitutive data.
  typename CONSTITUTIVE_TYPE::KernelWrapper const m_constitutiveUpdate;

  /// The finite element space/discretization object for the element type in
  /// the SUBREGION_TYPE.
  FE_TYPE const & m_finiteElementSpace;
};

} // namespace finiteElement
} // namespace geosx



#endif /* GEOSX_FINITEELEMENT_TEAMKERNELBASE_HPP_ */
