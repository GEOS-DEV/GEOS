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
#include "tensor/tensor_types.hpp"
#include "common.hpp"
// #include <cuda.h>

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
  template < typename KernelConfig >
  struct StackVariables: public KernelConfig
  {
    /**
     * @brief Constructor
     */
    GEOSX_HOST_DEVICE
    StackVariables( LaunchContext & ctx )
      : KernelConfig( ctx ), element_index( -1 )
    { }

    /// Index of the finite element
    localIndex element_index;
  };

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
