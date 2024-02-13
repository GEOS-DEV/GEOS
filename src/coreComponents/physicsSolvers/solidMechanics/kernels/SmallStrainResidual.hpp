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
 * @file SmallStrainResidual.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_KERNELS_SMALLSTRAINRESIDUAL_HPP_
#define GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_KERNELS_SMALLSTRAINRESIDUAL_HPP_

#include "finiteElement/kernelInterface/KernelBase.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsFields.hpp"


namespace geos
{

/// Namespace to contain the solid mechanics kernels.
namespace solidMechanicsLagrangianFEMKernels
{

/**
 * @brief Implements kernels for solving the equations of motion using the
 *   explicit Newmark method under the small strain assumption.
 * @copydoc geosx::finiteElement::KernelBase
 * @tparam SUBREGION_TYPE The type of subregion that the kernel will act on.
 *
 * ### Explicit Small Strain Description
 * Implements the KernelBase interface functions required for explicit time
 * integration of the equations of motion using the
 * "finite element kernel application" functions such as
 * geosx::finiteElement::RegionBasedKernelApplication.
 *
 * In this implementation, the interface for KernelBase is used, but
 * SmallStrainResidual only conforms to the interface set by KernelBase, and
 * does not inherit from KernelBase.
 * The number of degrees of freedom per support point for both
 * the test and trial spaces are specified as `3`.
 */
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class SmallStrainResidual : public finiteElement::KernelBase< SUBREGION_TYPE,
                                                              CONSTITUTIVE_TYPE,
                                                              FE_TYPE,
                                                              3,
                                                              3 >
{
public:

  /// Alias for the base class;
  using Base = finiteElement::KernelBase< SUBREGION_TYPE,
                                          CONSTITUTIVE_TYPE,
                                          FE_TYPE,
                                          3,
                                          3 >;

  /// Maximum number of nodes per element, which is equal to the maxNumTestSupportPointPerElem and
  /// maxNumTrialSupportPointPerElem by definition. When the FE_TYPE is not a Virtual Element, this
  /// will be the actual number of nodes per element.
  static constexpr int numNodesPerElem = Base::maxNumTestSupportPointsPerElem;

  using Base::numDofPerTestSupportPoint;
  using Base::numDofPerTrialSupportPoint;
  using Base::m_elemsToNodes;
  using Base::m_elemGhostRank;
  using Base::m_constitutiveUpdate;
  using Base::m_finiteElementSpace;

//*****************************************************************************
  /**
   * @brief Constructor
   * @copydoc geosx::finiteElement::KernelBase::KernelBase
   * @param nodeManager Reference to the NodeManager object.
   * @param edgeManager Reference to the EdgeManager object.
   * @param faceManager Reference to the FaceManager object.
   * @param targetRegionIndex Index of the region the subregion belongs to.
   * @param dt The time interval for the step.
   * @param elementListName The name of the entry that holds the list of
   *   elements to be processed during this kernel launch.
   */
  SmallStrainResidual( NodeManager & nodeManager,
                       EdgeManager const & edgeManager,
                       FaceManager const & faceManager,
                       localIndex const targetRegionIndex,
                       SUBREGION_TYPE const & elementSubRegion,
                       FE_TYPE const & finiteElementSpace,
                       CONSTITUTIVE_TYPE & inputConstitutiveType,
                       arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const inputSrc,
                       arrayView2d< real64, nodes::TOTAL_DISPLACEMENT_USD > const inputDst,
                       real64 const dt,
                       string const elementListName,
                       int const kernelOptimizationOption );

  //*****************************************************************************
  /**
   * @copydoc geosx::finiteElement::KernelBase::StackVariables
   *
   * ### SmallStrainResidual Description
   * Adds a stack arrays for the nodal force, primary displacement variable, etc.
   */
  struct StackVariables : Base::StackVariables
  {
public:
    GEOS_HOST_DEVICE
    StackVariables():
      fLocal{ { 0.0} },
      varLocal{ {0.0} },
      xLocal()
    {}

    /// C-array stack storage for the element local force
    real64 fLocal[ numNodesPerElem ][ numDofPerTrialSupportPoint ];

    /// C-array stack storage for element local primary variable values.
    real64 varLocal[ numNodesPerElem ][ numDofPerTestSupportPoint ];

// #if !defined(CALC_FEM_SHAPE_IN_KERNEL)
//     /// Dummy
//     int xLocal;
// #else
    /// C-array stack storage for element local the nodal positions.
    real64 xLocal[ numNodesPerElem ][ 3 ];
//#endif
  };
  //***************************************************************************


  /**
   * @copydoc geosx::finiteElement::KernelBase::setup
   *
   * Copies the primary variable, and position into the local stack array.
   */
  GEOS_HOST_DEVICE
  void setup( localIndex const k,
              StackVariables & stack ) const;

  /**
   * @copydoc geosx::finiteElement::KernelBase::quadraturePointKernel
   *
   * ### SmallStrainResidual Description
   * Calculates the shape function derivatives, and the strain tensor. Then
   * calls the constitutive update, and also performs the integration of
   * the stress divergence, rather than using the dedicated component function
   * to allow for some variable reuse.
   */
  GEOS_HOST_DEVICE
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack ) const;

  /**
   * @copydoc geosx::finiteElement::KernelBase::complete
   *
   * ### SmallStrainResidual Description
   * Performs the distribution of the nodal force out to the rank local arrays.
   */
  GEOS_HOST_DEVICE
  real64 complete( localIndex const k,
                   StackVariables const & stack ) const;

  /**
   * @copydoc geosx::finiteElement::KernelBase::complete
   *
   * ### SmallStrainResidual Description
   * Performs the distribution of the nodal force out to the rank local arrays.
   */
  GEOS_HOST_DEVICE
  real64 complete( localIndex const k,
                   real64 const (&fLocal) [ numNodesPerElem ][ numDofPerTestSupportPoint ] ) const;

  GEOS_HOST_DEVICE
  real64 complete( localIndex const ( &elemToNodeMap )[numNodesPerElem],
                   real64 const (&fLocal) [ numNodesPerElem ][ numDofPerTestSupportPoint ] ) const;


  GEOS_HOST_DEVICE
  void setup( localIndex const k,
              real64 ( &xLocal ) [ numNodesPerElem ][ numDofPerTrialSupportPoint ],
              real64 ( &varLocal ) [ numNodesPerElem ][ numDofPerTrialSupportPoint ] ) const;

  GEOS_HOST_DEVICE
  void setup( localIndex const k,
              localIndex ( &elemToNodeMap )[numNodesPerElem],
              real64 ( &xLocal ) [ numNodesPerElem ][ numDofPerTrialSupportPoint ],
              real64 ( &varLocal ) [ numNodesPerElem ][ numDofPerTrialSupportPoint ] ) const;



  /**
   * @copydoc geosx::finiteElement::KernelBase::quadraturePointKernel
   *
   * ### SmallStrainResidual Description
   * Calculates the shape function derivatives, and the strain tensor. Then
   * calls the constitutive update, and also performs the integration of
   * the stress divergence, rather than using the dedicated component function
   * to allow for some variable reuse.
   */
  GEOS_HOST_DEVICE
  void quadraturePointKernel( localIndex const k,
                              localIndex const qa,
                              localIndex const qb,
                              localIndex const qc,
                              real64 const (&xLocal)[ numNodesPerElem ][ numDofPerTrialSupportPoint ],
                              real64 const (&varLocal)[ numNodesPerElem ][ numDofPerTrialSupportPoint ],
                              real64 ( &fLocal ) [ numNodesPerElem ][ numDofPerTrialSupportPoint ] ) const;


  template< int qa, int qb, int qc >
  GEOS_HOST_DEVICE
  void quadraturePointKernel( localIndex const k,
                              real64 const (&xLocal)[ numNodesPerElem ][ numDofPerTrialSupportPoint ],
                              real64 const (&varLocal)[ numNodesPerElem ][ numDofPerTrialSupportPoint ],
                              real64 ( &fLocal ) [ numNodesPerElem ][ numDofPerTrialSupportPoint ] ) const;



  /**
   * @copydoc geosx::finiteElement::KernelBase::kernelLaunch
   *
   * ### SmallStrainResidual Description
   * Copy of the KernelBase::kernelLaunch function without the exclusion of ghost
   * elements.
   */
  template< typename POLICY,
            typename KERNEL_TYPE >
  static real64
  kernelLaunch( localIndex const numElems,
                KERNEL_TYPE const & kernelComponent );


protected:
  /// The array containing the nodal position array.
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const m_X;

  /// The array containing the nodal displacement array.
  arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const m_input;

  /// The array containing the nodal acceleration array, which is used to store
  /// the force.
  arrayView2d< real64, nodes::ACCELERATION_USD > const m_res;

  int const m_kernelOptimizationOption;

  /// The list of elements to process for the kernel launch.


};



/// The factory used to construct a SmallStrainResidual kernel.
using SmallStrainResidualFactory = finiteElement::KernelFactory< SmallStrainResidual,
                                                                 arrayView2d< real64 const, nodes::ACCELERATION_USD > const,
                                                                 arrayView2d< real64, nodes::ACCELERATION_USD > const,
                                                                 real64,
                                                                 string const,
                                                                 int const >;



} // namespace solidMechanicsLagrangianFEMKernels

} // namespace geos

#endif //GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_KERNELS_SMALLSTRAINRESIDUAL_HPP_
