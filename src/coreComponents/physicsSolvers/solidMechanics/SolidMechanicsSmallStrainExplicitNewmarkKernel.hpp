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
 * @file SolidMechanicsSmallStrainExplicitNewmarkKernel.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSSMALLSTRAINEXPLICITNEWMARK_HPP_
#define GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSSMALLSTRAINEXPLICITNEWMARK_HPP_

#include "finiteElement/kernelInterface/KernelBase.hpp"


namespace geosx
{

/// Namespace to contain the solid mechanics kernels.
namespace SolidMechanicsLagrangianFEMKernels
{

/// If UPDATE_STRESS is undef, uses total displacement and stress is not
/// updated at all.
/// If UPDATE_STRESS 1, uses total displacement to and adds material stress
/// state to integral for nodalforces.
/// If UPDATE_STRESS 2 then velocity*dt is used to update material stress state
#define UPDATE_STRESS 2

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
 * ExplicitSmallStrain only conforms to the interface set by KernelBase, and
 * does not inherit from KernelBase.
 * The number of degrees of freedom per support point for both
 * the test and trial spaces are specified as `3`.
 */
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class ExplicitSmallStrain : public finiteElement::KernelBase< SUBREGION_TYPE,
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
  ExplicitSmallStrain( NodeManager & nodeManager,
                       EdgeManager const & edgeManager,
                       FaceManager const & faceManager,
                       localIndex const targetRegionIndex,
                       SUBREGION_TYPE const & elementSubRegion,
                       FE_TYPE const & finiteElementSpace,
                       CONSTITUTIVE_TYPE & inputConstitutiveType,
                       real64 const dt,
                       string const & elementListName ):
    Base( elementSubRegion,
          finiteElementSpace,
          inputConstitutiveType ),
    m_X( nodeManager.referencePosition()),
    m_u( nodeManager.totalDisplacement()),
    m_vel( nodeManager.velocity()),
    m_acc( nodeManager.acceleration() ),
    m_dt( dt ),
    m_elementList( elementSubRegion.template getReference< SortedArray< localIndex > >( elementListName ).toViewConst() )
  {
    GEOSX_UNUSED_VAR( edgeManager );
    GEOSX_UNUSED_VAR( faceManager );
    GEOSX_UNUSED_VAR( targetRegionIndex );
  }

  //*****************************************************************************
  /**
   * @copydoc geosx::finiteElement::KernelBase::StackVariables
   *
   * ### ExplicitSmallStrain Description
   * Adds a stack arrays for the nodal force, primary displacement variable, etc.
   */
  struct StackVariables : Base::StackVariables
  {
public:
    GEOSX_HOST_DEVICE
    StackVariables():
      fLocal{ { 0.0} },
      varLocal{ {0.0} },
      xLocal()
    {}

    /// C-array stack storage for the element local force
    real64 fLocal[ numNodesPerElem ][ numDofPerTrialSupportPoint ];

    /// C-array stack storage for element local primary variable values.
    real64 varLocal[ numNodesPerElem ][ numDofPerTestSupportPoint ];

#if !defined(CALC_FEM_SHAPE_IN_KERNEL)
    /// Dummy
    int xLocal;
#else
    /// C-array stack storage for element local the nodal positions.
    real64 xLocal[ numNodesPerElem ][ 3 ];
#endif
  };
  //***************************************************************************


  /**
   * @copydoc geosx::finiteElement::KernelBase::setup
   *
   * Copies the primary variable, and position into the local stack array.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void setup( localIndex const k,
              StackVariables & stack ) const
  {
    for( localIndex a=0; a< numNodesPerElem; ++a )
    {
      localIndex const nodeIndex = m_elemsToNodes( k, a );
      for( int i=0; i<numDofPerTrialSupportPoint; ++i )
      {
#if defined(CALC_FEM_SHAPE_IN_KERNEL)
        stack.xLocal[ a ][ i ] = m_X[ nodeIndex ][ i ];
#endif

#if UPDATE_STRESS==2
        stack.varLocal[ a ][ i ] = m_vel[ nodeIndex ][ i ] * m_dt;
#else
        stack.varLocal[ a ][ i ] = m_u[ nodeIndex ][ i ];
#endif
      }
    }
  }

  /**
   * @copydoc geosx::finiteElement::KernelBase::quadraturePointKernel
   *
   * ### ExplicitSmallStrain Description
   * Calculates the shape function derivatives, and the strain tensor. Then
   * calls the constitutive update, and also performs the integration of
   * the stress divergence, rather than using the dedicated component function
   * to allow for some variable reuse.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack ) const
  {
//#define USE_JACOBIAN
#if !defined( USE_JACOBIAN )
    real64 dNdX[ numNodesPerElem ][ 3 ];
    real64 const detJ = m_finiteElementSpace.template getGradN< FE_TYPE >( k, q, stack.xLocal, dNdX );
    /// Macro to substitute in the shape function derivatives.
    real64 strain[6] = {0};
    FE_TYPE::symmetricGradient( dNdX, stack.varLocal, strain );

    real64 stressLocal[ 6 ] = {0};
#if UPDATE_STRESS == 2
    m_constitutiveUpdate.smallStrainUpdate_StressOnly( k, q, strain, stressLocal );
#else
    m_constitutiveUpdate.smallStrainNoStateUpdate_StressOnly( k, q, strain, stressLocal );
#endif

    for( localIndex c = 0; c < 6; ++c )
    {
#if UPDATE_STRESS == 2
      stressLocal[ c ] *= -detJ;
#elif UPDATE_STRESS == 1
      stressLocal[ c ] = -( stressLocal[ c ] + m_constitutiveUpdate.m_newStress( k, q, c ) ) * detJ; // TODO: decide on
                                                                                                     // initial stress
                                                                                                     // strategy
#else
      stressLocal[ c ] *= -detJ;
#endif
    }

    FE_TYPE::plusGradNajAij( dNdX, stressLocal, stack.fLocal );

#else
    real64 invJ[3][3];
    real64 const detJ = FE_TYPE::inverseJacobianTransformation( q, stack.xLocal, invJ );

    real64 strain[6] = {0};
    FE_TYPE::symmetricGradient( q, invJ, stack.varLocal, strain );

    real64 stressLocal[ 6 ] = {0};
#if UPDATE_STRESS == 2
    m_constitutiveUpdate.smallStrainUpdate_StressOnly( k, q, strain, stressLocal );
#else
    m_constitutiveUpdate.smallStrainNoStateUpdate_StressOnly( k, q, strain, stressLocal );
#endif

    for( localIndex c = 0; c < 6; ++c )
    {
#if UPDATE_STRESS == 2
      stressLocal[ c ] *= detJ;
#elif UPDATE_STRESS == 1
      stressLocal[ c ] = ( stressLocal[ c ] + m_constitutiveUpdate.m_newStress( k, q, c ) ) * DETJ; // TODO: decide on
                                                                                                    // initial stress
                                                                                                    // strategy
#else
      stressLocal[ c ] *= DETJ;
#endif
    }

    FE_TYPE::plusGradNajAij( q, invJ, stressLocal, stack.fLocal );
#endif
  }

  /**
   * @copydoc geosx::finiteElement::KernelBase::complete
   *
   * ### ExplicitSmallStrain Description
   * Performs the distribution of the nodal force out to the rank local arrays.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  real64 complete( localIndex const k,
                   StackVariables const & stack ) const
  {
    for( localIndex a = 0; a < numNodesPerElem; ++a )
    {
      localIndex const nodeIndex = m_elemsToNodes( k, a );
      for( int b = 0; b < numDofPerTestSupportPoint; ++b )
      {
        RAJA::atomicAdd< parallelDeviceAtomic >( &m_acc( nodeIndex, b ), stack.fLocal[ a ][ b ] );
      }
    }
    return 0;
  }

  /**
   * @copydoc geosx::finiteElement::KernelBase::kernelLaunch
   *
   * ### ExplicitSmallStrain Description
   * Copy of the KernelBase::kernelLaunch function without the exclusion of ghost
   * elements.
   */
  template< typename POLICY,
            typename KERNEL_TYPE >
  static real64
  kernelLaunch( localIndex const numElems,
                KERNEL_TYPE const & kernelComponent )
  {
    GEOSX_MARK_FUNCTION;

    GEOSX_UNUSED_VAR( numElems );

    localIndex const numProcElems = kernelComponent.m_elementList.size();
    forAll< POLICY >( numProcElems,
                      [=] GEOSX_DEVICE ( localIndex const index )
    {
      localIndex const k = kernelComponent.m_elementList[ index ];

      typename KERNEL_TYPE::StackVariables stack;

      kernelComponent.setup( k, stack );
      for( integer q=0; q<KERNEL_TYPE::numQuadraturePointsPerElem; ++q )
      {
        kernelComponent.quadraturePointKernel( k, q, stack );
      }
      kernelComponent.complete( k, stack );
    } );
    return 0;
  }


protected:
  /// The array containing the nodal position array.
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const m_X;

  /// The array containing the nodal displacement array.
  arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const m_u;

  /// The array containing the nodal velocity array.
  arrayView2d< real64 const, nodes::VELOCITY_USD > const m_vel;

  /// The array containing the nodal acceleration array, which is used to store
  /// the force.
  arrayView2d< real64, nodes::ACCELERATION_USD > const m_acc;

  /// The time increment for this time integration step.
  real64 const m_dt;

  /// The list of elements to process for the kernel launch.
  SortedArrayView< localIndex const > const m_elementList;


};
#undef UPDATE_STRESS

/// The factory used to construct a ExplicitSmallStrain kernel.
using ExplicitSmallStrainFactory = finiteElement::KernelFactory< ExplicitSmallStrain,
                                                                 real64,
                                                                 string const & >;

} // namespace SolidMechanicsLagrangianFEMKernels

} // namespace geosx

#endif //GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSSMALLSTRAINEXPLICITNEWMARK_HPP_
