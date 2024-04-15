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
 * @file Bis_AcousticVTIZhangWaveEquationSEMKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_BISACOUSTICVTIZHANGWAVEEQUATIONSEMKERNEL_HPP_
#define GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_BISACOUSTICVTIZHANGWAVEEQUATIONSEMKERNEL_HPP_

#include "finiteElement/kernelInterface/KernelBase.hpp"
#include "WaveSolverUtils.hpp"
#if !defined( GEOS_USE_HIP )
#include "finiteElement/elementFormulations/Qk_Hexahedron_Lagrange_GaussLobatto.hpp"
#endif
#include "AcousticFields.hpp"

namespace geos
{

using namespace fields;

/// Namespace to contain the acoustic wave kernels.
namespace Bis_acousticVTIZhangWaveEquationSEMKernels
{

/**
 * @brief Implements kernels for solving the acoustic wave equations
 *   explicit central FD method and SEM
 * @copydoc geos::finiteElement::KernelBase
 * @tparam SUBREGION_TYPE The type of subregion that the kernel will act on.
 *
 * ### AcousticWaveEquationSEMKernel Description
 * Implements the KernelBase interface functions required for solving
 * the acoustic wave equations using the
 * "finite element kernel application" functions such as
 * geos::finiteElement::RegionBasedKernelApplication.
 *
 * The number of degrees of freedom per support point for both
 * the test and trial spaces are specified as `1`.
 */


template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class ExplicitBISAcousticVTIZhangSEM : public finiteElement::KernelBase< SUBREGION_TYPE,
                                                              CONSTITUTIVE_TYPE,
                                                              FE_TYPE,
                                                              1,
                                                              1 >
{
public:

  /// Alias for the base class;
  using Base = finiteElement::KernelBase< SUBREGION_TYPE,
                                          CONSTITUTIVE_TYPE,
                                          FE_TYPE,
                                          1,
                                          1 >;

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
   * @copydoc geos::finiteElement::KernelBase::KernelBase
   * @param nodeManager Reference to the NodeManager object.
   * @param edgeManager Reference to the EdgeManager object.
   * @param faceManager Reference to the FaceManager object.
   * @param targetRegionIndex Index of the region the subregion belongs to.
   * @param dt The time interval for the step.
   *   elements to be processed during this kernel launch.
   */
  ExplicitBISAcousticVTIZhangSEM( NodeManager & nodeManager,
                       EdgeManager const & edgeManager,
                       FaceManager const & faceManager,
                       localIndex const targetRegionIndex,
                       SUBREGION_TYPE const & elementSubRegion,
                       FE_TYPE const & finiteElementSpace,
                       CONSTITUTIVE_TYPE & inputConstitutiveType,
                       real64 const dt ):
    Base( elementSubRegion,
          finiteElementSpace,
          inputConstitutiveType ),
    m_nodeCoords( nodeManager.getField< fields::referencePosition32 >() ),
    m_p_n( nodeManager.getField< acousticvtifields::Pressure_p_n >() ),
    m_q_n( nodeManager.getField< acousticvtifields::Pressure_q_n >() ),
    m_stiffnessVector_p( nodeManager.getField< acousticvtifields::StiffnessVector_p >() ),
    m_stiffnessVector_q( nodeManager.getField< acousticvtifields::StiffnessVector_q >() ),    
    m_density( elementSubRegion.template getField< acousticfields::AcousticDensity >() ),
    m_vti_epsilon( elementSubRegion.template getField< acousticvtifields::AcousticEpsilon >() ),
    m_vti_delta( elementSubRegion.template getField< acousticvtifields::AcousticDelta >() ),
    m_dt( dt )
  {
    GEOS_UNUSED_VAR( edgeManager );
    GEOS_UNUSED_VAR( faceManager );
    GEOS_UNUSED_VAR( targetRegionIndex );
  }

  //*****************************************************************************
  /**
   * @copydoc geos::finiteElement::KernelBase::StackVariables
   *
   * ### ExplicitBISAcousticVTIZhangSEM Description
   * Adds a stack arrays for the nodal force, primary displacement variable, etc.
   */
  struct StackVariables : Base::StackVariables
  {
public:
    GEOS_HOST_DEVICE
    StackVariables():
      xLocal(),
      stiffnessVectorLocal_p(),
      stiffnessVectorLocal_q()
    {}

    /// C-array stack storage for element local the nodal positions.
    real64 xLocal[ 8 ][ 3 ];
    real32 stiffnessVectorLocal_p[ numNodesPerElem ]{};
    real32 stiffnessVectorLocal_q[ numNodesPerElem ]{};
    real32 invDensity;
    real32 vti_epsi; // (1 + 2*epsilon)
    real32 vti_sqrtDelta; // sqrt(1 + 2*delta)
  };
  //***************************************************************************


  /**
   * @copydoc geos::finiteElement::KernelBase::setup
   *
   * Copies the primary variable, and position into the local stack array.
   */
  GEOS_HOST_DEVICE
  inline
  void setup( localIndex const k,
              StackVariables & stack ) const
  {
    real32 epsi = std::fabs( m_vti_epsilon[k] );
    real32 delt = std::fabs( m_vti_delta[k] );
    if( std::fabs( epsi ) < 1e-5 )
      epsi = 0;
    if( std::fabs( delt ) < 1e-5 )
      delt = 0;
    if( delt > epsi )
      delt = epsi;
    stack.vti_epsi = (1 + 2 * epsi);
    stack.vti_sqrtDelta = sqrt( 1 + 2 * delt );
    stack.invDensity = 1./m_density[k];
    for( localIndex a=0; a< 8; a++ )
    {
      localIndex const nodeIndex =  m_elemsToNodes( k, FE_TYPE::meshIndexToLinearIndex3D( a ) );
      for( int i=0; i< 3; ++i )
      {
        stack.xLocal[ a ][ i ] = m_nodeCoords[ nodeIndex ][ i ];
      }
    }
  }

  /**
   * @copydoc geos::finiteElement::KernelBase::complete
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  real64 complete( localIndex const k,
                   StackVariables & stack ) const
  {
    for( int i=0; i<numNodesPerElem; i++ )
    {
      RAJA::atomicAdd< parallelDeviceAtomic >( &m_stiffnessVector_p[m_elemsToNodes( k, i )], stack.stiffnessVectorLocal_p[i] );
      RAJA::atomicAdd< parallelDeviceAtomic >( &m_stiffnessVector_q[m_elemsToNodes( k, i )], stack.stiffnessVectorLocal_q[i] );
    }
    return 0;
  }


  /**
   * @copydoc geos::finiteElement::KernelBase::quadraturePointKernel
   *
   * ### ExplicitBISAcousticVTIZhangSEM Description
   * Calculates stiffness vector
   *
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack ) const
  {
    // Pseudo Stiffness xy
    m_finiteElementSpace.template computeStiffnessxyTerm( q, stack.xLocal, [&] ( int i, int j, real64 val )
    {
      real32 const localIncrement_p = -val * stack.invDensity * stack.vti_epsi * m_p_n[m_elemsToNodes( k, j )];
      stack.stiffnessVectorLocal_p[i] += localIncrement_p;
      real32 const localIncrement_q = -val * stack.invDensity * stack.vti_sqrtDelta * m_p_n[m_elemsToNodes( k, j )];
      stack.stiffnessVectorLocal_q[i] += localIncrement_q;
    } );

    // Pseudo-Stiffness z

    m_finiteElementSpace.template computeStiffnesszTerm( q, stack.xLocal, [&] ( int i, int j, real64 val )
    {
      real32 const localIncrement_p = -val * stack.invDensity * stack.vti_sqrtDelta* m_q_n[m_elemsToNodes( k, j )];
      stack.stiffnessVectorLocal_p[i] += localIncrement_p;
      real32 const localIncrement_q = -val * stack.invDensity * m_q_n[m_elemsToNodes( k, j )];
      stack.stiffnessVectorLocal_q[i] += localIncrement_q;
    } );
  }

protected:
  /// The array containing the nodal position array.
  arrayView2d< WaveSolverBase::wsCoordType const, nodes::REFERENCE_POSITION_USD > const m_nodeCoords;

  /// The array containing the nodal pressure array.
  arrayView1d< real32 const > const m_p_n;

  /// The array containing the nodal auxiliary variable array.
  arrayView1d< real32 const > const m_q_n;

  /// The array containing the product of the stiffness matrix and the nodal pressure for the equation in p.
  arrayView1d< real32 > const m_stiffnessVector_p;

  /// The array containing the product of the stiffness matrix and the nodal pressure for the equation in q.
  arrayView1d< real32 > const m_stiffnessVector_q;

  /// The array containing the medium density.
  arrayView1d< real32 const > const m_density;

  /// The array containing the epsilon Thomsen parameter.
  arrayView1d< real32 const > const m_vti_epsilon;

  /// The array containing the delta Thomsen parameter.
  arrayView1d< real32 const > const m_vti_delta;

  /// The time increment for this time integration step.
  real64 const m_dt;
};



/// The factory used to construct a ExplicitAcousticWaveEquation kernel.
using ExplicitBISAcousticVTIZhangSEMFactory = finiteElement::KernelFactory< ExplicitBISAcousticVTIZhangSEM,
                                                                 real64 >;


} // namespace Bis_acousticVTIZhangWaveEquationSEMKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICWAVEEQUATIONSEMKERNEL_HPP_
