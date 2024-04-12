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

/**
 * @file ElasticWaveEquationSEMKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ELASTICWAVEEQUATIONSEMATTENUATIONKERNEL_HPP_
#define GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ELASTICWAVEEQUATIONSEMATTENUATIONKERNEL_HPP_

#include "ElasticWaveEquationSEMKernel.hpp"
#include "WaveSolverUtils.hpp"
#include "ElasticFields.hpp"

namespace geos
{
using namespace fields;
/// Namespace to contain the elastic wave kernels.
namespace elasticWaveEquationSEMKernels
{

/**
 * @brief Implements kernels for solving the elastic wave equations
 *   explicit central FD method and SEM
 * @copydoc geos::finiteElement::KernelBase
 * @tparam SUBREGION_TYPE The type of subregion that the kernel will act on.
 *
 * ### ElasticWaveEquationSEMKernel Description
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
class ExplicitElasticSEMAttenuation : public finiteElement::KernelBase< SUBREGION_TYPE,
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

  /// Number of nodes per element...which is equal to the
  /// numTestSupportPointPerElem and numTrialSupportPointPerElem by definition.
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
   */
  ExplicitElasticSEMAttenuation( NodeManager & nodeManager,
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
    m_elasticKernel( nodeManager, edgeManager, faceManager, targetRegionIndex, elementSubRegion, finiteElementSpace, inputConstitutiveType, dt ), 
    m_anelasticKernel( nodeManager, edgeManager, faceManager, targetRegionIndex, elementSubRegion, finiteElementSpace, inputConstitutiveType, dt ), 
    m_qualityFactorP( elementSubRegion.template getField< elasticfields::ElasticQualityFactorP >() ),
    m_qualityFactorS( elementSubRegion.template getField< elasticfields::ElasticQualityFactorS >() )
  {
  }



  //*****************************************************************************
  /**
   * @copydoc geos::finiteElement::KernelBase::StackVariables
   *
   * ### ExplicitElasticSEMAttenuation Description
   * Stack arrays for both attenuated and standard elastic kernels
   */
  struct StackVariables : Base::StackVariables
  {
public:
    GEOS_HOST_DEVICE
    StackVariables():
      elastic(),
      anelastic()
    {}

    // sub-stacks for elastic and anelastic kernels
    typename ExplicitElasticSEM< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::StackVariables elastic;
    typename ExplicitElasticSEM< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE,
                      elasticfields::StiffnessVectorAx, elasticfields::StiffnessVectorAy, elasticfields::StiffnessVectorAz >::StackVariables anelastic;
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
    m_elasticKernel.setup( k, stack.elastic );
    m_anelasticKernel.setup( k, stack.anelastic );
    stack.anelastic.mu =  stack.elastic.mu / m_qualityFactorS[ k ];
    real64 lambdap2mua= (stack.elastic.lambda + 2.0 * stack.elastic.mu ) / m_qualityFactorP[ k ];
    stack.anelastic.lambda = (lambdap2mua - 2.0 * stack.anelastic.mu );
  }

  /**
   * @copydoc geos::finiteElement::KernelBase::complete
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  real64 complete( localIndex const k,
                   StackVariables & stack ) const
  {
    m_elasticKernel.complete( k, stack.elastic );
    m_anelasticKernel.complete( k, stack.anelastic );
    return 0;
  }

  /**
   * @copydoc geos::finiteElement::KernelBase::quadraturePointKernel
   *
   * ### ExplicitElasticSEM Description
   * Calculates stiffness vector
   *
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack ) const
  {
    m_elasticKernel.quadraturePointKernel( k, q, stack.elastic );
    m_anelasticKernel.quadraturePointKernel( k, q, stack.anelastic );
  }


protected:

  /// The kernels for the attenuated and standard stiffness vectors
  ExplicitElasticSEM< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE > m_elasticKernel; 
  ExplicitElasticSEM< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE,
                      elasticfields::StiffnessVectorAx, elasticfields::StiffnessVectorAy, elasticfields::StiffnessVectorAz > m_anelasticKernel; 

  /// The array containing the P-wave attenuation quality factor
  arrayView1d< real32 const > const m_qualityFactorP;

  /// The array containing the S-wave attenuation quality factor
  arrayView1d< real32 const > const m_qualityFactorS;

};


/// The factory used to construct a ExplicitAcousticWaveEquation kernel.
using ExplicitElasticSEMAttenuationFactory = finiteElement::KernelFactory< ExplicitElasticSEMAttenuation,
                                                                          real64 >;

} // namespace ElasticWaveEquationSEMKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ELASTICWAVEEQUATIONSEMATTENUATIONKERNEL_HPP_
