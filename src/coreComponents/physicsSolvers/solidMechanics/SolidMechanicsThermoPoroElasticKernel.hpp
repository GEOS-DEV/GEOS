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
 * @file SolidMechanicsThermoPoroElasticKernel.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSTHERMOPOROELASTIC_HPP_
#define GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSTHERMOPOROELASTIC_HPP_

#include "physicsSolvers/solidMechanics/SolidMechanicsSmallStrainQuasiStaticKernel.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseExtrinsicData.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBaseExtrinsicData.hpp"

namespace geosx
{

namespace solidMechanicsLagrangianFEMKernels
{

/**
 * @brief Implements kernels for solving quasi-static equilibrium.
 * @copydoc QuasiStatic
 * @tparam NUM_NODES_PER_ELEM The number of nodes per element for the
 *                            @p SUBREGION_TYPE.
 * @tparam UNUSED An unused parameter since we are assuming that the test and
 *                trial space have the same number of support points.
 *
 * ### ThermoPoroElastic Description
 * Implements the KernelBase interface functions required for using the
 * effective stress for the integration of the stress divergence. This is
 * templated on one of the existing solid mechanics kernel implementations
 * such as QuasiStatic.
 */
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE,
          template< typename,
                    typename,
                    typename > class BASE >
class ThermoPoroElastic : public BASE< SUBREGION_TYPE,
                                       CONSTITUTIVE_TYPE,
                                       FE_TYPE >
{
public:
  /// Alias for the base class;
  using Base = BASE< SUBREGION_TYPE,
                     CONSTITUTIVE_TYPE,
                     FE_TYPE >;

  using Base::m_constitutiveUpdate;
  using typename Base::StackVariables;

  /**
   * @brief Constructor
   * @copydoc geosx::finiteElement::ImplicitKernelBase::ImplicitKernelBase
   * @param inputGravityVector The gravity vector.
   */
  ThermoPoroElastic( NodeManager const & nodeManager,
                     EdgeManager const & edgeManager,
                     FaceManager const & faceManager,
                     localIndex const targetRegionIndex,
                     SUBREGION_TYPE const & elementSubRegion,
                     FE_TYPE const & finiteElementSpace,
                     CONSTITUTIVE_TYPE & inputConstitutiveType,
                     arrayView1d< globalIndex const > const inputDofNumber,
                     globalIndex const rankOffset,
                     CRSMatrixView< real64, globalIndex const > const inputMatrix,
                     arrayView1d< real64 > const inputRhs,
                     real64 const (&inputGravityVector)[3] ):
    Base( nodeManager,
          edgeManager,
          faceManager,
          targetRegionIndex,
          elementSubRegion,
          finiteElementSpace,
          inputConstitutiveType,
          inputDofNumber,
          rankOffset,
          inputMatrix,
          inputRhs,
          inputGravityVector ),
    m_fluidDensity( elementSubRegion.template getConstitutiveModel< constitutive::SingleFluidBase >( elementSubRegion.template getReference< string >( fluidModelKey ) ).density() ),
    m_fluidDensity_n( elementSubRegion.template getConstitutiveModel< constitutive::SingleFluidBase >( elementSubRegion.template getReference< string >( fluidModelKey ) ).density_n() ),
    m_initialFluidDensity( elementSubRegion.template getConstitutiveModel< constitutive::SingleFluidBase >( elementSubRegion.template getReference< string >( fluidModelKey ) ).initialDensity() ),
    m_initialFluidPressure( elementSubRegion.template getExtrinsicData< extrinsicMeshData::flow::initialPressure >() ),
    m_fluidPressure_n( elementSubRegion.template getExtrinsicData< extrinsicMeshData::flow::pressure_n >() ),
    m_fluidPressure( elementSubRegion.template getExtrinsicData< extrinsicMeshData::flow::pressure >() )  
  {}

  /**
   * @copydoc geosx::finiteElement::KernelBase::quadraturePointKernel
   * @tparam STRESS_MODIFIER Type of optional functor to allow for the
   * modification of stress prior to integration.
   * @param stressModifier An optional functor to allow for the modification
   *  of stress prior to integration.
   * For solid mechanics kernels, the strain increment is calculated, and the
   * constitutive update is called. In addition, the constitutive stiffness
   * stack variable is filled by the constitutive model.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack ) const
  {
    real64 dNdX[ numNodesPerElem ][ 3 ];
    real64 const detJ = m_finiteElementSpace.template getGradN< FE_TYPE >( k, q, stack.xLocal, dNdX );

    real64 strainInc[6] = {0};
    real64 stress[6] = {0};

    typename CONSTITUTIVE_TYPE::KernelWrapper::DiscretizationOps stiffness;

    FE_TYPE::symmetricGradient( dNdX, stack.uhat_local, strainInc );

    m_constitutiveUpdate.smallStrainUpdate( k, q, strainInc, stress, stiffness );

    stressModifier( stress );
    for( localIndex i=0; i<6; ++i )
    {
      stress[i] *= -detJ;
    }

    real64 const gravityForce[3] = { m_gravityVector[0] * m_density( k, q )* detJ,
                                     m_gravityVector[1] * m_density( k, q )* detJ,
                                     m_gravityVector[2] * m_density( k, q )* detJ };

    real64 N[numNodesPerElem];
    FE_TYPE::calcN( q, N );
    FE_TYPE::plusGradNajAijPlusNaFi( dNdX,
                                     stress,
                                     N,
                                     gravityForce,
                                     reinterpret_cast< real64 (&)[numNodesPerElem][3] >(stack.localResidual) );
    stiffness.template upperBTDB< numNodesPerElem >( dNdX, -detJ, stack.localJacobian );
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

    // TODO: Does this work if BTDB is non-symmetric?
    CONSTITUTIVE_TYPE::KernelWrapper::DiscretizationOps::template fillLowerBTDB< numNodesPerElem >( stack.localJacobian );

    for( int localNode = 0; localNode < numNodesPerElem; ++localNode )
    {
      for( int dim = 0; dim < numDofPerTestSupportPoint; ++dim )
      {
        localIndex const dof =
          LvArray::integerConversion< localIndex >( stack.localRowDofIndex[ numDofPerTestSupportPoint * localNode + dim ] - m_dofRankOffset );
        if( dof < 0 || dof >= m_matrix.numRows() ) continue;
        m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                                stack.localRowDofIndex,
                                                                                stack.localJacobian[ numDofPerTestSupportPoint * localNode + dim ],
                                                                                numNodesPerElem * numDofPerTrialSupportPoint );

        RAJA::atomicAdd< parallelDeviceAtomic >( &m_rhs[ dof ], stack.localResidual[ numDofPerTestSupportPoint * localNode + dim ] );
        maxForce = fmax( maxForce, fabs( stack.localResidual[ numDofPerTestSupportPoint * localNode + dim ] ) );
      }
    }


    return maxForce;
  }



protected:
  
  /// The rank global densities
  arrayView2d< real64 const > const m_fluidDensity;
  arrayView2d< real64 const > const m_fluidDensity_n;
  arrayView2d< real64 const > const m_initialFluidDensity;

  /// The rank-global initial fluid pressure array
  arrayView1d< real64 const > const m_initialFluidPressure;

  /// The rank-global fluid pressure arrays.
  arrayView1d< real64 const > const m_fluidPressure_n;
  arrayView1d< real64 const > const m_fluidPressure;

  /// The rank-global initial temperature array
  arrayView1d< real64 const > const m_initialTemperature;

  /// The rank-global temperature arrays.
  arrayView1d< real64 const > const m_temperature;

};

/// The factory used to construct a ThermoPoroElastic kernel.
using ThermoPoroElasticFactory = finiteElement::KernelFactory< ThermoPoroElastic,
                                                               arrayView1d< globalIndex const > const,
                                                               globalIndex,
                                                               CRSMatrixView< real64, globalIndex const > const,
                                                               arrayView1d< real64 > const,
                                                               real64 const (&)[3], 
                                                               string const >;

} // namespace solidMechanicsLagrangianFEMKernels

} // namespace geosx

#include "finiteElement/kernelInterface/SparsityKernelBase.hpp"

#endif // GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSSMALLSTRAINQUASISTATIC_HPP_
