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
 * @file ResidualNormKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_THERMALSINGLEPHASEBASEKERNELS_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_THERMALSINGLEPHASEBASEKERNELS_HPP

#include "physicsSolvers/SolverBaseKernels.hpp"
#include "constitutive/solid/CoupledSolidBase.hpp"
#include "constitutive/solid/SolidInternalEnergy.hpp"


namespace geos
{

namespace thermalSinglePhaseBaseKernels
{

/******************************** ResidualNormKernel ********************************/

/**
 * @class ResidualNormKernel
 */
class ResidualNormKernel : public solverBaseKernels::ResidualNormKernelBase< 2 >
{
public:

  using Base = solverBaseKernels::ResidualNormKernelBase< 2 >;
  using Base::minNormalizer;
  using Base::m_rankOffset;
  using Base::m_localResidual;
  using Base::m_dofNumber;

  ResidualNormKernel( globalIndex const rankOffset,
                      arrayView1d< real64 const > const & localResidual,
                      arrayView1d< globalIndex const > const & dofNumber,
                      arrayView1d< localIndex const > const & ghostRank,
                      ElementSubRegionBase const & subRegion,
                      constitutive::SingleFluidBase const & fluid,
                      constitutive::CoupledSolidBase const & solid,
                      constitutive::SolidInternalEnergy const & solidInternalEnergy )
    : Base( rankOffset,
            localResidual,
            dofNumber,
            ghostRank ),
    m_volume( subRegion.getElementVolume() ),
    m_porosity_n( solid.getPorosity_n() ),
    m_density_n( fluid.density_n() ),
    m_fluidInternalEnergy_n( fluid.internalEnergy_n() ),
    m_solidInternalEnergy_n( solidInternalEnergy.getInternalEnergy_n() )
  {}

  GEOS_HOST_DEVICE
  void computeMassEnergyNormalizers( localIndex const ei,
                                     real64 & massNormalizer,
                                     real64 & energyNormalizer ) const
  {
    massNormalizer = LvArray::math::max( minNormalizer, m_density_n[ei][0] * m_porosity_n[ei][0] * m_volume[ei] );
    energyNormalizer =
      LvArray::math::max( minNormalizer,
                          LvArray::math::abs( m_solidInternalEnergy_n[ei][0] * ( 1.0 - m_porosity_n[ei][0] ) * m_volume[ei]
                                              + m_fluidInternalEnergy_n[ei][0] * m_density_n[ei][0] * m_porosity_n[ei][0] * m_volume[ei] ) );
  }

  GEOS_HOST_DEVICE
  virtual void computeLinf( localIndex const ei,
                            LinfStackVariables & stack ) const override
  {
    real64 massNormalizer = 0.0, energyNormalizer = 0.0;
    computeMassEnergyNormalizers( ei, massNormalizer, energyNormalizer );

    // step 1: mass residual

    real64 const valMass = LvArray::math::abs( m_localResidual[stack.localRow] ) / massNormalizer;
    if( valMass > stack.localValue[0] )
    {
      stack.localValue[0] = valMass;
    }

    // step 2: energy residual
    real64 const valEnergy = LvArray::math::abs( m_localResidual[stack.localRow + 1] ) / energyNormalizer;
    if( valEnergy > stack.localValue[1] )
    {
      stack.localValue[1] = valEnergy;
    }
  }

  GEOS_HOST_DEVICE
  virtual void computeL2( localIndex const ei,
                          L2StackVariables & stack ) const override
  {
    real64 massNormalizer = 0.0, energyNormalizer = 0.0;
    computeMassEnergyNormalizers( ei, massNormalizer, energyNormalizer );

    // step 1: mass residual

    stack.localValue[0] += m_localResidual[stack.localRow] * m_localResidual[stack.localRow];
    stack.localNormalizer[0] += massNormalizer;

    // step 2: energy residual

    stack.localValue[1] += m_localResidual[stack.localRow + 1] * m_localResidual[stack.localRow + 1];
    stack.localNormalizer[1] += energyNormalizer;
  }


protected:

  /// View on the volume
  arrayView1d< real64 const > const m_volume;

  /// View on porosity at the previous converged time step
  arrayView2d< real64 const > const m_porosity_n;

  /// View on total mass/molar density at the previous converged time step
  arrayView2d< real64 const > const m_density_n;
  arrayView2d< real64 const > const m_fluidInternalEnergy_n;

  /// View on solid internal energy at the previous converged time step
  arrayView2d< real64 const > const m_solidInternalEnergy_n;

};

/**
 * @class ResidualNormKernelFactory
 */
class ResidualNormKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] normType the type of norm used (Linf or L2)
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey the string key to retrieve the degress of freedom numbers
   * @param[in] localResidual the residual vector on my MPI rank
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   * @param[in] solid the solid model
   * @param[in] solidInternalEnergy the solid internal energy model
   * @param[out] residualNorm the residual norm on the subRegion
   * @param[out] residualNormalizer the residual normalizer on the subRegion
   */
  template< typename POLICY >
  static void
  createAndLaunch( solverBaseKernels::NormType const normType,
                   globalIndex const rankOffset,
                   string const dofKey,
                   arrayView1d< real64 const > const & localResidual,
                   ElementSubRegionBase const & subRegion,
                   constitutive::SingleFluidBase const & fluid,
                   constitutive::CoupledSolidBase const & solid,
                   constitutive::SolidInternalEnergy const & solidInternalEnergy,
                   real64 (& residualNorm)[2],
                   real64 (& residualNormalizer)[2] )
  {
    arrayView1d< globalIndex const > const dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
    arrayView1d< integer const > const ghostRank = subRegion.ghostRank();

    ResidualNormKernel kernel( rankOffset, localResidual, dofNumber, ghostRank, subRegion, fluid, solid, solidInternalEnergy );
    if( normType == solverBaseKernels::NormType::Linf )
    {
      ResidualNormKernel::launchLinf< POLICY >( subRegion.size(), kernel, residualNorm );
    }
    else // L2 norm
    {
      ResidualNormKernel::launchL2< POLICY >( subRegion.size(), kernel, residualNorm, residualNormalizer );
    }
  }

};

} // namespace thermalSinglePhaseBaseKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_THERMALSINGLEPHASEBASEKERNELS_HPP
