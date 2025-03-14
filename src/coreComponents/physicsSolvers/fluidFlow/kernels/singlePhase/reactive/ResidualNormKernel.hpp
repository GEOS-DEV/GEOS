/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file ResidualNormKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_REACTIVE_RESIDUALNORMKERNEL_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_REACTIVE_RESIDUALNORMKERNEL_HPP

#include "physicsSolvers/PhysicsSolverBaseKernels.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseReactiveTransportFields.hpp"

namespace geos
{

namespace singlePhaseReactiveBaseKernels
{

/******************************** ResidualNormKernel ********************************/

/**
 * @class IsothermalResidualNormKernel
 */
class IsothermalResidualNormKernel : public physicsSolverBaseKernels::ResidualNormKernelBase< 2 >
{
public:

  using Base = physicsSolverBaseKernels::ResidualNormKernelBase< 2 >;
  using Base::m_minNormalizer;
  using Base::m_rankOffset;
  using Base::m_localResidual;
  using Base::m_dofNumber;

  IsothermalResidualNormKernel( globalIndex const rankOffset,
                                arrayView1d< real64 const > const & localResidual,
                                arrayView1d< globalIndex const > const & dofNumber,
                                arrayView1d< localIndex const > const & ghostRank,
                                integer const numPrimarySpecies,
                                ElementSubRegionBase const & subRegion,
                                real64 const minNormalizer )
    : Base( rankOffset,
            localResidual,
            dofNumber,
            ghostRank,
            minNormalizer ),
    m_numPrimarySpecies( numPrimarySpecies ),
    m_mass_n( subRegion.template getField< fields::flow::mass_n >() ),
    m_totalPrimarySpeciesAmount_n( subRegion.getField< fields::flow::totalPrimarySpeciesAmount_n >() )
  {}

  GEOS_HOST_DEVICE
  virtual void computeLinf( localIndex const ei,
                            LinfStackVariables & stack ) const override
  {
    real64 const totalMassNormalizer = LvArray::math::max( m_minNormalizer, m_mass_n[ei] );

    // step 1: total mass residuals
    real64 const valMass = LvArray::math::abs( m_localResidual[stack.localRow] ) / totalMassNormalizer;
    if( valMass > stack.localValue[0] )
    {
      stack.localValue[0] = valMass;
    }
    
    // step 2: species amount residuals
    for( integer idof = 0; idof < m_numPrimarySpecies; ++idof )
    {
      real64 const speciesAmountNormalizer = LvArray::math::max( m_minNormalizer, m_totalPrimarySpeciesAmount_n[ei][idof] );
      real64 const valAmount = LvArray::math::abs( m_localResidual[stack.localRow + idof + 1] ) / speciesAmountNormalizer;
      if( valAmount > stack.localValue[1] )
      {
        stack.localValue[1] = valAmount;
      }
    }
  }

  GEOS_HOST_DEVICE
  virtual void computeL2( localIndex const ei,
                          L2StackVariables & stack ) const override
  {
    real64 const totalMassNormalizer = LvArray::math::max( m_minNormalizer, m_mass_n[ei] );

    // step 1: total mass residuals
    stack.localValue[0] += m_localResidual[stack.localRow] * m_localResidual[stack.localRow];
    stack.localNormalizer[0] += totalMassNormalizer;

    // step 2: species amount residuals
    for( integer idof = 0; idof < m_numPrimarySpecies; ++idof )
    {
      real64 const speciesAmountNormalizer = LvArray::math::max( m_minNormalizer, m_totalPrimarySpeciesAmount_n[ei][idof] );

      stack.localValue[1] += m_localResidual[stack.localRow + idof + 1] * m_localResidual[stack.localRow + idof + 1];
      stack.localNormalizer[1] += speciesAmountNormalizer;
    }
  }


protected:

  /// Number of primary species
  integer const m_numPrimarySpecies;

  /// View on mass at the previous converged time step
  arrayView1d< real64 const > const m_mass_n;

  // View on primary species amount (moles) from previous time step
  arrayView2d< real64 const, compflow::USD_COMP > m_totalPrimarySpeciesAmount_n;

};

/**
 * @class ThermalResidualNormKernel
 */
class ThermalResidualNormKernel : public physicsSolverBaseKernels::ResidualNormKernelBase< 3 >
{
public:

  using Base = physicsSolverBaseKernels::ResidualNormKernelBase< 3 >;
  using Base::m_minNormalizer;
  using Base::m_rankOffset;
  using Base::m_localResidual;
  using Base::m_dofNumber;

  ThermalResidualNormKernel( globalIndex const rankOffset,
                             arrayView1d< real64 const > const & localResidual,
                             arrayView1d< globalIndex const > const & dofNumber,
                             arrayView1d< localIndex const > const & ghostRank,
                             integer const numPrimarySpecies,
                             ElementSubRegionBase const & subRegion,
                             real64 const minNormalizer )
    : Base( rankOffset,
            localResidual,
            dofNumber,
            ghostRank,
            minNormalizer ),
    m_numPrimarySpecies( numPrimarySpecies ),
    m_mass_n( subRegion.template getField< fields::flow::mass_n >() ),
    m_totalPrimarySpeciesAmount_n( subRegion.getField< fields::flow::totalPrimarySpeciesAmount_n >() ),
    m_energy_n( subRegion.template getField< fields::flow::energy_n >() )
  {}

  GEOS_HOST_DEVICE
  void computeMassEnergyNormalizers( localIndex const ei,
                                     real64 & totalMassNormalizer,
                                     real64 & energyNormalizer ) const
  {
    totalMassNormalizer = LvArray::math::max( m_minNormalizer, m_mass_n[ei] );
    energyNormalizer = LvArray::math::max( m_minNormalizer, LvArray::math::abs( m_energy_n[ei] ) ); // energy can be negative
  }

  GEOS_HOST_DEVICE
  virtual void computeLinf( localIndex const ei,
                            LinfStackVariables & stack ) const override
  {
    real64 totalMassNormalizer = 0.0, energyNormalizer = 0.0;
    computeMassEnergyNormalizers( ei, totalMassNormalizer, energyNormalizer );

    // step 1: mass residual

    real64 const valMass = LvArray::math::abs( m_localResidual[stack.localRow] ) / totalMassNormalizer;
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

    // step 3: species amount residuals
    for( integer idof = 0; idof < m_numPrimarySpecies; ++idof )
    {
      real64 const speciesAmountNormalizer = LvArray::math::max( m_minNormalizer, m_totalPrimarySpeciesAmount_n[ei][idof] );
      real64 const valAmount = LvArray::math::abs( m_localResidual[stack.localRow + idof + 2] ) / speciesAmountNormalizer;
      if( valAmount > stack.localValue[2] )
      {
        stack.localValue[2] = valAmount;
      }
    }
  }

  GEOS_HOST_DEVICE
  virtual void computeL2( localIndex const ei,
                          L2StackVariables & stack ) const override
  {
    real64 totalMassNormalizer = 0.0, energyNormalizer = 0.0;
    computeMassEnergyNormalizers( ei, totalMassNormalizer, energyNormalizer );

    // step 1: mass residual

    stack.localValue[0] += m_localResidual[stack.localRow] * m_localResidual[stack.localRow];
    stack.localNormalizer[0] += totalMassNormalizer;

    // step 2: energy residual

    stack.localValue[1] += m_localResidual[stack.localRow + 1] * m_localResidual[stack.localRow + 1];
    stack.localNormalizer[1] += energyNormalizer;

    // step 3: species amount residuals
    for( integer idof = 0; idof < m_numPrimarySpecies; ++idof )
    {
      real64 const speciesAmountNormalizer = LvArray::math::max( m_minNormalizer, m_totalPrimarySpeciesAmount_n[ei][idof] );

      stack.localValue[2] += m_localResidual[stack.localRow + idof + 2] * m_localResidual[stack.localRow + idof + 2];
      stack.localNormalizer[2] += speciesAmountNormalizer;
    }
  }


protected:

  /// Number of primary species
  integer const m_numPrimarySpecies;

  /// View on mass at the previous converged time step
  arrayView1d< real64 const > const m_mass_n;

  // View on primary species amount (moles) from previous time step
  arrayView2d< real64 const, compflow::USD_COMP > m_totalPrimarySpeciesAmount_n;

  /// View on energy at the previous converged time step
  arrayView1d< real64 const > const m_energy_n;

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
   * @param[in] numPrimarySpecies the number of primary species
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey the string key to retrieve the degress of freedom numbers
   * @param[in] localResidual the residual vector on my MPI rank
   * @param[in] subRegion the element subregion
   * @param[out] residualNorm the residual norm on the subRegion
   * @param[out] residualNormalizer the residual normalizer on the subRegion
   */
  template< typename POLICY >
  static void
  createAndLaunch( physicsSolverBaseKernels::NormType const normType,
                   integer const numPrimarySpecies,
                   globalIndex const rankOffset,
                   string const dofKey,
                   arrayView1d< real64 const > const & localResidual,
                   ElementSubRegionBase const & subRegion,
                   real64 const minNormalizer,
                   real64 (& residualNorm)[2],
                   real64 (& residualNormalizer)[2] )
  {
    arrayView1d< globalIndex const > const dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
    arrayView1d< integer const > const ghostRank = subRegion.ghostRank();

    IsothermalResidualNormKernel kernel( rankOffset, localResidual, dofNumber, ghostRank, numPrimarySpecies, subRegion, minNormalizer );
    if( normType == physicsSolverBaseKernels::NormType::Linf )
    {
      IsothermalResidualNormKernel::launchLinf< POLICY >( subRegion.size(), kernel, residualNorm );
    }
    else // L2 norm
    {
      IsothermalResidualNormKernel::launchL2< POLICY >( subRegion.size(), kernel, residualNorm, residualNormalizer );
    }
  }

  /**
   * @brief Create a new kernel and launch (thermal version)
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] normType the type of norm used (Linf or L2)
   * @param[in] numPrimarySpecies the number of primary species
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey the string key to retrieve the degress of freedom numbers
   * @param[in] localResidual the residual vector on my MPI rank
   * @param[in] subRegion the element subregion
   * @param[out] residualNorm the residual norm on the subRegion
   * @param[out] residualNormalizer the residual normalizer on the subRegion
   */
  template< typename POLICY >
  static void
  createAndLaunch( physicsSolverBaseKernels::NormType const normType,
                   integer const numPrimarySpecies,
                   globalIndex const rankOffset,
                   string const & dofKey,
                   arrayView1d< real64 const > const & localResidual,
                   ElementSubRegionBase const & subRegion,
                   real64 const minNormalizer,
                   real64 (& residualNorm)[3],
                   real64 (& residualNormalizer)[3] )
  {
    arrayView1d< globalIndex const > const dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
    arrayView1d< integer const > const ghostRank = subRegion.ghostRank();

    ThermalResidualNormKernel kernel( rankOffset, localResidual, dofNumber, ghostRank, numPrimarySpecies, subRegion, minNormalizer );
    if( normType == physicsSolverBaseKernels::NormType::Linf )
    {
      ThermalResidualNormKernel::launchLinf< POLICY >( subRegion.size(), kernel, residualNorm );
    }
    else // L2 norm
    {
      ThermalResidualNormKernel::launchL2< POLICY >( subRegion.size(), kernel, residualNorm, residualNormalizer );
    }
  }

};

} // namespace singlePhaseReactiveBaseKernels

} // namespace geos


#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_REACTIVE_RESIDUALNORMKERNEL_HPP
