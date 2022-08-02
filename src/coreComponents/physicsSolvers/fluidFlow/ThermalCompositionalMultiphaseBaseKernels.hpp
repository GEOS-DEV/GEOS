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
 * @file ThermalCompositionalMultiphaseBaseKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_THERMALCOMPOSITIONALMULTIPHASEBASEKERNELS_HPP
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_THERMALCOMPOSITIONALMULTIPHASEBASEKERNELS_HPP

#include "physicsSolvers/fluidFlow/IsothermalCompositionalMultiphaseBaseKernels.hpp"

namespace geosx
{

namespace thermalCompositionalMultiphaseBaseKernels
{

using namespace constitutive;


/******************************** PhaseVolumeFractionKernel ********************************/

/**
 * @class PhaseVolumeFractionKernel
 * @tparam NUM_COMP number of fluid components
 * @tparam NUM_PHASE number of fluid phases
 * @brief Define the interface for the property kernel in charge of computing the phase volume fractions
 */
template< integer NUM_COMP, integer NUM_PHASE >
class PhaseVolumeFractionKernel : public isothermalCompositionalMultiphaseBaseKernels::PhaseVolumeFractionKernel< NUM_COMP, NUM_PHASE >
{
public:

  using Base = isothermalCompositionalMultiphaseBaseKernels::PhaseVolumeFractionKernel< NUM_COMP, NUM_PHASE >;
  using Base::m_dPhaseDens;
  using Base::m_dPhaseFrac;

  /**
   * @brief Constructor
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   */
  PhaseVolumeFractionKernel( ObjectManagerBase & subRegion,
                             MultiFluidBase const & fluid )
    : Base( subRegion, fluid ),
    m_dPhaseVolFrac_dTemp( subRegion.getExtrinsicData< extrinsicMeshData::flow::dPhaseVolumeFraction_dTemperature >() )
  {}

  /**
   * @brief Compute the phase volume fractions in an element
   * @param[in] ei the element index
   */
  GEOSX_HOST_DEVICE
  void compute( localIndex const ei ) const
  {
    using Deriv = multifluid::DerivativeOffset;

    arraySlice2d< real64 const, multifluid::USD_PHASE_DC - 2 > const dPhaseDens = m_dPhaseDens[ei][0];
    arraySlice2d< real64 const, multifluid::USD_PHASE_DC - 2 > const dPhaseFrac = m_dPhaseFrac[ei][0];

    arraySlice1d< real64, compflow::USD_PHASE - 1 > const dPhaseVolFrac_dTemp = m_dPhaseVolFrac_dTemp[ei];
    LvArray::forValuesInSlice( dPhaseVolFrac_dTemp, []( real64 & val ){ val = 0.0; } );

    // Call the base compute the compute the phase volume fraction
    Base::compute( ei, [&] ( localIndex const ip,
                             real64 const & phaseVolFrac,
                             real64 const & phaseDensInv,
                             real64 const & totalDensity )
    {
      // when this lambda is called, we are in the phase loop
      // for each phase ip, compute the derivative of phase volume fraction wrt temperature
      dPhaseVolFrac_dTemp[ip] = (dPhaseFrac[ip][Deriv::dT] - phaseVolFrac * dPhaseDens[ip][Deriv::dT]) * phaseDensInv;
      dPhaseVolFrac_dTemp[ip] *= totalDensity;
    } );
  }

protected:

  // outputs

  /// Views on thermal derivatives of phase volume fractions
  arrayView2d< real64, compflow::USD_PHASE > m_dPhaseVolFrac_dTemp;

};

/**
 * @class PhaseVolumeFractionKernelFactory
 */
class PhaseVolumeFractionKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] numComp the number of fluid components
   * @param[in] numPhase the number of fluid phases
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   */
  template< typename POLICY >
  static void
  createAndLaunch( integer const numComp,
                   integer const numPhase,
                   ObjectManagerBase & subRegion,
                   MultiFluidBase const & fluid )
  {
    if( numPhase == 2 )
    {
      isothermalCompositionalMultiphaseBaseKernels::
        internal::kernelLaunchSelectorCompSwitch( numComp, [&] ( auto NC )
      {
        integer constexpr NUM_COMP = NC();
        PhaseVolumeFractionKernel< NUM_COMP, 2 > kernel( subRegion, fluid );
        PhaseVolumeFractionKernel< NUM_COMP, 2 >::template launch< POLICY >( subRegion.size(), kernel );
      } );
    }
    else if( numPhase == 3 )
    {
      isothermalCompositionalMultiphaseBaseKernels::
        internal::kernelLaunchSelectorCompSwitch( numComp, [&] ( auto NC )
      {
        integer constexpr NUM_COMP = NC();
        PhaseVolumeFractionKernel< NUM_COMP, 3 > kernel( subRegion, fluid );
        PhaseVolumeFractionKernel< NUM_COMP, 3 >::template launch< POLICY >( subRegion.size(), kernel );
      } );
    }
  }
};


/******************************** ElementBasedAssemblyKernel ********************************/

/**
 * @class ElementBasedAssemblyKernel
 * @tparam NUM_COMP number of fluid components
 * @tparam NUM_DOF number of degrees of freedom
 * @brief Define the interface for the assembly kernel in charge of thermal accumulation and volume balance
 */
template< localIndex NUM_COMP, localIndex NUM_DOF >
class ElementBasedAssemblyKernel : public isothermalCompositionalMultiphaseBaseKernels::ElementBasedAssemblyKernel< NUM_COMP, NUM_DOF >
{
public:

  using Base = isothermalCompositionalMultiphaseBaseKernels::ElementBasedAssemblyKernel< NUM_COMP, NUM_DOF >;
  using Base::numComp;
  using Base::numDof;
  using Base::numEqn;
  using Base::m_numPhases;
  using Base::m_rankOffset;
  using Base::m_dofNumber;
  using Base::m_elemGhostRank;
  using Base::m_volume;
  using Base::m_porosity_n;
  using Base::m_porosity;
  using Base::m_dPoro_dPres;
  using Base::m_dCompFrac_dCompDens;
  using Base::m_phaseVolFrac_n;
  using Base::m_phaseVolFrac;
  using Base::m_dPhaseVolFrac_dPres;
  using Base::m_dPhaseVolFrac_dCompDens;
  using Base::m_phaseDens_n;
  using Base::m_phaseDens;
  using Base::m_dPhaseDens;
  using Base::m_phaseCompFrac_n;
  using Base::m_phaseCompFrac;
  using Base::m_dPhaseCompFrac;
  using Base::m_localMatrix;
  using Base::m_localRhs;

  /**
   * @brief Constructor
   * @param[in] numPhases the number of fluid phases
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey the string key to retrieve the degress of freedom numbers
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   * @param[in] solid the solid model
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  ElementBasedAssemblyKernel( localIndex const numPhases,
                              globalIndex const rankOffset,
                              string const dofKey,
                              ElementSubRegionBase const & subRegion,
                              MultiFluidBase const & fluid,
                              CoupledSolidBase const & solid,
                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                              arrayView1d< real64 > const & localRhs )
    : Base( numPhases, rankOffset, dofKey, subRegion, fluid, solid, localMatrix, localRhs ),
    m_dPhaseVolFrac_dTemp( subRegion.getExtrinsicData< extrinsicMeshData::flow::dPhaseVolumeFraction_dTemperature >() ),
    m_phaseInternalEnergy_n( fluid.phaseInternalEnergy_n() ),
    m_phaseInternalEnergy( fluid.phaseInternalEnergy() ),
    m_dPhaseInternalEnergy( fluid.dPhaseInternalEnergy() ),
    m_rockInternalEnergy_n( solid.getInternalEnergy_n() ),
    m_rockInternalEnergy( solid.getInternalEnergy() ),
    m_dRockInternalEnergy_dTemp( solid.getDinternalEnergy_dTemperature() )
  {}

  struct StackVariables : public Base::StackVariables
  {
public:

    GEOSX_HOST_DEVICE
    StackVariables()
      : Base::StackVariables()
    {}

    using Base::StackVariables::poreVolume;
    using Base::StackVariables::poreVolume_n;
    using Base::StackVariables::dPoreVolume_dPres;
    using Base::StackVariables::localRow;
    using Base::StackVariables::dofIndices;
    using Base::StackVariables::localResidual;
    using Base::StackVariables::localJacobian;

    // Solid energy

    /// Solid energy at time n+1
    real64 solidEnergy = 0.0;

    /// Solid energy at the previous converged time step
    real64 solidEnergy_n = 0.0;

    /// Derivative of solid internal energy with respect to pressure
    real64 dSolidEnergy_dPres = 0.0;

    /// Derivative of solid internal energy with respect to temperature
    real64 dSolidEnergy_dTemp = 0.0;

  };

  /**
   * @brief Performs the setup phase for the kernel.
   * @param[in] ei the element index
   * @param[in] stack the stack variables
   */
  GEOSX_HOST_DEVICE
  void setup( localIndex const ei,
              StackVariables & stack ) const
  {
    Base::setup( ei, stack );

    // initialize the solid volume
    real64 const solidVolume = m_volume[ei] * ( 1.0 - m_porosity[ei][0] );
    real64 const solidVolume_n = m_volume[ei] * ( 1.0 - m_porosity_n[ei][0] );
    real64 const dSolidVolume_dPres = -m_volume[ei] * m_dPoro_dPres[ei][0];

    // initialize the solid internal energy
    stack.solidEnergy = solidVolume * m_rockInternalEnergy[ei][0];
    stack.solidEnergy_n = solidVolume_n * m_rockInternalEnergy_n[ei][0];
    stack.dSolidEnergy_dPres = dSolidVolume_dPres * m_rockInternalEnergy[ei][0];
    stack.dSolidEnergy_dTemp = solidVolume * m_dRockInternalEnergy_dTemp[ei][0];
  }

  /**
   * @brief Compute the local accumulation contributions to the residual and Jacobian
   * @tparam FUNC the type of the function that can be used to customize the kernel
   * @param[in] ei the element index
   * @param[inout] stack the stack variables
   */
  GEOSX_HOST_DEVICE
  void computeAccumulation( localIndex const ei,
                            StackVariables & stack ) const
  {
    using Deriv = multifluid::DerivativeOffset;

    Base::computeAccumulation( ei, stack, [&] ( integer const ip,
                                                real64 const & phaseAmount,
                                                real64 const & phaseAmount_n,
                                                real64 const & dPhaseAmount_dP,
                                                real64 const (&dPhaseAmount_dC)[numComp] )
    {
      // We are in the loop over phases, ip provides the current phase index.
      // We have to do two things:
      //   1- Assemble the derivatives of the component mass balance equations with respect to temperature
      //   2- Assemble the phase-dependent part of the accumulation term of the energy equation

      real64 dPhaseInternalEnergy_dC[numComp]{};

      // construct the slices
      arraySlice2d< real64 const, compflow::USD_COMP_DC - 1 > dCompFrac_dCompDens = m_dCompFrac_dCompDens[ei];
      arraySlice1d< real64 const, compflow::USD_PHASE - 1 > phaseVolFrac = m_phaseVolFrac[ei];
      arraySlice1d< real64 const, compflow::USD_PHASE - 1 > dPhaseVolFrac_dTemp = m_dPhaseVolFrac_dTemp[ei];
      arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > phaseDens = m_phaseDens[ei][0];
      arraySlice2d< real64 const, multifluid::USD_PHASE_DC - 2 > dPhaseDens = m_dPhaseDens[ei][0];
      arraySlice2d< real64 const, multifluid::USD_PHASE_COMP - 2 > phaseCompFrac = m_phaseCompFrac[ei][0];
      arraySlice3d< real64 const, multifluid::USD_PHASE_COMP_DC - 2 > dPhaseCompFrac = m_dPhaseCompFrac[ei][0];
      arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > phaseInternalEnergy_n = m_phaseInternalEnergy_n[ei][0];
      arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > phaseInternalEnergy = m_phaseInternalEnergy[ei][0];
      arraySlice2d< real64 const, multifluid::USD_PHASE_DC - 2 > dPhaseInternalEnergy = m_dPhaseInternalEnergy[ei][0];

      // Step 1: assemble the derivatives of the component mass balance equations with respect to temperature

      real64 const dPhaseAmount_dT = stack.poreVolume
                                     * (dPhaseVolFrac_dTemp[ip] * phaseDens[ip] + phaseVolFrac[ip] * dPhaseDens[ip][Deriv::dT] );
      for( integer ic = 0; ic < numComp; ++ic )
      {
        stack.localJacobian[ic][numDof-1] += dPhaseAmount_dT * phaseCompFrac[ip][ic]
                                             + phaseAmount * dPhaseCompFrac[ip][ic][Deriv::dT];
      }

      // Step 2: assemble the phase-dependent part of the accumulation term of the energy equation

      real64 const phaseEnergy = phaseAmount * phaseInternalEnergy[ip];
      real64 const phaseEnergy_n = phaseAmount_n * phaseInternalEnergy_n[ip];
      real64 const dPhaseEnergy_dP = dPhaseAmount_dP * phaseInternalEnergy[ip]
                                     + phaseAmount * dPhaseInternalEnergy[ip][Deriv::dP];
      real64 const dPhaseEnergy_dT = dPhaseAmount_dT * phaseInternalEnergy[ip]
                                     + phaseAmount * dPhaseInternalEnergy[ip][Deriv::dT];

      // local accumulation
      stack.localResidual[numEqn-1] += phaseEnergy - phaseEnergy_n;

      // derivatives w.r.t. pressure and temperature
      stack.localJacobian[numEqn-1][0]        += dPhaseEnergy_dP;
      stack.localJacobian[numEqn-1][numDof-1] += dPhaseEnergy_dT;

      // derivatives w.r.t. component densities
      applyChainRule( numComp, dCompFrac_dCompDens, dPhaseInternalEnergy[ip], dPhaseInternalEnergy_dC, Deriv::dC );
      for( integer jc = 0; jc < numComp; ++jc )
      {
        stack.localJacobian[numEqn-1][jc + 1] += phaseInternalEnergy[ip] * dPhaseAmount_dC[jc]
                                                 + dPhaseInternalEnergy_dC[jc] * phaseAmount;
      }
    } );

    // Step 3: assemble the solid part of the accumulation term

    // local accumulation and derivatives w.r.t. pressure and temperature
    stack.localResidual[numEqn-1] += stack.solidEnergy - stack.solidEnergy_n;
    stack.localJacobian[numEqn-1][0] += stack.dSolidEnergy_dPres;
    stack.localJacobian[numEqn-1][numDof-1] += stack.dSolidEnergy_dTemp;

  }

  /**
   * @brief Compute the local volume balance contributions to the residual and Jacobian
   * @tparam FUNC the type of the function that can be used to customize the kernel
   * @param[in] ei the element index
   * @param[inout] stack the stack variables
   */
  GEOSX_HOST_DEVICE
  void computeVolumeBalance( localIndex const ei,
                             StackVariables & stack ) const
  {
    Base::computeVolumeBalance( ei, stack, [&] ( real64 const & oneMinusPhaseVolFraction )
    {
      GEOSX_UNUSED_VAR( oneMinusPhaseVolFraction );

      arraySlice1d< real64 const, compflow::USD_PHASE - 1 > dPhaseVolFrac_dTemp = m_dPhaseVolFrac_dTemp[ei];

      for( integer ip = 0; ip < m_numPhases; ++ip )
      {
        stack.localJacobian[numEqn-2][numDof-1] -= dPhaseVolFrac_dTemp[ip];
      }
    } );
  }

  GEOSX_HOST_DEVICE
  void complete( localIndex const ei,
                 StackVariables & stack ) const
  {
    // Step 1: assemble the component mass balance equations and volume balance equations
    Base::complete( ei, stack );

    // Step 2: assemble the energy equation
    m_localRhs[stack.localRow + numEqn-1] += stack.localResidual[numEqn-1];
    m_localMatrix.template addToRow< serialAtomic >( stack.localRow + numEqn-1,
                                                     stack.dofIndices,
                                                     stack.localJacobian[numEqn-1],
                                                     numDof );
  }

protected:

  /// Views on derivatives wrt to temperature for phase volume fraction
  arrayView2d< real64 const, compflow::USD_PHASE > m_dPhaseVolFrac_dTemp;

  /// Views on phase internal energy
  arrayView3d< real64 const, multifluid::USD_PHASE > m_phaseInternalEnergy_n;
  arrayView3d< real64 const, multifluid::USD_PHASE > m_phaseInternalEnergy;
  arrayView4d< real64 const, multifluid::USD_PHASE_DC > m_dPhaseInternalEnergy;

  /// Views on rock internal energy
  arrayView2d< real64 const > m_rockInternalEnergy_n;
  arrayView2d< real64 const > m_rockInternalEnergy;
  arrayView2d< real64 const > m_dRockInternalEnergy_dTemp;

};

/**
 * @class ElementBasedAssemblyKernelFactory
 */
class ElementBasedAssemblyKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] numComps the number of fluid components
   * @param[in] numPhases the number of fluid phases
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey the string key to retrieve the degress of freedom numbers
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   * @param[in] solid the solid model
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  template< typename POLICY >
  static void
  createAndLaunch( localIndex const numComps,
                   localIndex const numPhases,
                   globalIndex const rankOffset,
                   string const dofKey,
                   ElementSubRegionBase const & subRegion,
                   MultiFluidBase const & fluid,
                   CoupledSolidBase const & solid,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs )
  {
    isothermalCompositionalMultiphaseBaseKernels::
      internal::kernelLaunchSelectorCompSwitch( numComps, [&] ( auto NC )
    {
      localIndex constexpr NUM_COMP = NC();
      localIndex constexpr NUM_DOF = NC()+2;
      ElementBasedAssemblyKernel< NUM_COMP, NUM_DOF >
      kernel( numPhases, rankOffset, dofKey, subRegion, fluid, solid, localMatrix, localRhs );
      ElementBasedAssemblyKernel< NUM_COMP, NUM_DOF >::template
      launch< POLICY, ElementBasedAssemblyKernel< NUM_COMP, NUM_DOF > >( subRegion.size(), kernel );
    } );
  }

};

/******************************** FluidUpdateKernel ********************************/

struct FluidUpdateKernel
{
  template< typename POLICY, typename FLUID_WRAPPER >
  static void
  launch( localIndex const size,
          FLUID_WRAPPER const & fluidWrapper,
          arrayView1d< real64 const > const & pres,
          arrayView1d< real64 const > const & temp,
          arrayView2d< real64 const, compflow::USD_COMP > const & compFrac )
  {
    forAll< POLICY >( size, [=] GEOSX_HOST_DEVICE ( localIndex const k )
    {
      for( localIndex q = 0; q < fluidWrapper.numGauss(); ++q )
      {
        fluidWrapper.update( k, q, pres[k], temp[k], compFrac[k] );
      }
    } );
  }

  template< typename POLICY, typename FLUID_WRAPPER >
  static void
  launch( SortedArrayView< localIndex const > const & targetSet,
          FLUID_WRAPPER const & fluidWrapper,
          arrayView1d< real64 const > const & pres,
          arrayView1d< real64 const > const & temp,
          arrayView2d< real64 const, compflow::USD_COMP > const & compFrac )
  {
    forAll< POLICY >( targetSet.size(), [=] GEOSX_HOST_DEVICE ( localIndex const a )
    {
      localIndex const k = targetSet[a];
      for( localIndex q = 0; q < fluidWrapper.numGauss(); ++q )
      {
        fluidWrapper.update( k, q, pres[k], temp[k], compFrac[k] );
      }
    } );
  }
};

/******************************** SolidInternalEnergyUpdateKernel ********************************/

struct SolidInternalEnergyUpdateKernel
{

  template< typename POLICY, typename SOLID_INTERNAL_ENERGY_WRAPPER >
  static void
  launch( localIndex const size,
          SOLID_INTERNAL_ENERGY_WRAPPER const & solidInternalEnergyWrapper,
          arrayView1d< real64 const > const & temp )
  {
    forAll< POLICY >( size, [=] GEOSX_HOST_DEVICE ( localIndex const k )
    {
      solidInternalEnergyWrapper.update( k, temp[k] );
    } );
  }
};

/******************************** ResidualNormKernel ********************************/

/**
 * @class ResidualNormKernel
 */
class ResidualNormKernel : public solverBaseKernels::ResidualNormKernelBase< 2 >
{
public:

  using Base = ResidualNormKernelBase< 2 >;
  using Base::minNormalizer;
  using Base::m_rankOffset;
  using Base::m_localResidual;
  using Base::m_dofNumber;

  ResidualNormKernel( globalIndex const rankOffset,
                      arrayView1d< real64 const > const & localResidual,
                      arrayView1d< globalIndex const > const & dofNumber,
                      arrayView1d< localIndex const > const & ghostRank,
                      integer const numComponents,
                      integer const numPhases,
                      ElementSubRegionBase const & subRegion,
                      MultiFluidBase const & fluid,
                      CoupledSolidBase const & solid,
                      SolidInternalEnergy const & solidInternalEnergy )
    : Base( rankOffset,
            localResidual,
            dofNumber,
            ghostRank ),
    m_numComponents( numComponents ),
    m_numPhases( numPhases ),
    m_volume( subRegion.getElementVolume() ),
    m_porosity_n( solid.getPorosity_n() ),
    m_phaseVolFrac_n( subRegion.getExtrinsicData< extrinsicMeshData::flow::phaseVolumeFraction_n >() ),
    m_totalDens_n( fluid.totalDensity_n() ),
    m_phaseDens_n( fluid.phaseDensity_n() ),
    m_phaseInternalEnergy_n( fluid.phaseInternalEnergy_n() ),
    m_solidInternalEnergy_n( solidInternalEnergy.getInternalEnergy_n() )
  {}

  GEOSX_HOST_DEVICE
  void computeMassEnergyNormalizers( localIndex const ei,
                                     real64 & massNormalizer,
                                     real64 & energyNormalizer ) const
  {
    massNormalizer = LvArray::math::max( minNormalizer, m_totalDens_n[ei][0] * m_porosity_n[ei][0] * m_volume[ei] );
    real64 const poreVolume = m_porosity_n[ei][0] * m_volume[ei];
    energyNormalizer = m_solidInternalEnergy_n[ei][0] * ( 1.0 - m_porosity_n[ei][0] ) * m_volume[ei];
    for( integer ip = 0; ip < m_numPhases; ++ip )
    {
      energyNormalizer += m_phaseInternalEnergy_n[ei][0][ip] * m_phaseDens_n[ei][0][ip] * m_phaseVolFrac_n[ei][ip] * poreVolume;
    }
    energyNormalizer = LvArray::math::max( minNormalizer, energyNormalizer );
  }

  GEOSX_HOST_DEVICE
  virtual void computeLinf( localIndex const ei,
                            LinfStackVariables & stack ) const override
  {
    real64 massNormalizer = 0.0, energyNormalizer = 0.0;
    computeMassEnergyNormalizers( ei, massNormalizer, energyNormalizer );
    real64 const volumeNormalizer = LvArray::math::max( minNormalizer, m_porosity_n[ei][0] * m_volume[ei] );

    // step 1: mass residual

    for( integer idof = 0; idof < m_numComponents; ++idof )
    {
      real64 const valMass = LvArray::math::abs( m_localResidual[stack.localRow + idof] ) / massNormalizer;
      if( valMass > stack.localValue[0] )
      {
        stack.localValue[0] = valMass;
      }
    }

    // step 2: volume residual

    real64 const valVol = LvArray::math::abs( m_localResidual[stack.localRow + m_numComponents] ) / volumeNormalizer;
    if( valVol > stack.localValue[0] )
    {
      stack.localValue[0] = valVol;
    }

    // step 3: energy residual

    real64 const valEnergy = LvArray::math::abs( m_localResidual[stack.localRow + m_numComponents + 1] ) / energyNormalizer;
    if( valEnergy > stack.localValue[1] )
    {
      stack.localValue[1] = valEnergy;
    }
  }

  GEOSX_HOST_DEVICE
  virtual void computeL2( localIndex const ei,
                          L2StackVariables & stack ) const override
  {
    // note: for the L2 norm, we bundle the volume and mass residuals/normalizers
    real64 massNormalizer = 0.0, energyNormalizer = 0.0;
    computeMassEnergyNormalizers( ei, massNormalizer, energyNormalizer );

    // step 1: mass residual

    for( integer idof = 0; idof < m_numComponents; ++idof )
    {
      stack.localValue[0] += m_localResidual[stack.localRow + idof] * m_localResidual[stack.localRow + idof];
      stack.localNormalizer[0] += massNormalizer;
    }

    // step 2: volume residual

    real64 const valVol = m_localResidual[stack.localRow + m_numComponents] * m_totalDens_n[ei][0]; // we need a mass here, hence the
                                                                                                    // multiplication
    stack.localValue[0] += valVol * valVol;
    stack.localNormalizer[0] += massNormalizer;

    // step 3: energy residual

    stack.localValue[1] += m_localResidual[stack.localRow + m_numComponents + 1] * m_localResidual[stack.localRow + m_numComponents + 1];
    stack.localNormalizer[1] += energyNormalizer;
  }

protected:

  /// Number of fluid components
  integer const m_numComponents;

  /// Number of fluid phases
  integer const m_numPhases;

  /// View on the volume
  arrayView1d< real64 const > const m_volume;

  /// View on porosity at the previous converged time step
  arrayView2d< real64 const > const m_porosity_n;

  /// View on phase properties at the previous converged time step
  arrayView2d< real64 const, compflow::USD_PHASE > const m_phaseVolFrac_n;
  arrayView2d< real64 const, multifluid::USD_FLUID > const m_totalDens_n;
  arrayView3d< real64 const, multifluid::USD_PHASE > const m_phaseDens_n;
  arrayView3d< real64 const, multifluid::USD_PHASE > const m_phaseInternalEnergy_n;

  /// View on solid properties at the previous converged time step
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
   * @param[in] numComps the number of fluid components
   * @param[in] numPhases the number of fluid phases
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
                   integer const numComps,
                   integer const numPhases,
                   globalIndex const rankOffset,
                   string const dofKey,
                   arrayView1d< real64 const > const & localResidual,
                   ElementSubRegionBase const & subRegion,
                   MultiFluidBase const & fluid,
                   CoupledSolidBase const & solid,
                   SolidInternalEnergy const & solidInternalEnergy,
                   real64 (& residualNorm)[2],
                   real64 (& residualNormalizer)[2] )
  {
    arrayView1d< globalIndex const > const dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
    arrayView1d< integer const > const ghostRank = subRegion.ghostRank();

    ResidualNormKernel kernel( rankOffset, localResidual, dofNumber, ghostRank,
                               numComps, numPhases, subRegion, fluid, solid, solidInternalEnergy );
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


} // namespace thermalCompositionalMultiphaseBaseKernels

} // namespace geosx


#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_THERMALCOMPOSITIONALMULTIPHASEBASEKERNELS_HPP
