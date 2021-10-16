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
 * @file CompositionalMultiphaseFluid.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_COMPOSITIONALMULTIPHASEFLUID_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_COMPOSITIONALMULTIPHASEFLUID_HPP_

#include "constitutive/fluid/MultiFluidBase.hpp"
#include "constitutive/fluid/MultiFluidUtils.hpp"

#include "constitutive/PVTPackage/PVTPackage/source/pvt/pvt.hpp"

namespace geosx
{
namespace constitutive
{

class CompositionalMultiphaseFluid : public MultiFluidBase
{
public:

  using exec_policy = serialPolicy;

  CompositionalMultiphaseFluid( string const & name, Group * const parent );

  virtual std::unique_ptr< ConstitutiveBase >
  deliverClone( string const & name,
                Group * const parent ) const override;

  static string catalogName() { return "CompositionalMultiphaseFluid"; }

  virtual string getCatalogName() const override { return catalogName(); }

  virtual integer getWaterPhaseIndex() const override final;

  struct viewKeyStruct : MultiFluidBase::viewKeyStruct
  {
    static constexpr char const * equationsOfStateString() { return "equationsOfState"; }
    static constexpr char const * componentCriticalPressureString() { return "componentCriticalPressure"; }
    static constexpr char const * componentCriticalTemperatureString() { return "componentCriticalTemperature"; }
    static constexpr char const * componentAcentricFactorString() { return "componentAcentricFactor"; }
    static constexpr char const * componentVolumeShiftString() { return "componentVolumeShift"; }
    static constexpr char const * componentBinaryCoeffString() { return "componentBinaryCoeff"; }
  };

  /**
   * @brief Kernel wrapper class for CompositionalMultiphaseFluid.
   */
  class KernelWrapper final : public MultiFluidBase::KernelWrapper
  {
public:

    GEOSX_HOST_DEVICE
    virtual void compute( real64 const pressure,
                          real64 const temperature,
                          arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition,
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseFraction,
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseDensity,
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseMassDensity,
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseEnthalpy,
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseInternalEnergy,
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseViscosity,
                          arraySlice2d< real64, multifluid::USD_PHASE_COMP-2 > const & phaseCompFraction,
                          real64 & totalDensity ) const override;

    GEOSX_HOST_DEVICE
    virtual void compute( real64 const pressure,
                          real64 const temperature,
                          arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition,
                          PhaseProp::SliceType const phaseFraction,
                          PhaseProp::SliceType const phaseDensity,
                          PhaseProp::SliceType const phaseMassDensity,
                          PhaseProp::SliceType const phaseViscosity,
                          PhaseProp::SliceType const phaseEnthalpy,
                          PhaseProp::SliceType const phaseInternalEnergy,
                          PhaseComp::SliceType const phaseCompFraction,
                          FluidProp::SliceType const totalDensity ) const override;

    GEOSX_HOST_DEVICE
    virtual void update( localIndex const k,
                         localIndex const q,
                         real64 const pressure,
                         real64 const temperature,
                         arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition ) const override;

private:

    friend class CompositionalMultiphaseFluid;

    KernelWrapper( pvt::MultiphaseSystem & fluid,
                   arrayView1d< pvt::PHASE_TYPE > const & phaseTypes,
                   arrayView1d< real64 const > const & componentMolarWeight,
                   bool const useMass,
                   PhaseProp::ViewType phaseFraction,
                   PhaseProp::ViewType phaseDensity,
                   PhaseProp::ViewType phaseMassDensity,
                   PhaseProp::ViewType phaseViscosity,
                   PhaseProp::ViewType phaseEnthalpy,
                   PhaseProp::ViewType phaseInternalEnergy,
                   PhaseComp::ViewType phaseCompFraction,
                   FluidProp::ViewType totalDensity );

    pvt::MultiphaseSystem & m_fluid;

    arrayView1d< pvt::PHASE_TYPE > m_phaseTypes;
  };

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper();

protected:

  virtual void postProcessInput() override;

  virtual void initializePostSubGroups() override;

private:

  void createFluid();

  /// PVTPackage fluid object
  std::unique_ptr< pvt::MultiphaseSystem > m_fluid{};

  /// PVTPackage phase labels
  array1d< pvt::PHASE_TYPE > m_phaseTypes;

  // names of equations of state to use for each phase
  string_array m_equationsOfState;

  // standard EOS component input
  array1d< real64 > m_componentCriticalPressure;
  array1d< real64 > m_componentCriticalTemperature;
  array1d< real64 > m_componentAcentricFactor;
  array1d< real64 > m_componentVolumeShift;
  array2d< real64 > m_componentBinaryCoeff;

};

GEOSX_HOST_DEVICE
inline void
CompositionalMultiphaseFluid::KernelWrapper::
  compute( real64 const pressure,
           real64 const temperature,
           arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition,
           arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseFrac,
           arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseDens,
           arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseMassDens,
           arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseVisc,
           arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseEnthalpy,
           arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseInternalEnergy,
           arraySlice2d< real64, multifluid::USD_PHASE_COMP - 2 > const & phaseCompFrac,
           real64 & totalDens ) const
{
  GEOSX_UNUSED_VAR( phaseEnthalpy, phaseInternalEnergy );
#if defined(__CUDA_ARCH__)
  GEOSX_ERROR( "This function cannot be used on GPU" );
#else
  localIndex const NC = m_componentMolarWeight.size();
  localIndex const NP = m_phaseTypes.size();

  // 1. Convert input mass fractions to mole fractions and keep derivatives
  std::vector< double > compMoleFrac( NC );

  if( m_useMass )
  {
    real64 totalMolality = 0.0;
    for( localIndex ic = 0; ic < NC; ++ic )
    {
      compMoleFrac[ic] = composition[ic] / m_componentMolarWeight[ic]; // this is molality (units of mole/mass)
      totalMolality += compMoleFrac[ic];
    }

    real64 const totalMolalityInv = 1.0 / totalMolality;
    for( localIndex ic = 0; ic < NC; ++ic )
    {
      compMoleFrac[ic] *= totalMolalityInv;
    }
  }
  else
  {
    for( localIndex ic = 0; ic < NC; ++ic )
    {
      compMoleFrac[ic] = composition[ic];
    }
  }

  // 2. Trigger PVTPackage compute and get back phase split
  m_fluid.Update( pressure, temperature, compMoleFrac );

  GEOSX_WARNING_IF( !m_fluid.hasSucceeded(),
                    "Phase equilibrium calculations not converged" );

  pvt::MultiphaseSystemProperties const & props = m_fluid.getMultiphaseSystemProperties();

  // 3. Extract phase split and phase properties from PVTPackage
  for( localIndex ip = 0; ip < NP; ++ip )
  {
    pvt::PHASE_TYPE const & phaseType = m_phaseTypes[ip];

    auto const & frac = props.getPhaseMoleFraction( phaseType );
    auto const & comp = props.getMoleComposition( phaseType );
    auto const & dens = m_useMass ? props.getMassDensity( phaseType ) : props.getMoleDensity( phaseType );
    auto const & massDens = props.getMassDensity( phaseType );

    phaseFrac[ip] = frac.value;
    phaseDens[ip] = dens.value;
    phaseMassDens[ip] = massDens.value;
    phaseVisc[ip] = 0.001;   // TODO
    for( localIndex jc = 0; jc < NC; ++jc )
    {
      phaseCompFrac[ip][jc] = comp.value[jc];
    }
  }

  // 4. if mass variables used instead of molar, perform the conversion
  if( m_useMass )
  {
    // 4.1. Convert phase fractions (requires two passes)
    real64 totalMass{};

    // 4.1.1. Compute mass of each phase and total mass (on a 1-mole basis)
    for( localIndex ip = 0; ip < NP; ++ip )
    {
      auto const & phaseMW = props.getMolecularWeight( m_phaseTypes[ip] );
      phaseFrac[ip] *= phaseMW.value;
      totalMass += phaseFrac[ip];
    }

    // 4.1.2. Normalize to get mass fractions
    real64 const totalMassInv = 1.0 / totalMass;
    for( localIndex ip = 0; ip < NP; ++ip )
    {
      phaseFrac[ip] *= totalMassInv;
    }

    // 4.2. Convert phase compositions
    for( localIndex ip = 0; ip < NP; ++ip )
    {
      auto const & phaseMW = props.getMolecularWeight( m_phaseTypes[ip] );
      real64 const phaseMWInv = 1.0 / phaseMW.value;

      for( localIndex ic = 0; ic < NC; ++ic )
      {
        phaseCompFrac[ip][ic] = phaseCompFrac[ip][ic] * m_componentMolarWeight[ic] * phaseMWInv;
      }
    }
  }

  // 5. Compute total fluid mass/molar density
  {
    totalDens = 0.0;

    // 5.1. Sum mass/molar fraction/density ratio over all phases to get the inverse of density
    for( localIndex ip = 0; ip < NP; ++ip )
    {
      totalDens += phaseFrac[ip] / phaseDens[ip];
    }

    // 5.2. Invert the previous quantity to get actual density
    totalDens = 1.0 / totalDens;
  }
#endif
}

GEOSX_HOST_DEVICE
inline void
CompositionalMultiphaseFluid::KernelWrapper::
  update( localIndex const k,
          localIndex const q,
          real64 const pressure,
          real64 const temperature,
          arraySlice1d< geosx::real64 const, compflow::USD_COMP - 1 > const & composition ) const
{
  compute( pressure,
           temperature,
           composition,
           m_phaseFraction( k, q ),
           m_phaseDensity( k, q ),
           m_phaseMassDensity( k, q ),
           m_phaseViscosity( k, q ),
           m_phaseEnthalpy( k, q ),
           m_phaseInternalEnergy( k, q ),
           m_phaseCompFraction( k, q ),
           m_totalDensity( k, q ) );
}

GEOSX_HOST_DEVICE
inline void
CompositionalMultiphaseFluid::KernelWrapper::
  compute( real64 const pressure,
           real64 const temperature,
           arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition,
           PhaseProp::SliceType const phaseFraction,
           PhaseProp::SliceType const phaseDensity,
           PhaseProp::SliceType const phaseMassDensity,
           PhaseProp::SliceType const phaseViscosity,
           PhaseProp::SliceType const phaseEnthalpy,
           PhaseProp::SliceType const phaseInternalEnergy,
           PhaseComp::SliceType const phaseCompFraction,
           FluidProp::SliceType const totalDensity ) const
{
  GEOSX_UNUSED_VAR( phaseEnthalpy, phaseInternalEnergy );
#if defined(__CUDA_ARCH__)
  GEOSX_ERROR( "This function cannot be used on GPU" );
#else

  localIndex constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS;
  localIndex const NC = numComponents();
  localIndex const NP = numPhases();

  // 1. Convert input mass fractions to mole fractions and keep derivatives

  std::vector< double > compMoleFrac( NC );
  stackArray2d< real64, maxNumComp * maxNumComp > dCompMoleFrac_dCompMassFrac( NC, NC );

  if( m_useMass )
  {
    dCompMoleFrac_dCompMassFrac.resize( NC, NC );
    dCompMoleFrac_dCompMassFrac.zero();

    real64 totalMolality = 0.0;
    for( localIndex ic = 0; ic < NC; ++ic )
    {
      real64 const mwInv = 1.0 / m_componentMolarWeight[ic];
      compMoleFrac[ic] = composition[ic] * mwInv; // this is molality (units of mole/mass)
      dCompMoleFrac_dCompMassFrac[ic][ic] = mwInv;
      totalMolality += compMoleFrac[ic];
    }

    real64 const totalMolalityInv = 1.0 / totalMolality;
    for( localIndex ic = 0; ic < NC; ++ic )
    {
      compMoleFrac[ic] *= totalMolalityInv;

      for( localIndex jc = 0; jc < NC; ++jc )
      {
        dCompMoleFrac_dCompMassFrac[ic][jc] -= compMoleFrac[ic] / m_componentMolarWeight[jc];
        dCompMoleFrac_dCompMassFrac[ic][jc] *= totalMolalityInv;
      }
    }
  }
  else
  {
    for( localIndex ic = 0; ic < NC; ++ic )
    {
      compMoleFrac[ic] = composition[ic];
    }
  }

  // 2. Trigger PVTPackage compute and get back phase split
  m_fluid.Update( pressure, temperature, compMoleFrac );

  GEOSX_WARNING_IF( !m_fluid.hasSucceeded(),
                    "Phase equilibrium calculations not converged" );

  pvt::MultiphaseSystemProperties const & props = m_fluid.getMultiphaseSystemProperties();

  // 3. Extract phase split, phase properties and derivatives from PVTPackage
  for( localIndex ip = 0; ip < NP; ++ip )
  {
    pvt::PHASE_TYPE const & phaseType = m_phaseTypes[ip];

    auto const & frac = props.getPhaseMoleFraction( phaseType );
    auto const & comp = props.getMoleComposition( phaseType );
    auto const & dens = m_useMass ? props.getMassDensity( phaseType ) : props.getMoleDensity( phaseType );
    auto const & massDens = props.getMassDensity( phaseType );

    phaseFraction.value[ip] = frac.value;
    phaseFraction.dPres[ip] = frac.dP;
    phaseFraction.dTemp[ip] = frac.dT;

    phaseDensity.value[ip] = dens.value;
    phaseDensity.dPres[ip] = dens.dP;
    phaseDensity.dTemp[ip] = dens.dT;

    phaseMassDensity.value[ip] = massDens.value;
    phaseMassDensity.dPres[ip] = massDens.dP;
    phaseMassDensity.dTemp[ip] = massDens.dT;

    // TODO
    phaseViscosity.value[ip] = 0.001;
    phaseViscosity.dPres[ip] = 0.0;
    phaseViscosity.dTemp[ip] = 0.0;

    for( localIndex jc = 0; jc < NC; ++jc )
    {
      phaseFraction.dComp[ip][jc] = frac.dz[jc];
      phaseDensity.dComp[ip][jc] = dens.dz[jc];
      phaseMassDensity.dComp[ip][ip] = massDens.dz[jc];
      phaseViscosity.dComp[ip][jc] = 0.0; // TODO

      phaseCompFraction.value[ip][jc] = comp.value[jc];
      phaseCompFraction.dPres[ip][jc] = comp.dP[jc];
      phaseCompFraction.dTemp[ip][jc] = comp.dT[jc];

      for( localIndex ic = 0; ic < NC; ++ic )
      {
        phaseCompFraction.dComp[ip][ic][jc] = comp.dz[ic][jc];
      }
    }
  }

  // 4. if mass variables used instead of molar, perform the conversion
  if( m_useMass )
  {
    // 4.1. Convert phase fractions (requires two passes)
    real64 totalMass{};
    real64 dTotalMass_dP{};
    real64 dTotalMass_dT{};
    real64 dTotalMass_dC[maxNumComp]{};

    // 4.1.1. Compute mass of each phase and total mass (on a 1-mole basis)
    for( localIndex ip = 0; ip < NP; ++ip )
    {
      auto const & phaseMW = props.getMolecularWeight( m_phaseTypes[ip] );
      real64 const nu = phaseFraction.value[ip];

      phaseFraction.value[ip] *= phaseMW.value;
      phaseFraction.dPres[ip] = phaseFraction.dPres[ip] * phaseMW.value + nu * phaseMW.dP;
      phaseFraction.dTemp[ip] = phaseFraction.dTemp[ip] * phaseMW.value + nu * phaseMW.dT;

      totalMass += phaseFraction.value[ip];
      dTotalMass_dP += phaseFraction.dPres[ip];
      dTotalMass_dT += phaseFraction.dTemp[ip];

      for( localIndex jc = 0; jc < NC; ++jc )
      {
        phaseFraction.dComp[ip][jc] = phaseFraction.dComp[ip][jc] * phaseMW.value + nu * phaseMW.dz[jc];
        dTotalMass_dC[jc] += phaseFraction.dComp[ip][jc];
      }
    }

    // 4.1.2. Normalize to get mass fractions
    real64 const totalMassInv = 1.0 / totalMass;
    for( localIndex ip = 0; ip < NP; ++ip )
    {
      phaseFraction.value[ip] *= totalMassInv;
      phaseFraction.dPres[ip] = ( phaseFraction.dPres[ip] - phaseFraction.value[ip] * dTotalMass_dP ) * totalMassInv;
      phaseFraction.dTemp[ip] = ( phaseFraction.dTemp[ip] - phaseFraction.value[ip] * dTotalMass_dT ) * totalMassInv;

      for( localIndex jc = 0; jc < NC; ++jc )
      {
        phaseFraction.dComp[ip][jc] = ( phaseFraction.dComp[ip][jc] - phaseFraction.value[ip] * dTotalMass_dC[jc] ) * totalMassInv;
      }
    }

    // 4.2. Convert phase compositions
    for( localIndex ip = 0; ip < NP; ++ip )
    {
      auto const & phaseMW = props.getMolecularWeight( m_phaseTypes[ip] );
      real64 const phaseMWInv = 1.0 / phaseMW.value;

      for( localIndex ic = 0; ic < NC; ++ic )
      {

        real64 const compMW = m_componentMolarWeight[ic];

        phaseCompFraction.value[ip][ic] = phaseCompFraction.value[ip][ic] * compMW * phaseMWInv;
        phaseCompFraction.dPres[ip][ic] =
          ( phaseCompFraction.dPres[ip][ic] * compMW - phaseCompFraction.value[ip][ic] * phaseMW.dP ) * phaseMWInv;
        phaseCompFraction.dTemp[ip][ic] =
          ( phaseCompFraction.dTemp[ip][ic] * compMW - phaseCompFraction.value[ip][ic] * phaseMW.dT ) * phaseMWInv;

        for( localIndex jc = 0; jc < NC; ++jc )
        {
          phaseCompFraction.dComp[ip][ic][jc] =
            ( phaseCompFraction.dComp[ip][ic][jc] * compMW - phaseCompFraction.value[ip][ic] * phaseMW.dz[jc] ) * phaseMWInv;
        }
      }
    }

    // 4.3. Update derivatives w.r.t. mole fractions to derivatives w.r.t mass fractions
    array1d< real64 > work( NC );
    for( localIndex ip = 0; ip < NP; ++ip )
    {
      applyChainRuleInPlace( NC, dCompMoleFrac_dCompMassFrac, phaseFraction.dComp[ip], work );
      applyChainRuleInPlace( NC, dCompMoleFrac_dCompMassFrac, phaseDensity.dComp[ip], work );

      for( localIndex ic = 0; ic < NC; ++ic )
      {
        applyChainRuleInPlace( NC, dCompMoleFrac_dCompMassFrac, phaseCompFraction.dComp[ip][ic], work );
      }
    }
  }

  // 5. Compute total fluid mass/molar density and derivatives
  {
    totalDensity.value = 0.0;
    totalDensity.dPres = 0.0;
    totalDensity.dTemp = 0.0;
    for( localIndex jc = 0; jc < NC; ++jc )
    {
      totalDensity.dComp[jc] = 0.0;
    }

    // 5.1. Sum mass/molar fraction/density ratio over all phases to get the inverse of density
    for( localIndex ip = 0; ip < NP; ++ip )
    {
      real64 const densInv = 1.0 / phaseDensity.value[ip];
      real64 const value = phaseFraction.value[ip] * densInv;

      totalDensity.value += value;
      totalDensity.dPres += ( phaseFraction.dPres[ip] - value * phaseDensity.dPres[ip] ) * densInv;
      totalDensity.dTemp += ( phaseFraction.dTemp[ip] - value * phaseDensity.dTemp[ip] ) * densInv;

      for( localIndex jc = 0; jc < NC; ++jc )
      {
        totalDensity.dComp[jc] += ( phaseFraction.dComp[ip][jc] - value * phaseDensity.dComp[ip][jc] ) * densInv;
      }
    }

    // 5.2. Invert the previous quantity to get actual density
    totalDensity.value = 1.0 / totalDensity.value;
    real64 const minusDens2 = -totalDensity.value * totalDensity.value;
    totalDensity.dPres *= minusDens2;
    totalDensity.dTemp *= minusDens2;

    for( localIndex jc = 0; jc < NC; ++jc )
    {
      totalDensity.dComp[jc] *= minusDens2;
    }
  }
#endif
}

} /* namespace constitutive */

} /* namespace geosx */

#endif //GEOSX_CONSTITUTIVE_FLUID_COMPOSITIONALMULTIPHASEFLUID_HPP_
