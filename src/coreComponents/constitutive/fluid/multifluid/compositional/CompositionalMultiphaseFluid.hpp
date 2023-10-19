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

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_COMPOSITIONALMULTIPHASEFLUID_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_COMPOSITIONALMULTIPHASEFLUID_HPP_

#include "constitutive/fluid/multifluid/compositional/CompositionalMultiphaseFluidParameters.hpp"

#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "constitutive/fluid/multifluid/MultiFluidUtils.hpp"

namespace geos
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

  static string catalogName() { return "CompositionalMultiphaseFluidGeos"; }

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

private:
  class IFluid
  {
public:
    virtual ~IFluid() = default;
  };

public:

  /**
   * @brief Kernel wrapper class for CompositionalMultiphaseFluid.
   */
  class KernelWrapper final : public MultiFluidBase::KernelWrapper
  {
public:

    GEOS_HOST_DEVICE
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

    GEOS_HOST_DEVICE
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

    GEOS_HOST_DEVICE
    virtual void update( localIndex const k,
                         localIndex const q,
                         real64 const pressure,
                         real64 const temperature,
                         arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition ) const override;

private:

    friend class CompositionalMultiphaseFluid;

    KernelWrapper( CompositionalMultiphaseFluid::IFluid & fluid,
                   arrayView1d< PhaseType > const & phaseTypes,
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

    CompositionalMultiphaseFluid::IFluid & m_fluid;
    arrayView1d< PhaseType > m_phaseTypes;
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

  PhaseType getPhaseType( string const & name ) const;

  /// Fluid object
  std::unique_ptr< IFluid > m_fluid{};

  /// Phase labels
  array1d< PhaseType > m_phaseTypes;

  /// Index of the water phase
  integer m_aqueousPhaseIndex{-1};

  /// Equation of state labels
  array1d< EquationOfStateType > m_equationsOfState;

  // standard EOS component input
  array1d< real64 > m_componentCriticalPressure;
  array1d< real64 > m_componentCriticalTemperature;
  array1d< real64 > m_componentAcentricFactor;
  array1d< real64 > m_componentVolumeShift;
  array2d< real64 > m_componentBinaryCoeff;

  // Currently not implemented so uses constant values
  static constexpr real64 MOLECULAR_WEIGHT[] = { 22.0, 33.0 };
  static constexpr real64 DENSITY[] = { 1.0, 1.0 };
  static constexpr real64 VISCOSITY[] = { 1.0, 1.0  };
};

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
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
  GEOS_UNUSED_VAR( phaseEnthalpy, phaseInternalEnergy );
  GEOS_UNUSED_VAR( pressure, temperature );

  integer constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS;
  integer constexpr maxNumPhase = MultiFluidBase::MAX_NUM_PHASES;
  integer const numComp = numComponents();
  integer const numPhase = numPhases();

  // 1. Convert input mass fractions to mole fractions and keep derivatives

  stackArray1d< real64, maxNumComp > compMoleFrac( numComp );
  if( m_useMass )
  {
    convertToMoleFractions< maxNumComp >( composition, compMoleFrac );
  }
  else
  {
    for( integer ic = 0; ic < numComp; ++ic )
    {
      compMoleFrac[ic] = composition[ic];
    }
  }

  // 2. Calculate values

  for( integer ip = 0; ip < numPhase; ++ip )
  {
    phaseFrac[ip] = 0.0;        // TODO
    phaseDens[ip] = 0.0;        // TODO
    phaseMassDens[ip] = 0.0;    // TODO
    phaseVisc[ip] = 0.0;        // TODO
    for( integer jc = 0; jc < numComp; ++jc )
    {
      phaseCompFrac[ip][jc] = composition[jc];  // TODO
    }
  }

  // 4. if mass variables used instead of molar, perform the conversion

  if( m_useMass )
  {

    // unfortunately here, we have to copy the molecular weight coming from PVT package...
    real64 phaseMolecularWeight[maxNumPhase]{};
    for( integer ip = 0; ip < numPhase; ++ip )
    {
      phaseMolecularWeight[ip] = 0.0;
    }

    // convert mole fractions to mass fractions
    convertToMassFractions< maxNumComp >( phaseMolecularWeight,
                                          phaseFrac,
                                          phaseCompFrac );

  }

  // 5. Compute total fluid mass/molar density

  computeTotalDensity< maxNumComp, maxNumPhase >( phaseFrac,
                                                  phaseDens,
                                                  totalDens );
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
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
  GEOS_UNUSED_VAR( phaseEnthalpy, phaseInternalEnergy );
  GEOS_UNUSED_VAR( pressure, temperature );

  using Deriv = multifluid::DerivativeOffset;

  integer constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS;
  integer constexpr maxNumPhase = MultiFluidBase::MAX_NUM_PHASES;
  integer const numComp = numComponents();
  integer const numPhase = numPhases();

  // 1. Convert input mass fractions to mole fractions and keep derivatives

  stackArray1d< real64, maxNumComp > compMoleFrac( numComp );
  real64 dCompMoleFrac_dCompMassFrac[maxNumComp][maxNumComp]{};

  if( m_useMass )
  {
    // convert mass fractions to mole fractions
    convertToMoleFractions( composition,
                            compMoleFrac,
                            dCompMoleFrac_dCompMassFrac );
  }
  else
  {
    for( integer ic = 0; ic < numComp; ++ic )
    {
      compMoleFrac[ic] = composition[ic];
    }
  }

  // 2. Calculate values


  // 3. Extract phase split, phase properties and derivatives from result

  for( integer ip = 0; ip < numPhase; ++ip )
  {
    phaseFraction.value[ip] = 0.0;
    phaseFraction.derivs[ip][Deriv::dP] = 0.0;
    phaseFraction.derivs[ip][Deriv::dT] = 0.0;

    phaseDensity.value[ip] = 0.0;
    phaseDensity.derivs[ip][Deriv::dP] = 0.0;
    phaseDensity.derivs[ip][Deriv::dT] = 0.0;

    phaseMassDensity.value[ip] = 0.0;
    phaseMassDensity.derivs[ip][Deriv::dP] = 0.0;
    phaseMassDensity.derivs[ip][Deriv::dT] = 0.0;

    phaseViscosity.value[ip] = 0.0;
    phaseViscosity.derivs[ip][Deriv::dP] = 0.0;
    phaseViscosity.derivs[ip][Deriv::dT] = 0.0;

    for( integer jc = 0; jc < numComp; ++jc )
    {
      phaseFraction.derivs[ip][Deriv::dC+jc] = 0.0;
      phaseDensity.derivs[ip][Deriv::dC+jc] = 0.0;
      phaseMassDensity.derivs[ip][Deriv::dC+jc] = 0.0;
      phaseViscosity.derivs[ip][Deriv::dC+jc] = 0.0;

      phaseCompFraction.value[ip][jc] = 0.0;
      phaseCompFraction.derivs[ip][jc][Deriv::dP] = 0.0;
      phaseCompFraction.derivs[ip][jc][Deriv::dT] = 0.0;

      for( integer ic = 0; ic < numComp; ++ic )
      {
        phaseCompFraction.derivs[ip][ic][Deriv::dC+jc] = 0.0;
      }
    }
  }

  // 4. if mass variables used instead of molar, perform the conversion
  if( m_useMass )
  {

    // unfortunately here, we have to copy the molecular weight coming from PVT package...
    real64 phaseMolecularWeight[maxNumPhase]{};
    real64 dPhaseMolecularWeight[maxNumPhase][maxNumComp+2]{};

    for( integer ip = 0; ip < numPhase; ++ip )
    {
      phaseMolecularWeight[ip] = 1.0;
      dPhaseMolecularWeight[ip][Deriv::dP] = 0.0;
      dPhaseMolecularWeight[ip][Deriv::dT] = 0.0;
      for( integer ic = 0; ic < numComp; ++ic )
      {
        dPhaseMolecularWeight[ip][Deriv::dC+ic] = 0.0;
      }
    }

    convertToMassFractions( dCompMoleFrac_dCompMassFrac,
                            phaseMolecularWeight,
                            dPhaseMolecularWeight,
                            phaseFraction,
                            phaseCompFraction,
                            phaseDensity.derivs,
                            phaseViscosity.derivs,
                            phaseEnthalpy.derivs,
                            phaseInternalEnergy.derivs );
  }

  // 5. Compute total fluid mass/molar density and derivatives

  computeTotalDensity( phaseFraction,
                       phaseDensity,
                       totalDensity );
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
CompositionalMultiphaseFluid::KernelWrapper::
  update( localIndex const k,
          localIndex const q,
          real64 const pressure,
          real64 const temperature,
          arraySlice1d< geos::real64 const, compflow::USD_COMP - 1 > const & composition ) const
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

} /* namespace constitutive */

} /* namespace geos */

#endif //GEOS_CONSTITUTIVE_FLUID_COMPOSITIONALMULTIPHASEFLUID_HPP_
