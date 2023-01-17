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
 * @file ExampleFluidModel.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_EXAMPLEFLUIDMODEL_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_EXAMPLEFLUIDMODEL_HPP_

#include "constitutive/fluid/MultiFluidBase.hpp"
#include "constitutive/fluid/MultiFluidUtils.hpp"

namespace geosx
{
namespace constitutive
{

class ExampleFluidModel : public MultiFluidBase
{
public:

  using exec_policy = serialPolicy;

  ExampleFluidModel( string const & name, Group * const parent );

  virtual std::unique_ptr< ConstitutiveBase >
  deliverClone( string const & name,
                Group * const parent ) const override;

  static string catalogName() { return "ExampleFluidModel"; }

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
   * @brief Kernel wrapper class for ExampleFluidModel.
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

    friend class ExampleFluidModel;

    KernelWrapper( arrayView1d< real64 const > const & componentMolarWeight,
                   bool const useMass,
                   PhaseProp::ViewType phaseFraction,
                   PhaseProp::ViewType phaseDensity,
                   PhaseProp::ViewType phaseMassDensity,
                   PhaseProp::ViewType phaseViscosity,
                   PhaseProp::ViewType phaseEnthalpy,
                   PhaseProp::ViewType phaseInternalEnergy,
                   PhaseComp::ViewType phaseCompFraction,
                   FluidProp::ViewType totalDensity );

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
ExampleFluidModel::KernelWrapper::
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

  integer constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS;
  integer constexpr maxNumPhase = MultiFluidBase::MAX_NUM_PHASES;
  integer const numComp = numComponents();
  integer const numPhase = numPhases();

  // 1. Convert input mass fractions to mole fractions and keep derivatives

  stackArray1d< real64, maxNumComp > compMoleFrac( numComp );

  if( m_useMass )
  {
    convertToMoleFractions< maxNumComp >( composition,
                                          compMoleFrac );
  }
  else
  {
    for( integer ic = 0; ic < numComp; ++ic )
    {
      compMoleFrac[ic] = composition[ic];
    }
  }

  // 2. Perform PVT calculations here

  GEOSX_UNUSED_VAR( pressure, temperature, composition, phaseFrac, phaseDens, phaseMassDens, phaseVisc, phaseEnthalpy, phaseInternalEnergy, phaseCompFrac, totalDens );

  // Compute phase fractions and phase component fractions

  // Compute phase densities, total density and phase viscosities

  // 3. if mass variables used instead of molar, perform the conversion

  if( m_useMass )
  {

    real64 phaseMolecularWeight[maxNumPhase]{};
    for( integer ip = 0; ip < numPhase; ++ip )
    {
      phaseMolecularWeight[ip] = 0.0; // get phase molecular weight here
    }

    // convert mole fractions to mass fractions
    convertToMassFractions< maxNumComp >( phaseMolecularWeight,
                                          phaseFrac,
                                          phaseCompFrac );

  }

  // 4. Compute total fluid mass/molar density

  computeTotalDensity< maxNumComp, maxNumPhase >( phaseFrac,
                                                  phaseDens,
                                                  totalDens );

#endif
}

GEOSX_HOST_DEVICE
inline void
ExampleFluidModel::KernelWrapper::
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

  using Deriv = multifluid::DerivativeOffset;

  integer constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS;
  integer constexpr maxNumPhase = MultiFluidBase::MAX_NUM_PHASES;
  integer const numComp = numComponents();
  integer const numPhase = numPhases();

  // 1. Convert input mass fractions to mole fractions and keep derivatives

  stackArray1d< real64, maxNumComp  > compMoleFrac( numComp );
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

  // 2. Perform PVT calculations here

  GEOSX_UNUSED_VAR( pressure, temperature, composition, phaseFraction, phaseDensity, phaseMassDensity, phaseViscosity, phaseEnthalpy, phaseInternalEnergy, phaseCompFraction, totalDensity );

  // Compute phase fractions and phase component fractions (+ derivatives)

  // Compute phase densities, total density and phase viscosities (+ derivatives)


  // 3. if mass variables used instead of molar, perform the conversion
  if( m_useMass )
  {

    real64 phaseMolecularWeight[maxNumPhase]{};
    real64 dPhaseMolecularWeight[maxNumPhase][maxNumComp+2]{};

    for( integer ip = 0; ip < numPhase; ++ip )
    {
      phaseMolecularWeight[ip] = 0.0; // get phase molecular weight here (+ derivatives)
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

  // 4. Compute total fluid mass/molar density and derivatives

  computeTotalDensity( phaseFraction,
                       phaseDensity,
                       totalDensity );

#endif
}

GEOSX_HOST_DEVICE
inline void
ExampleFluidModel::KernelWrapper::
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

} /* namespace constitutive */

} /* namespace geosx */

#endif //GEOSX_CONSTITUTIVE_FLUID_EXAMPLEFLUIDMODEL_HPP_
