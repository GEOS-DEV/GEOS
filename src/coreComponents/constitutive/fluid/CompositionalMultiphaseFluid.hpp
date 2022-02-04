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
           arraySlice2d< real64, multifluid::USD_PHASE_COMP - 2 > const & phaseCompFrac,
           real64 & totalDens ) const
{
#if defined(__CUDA_ARCH__)
  GEOSX_ERROR( "This function cannot be used on GPU" );
#else

  integer constexpr maxNumPhase = MultiFluidBase::MAX_NUM_PHASES;
  integer const numComp = numComponents();
  integer const numPhase = numPhases();

  // 1. Convert input mass fractions to mole fractions and keep derivatives

  std::vector< double > compMoleFrac( numComp );

  if( m_useMass )
  {
    convertToMoleFractions( composition,
                            m_componentMolarWeight,
                            compMoleFrac );
  }
  else
  {
    for( integer ic = 0; ic < numComp; ++ic )
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

  for( integer ip = 0; ip < numPhase; ++ip )
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
    for( integer jc = 0; jc < numComp; ++jc )
    {
      phaseCompFrac[ip][jc] = comp.value[jc];
    }
  }

  // 4. if mass variables used instead of molar, perform the conversion

  if( m_useMass )
  {

    // unfortunately here, we have to copy the molecular weight coming from PVT package...
    real64 phaseMolecularWeight[maxNumPhase]{};
    for( integer ip = 0; ip < numPhase; ++ip )
    {
      phaseMolecularWeight[ip] = props.getMolecularWeight( m_phaseTypes[ip] ).value;
    }

    // convert mole fractions to mass fractions
    convertToMassFractions< maxNumPhase >( phaseMolecularWeight,
                                           phaseFrac,
                                           phaseCompFrac );
  }

  // 5. Compute total fluid mass/molar density

  computeTotalDensity( phaseFrac,
                       phaseDens,
                       totalDens );

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
           PhaseComp::SliceType const phaseCompFraction,
           FluidProp::SliceType const totalDensity ) const
{
#if defined(__CUDA_ARCH__)
  GEOSX_ERROR( "This function cannot be used on GPU" );
#else

  integer constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS;
  integer constexpr maxNumPhase = MultiFluidBase::MAX_NUM_PHASES;
  integer const numComp = numComponents();
  integer const numPhase = numPhases();

  // 1. Convert input mass fractions to mole fractions and keep derivatives

  std::vector< double > compMoleFrac( numComp );
  real64 dCompMoleFrac_dCompMassFrac[maxNumComp][maxNumComp]{};

  if( m_useMass )
  {
    // convert mass fractions to mole fractions
    convertToMoleFractions( composition,
                            m_componentMolarWeight,
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

  // 2. Trigger PVTPackage compute and get back phase split

  m_fluid.Update( pressure, temperature, compMoleFrac );

  GEOSX_WARNING_IF( !m_fluid.hasSucceeded(),
                    "Phase equilibrium calculations not converged" );

  pvt::MultiphaseSystemProperties const & props = m_fluid.getMultiphaseSystemProperties();

  // 3. Extract phase split, phase properties and derivatives from PVTPackage

  for( integer ip = 0; ip < numPhase; ++ip )
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

    for( integer jc = 0; jc < numComp; ++jc )
    {
      phaseFraction.dComp[ip][jc] = frac.dz[jc];
      phaseDensity.dComp[ip][jc] = dens.dz[jc];
      phaseMassDensity.dComp[ip][ip] = massDens.dz[jc];
      phaseViscosity.dComp[ip][jc] = 0.0; // TODO

      phaseCompFraction.value[ip][jc] = comp.value[jc];
      phaseCompFraction.dPres[ip][jc] = comp.dP[jc];
      phaseCompFraction.dTemp[ip][jc] = comp.dT[jc];

      for( integer ic = 0; ic < numComp; ++ic )
      {
        phaseCompFraction.dComp[ip][ic][jc] = comp.dz[ic][jc];
      }
    }
  }

  // 4. if mass variables used instead of molar, perform the conversion
  if( m_useMass )
  {

    // unfortunately here, we have to copy the molecular weight coming from PVT package...
    real64 phaseMolecularWeight[maxNumPhase]{};
    real64 dPhaseMolecularWeight_dPres[maxNumPhase]{};
    real64 dPhaseMolecularWeight_dTemp[maxNumPhase]{};
    real64 dPhaseMolecularWeight_dComp[maxNumPhase][maxNumComp]{};

    for( integer ip = 0; ip < numPhase; ++ip )
    {
      phaseMolecularWeight[ip] = props.getMolecularWeight( m_phaseTypes[ip] ).value;
      dPhaseMolecularWeight_dPres[ip] = props.getMolecularWeight( m_phaseTypes[ip] ).dP;
      dPhaseMolecularWeight_dTemp[ip] = props.getMolecularWeight( m_phaseTypes[ip] ).dT;
      for( integer ic = 0; ic < numComp; ++ic )
      {
        dPhaseMolecularWeight_dComp[ip][ic] = props.getMolecularWeight( m_phaseTypes[ip] ).dz[ic];
      }
    }

    convertToMassFractions< maxNumComp, maxNumPhase >( dCompMoleFrac_dCompMassFrac,
                                                       phaseMolecularWeight,
                                                       dPhaseMolecularWeight_dPres,
                                                       dPhaseMolecularWeight_dTemp,
                                                       dPhaseMolecularWeight_dComp,
                                                       phaseFraction,
                                                       phaseCompFraction,
                                                       phaseDensity.dComp,
                                                       phaseViscosity.dComp );
  }

  // 5. Compute total fluid mass/molar density and derivatives

  computeTotalDensity( phaseFraction,
                       phaseDensity,
                       totalDensity );

#endif
}

} /* namespace constitutive */

} /* namespace geosx */

#endif //GEOSX_CONSTITUTIVE_FLUID_COMPOSITIONALMULTIPHASEFLUID_HPP_
