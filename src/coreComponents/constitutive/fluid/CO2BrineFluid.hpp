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
 * @file CO2BrineFluid.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_CO2BRINEFLUID_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_CO2BRINEFLUID_HPP_

#include "constitutive/fluid/MultiFluidBase.hpp"
#include "constitutive/fluid/MultiFluidUtils.hpp"
#include "constitutive/fluid/PVTFunctions/BrineEnthalpy.hpp"
#include "constitutive/fluid/PVTFunctions/BrineInternalEnergy.hpp"
#include "constitutive/fluid/PVTFunctions/CO2Enthalpy.hpp"
#include "constitutive/fluid/PVTFunctions/CO2InternalEnergy.hpp"
#include "constitutive/fluid/PVTFunctions/CO2Solubility.hpp"
#include "constitutive/fluid/PVTFunctions/EzrokhiBrineDensity.hpp"
#include "constitutive/fluid/PVTFunctions/EzrokhiBrineViscosity.hpp"
#include "constitutive/fluid/PVTFunctions/FenghourCO2Viscosity.hpp"
#include "constitutive/fluid/PVTFunctions/NoOpPVTFunction.hpp"
#include "constitutive/fluid/PVTFunctions/PhillipsBrineDensity.hpp"
#include "constitutive/fluid/PVTFunctions/PhillipsBrineViscosity.hpp"
#include "constitutive/fluid/PVTFunctions/SpanWagnerCO2Density.hpp"


#include <memory>

namespace geosx
{

namespace constitutive
{

template< typename P1DENS, typename P1VISC, typename P1ENTH, typename P1INTENERGY,
          typename P2DENS, typename P2VISC, typename P2ENTH, typename P2INTENERGY,
          typename FLASH >
class CO2BrineFluid : public MultiFluidBase
{
public:

  using exec_policy = parallelDevicePolicy<>;

  CO2BrineFluid( string const & name,
                 Group * const parent );

  virtual std::unique_ptr< ConstitutiveBase >
  deliverClone( string const & name,
                Group * const parent ) const override;

  static string catalogName();

  virtual string getCatalogName() const override { return catalogName(); }

  /**
   * @brief Kernel wrapper class for CO2BrineFluid.
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
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseEnthalpy,
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseInternalEnergy,
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

    friend class CO2BrineFluid;

    KernelWrapper( integer p1Index,
                   integer p2Index,
                   P1DENS const & p1Density,
                   P1VISC const & p1Viscosity,
                   P1ENTH const & p1Enthalpy,
                   P1INTENERGY const & p1IntEnergy,
                   P2DENS const & p2Density,
                   P2VISC const & p2Viscosity,
                   P2ENTH const & p2Enthalpy,
                   P2INTENERGY const & p2IntEnergy,
                   FLASH const & flash,
                   arrayView1d< real64 const > componentMolarWeight,
                   bool const useMass,
                   PhaseProp::ViewType phaseFraction,
                   PhaseProp::ViewType phaseDensity,
                   PhaseProp::ViewType phaseMassDensity,
                   PhaseProp::ViewType phaseViscosity,
                   PhaseProp::ViewType phaseEnthalpy,
                   PhaseProp::ViewType phaseInternalEnergy,
                   PhaseComp::ViewType phaseCompFraction,
                   FluidProp::ViewType totalDensity );

    /// Index of the liquid phase
    integer m_p1Index;

    /// Index of the gas phase
    integer m_p2Index;


    // Brine constitutive kernel wrappers

    /// Kernel wrapper for brine density updates
    typename P1DENS::KernelWrapper m_p1Density;

    /// Kernel wrapper for brine viscosity updates
    typename P1VISC::KernelWrapper m_p1Viscosity;

    /// Kernel wrapper for brine enthalpy updates
    typename P1ENTH::KernelWrapper m_p1Enthalpy;

    /// Kernel wrapper for brine internal energy updates
    typename P1INTENERGY::KernelWrapper m_p1IntEnergy;


    // CO2 constitutive kernel wrapper

    /// Kernel wrapper for CO2 density updates
    typename P2DENS::KernelWrapper m_p2Density;

    /// Kernel wrapper for CO2 viscosity updates
    typename P2VISC::KernelWrapper m_p2Viscosity;

    /// Kernel wrapper for CO2 enthalpy updates
    typename P2ENTH::KernelWrapper m_p2Enthalpy;

    /// Kernel wrapper for CO2 internal energy updates
    typename P2INTENERGY::KernelWrapper m_p2IntEnergy;


    // Flash kernel wrapper

    /// Kernel wrapper for phase fraction and phase component fraction updates
    typename FLASH::KernelWrapper m_flash;
  };

  virtual integer getWaterPhaseIndex() const override final;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper();

  struct viewKeyStruct : MultiFluidBase::viewKeyStruct
  {
    static constexpr char const * flashModelParaFileString() { return "flashModelParaFile"; }
    static constexpr char const * phasePVTParaFilesString() { return "phasePVTParaFiles"; }
  };

protected:

  virtual void postProcessInput() override;

private:

  void createPVTModels();

  /// Names of the files defining the viscosity and density models
  path_array m_phasePVTParaFiles;

  /// Name of the file defining the flash model
  Path m_flashModelParaFile;

  /// Index of the liquid phase
  integer m_p1Index;

  /// Index of the gas phase
  integer m_p2Index;


  // Brine constitutive models

  /// Pointer to the brine density model
  std::unique_ptr< P1DENS > m_p1Density;

  /// Pointer to the brine viscosity model
  std::unique_ptr< P1VISC > m_p1Viscosity;

  /// Pointer to the brine enthalpy model
  std::unique_ptr< P1ENTH > m_p1Enthalpy;

  /// Pointer to the brine internal energy model
  std::unique_ptr< P1INTENERGY > m_p1IntEnergy;


  // CO2 constitutive models

  /// Pointer to the CO2 density model
  std::unique_ptr< P2DENS > m_p2Density;

  /// Pointer to the CO2 viscosity model
  std::unique_ptr< P2VISC > m_p2Viscosity;

  /// Pointer to the CO2 enthalpy model
  std::unique_ptr< P2ENTH > m_p2Enthalpy;

  /// Pointer to the CO2 internal energy model
  std::unique_ptr< P2INTENERGY > m_p2IntEnergy;


  // Flash model

  /// Pointer to the flash model
  std::unique_ptr< FLASH > m_flash;

};

// these aliases are useful in constitutive dispatch
using CO2BrinePhillipsFluid =
  CO2BrineFluid< PVTProps::PhillipsBrineDensity, PVTProps::PhillipsBrineViscosity, PVTProps::NoOpPVTFunction, PVTProps::NoOpPVTFunction,
                 PVTProps::SpanWagnerCO2Density, PVTProps::FenghourCO2Viscosity, PVTProps::NoOpPVTFunction, PVTProps::NoOpPVTFunction,
                 PVTProps::CO2Solubility >;
using CO2BrinePhillipsThermalFluid =
  CO2BrineFluid< PVTProps::PhillipsBrineDensity, PVTProps::PhillipsBrineViscosity, PVTProps::BrineEnthalpy, PVTProps::BrineInternalEnergy,
                 PVTProps::SpanWagnerCO2Density, PVTProps::FenghourCO2Viscosity, PVTProps::CO2Enthalpy, PVTProps::CO2InternalEnergy,
                 PVTProps::CO2Solubility >;

using CO2BrineEzrokhiFluid =
  CO2BrineFluid< PVTProps::EzrokhiBrineDensity, PVTProps::EzrokhiBrineViscosity, PVTProps::NoOpPVTFunction, PVTProps::NoOpPVTFunction,
                 PVTProps::SpanWagnerCO2Density, PVTProps::FenghourCO2Viscosity, PVTProps::NoOpPVTFunction, PVTProps::NoOpPVTFunction,
                 PVTProps::CO2Solubility >;
using CO2BrineEzrokhiThermalFluid =
  CO2BrineFluid< PVTProps::EzrokhiBrineDensity, PVTProps::EzrokhiBrineViscosity, PVTProps::BrineEnthalpy, PVTProps::BrineInternalEnergy,
                 PVTProps::SpanWagnerCO2Density, PVTProps::FenghourCO2Viscosity, PVTProps::CO2Enthalpy, PVTProps::CO2InternalEnergy,
                 PVTProps::CO2Solubility >;

template< typename P1DENS, typename P1VISC, typename P1ENTH, typename P1INTENERGY,
          typename P2DENS, typename P2VISC, typename P2ENTH, typename P2INTENERGY,
          typename FLASH >
GEOSX_HOST_DEVICE
inline void
CO2BrineFluid< P1DENS, P1VISC, P1ENTH, P1INTENERGY,
               P2DENS, P2VISC, P2ENTH, P2INTENERGY,
               FLASH >::KernelWrapper::
  compute( real64 pressure,
           real64 temperature,
           arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition,
           arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseFraction,
           arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseDensity,
           arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseMassDensity,
           arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseViscosity,
           arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseEnthalpy,
           arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseInternalEnergy,
           arraySlice2d< real64, multifluid::USD_PHASE_COMP-2 > const & phaseCompFraction,
           real64 & totalDensity ) const
{
  integer constexpr numComps = 2;
  integer constexpr numPhases = 2;
  integer const ip1 = m_p1Index;
  integer const ip2 = m_p2Index;

  // 1. Convert input mass fractions to mole fractions and keep derivatives

  stackArray1d< real64, numComps > compMoleFrac( numComps );
  if( m_useMass )
  {
    real64 totalMolality = 0.0;
    for( integer ic = 0; ic < numComps; ++ic )
    {
      real64 const mwInv = 1.0 / m_componentMolarWeight[ic];
      compMoleFrac[ic] = composition[ic] * mwInv; // this is molality (units of mole/mass)
      totalMolality += compMoleFrac[ic];
    }

    real64 const totalMolalityInv = 1.0 / totalMolality;
    for( integer ic = 0; ic < numComps; ++ic )
    {
      compMoleFrac[ic] *= totalMolalityInv;
    }
  }
  else
  {
    for( integer ic = 0; ic < numComps; ++ic )
    {
      compMoleFrac[ic] = composition[ic];
    }
  }

  // 2. Compute phase fractions and phase component fractions

  real64 const temperatureInCelsius = temperature - 273.15;
  m_flash.compute( pressure,
                   temperatureInCelsius,
                   compMoleFrac.toSliceConst(),
                   phaseFraction,
                   phaseCompFraction );

  // 3. Compute phase densities and phase viscosities

  m_p1Density.compute( pressure,
                       temperatureInCelsius,
                       phaseCompFraction[ip1].toSliceConst(),
                       phaseDensity[ip1],
                       m_useMass );
  m_p1Viscosity.compute( pressure,
                         temperatureInCelsius,
                         phaseCompFraction[ip1].toSliceConst(),
                         phaseViscosity[ip1],
                         m_useMass );

  m_p2Density.compute( pressure,
                       temperatureInCelsius,
                       phaseCompFraction[ip2].toSliceConst(),
                       phaseDensity[ip2],
                       m_useMass );
  m_p2Viscosity.compute( pressure,
                         temperatureInCelsius,
                         phaseCompFraction[ip2].toSliceConst(),
                         phaseViscosity[ip2],
                         m_useMass );

  // 4. Compute enthalpy and internal energy

  m_p1Enthalpy.compute( pressure,
                        temperatureInCelsius,
                        phaseCompFraction[ip1].toSliceConst(),
                        phaseEnthalpy[ip1],
                        m_useMass );
  m_p1IntEnergy.compute( pressure,
                         temperatureInCelsius,
                         phaseCompFraction[ip1].toSliceConst(),
                         phaseInternalEnergy[ip1],
                         m_useMass );

  m_p2Enthalpy.compute( pressure,
                        temperatureInCelsius,
                        phaseCompFraction[ip2].toSliceConst(),
                        phaseEnthalpy[ip2],
                        m_useMass );
  m_p2IntEnergy.compute( pressure,
                         temperatureInCelsius,
                         phaseCompFraction[ip2].toSliceConst(),
                         phaseInternalEnergy[ip2],
                         m_useMass );

  // 5. Depending on the m_useMass flag, convert to mass variables or simply compute mass density

  // TODO: for now the treatment of molar/mass density requires too many interpolations in the tables, it needs to be fixed
  //       we should modify the PVT functions so that they can return phaseMassDens, phaseDens, and phaseMW in one call
  // TODO: extract the following piece of code and write a function that can be used here and in MultiFluidPVTPackageWrapper

  if( m_useMass )
  {
    // 5.1. Convert phase fractions (requires two passes)
    real64 totalMass{};

    // 5.1.0. Compute the phase molecular weights (ultimately, get that from the PVT function)
    real64 phaseMW[2]{};
    real64 phaseMolarDens = 0.0;
    m_p1Density.compute( pressure,
                         temperatureInCelsius,
                         phaseCompFraction[ip1].toSliceConst(),
                         phaseMolarDens,
                         false );
    phaseMW[ip1] = phaseDensity[ip1] / phaseMolarDens;

    m_p2Density.compute( pressure,
                         temperatureInCelsius,
                         phaseCompFraction[ip2].toSliceConst(),
                         phaseMolarDens,
                         false );
    phaseMW[ip2] = phaseDensity[ip2] / phaseMolarDens;

    // 5.1.1. Compute mass of each phase and total mass (on a 1-mole basis)
    for( integer ip = 0; ip < numPhases; ++ip )
    {
      phaseFraction[ip] *= phaseMW[ip];
      totalMass += phaseFraction[ip];
    }
    // 5.1.2. Normalize to get mass fractions
    real64 const totalMassInv = 1.0 / totalMass;
    for( integer ip = 0; ip < numPhases; ++ip )
    {
      phaseFraction[ip] *= totalMassInv;
    }
    // 5.2. Convert phase compositions
    for( integer ip = 0; ip < numPhases; ++ip )
    {
      real64 const phaseMWInv = 1.0 / phaseMW[ip];

      for( integer ic = 0; ic < numComps; ++ic )
      {
        real64 const compMW = m_componentMolarWeight[ic];
        phaseCompFraction[ip][ic] = phaseCompFraction[ip][ic] * compMW * phaseMWInv;
      }
    }
    // 5.3 Copy the phase densities into the phase mass densities
    for( integer ip = 0; ip < numPhases; ++ip )
    {
      phaseMassDensity[ip] = phaseDensity[ip];
    }
  }
  else
  {
    // for now, we have to compute the phase mass density here
    m_p1Density.compute( pressure,
                         temperatureInCelsius,
                         phaseCompFraction[ip1].toSliceConst(),
                         phaseMassDensity[ip1],
                         true );
    m_p2Density.compute( pressure,
                         temperatureInCelsius,
                         phaseCompFraction[ip2].toSliceConst(),
                         phaseMassDensity[ip2],
                         true );
  }

  // TODO: extract the following piece of code and write a function that can be used here and in MultiFluidPVTPackageWrapper

  // 6. Compute total fluid mass/molar density and derivatives
  {
    totalDensity = 0.0;
    // 6.1. Sum mass/molar fraction/density ratio over all phases to get the inverse of density
    for( integer ip = 0; ip < numPhases; ++ip )
    {
      real64 const densInv = 1.0 / phaseDensity[ip];
      real64 const value = phaseFraction[ip] * densInv;
      totalDensity += value;
    }
    // 6.2. Invert the previous quantity to get actual density
    totalDensity = 1.0 / totalDensity;
  }
}

template< typename P1DENS, typename P1VISC, typename P1ENTH, typename P1INTENERGY,
          typename P2DENS, typename P2VISC, typename P2ENTH, typename P2INTENERGY,
          typename FLASH >
GEOSX_HOST_DEVICE
inline void
CO2BrineFluid< P1DENS, P1VISC, P1ENTH, P1INTENERGY,
               P2DENS, P2VISC, P2ENTH, P2INTENERGY,
               FLASH >::KernelWrapper::
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
  integer constexpr numComps = 2;
  integer constexpr numPhases = 2;
  integer const ip1 = m_p1Index;
  integer const ip2 = m_p2Index;

  // 1. Convert input mass fractions to mole fractions and keep derivatives

  stackArray1d< real64, numComps > compMoleFrac( numComps );
  stackArray2d< real64, numComps * numComps > dCompMoleFrac_dCompMassFrac( numComps, numComps );

  if( m_useMass )
  {
    real64 totalMolality = 0.0;
    for( integer ic = 0; ic < numComps; ++ic )
    {
      real64 const mwInv = 1.0 / m_componentMolarWeight[ic];
      compMoleFrac[ic] = composition[ic] * mwInv; // this is molality (units of mole/mass)
      dCompMoleFrac_dCompMassFrac[ic][ic] = mwInv;
      totalMolality += compMoleFrac[ic];
    }

    real64 const totalMolalityInv = 1.0 / totalMolality;
    for( integer ic = 0; ic < numComps; ++ic )
    {
      compMoleFrac[ic] *= totalMolalityInv;

      for( integer jc = 0; jc < numComps; ++jc )
      {
        dCompMoleFrac_dCompMassFrac[ic][jc] -= compMoleFrac[ic] / m_componentMolarWeight[jc];
        dCompMoleFrac_dCompMassFrac[ic][jc] *= totalMolalityInv;
      }
    }
  }
  else
  {
    for( integer ic = 0; ic < numComps; ++ic )
    {
      compMoleFrac[ic] = composition[ic];
    }
  }

  // 2. Compute phase fractions and phase component fractions

  real64 const temperatureInCelsius = temperature - 273.15;
  m_flash.compute( pressure,
                   temperatureInCelsius,
                   compMoleFrac.toSliceConst(),
                   phaseFraction.value, phaseFraction.dPres, phaseFraction.dTemp, phaseFraction.dComp,
                   phaseCompFraction.value, phaseCompFraction.dPres, phaseCompFraction.dTemp, phaseCompFraction.dComp );

  // 3. Compute phase densities and phase viscosities

  m_p1Density.compute( pressure,
                       temperatureInCelsius,
                       phaseCompFraction.value[ip1].toSliceConst(), phaseCompFraction.dPres[ip1].toSliceConst(),
                       phaseCompFraction.dTemp[ip1].toSliceConst(), phaseCompFraction.dComp[ip1].toSliceConst(),
                       phaseDensity.value[ip1], phaseDensity.dPres[ip1],
                       phaseDensity.dTemp[ip1], phaseDensity.dComp[ip1],
                       m_useMass );
  m_p1Viscosity.compute( pressure,
                         temperatureInCelsius,
                         phaseCompFraction.value[ip1].toSliceConst(), phaseCompFraction.dPres[ip1].toSliceConst(),
                         phaseCompFraction.dTemp[ip1].toSliceConst(), phaseCompFraction.dComp[ip1].toSliceConst(),
                         phaseViscosity.value[ip1], phaseViscosity.dPres[ip1],
                         phaseViscosity.dTemp[ip1], phaseViscosity.dComp[ip1],
                         m_useMass );
  m_p2Density.compute( pressure,
                       temperatureInCelsius,
                       phaseCompFraction.value[ip2].toSliceConst(), phaseCompFraction.dPres[ip2].toSliceConst(),
                       phaseCompFraction.dTemp[ip2].toSliceConst(), phaseCompFraction.dComp[ip2].toSliceConst(),
                       phaseDensity.value[ip2], phaseDensity.dPres[ip2],
                       phaseDensity.dTemp[ip2], phaseDensity.dComp[ip2],
                       m_useMass );
  m_p2Viscosity.compute( pressure,
                         temperatureInCelsius,
                         phaseCompFraction.value[ip2].toSliceConst(), phaseCompFraction.dPres[ip2].toSliceConst(),
                         phaseCompFraction.dTemp[ip2].toSliceConst(), phaseCompFraction.dComp[ip2].toSliceConst(),
                         phaseViscosity.value[ip2], phaseViscosity.dPres[ip2],
                         phaseViscosity.dTemp[ip2], phaseViscosity.dComp[ip2],
                         m_useMass );


  // 4. Compute enthalpy and internal energy

  m_p1Enthalpy.compute( pressure,
                        temperatureInCelsius,
                        phaseCompFraction.value[ip1].toSliceConst(), phaseCompFraction.dPres[ip1].toSliceConst(),
                        phaseCompFraction.dTemp[ip1].toSliceConst(), phaseCompFraction.dComp[ip1].toSliceConst(),
                        phaseEnthalpy.value[ip1], phaseEnthalpy.dPres[ip1],
                        phaseEnthalpy.dTemp[ip1], phaseEnthalpy.dComp[ip1],
                        m_useMass );
  m_p1IntEnergy.compute( pressure,
                         temperatureInCelsius,
                         phaseCompFraction.value[ip1].toSliceConst(), phaseCompFraction.dPres[ip1].toSliceConst(),
                         phaseCompFraction.dTemp[ip1].toSliceConst(), phaseCompFraction.dComp[ip1].toSliceConst(),
                         phaseInternalEnergy.value[ip1], phaseInternalEnergy.dPres[ip1],
                         phaseInternalEnergy.dTemp[ip1], phaseInternalEnergy.dComp[ip1],
                         m_useMass );

  m_p2Enthalpy.compute( pressure,
                        temperatureInCelsius,
                        phaseCompFraction.value[ip2].toSliceConst(), phaseCompFraction.dPres[ip2].toSliceConst(),
                        phaseCompFraction.dTemp[ip2].toSliceConst(), phaseCompFraction.dComp[ip2].toSliceConst(),
                        phaseEnthalpy.value[ip2], phaseEnthalpy.dPres[ip2],
                        phaseEnthalpy.dTemp[ip2], phaseEnthalpy.dComp[ip2],
                        m_useMass );
  m_p2IntEnergy.compute( pressure,
                         temperatureInCelsius,
                         phaseCompFraction.value[ip2].toSliceConst(), phaseCompFraction.dPres[ip2].toSliceConst(),
                         phaseCompFraction.dTemp[ip2].toSliceConst(), phaseCompFraction.dComp[ip2].toSliceConst(),
                         phaseInternalEnergy.value[ip2], phaseInternalEnergy.dPres[ip2],
                         phaseInternalEnergy.dTemp[ip2], phaseInternalEnergy.dComp[ip2],
                         m_useMass );

  // 5. Depending on the m_useMass flag, convert to mass variables or simply compute mass density

  // TODO: for now the treatment of molar/mass density requires too many interpolations in the tables, it needs to be fixed
  //       we should modify the PVT functions so that they can return phaseMassDens, phaseDens, and phaseMW in one call
  // TODO: extract the following piece of code and write a function that can be used here and in MultiFluidPVTPackageWrapper

  if( m_useMass )
  {
    // 5.1. Convert phase fractions (requires two passes)
    real64 totalMass{};
    real64 dTotalMass_dP{};
    real64 dTotalMass_dT{};
    real64 dTotalMass_dC[2]{};

    // 5.1.0. Compute the phase molecular weights (ultimately, get that from the PVT function)
    real64 phaseMW[2]{};
    real64 dPhaseMW_dPres[2]{};
    real64 dPhaseMW_dTemp[2]{};
    real64 dPhaseMW_dComp[2][2]{};
    real64 phaseMolarDens = 0.0;
    real64 dPhaseMolarDens_dPres = 0.0;
    real64 dPhaseMolarDens_dTemp = 0.0;
    stackArray1d< real64, numComps > dPhaseMolarDens_dComp( 2 );
    m_p2Density.compute( pressure,
                         temperatureInCelsius,
                         phaseCompFraction.value[ip2].toSliceConst(), phaseCompFraction.dPres[ip2].toSliceConst(),
                         phaseCompFraction.dTemp[ip2].toSliceConst(), phaseCompFraction.dComp[ip2].toSliceConst(),
                         phaseMolarDens, dPhaseMolarDens_dPres,
                         dPhaseMolarDens_dTemp, dPhaseMolarDens_dComp.toSlice(),
                         false );
    phaseMW[ip2] = phaseDensity.value[ip2] / phaseMolarDens;
    dPhaseMW_dPres[ip2] = phaseDensity.dPres[ip2] / phaseMolarDens - phaseMW[ip2] * dPhaseMolarDens_dPres / phaseMolarDens;
    dPhaseMW_dTemp[ip2] = phaseDensity.dTemp[ip2] / phaseMolarDens - phaseMW[ip2] * dPhaseMolarDens_dTemp / phaseMolarDens;
    for( integer ic = 0; ic < numComps; ++ic )
    {
      dPhaseMW_dComp[ip2][ic] = phaseDensity.dComp[ip2][ic] / phaseMolarDens - phaseMW[ip2] * dPhaseMolarDens_dComp[ic] / phaseMolarDens;
    }
    m_p1Density.compute( pressure,
                         temperatureInCelsius,
                         phaseCompFraction.value[ip1].toSliceConst(), phaseCompFraction.dPres[ip1].toSliceConst(),
                         phaseCompFraction.dTemp[ip1].toSliceConst(), phaseCompFraction.dComp[ip1].toSliceConst(),
                         phaseMolarDens, dPhaseMolarDens_dPres,
                         dPhaseMolarDens_dTemp, dPhaseMolarDens_dComp.toSlice(),
                         false );
    phaseMW[ip1] = phaseDensity.value[ip1] / phaseMolarDens;
    dPhaseMW_dPres[ip1] = phaseDensity.dPres[ip1] / phaseMolarDens - phaseMW[ip1] * dPhaseMolarDens_dPres / phaseMolarDens;
    dPhaseMW_dTemp[ip1] = phaseDensity.dTemp[ip1] / phaseMolarDens - phaseMW[ip1] * dPhaseMolarDens_dTemp / phaseMolarDens;
    for( integer ic = 0; ic < numComps; ++ic )
    {
      dPhaseMW_dComp[ip1][ic] = phaseDensity.dComp[ip1][ic] / phaseMolarDens - phaseMW[ip1] * dPhaseMolarDens_dComp[ic] / phaseMolarDens;
    }


    // 5.1.1. Compute mass of each phase and total mass (on a 1-mole basis)
    for( integer ip = 0; ip < numPhases; ++ip )
    {
      real64 const nu = phaseFraction.value[ip];

      phaseFraction.value[ip] *= phaseMW[ip];
      phaseFraction.dPres[ip] = phaseFraction.dPres[ip] * phaseMW[ip] + nu * dPhaseMW_dPres[ip];
      phaseFraction.dTemp[ip] = phaseFraction.dTemp[ip] * phaseMW[ip] + nu * dPhaseMW_dTemp[ip];

      totalMass += phaseFraction.value[ip];
      dTotalMass_dP += phaseFraction.dPres[ip];
      dTotalMass_dT += phaseFraction.dTemp[ip];

      for( integer jc = 0; jc < numComps; ++jc )
      {
        phaseFraction.dComp[ip][jc] = phaseFraction.dComp[ip][jc] * phaseMW[ip] + nu * dPhaseMW_dComp[ip][jc];
        dTotalMass_dC[jc] += phaseFraction.dComp[ip][jc];
      }
    }

    // 5.1.2. Normalize to get mass fractions
    real64 const totalMassInv = 1.0 / totalMass;
    for( integer ip = 0; ip < numPhases; ++ip )
    {
      phaseFraction.value[ip] *= totalMassInv;
      phaseFraction.dPres[ip] = ( phaseFraction.dPres[ip] - phaseFraction.value[ip] * dTotalMass_dP ) * totalMassInv;
      phaseFraction.dTemp[ip] = ( phaseFraction.dTemp[ip] - phaseFraction.value[ip] * dTotalMass_dT ) * totalMassInv;

      for( integer jc = 0; jc < numComps; ++jc )
      {
        phaseFraction.dComp[ip][jc] = ( phaseFraction.dComp[ip][jc] - phaseFraction.value[ip] * dTotalMass_dC[jc] ) * totalMassInv;
      }
    }

    // 5.2. Convert phase compositions
    for( integer ip = 0; ip < numPhases; ++ip )
    {
      real64 const phaseMWInv = 1.0 / phaseMW[ip];

      for( integer ic = 0; ic < numComps; ++ic )
      {

        real64 const compMW = m_componentMolarWeight[ic];

        phaseCompFraction.value[ip][ic] = phaseCompFraction.value[ip][ic] * compMW * phaseMWInv;
        phaseCompFraction.dPres[ip][ic] =
          ( phaseCompFraction.dPres[ip][ic] * compMW - phaseCompFraction.value[ip][ic] * dPhaseMW_dPres[ip] ) * phaseMWInv;
        phaseCompFraction.dTemp[ip][ic] =
          ( phaseCompFraction.dTemp[ip][ic] * compMW - phaseCompFraction.value[ip][ic] * dPhaseMW_dTemp[ip] ) * phaseMWInv;

        for( integer jc = 0; jc < numComps; ++jc )
        {
          phaseCompFraction.dComp[ip][ic][jc] =
            ( phaseCompFraction.dComp[ip][ic][jc] * compMW - phaseCompFraction.value[ip][ic] * dPhaseMW_dComp[ip][jc] ) * phaseMWInv;
        }
      }
    }

    // 5.3. Update derivatives w.r.t. mole fractions to derivatives w.r.t mass fractions
    real64 work[numComps]{};
    for( integer ip = 0; ip < numPhases; ++ip )
    {
      applyChainRuleInPlace( numComps, dCompMoleFrac_dCompMassFrac, phaseFraction.dComp[ip], work );
      applyChainRuleInPlace( numComps, dCompMoleFrac_dCompMassFrac, phaseDensity.dComp[ip], work );

      for( integer ic = 0; ic < numComps; ++ic )
      {
        applyChainRuleInPlace( numComps, dCompMoleFrac_dCompMassFrac, phaseCompFraction.dComp[ip][ic], work );
      }
    }

    // 5.4 Copy the phase densities into the phase mass densities
    for( integer ip = 0; ip < numPhases; ++ip )
    {
      phaseMassDensity.value[ip] = phaseDensity.value[ip];
      phaseMassDensity.dPres[ip] = phaseDensity.dPres[ip];
      phaseMassDensity.dTemp[ip] = phaseDensity.dTemp[ip];
      for( integer ic = 0; ic < numComps; ++ic )
      {
        phaseMassDensity.dComp[ip][ic] = phaseDensity.dComp[ip][ic];
      }
    }
  }
  else
  {
    // for now, we have to compute the phase mass density here
    m_p1Density.compute( pressure,
                         temperatureInCelsius,
                         phaseCompFraction.value[ip1].toSliceConst(), phaseCompFraction.dPres[ip1].toSliceConst(),
                         phaseCompFraction.dTemp[ip1].toSliceConst(), phaseCompFraction.dComp[ip1].toSliceConst(),
                         phaseMassDensity.value[ip1], phaseMassDensity.dPres[ip1],
                         phaseMassDensity.dTemp[ip1], phaseMassDensity.dComp[ip1],
                         true );
    m_p2Density.compute( pressure,
                         temperatureInCelsius,
                         phaseCompFraction.value[ip2].toSliceConst(), phaseCompFraction.dPres[ip2].toSliceConst(),
                         phaseCompFraction.dTemp[ip2].toSliceConst(), phaseCompFraction.dComp[ip2].toSliceConst(),
                         phaseMassDensity.value[ip2], phaseMassDensity.dPres[ip2],
                         phaseMassDensity.dTemp[ip2], phaseMassDensity.dComp[ip2],
                         true );
  }

  // TODO: extract the following piece of code and write a function that can be used here and in MultiFluidPVTPackageWrapper

  // 6. Compute total fluid mass/molar density and derivatives
  {
    totalDensity.value = 0.0;
    totalDensity.dPres = 0.0;
    totalDensity.dTemp = 0.0;
    for( integer jc = 0; jc < numComps; ++jc )
    {
      totalDensity.dComp[jc] = 0.0;
    }

    // 6.1. Sum mass/molar fraction/density ratio over all phases to get the inverse of density
    for( integer ip = 0; ip < numPhases; ++ip )
    {
      real64 const densInv = 1.0 / phaseDensity.value[ip];
      real64 const value = phaseFraction.value[ip] * densInv;

      totalDensity.value += value;
      totalDensity.dPres += ( phaseFraction.dPres[ip] - value * phaseDensity.dPres[ip] ) * densInv;
      totalDensity.dTemp += ( phaseFraction.dTemp[ip] - value * phaseDensity.dTemp[ip] ) * densInv;

      for( integer jc = 0; jc < numComps; ++jc )
      {
        totalDensity.dComp[jc] += ( phaseFraction.dComp[ip][jc] - value * phaseDensity.dComp[ip][jc] ) * densInv;
      }
    }

    // 6.2. Invert the previous quantity to get actual density
    totalDensity.value = 1.0 / totalDensity.value;
    real64 const minusDens2 = -totalDensity.value * totalDensity.value;
    totalDensity.dPres *= minusDens2;
    totalDensity.dTemp *= minusDens2;

    for( integer jc = 0; jc < numComps; ++jc )
    {
      totalDensity.dComp[jc] *= minusDens2;
    }
  }
}

template< typename P1DENS, typename P1VISC, typename P1ENTH, typename P1INTENERGY,
          typename P2DENS, typename P2VISC, typename P2ENTH, typename P2INTENERGY,
          typename FLASH >
GEOSX_HOST_DEVICE inline void
CO2BrineFluid< P1DENS, P1VISC, P1ENTH, P1INTENERGY,
               P2DENS, P2VISC, P2ENTH, P2INTENERGY,
               FLASH >::KernelWrapper::
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

} // namespace constitutive

} // namespace geosx

#endif //GEOSX_CONSTITUTIVE_FLUID_MULTIPHASEMULTICOMPONENTFLUID_HPP_
