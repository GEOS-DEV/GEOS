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
 * @file MultiPhaseMultiComponentFluid.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_MULTIPHASEMULTICOMPONENTFLUID_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_MULTIPHASEMULTICOMPONENTFLUID_HPP_

#include "constitutive/fluid/MultiFluidBase.hpp"
#include "constitutive/fluid/MultiFluidUtils.hpp"
#include "constitutive/fluid/PVTFunctions/BrineCO2Density.hpp"
#include "constitutive/fluid/PVTFunctions/BrineViscosity.hpp"
#include "constitutive/fluid/PVTFunctions/CO2Solubility.hpp"
#include "constitutive/fluid/PVTFunctions/FenghourCO2Viscosity.hpp"
#include "constitutive/fluid/PVTFunctions/SpanWagnerCO2Density.hpp"

#include <memory>

namespace geosx
{

namespace constitutive
{

template< typename P1DENS, typename P1VISC, typename P2DENS, typename P2VISC, typename FLASH >
class MultiPhaseMultiComponentFluid : public MultiFluidBase
{
public:

  using exec_policy = parallelDevicePolicy<>;

  MultiPhaseMultiComponentFluid( string const & name,
                                 Group * const parent );

  virtual ~MultiPhaseMultiComponentFluid() override = default;

  virtual std::unique_ptr< ConstitutiveBase >
  deliverClone( string const & name,
                Group * const parent ) const override;


  static string catalogName();

  virtual string getCatalogName() const override { return catalogName(); }

  /**
   * @brief Kernel wrapper class for MultiPhaseMultiComponentFluid.
   */
  class KernelWrapper final : public MultiFluidBase::KernelWrapper
  {
public:

    GEOSX_HOST_DEVICE
    virtual void compute( real64 pressure,
                          real64 temperature,
                          arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition,
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseFraction,
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseDensity,
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseMassDensity,
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseViscosity,
                          arraySlice2d< real64, multifluid::USD_PHASE_COMP-2 > const & phaseCompFraction,
                          real64 & totalDensity ) const override;

    GEOSX_HOST_DEVICE
    virtual void compute( real64 pressure,
                          real64 temperature,
                          arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition,
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseFraction,
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & dPhaseFraction_dPressure,
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & dPhaseFraction_dTemperature,
                          arraySlice2d< real64, multifluid::USD_PHASE_DC - 2 > const & dPhaseFraction_dGlobalCompFraction,
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseDensity,
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & dPhaseDensity_dPressure,
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & dPhaseDensity_dTemperature,
                          arraySlice2d< real64, multifluid::USD_PHASE_DC - 2 > const & dPhaseDensity_dGlobalCompFraction,
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseMassDensity,
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & dPhaseMassDensity_dPressure,
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & dPhaseMassDensity_dTemperature,
                          arraySlice2d< real64, multifluid::USD_PHASE_DC - 2 > const & dPhaseMassDensity_dGlobalCompFraction,
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseViscosity,
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & dPhaseViscosity_dPressure,
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & dPhaseViscosity_dTemperature,
                          arraySlice2d< real64, multifluid::USD_PHASE_DC - 2 > const & dPhaseViscosity_dGlobalCompFraction,
                          arraySlice2d< real64, multifluid::USD_PHASE_COMP-2 > const & phaseCompFraction,
                          arraySlice2d< real64, multifluid::USD_PHASE_COMP-2 > const & dPhaseCompFraction_dPressure,
                          arraySlice2d< real64, multifluid::USD_PHASE_COMP-2 > const & dPhaseCompFraction_dTemperature,
                          arraySlice3d< real64, multifluid::USD_PHASE_COMP_DC-2 > const & dPhaseCompFraction_dGlobalCompFraction,
                          real64 & totalDensity,
                          real64 & dTotalDensity_dPressure,
                          real64 & dTotalDensity_dTemperature,
                          arraySlice1d< real64, multifluid::USD_FLUID_DC - 2 > const & dTotalDensity_dGlobalCompFraction ) const override;

    GEOSX_HOST_DEVICE
    virtual void update( localIndex k,
                         localIndex q,
                         real64 pressure,
                         real64 temperature,
                         arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition ) const override;

private:

    friend class MultiPhaseMultiComponentFluid;

    KernelWrapper( localIndex p1Index,
                   localIndex p2Index,
                   P1DENS const & p1DensityWrapper,
                   P1VISC const & p1ViscosityWrapper,
                   P2DENS const & p2DensityWrapper,
                   P2VISC const & p2ViscosityWrapper,
                   FLASH const & flashWrapper,
                   arrayView1d< real64 const > const & componentMolarWeight,
                   bool useMass,
                   MultiFluidBase::KernelWrapper::PhasePropViews const & phaseFraction,
                   MultiFluidBase::KernelWrapper::PhasePropViews const & phaseDensity,
                   MultiFluidBase::KernelWrapper::PhasePropViews const & phaseMassDensity,
                   MultiFluidBase::KernelWrapper::PhasePropViews const & phaseViscosity,
                   MultiFluidBase::KernelWrapper::PhaseCompViews const & phaseCompFraction,
                   MultiFluidBase::KernelWrapper::FluidPropViews const & totalDensity );

    /// Index of the liquid phase
    localIndex m_p1Index;

    /// Index of the gas phase
    localIndex m_p2Index;

    /// Kernel wrapper for brine density updates
    typename P1DENS::KernelWrapper m_p1Density;

    /// Kernel wrapper for brine viscosity updates
    typename P1VISC::KernelWrapper m_p1Viscosity;

    /// Kernel wrapper for CO2 density updates
    typename P2DENS::KernelWrapper m_p2Density;

    /// Kernel wrapper for CO2 viscosity updates
    typename P2VISC::KernelWrapper m_p2Viscosity;

    /// Kernel wrapper for phase fraction and phase component fraction updates
    typename FLASH::KernelWrapper m_flash;
  };

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const;

  struct viewKeyStruct : MultiFluidBase::viewKeyStruct
  {
    static constexpr char const * flashModelParaFileString() { return "flashModelParaFile"; }
    static constexpr char const * phasePVTParaFilesString() { return "phasePVTParaFiles"; }
  } viewKeysMultiPhaseMultiComponentFluid;

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

  /// Pointer to the brine density model
  std::unique_ptr< P1DENS > m_p1Density;

  /// Pointer to the brine viscosity model
  std::unique_ptr< P1VISC > m_p1Viscosity;

  /// Pointer to the CO2 density model
  std::unique_ptr< P2DENS > m_p2Density;

  /// Pointer to the CO2 viscosity model
  std::unique_ptr< P2VISC > m_p2Viscosity;

  /// Pointer to the flash model
  std::unique_ptr< FLASH > m_flash;
};

// this alias will be useful in constitutive dispatch
using CO2BrineFluid = MultiPhaseMultiComponentFluid< PVTProps::BrineCO2Density,
                                                     PVTProps::BrineViscosity,
                                                     PVTProps::SpanWagnerCO2Density,
                                                     PVTProps::FenghourCO2Viscosity,
                                                     PVTProps::CO2Solubility >;

template< typename P1DENSWRAPPER, typename P1VISCWRAPPER, typename P2DENSWRAPPER, typename P2VISCWRAPPER, typename FLASHWRAPPER >
GEOSX_HOST_DEVICE
inline void
MultiPhaseMultiComponentFluid< P1DENSWRAPPER, P1VISCWRAPPER, P2DENSWRAPPER, P2VISCWRAPPER, FLASHWRAPPER >::KernelWrapper::
  compute( real64 pressure,
           real64 temperature,
           arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition,
           arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseFraction,
           arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseDensity,
           arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseMassDensity,
           arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseViscosity,
           arraySlice2d< real64, multifluid::USD_PHASE_COMP-2 > const & phaseCompFraction,
           real64 & totalDensity ) const
{
  constexpr localIndex numComps = 2;
  constexpr localIndex numPhases = 2;
  localIndex const ip1 = m_p1Index;
  localIndex const ip2 = m_p2Index;

  // 1. Convert input mass fractions to mole fractions and keep derivatives

  stackArray1d< real64, numComps > compMoleFrac( numComps );
  if( m_useMass )
  {
    real64 totalMolality = 0.0;
    for( localIndex ic = 0; ic < numComps; ++ic )
    {
      real64 const mwInv = 1.0 / m_componentMolarWeight[ic];
      compMoleFrac[ic] = composition[ic] * mwInv; // this is molality (units of mole/mass)
      totalMolality += compMoleFrac[ic];
    }

    real64 const totalMolalityInv = 1.0 / totalMolality;
    for( localIndex ic = 0; ic < numComps; ++ic )
    {
      compMoleFrac[ic] *= totalMolalityInv;
    }
  }
  else
  {
    for( localIndex ic = 0; ic < numComps; ++ic )
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

  // 4. Depending on the m_useMass flag, convert to mass variables or simply compute mass density

  // TODO: for now the treatment of molar/mass density requires too many interpolations in the tables, it needs to be fixed
  //       we should modify the PVT functions so that they can return phaseMassDens, phaseDens, and phaseMW in one call
  // TODO: extract the following piece of code and write a function that can be used here and in MultiFluidPVTPackageWrapper

  if( m_useMass )
  {
    // 4.1. Convert phase fractions (requires two passes)
    real64 totalMass{};

    // 4.1.0. Compute the phase molecular weights (ultimately, get that from the PVT function)
    real64 phaseMW[2];
    real64 phaseMolarDens;
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

    // 4.1.1. Compute mass of each phase and total mass (on a 1-mole basis)
    for( localIndex ip = 0; ip < numPhases; ++ip )
    {
      phaseFraction[ip] *= phaseMW[ip];
      totalMass += phaseFraction[ip];
    }
    // 4.1.2. Normalize to get mass fractions
    real64 const totalMassInv = 1.0 / totalMass;
    for( localIndex ip = 0; ip < numPhases; ++ip )
    {
      phaseFraction[ip] *= totalMassInv;
    }
    // 4.2. Convert phase compositions
    for( localIndex ip = 0; ip < numPhases; ++ip )
    {
      real64 const phaseMWInv = 1.0 / phaseMW[ip];

      for( localIndex ic = 0; ic < numComps; ++ic )
      {
        real64 const compMW = m_componentMolarWeight[ic];
        phaseCompFraction[ip][ic] = phaseCompFraction[ip][ic] * compMW * phaseMWInv;
      }
    }
    // 4.3 Copy the phase densities into the phase mass densities
    for( localIndex ip = 0; ip < numPhases; ++ip )
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

  // 5. Compute total fluid mass/molar density and derivatives
  {
    totalDensity = 0.0;
    // 5.1. Sum mass/molar fraction/density ratio over all phases to get the inverse of density
    for( localIndex ip = 0; ip < numPhases; ++ip )
    {
      real64 const densInv = 1.0 / phaseDensity[ip];
      real64 const value = phaseFraction[ip] * densInv;
      totalDensity += value;
    }
    // 5.2. Invert the previous quantity to get actual density
    totalDensity = 1.0 / totalDensity;
  }
}

template< typename P1DENS, typename P1VISC, typename P2DENS, typename P2VISC, typename FLASH >
GEOSX_HOST_DEVICE
inline void
MultiPhaseMultiComponentFluid< P1DENS, P1VISC, P2DENS, P2VISC, FLASH >::KernelWrapper::
  compute( real64 pressure,
           real64 temperature,
           arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition,
           arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseFraction,
           arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & dPhaseFraction_dPressure,
           arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & dPhaseFraction_dTemperature,
           arraySlice2d< real64, multifluid::USD_PHASE_DC - 2 > const & dPhaseFraction_dGlobalCompFraction,
           arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseDensity,
           arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & dPhaseDensity_dPressure,
           arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & dPhaseDensity_dTemperature,
           arraySlice2d< real64, multifluid::USD_PHASE_DC - 2 > const & dPhaseDensity_dGlobalCompFraction,
           arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseMassDensity,
           arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & dPhaseMassDensity_dPressure,
           arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & dPhaseMassDensity_dTemperature,
           arraySlice2d< real64, multifluid::USD_PHASE_DC - 2 > const & dPhaseMassDensity_dGlobalCompFraction,
           arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseViscosity,
           arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & dPhaseViscosity_dPressure,
           arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & dPhaseViscosity_dTemperature,
           arraySlice2d< real64, multifluid::USD_PHASE_DC - 2 > const & dPhaseViscosity_dGlobalCompFraction,
           arraySlice2d< real64, multifluid::USD_PHASE_COMP-2 > const & phaseCompFraction,
           arraySlice2d< real64, multifluid::USD_PHASE_COMP-2 > const & dPhaseCompFraction_dPressure,
           arraySlice2d< real64, multifluid::USD_PHASE_COMP-2 > const & dPhaseCompFraction_dTemperature,
           arraySlice3d< real64, multifluid::USD_PHASE_COMP_DC-2 > const & dPhaseCompFraction_dGlobalCompFraction,
           real64 & totalDensity,
           real64 & dTotalDensity_dPressure,
           real64 & dTotalDensity_dTemperature,
           arraySlice1d< real64, multifluid::USD_FLUID_DC - 2 > const & dTotalDensity_dGlobalCompFraction ) const
{
  // 0. make shortcut structs to avoid long names (TODO maybe remove)
  CompositionalVarContainer< 1, multifluid::USD_PHASE - 2, multifluid::USD_PHASE_DC - 2 > phaseFrac {
    phaseFraction,
    dPhaseFraction_dPressure,
    dPhaseFraction_dTemperature,
    dPhaseFraction_dGlobalCompFraction
  };

  CompositionalVarContainer< 1, multifluid::USD_PHASE - 2, multifluid::USD_PHASE_DC - 2 > phaseDens {
    phaseDensity,
    dPhaseDensity_dPressure,
    dPhaseDensity_dTemperature,
    dPhaseDensity_dGlobalCompFraction
  };

  CompositionalVarContainer< 1, multifluid::USD_PHASE - 2, multifluid::USD_PHASE_DC - 2 > phaseMassDens {
    phaseMassDensity,
    dPhaseMassDensity_dPressure,
    dPhaseMassDensity_dTemperature,
    dPhaseMassDensity_dGlobalCompFraction
  };

  CompositionalVarContainer< 1, multifluid::USD_PHASE - 2, multifluid::USD_PHASE_DC - 2 > phaseVisc {
    phaseViscosity,
    dPhaseViscosity_dPressure,
    dPhaseViscosity_dTemperature,
    dPhaseViscosity_dGlobalCompFraction
  };

  CompositionalVarContainer< 2, multifluid::USD_PHASE_COMP-2, multifluid::USD_PHASE_COMP_DC-2 > phaseCompFrac {
    phaseCompFraction,
    dPhaseCompFraction_dPressure,
    dPhaseCompFraction_dTemperature,
    dPhaseCompFraction_dGlobalCompFraction
  };

  CompositionalVarContainer< 0, multifluid::USD_FLUID - 2, multifluid::USD_FLUID_DC - 2 > totalDens {
    totalDensity,
    dTotalDensity_dPressure,
    dTotalDensity_dTemperature,
    dTotalDensity_dGlobalCompFraction
  };

#if defined(__CUDACC__)
  // For some reason nvcc thinks these aren't used.
  GEOSX_UNUSED_VAR( phaseFrac, phaseDens, phaseMassDens, phaseVisc, phaseCompFrac, totalDens );
#endif

  constexpr localIndex numComps = 2;
  constexpr localIndex numPhases = 2;
  localIndex const ip1 = m_p1Index;
  localIndex const ip2 = m_p2Index;

  // 1. Convert input mass fractions to mole fractions and keep derivatives

  stackArray1d< real64, numComps > compMoleFrac( numComps );
  stackArray2d< real64, numComps * numComps > dCompMoleFrac_dCompMassFrac( numComps, numComps );

  if( m_useMass )
  {
    real64 totalMolality = 0.0;
    for( localIndex ic = 0; ic < numComps; ++ic )
    {
      real64 const mwInv = 1.0 / m_componentMolarWeight[ic];
      compMoleFrac[ic] = composition[ic] * mwInv; // this is molality (units of mole/mass)
      dCompMoleFrac_dCompMassFrac[ic][ic] = mwInv;
      totalMolality += compMoleFrac[ic];
    }

    real64 const totalMolalityInv = 1.0 / totalMolality;
    for( localIndex ic = 0; ic < numComps; ++ic )
    {
      compMoleFrac[ic] *= totalMolalityInv;

      for( localIndex jc = 0; jc < numComps; ++jc )
      {
        dCompMoleFrac_dCompMassFrac[ic][jc] -= compMoleFrac[ic] / m_componentMolarWeight[jc];
        dCompMoleFrac_dCompMassFrac[ic][jc] *= totalMolalityInv;
      }
    }
  }
  else
  {
    for( localIndex ic = 0; ic < numComps; ++ic )
    {
      compMoleFrac[ic] = composition[ic];
    }
  }

  // 2. Compute phase fractions and phase component fractions

  real64 const temperatureInCelsius = temperature - 273.15;
  m_flash.compute( pressure,
                   temperatureInCelsius,
                   compMoleFrac.toSliceConst(),
                   phaseFrac.value, phaseFrac.dPres, phaseFrac.dTemp, phaseFrac.dComp,
                   phaseCompFrac.value, phaseCompFrac.dPres, phaseCompFrac.dTemp, phaseCompFrac.dComp );

  // 3. Compute phase densities and phase viscosities

  m_p1Density.compute( pressure,
                       temperatureInCelsius,
                       phaseCompFrac.value[ip1].toSliceConst(), phaseCompFrac.dPres[ip1].toSliceConst(),
                       phaseCompFrac.dTemp[ip1].toSliceConst(), phaseCompFrac.dComp[ip1].toSliceConst(),
                       phaseDens.value[ip1], phaseDens.dPres[ip1],
                       phaseDens.dTemp[ip1], phaseDens.dComp[ip1],
                       m_useMass );
  m_p1Viscosity.compute( pressure,
                         temperatureInCelsius,
                         phaseCompFrac.value[ip1].toSliceConst(), phaseCompFrac.dPres[ip1].toSliceConst(),
                         phaseCompFrac.dTemp[ip1].toSliceConst(), phaseCompFrac.dComp[ip1].toSliceConst(),
                         phaseVisc.value[ip1], phaseVisc.dPres[ip1],
                         phaseVisc.dTemp[ip1], phaseVisc.dComp[ip1],
                         m_useMass );
  m_p2Density.compute( pressure,
                       temperatureInCelsius,
                       phaseCompFrac.value[ip2].toSliceConst(), phaseCompFrac.dPres[ip2].toSliceConst(),
                       phaseCompFrac.dTemp[ip2].toSliceConst(), phaseCompFrac.dComp[ip2].toSliceConst(),
                       phaseDens.value[ip2], phaseDens.dPres[ip2],
                       phaseDens.dTemp[ip2], phaseDens.dComp[ip2],
                       m_useMass );
  m_p2Viscosity.compute( pressure,
                         temperatureInCelsius,
                         phaseCompFrac.value[ip2].toSliceConst(), phaseCompFrac.dPres[ip2].toSliceConst(),
                         phaseCompFrac.dTemp[ip2].toSliceConst(), phaseCompFrac.dComp[ip2].toSliceConst(),
                         phaseVisc.value[ip2], phaseVisc.dPres[ip2],
                         phaseVisc.dTemp[ip2], phaseVisc.dComp[ip2],
                         m_useMass );

  // 4. Depending on the m_useMass flag, convert to mass variables or simply compute mass density

  // TODO: for now the treatment of molar/mass density requires too many interpolations in the tables, it needs to be fixed
  //       we should modify the PVT functions so that they can return phaseMassDens, phaseDens, and phaseMW in one call
  // TODO: extract the following piece of code and write a function that can be used here and in MultiFluidPVTPackageWrapper

  if( m_useMass )
  {
    // 4.1. Convert phase fractions (requires two passes)
    real64 totalMass{};
    real64 dTotalMass_dP{};
    real64 dTotalMass_dT{};
    real64 dTotalMass_dC[2]{};

    // 4.1.0. Compute the phase molecular weights (ultimately, get that from the PVT function)
    real64 phaseMW[2];
    real64 dPhaseMW_dPres[2];
    real64 dPhaseMW_dTemp[2];
    real64 dPhaseMW_dComp[2][2];
    real64 phaseMolarDens;
    real64 dPhaseMolarDens_dPres;
    real64 dPhaseMolarDens_dTemp;
    stackArray1d< real64, numComps > dPhaseMolarDens_dComp( 2 );
    m_p2Density.compute( pressure,
                         temperatureInCelsius,
                         phaseCompFrac.value[ip2].toSliceConst(), phaseCompFrac.dPres[ip2].toSliceConst(),
                         phaseCompFrac.dTemp[ip2].toSliceConst(), phaseCompFrac.dComp[ip2].toSliceConst(),
                         phaseMolarDens, dPhaseMolarDens_dPres,
                         dPhaseMolarDens_dTemp, dPhaseMolarDens_dComp.toSlice(),
                         false );
    phaseMW[ip2] = phaseDens.value[ip2] / phaseMolarDens;
    dPhaseMW_dPres[ip2] = phaseDens.dPres[ip2] / phaseMolarDens - phaseMW[ip2] * dPhaseMolarDens_dPres / phaseMolarDens;
    dPhaseMW_dTemp[ip2] = phaseDens.dTemp[ip2] / phaseMolarDens - phaseMW[ip2] * dPhaseMolarDens_dTemp / phaseMolarDens;
    for( localIndex ic = 0; ic < numComps; ++ic )
    {
      dPhaseMW_dComp[ip2][ic] = phaseDens.dComp[ip2][ic] / phaseMolarDens - phaseMW[ip2] * dPhaseMolarDens_dComp[ic] / phaseMolarDens;
    }
    m_p1Density.compute( pressure,
                         temperatureInCelsius,
                         phaseCompFrac.value[ip1].toSliceConst(), phaseCompFrac.dPres[ip1].toSliceConst(),
                         phaseCompFrac.dTemp[ip1].toSliceConst(), phaseCompFrac.dComp[ip1].toSliceConst(),
                         phaseMolarDens, dPhaseMolarDens_dPres,
                         dPhaseMolarDens_dTemp, dPhaseMolarDens_dComp.toSlice(),
                         false );
    phaseMW[ip1] = phaseDens.value[ip1] / phaseMolarDens;
    dPhaseMW_dPres[ip1] = phaseDens.dPres[ip1] / phaseMolarDens - phaseMW[ip1] * dPhaseMolarDens_dPres / phaseMolarDens;
    dPhaseMW_dTemp[ip1] = phaseDens.dTemp[ip1] / phaseMolarDens - phaseMW[ip1] * dPhaseMolarDens_dTemp / phaseMolarDens;
    for( localIndex ic = 0; ic < numComps; ++ic )
    {
      dPhaseMW_dComp[ip1][ic] = phaseDens.dComp[ip1][ic] / phaseMolarDens - phaseMW[ip1] * dPhaseMolarDens_dComp[ic] / phaseMolarDens;
    }


    // 4.1.1. Compute mass of each phase and total mass (on a 1-mole basis)
    for( localIndex ip = 0; ip < numPhases; ++ip )
    {
      real64 const nu = phaseFrac.value[ip];

      phaseFrac.value[ip] *= phaseMW[ip];
      phaseFrac.dPres[ip] = phaseFrac.dPres[ip] * phaseMW[ip] + nu * dPhaseMW_dPres[ip];
      phaseFrac.dTemp[ip] = phaseFrac.dTemp[ip] * phaseMW[ip] + nu * dPhaseMW_dTemp[ip];

      totalMass += phaseFrac.value[ip];
      dTotalMass_dP += phaseFrac.dPres[ip];
      dTotalMass_dT += phaseFrac.dTemp[ip];

      for( localIndex jc = 0; jc < numComps; ++jc )
      {
        phaseFrac.dComp[ip][jc] = phaseFrac.dComp[ip][jc] * phaseMW[ip] + nu * dPhaseMW_dComp[ip][jc];
        dTotalMass_dC[jc] += phaseFrac.dComp[ip][jc];
      }
    }

    // 4.1.2. Normalize to get mass fractions
    real64 const totalMassInv = 1.0 / totalMass;
    for( localIndex ip = 0; ip < numPhases; ++ip )
    {
      phaseFrac.value[ip] *= totalMassInv;
      phaseFrac.dPres[ip] = ( phaseFrac.dPres[ip] - phaseFrac.value[ip] * dTotalMass_dP ) * totalMassInv;
      phaseFrac.dTemp[ip] = ( phaseFrac.dTemp[ip] - phaseFrac.value[ip] * dTotalMass_dT ) * totalMassInv;

      for( localIndex jc = 0; jc < numComps; ++jc )
      {
        phaseFrac.dComp[ip][jc] = ( phaseFrac.dComp[ip][jc] - phaseFrac.value[ip] * dTotalMass_dC[jc] ) * totalMassInv;
      }
    }

    // 4.2. Convert phase compositions
    for( localIndex ip = 0; ip < numPhases; ++ip )
    {
      real64 const phaseMWInv = 1.0 / phaseMW[ip];

      for( localIndex ic = 0; ic < numComps; ++ic )
      {

        real64 const compMW = m_componentMolarWeight[ic];

        phaseCompFrac.value[ip][ic] = phaseCompFrac.value[ip][ic] * compMW * phaseMWInv;
        phaseCompFrac.dPres[ip][ic] =
          ( phaseCompFrac.dPres[ip][ic] * compMW - phaseCompFrac.value[ip][ic] * dPhaseMW_dPres[ip] ) * phaseMWInv;
        phaseCompFrac.dTemp[ip][ic] =
          ( phaseCompFrac.dTemp[ip][ic] * compMW - phaseCompFrac.value[ip][ic] * dPhaseMW_dTemp[ip] ) * phaseMWInv;

        for( localIndex jc = 0; jc < numComps; ++jc )
        {
          phaseCompFrac.dComp[ip][ic][jc] =
            ( phaseCompFrac.dComp[ip][ic][jc] * compMW - phaseCompFrac.value[ip][ic] * dPhaseMW_dComp[ip][jc] ) * phaseMWInv;
        }
      }
    }

    // 4.3. Update derivatives w.r.t. mole fractions to derivatives w.r.t mass fractions
    real64 work[numComps];
    for( localIndex ip = 0; ip < numPhases; ++ip )
    {
      applyChainRuleInPlace( numComps, dCompMoleFrac_dCompMassFrac, phaseFrac.dComp[ip], work );
      applyChainRuleInPlace( numComps, dCompMoleFrac_dCompMassFrac, phaseDens.dComp[ip], work );

      for( localIndex ic = 0; ic < numComps; ++ic )
      {
        applyChainRuleInPlace( numComps, dCompMoleFrac_dCompMassFrac, phaseCompFrac.dComp[ip][ic], work );
      }
    }

    // 4.4 Copy the phase densities into the phase mass densities
    for( localIndex ip = 0; ip < numPhases; ++ip )
    {
      phaseMassDens.value[ip] = phaseDens.value[ip];
      phaseMassDens.dPres[ip] = phaseDens.dPres[ip];
      phaseMassDens.dTemp[ip] = phaseDens.dTemp[ip];
      for( localIndex ic = 0; ic < numComps; ++ic )
      {
        phaseMassDens.dComp[ip][ic] = phaseDens.dComp[ip][ic];
      }
    }
  }
  else
  {
    // for now, we have to compute the phase mass density here
    m_p1Density.compute( pressure,
                         temperatureInCelsius,
                         phaseCompFrac.value[ip1].toSliceConst(), phaseCompFrac.dPres[ip1].toSliceConst(),
                         phaseCompFrac.dTemp[ip1].toSliceConst(), phaseCompFrac.dComp[ip1].toSliceConst(),
                         phaseMassDens.value[ip1], phaseMassDens.dPres[ip1],
                         phaseMassDens.dTemp[ip1], phaseMassDens.dComp[ip1],
                         true );
    m_p2Density.compute( pressure,
                         temperatureInCelsius,
                         phaseCompFrac.value[ip2].toSliceConst(), phaseCompFrac.dPres[ip2].toSliceConst(),
                         phaseCompFrac.dTemp[ip2].toSliceConst(), phaseCompFrac.dComp[ip2].toSliceConst(),
                         phaseMassDens.value[ip2], phaseMassDens.dPres[ip2],
                         phaseMassDens.dTemp[ip2], phaseMassDens.dComp[ip2],
                         true );
  }

  // TODO: extract the following piece of code and write a function that can be used here and in MultiFluidPVTPackageWrapper

  // 5. Compute total fluid mass/molar density and derivatives
  {
    totalDens.value = 0.0;
    totalDens.dPres = 0.0;
    totalDens.dTemp = 0.0;
    for( localIndex jc = 0; jc < numComps; ++jc )
    {
      totalDens.dComp[jc] = 0.0;
    }

    // 5.1. Sum mass/molar fraction/density ratio over all phases to get the inverse of density
    for( localIndex ip = 0; ip < numPhases; ++ip )
    {
      real64 const densInv = 1.0 / phaseDens.value[ip];
      real64 const value = phaseFrac.value[ip] * densInv;

      totalDens.value += value;
      totalDens.dPres += ( phaseFrac.dPres[ip] - value * phaseDens.dPres[ip] ) * densInv;
      totalDens.dTemp += ( phaseFrac.dTemp[ip] - value * phaseDens.dTemp[ip] ) * densInv;

      for( localIndex jc = 0; jc < numComps; ++jc )
      {
        totalDens.dComp[jc] += ( phaseFrac.dComp[ip][jc] - value * phaseDens.dComp[ip][jc] ) * densInv;
      }
    }

    // 5.2. Invert the previous quantity to get actual density
    totalDens.value = 1.0 / totalDens.value;
    real64 const minusDens2 = -totalDens.value * totalDens.value;
    totalDens.dPres *= minusDens2;
    totalDens.dTemp *= minusDens2;

    for( localIndex jc = 0; jc < numComps; ++jc )
    {
      totalDens.dComp[jc] *= minusDens2;
    }
  }
}

template< typename P1DENS, typename P1VISC, typename P2DENS, typename P2VISC, typename FLASH >
GEOSX_HOST_DEVICE inline void
MultiPhaseMultiComponentFluid< P1DENS, P1VISC, P2DENS, P2VISC, FLASH >::KernelWrapper::
  update( localIndex const k,
          localIndex const q,
          real64 const pressure,
          real64 const temperature,
          arraySlice1d< geosx::real64 const, compflow::USD_COMP - 1 > const & composition ) const
{
  compute( pressure,
           temperature,
           composition,
           m_phaseFraction.value[k][q],
           m_phaseFraction.dPressure[k][q],
           m_phaseFraction.dTemperature[k][q],
           m_phaseFraction.dGlobalCompFraction[k][q],
           m_phaseDensity.value[k][q],
           m_phaseDensity.dPressure[k][q],
           m_phaseDensity.dTemperature[k][q],
           m_phaseDensity.dGlobalCompFraction[k][q],
           m_phaseMassDensity.value[k][q],
           m_phaseMassDensity.dPressure[k][q],
           m_phaseMassDensity.dTemperature[k][q],
           m_phaseMassDensity.dGlobalCompFraction[k][q],
           m_phaseViscosity.value[k][q],
           m_phaseViscosity.dPressure[k][q],
           m_phaseViscosity.dTemperature[k][q],
           m_phaseViscosity.dGlobalCompFraction[k][q],
           m_phaseCompFraction.value[k][q],
           m_phaseCompFraction.dPressure[k][q],
           m_phaseCompFraction.dTemperature[k][q],
           m_phaseCompFraction.dGlobalCompFraction[k][q],
           m_totalDensity.value[k][q],
           m_totalDensity.dPressure[k][q],
           m_totalDensity.dTemperature[k][q],
           m_totalDensity.dGlobalCompFraction[k][q] );
}

} //namespace constitutive

} //namespace geosx

#endif //GEOSX_CONSTITUTIVE_FLUID_MULTIPHASEMULTICOMPONENTFLUID_HPP_
