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
 * @file MultiPhaseMultiComponentFluid.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_MULTIPHASEMULTICOMPONENTFLUID_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_MULTIPHASEMULTICOMPONENTFLUID_HPP_

#include "constitutive/fluid/MultiFluidBase.hpp"
#include "constitutive/fluid/MultiFluidUtils.hpp"
#include "constitutive/fluid/PVTFunctions/PhillipsBrineDensity.hpp"
#include "constitutive/fluid/PVTFunctions/PhillipsBrineViscosity.hpp"
#include "constitutive/fluid/PVTFunctions/EzrokhiBrineDensity.hpp"
#include "constitutive/fluid/PVTFunctions/EzrokhiBrineViscosity.hpp"
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

    friend class MultiPhaseMultiComponentFluid;

    KernelWrapper( integer p1Index,
                   integer p2Index,
                   P1DENS const & p1DensityWrapper,
                   P1VISC const & p1ViscosityWrapper,
                   P2DENS const & p2DensityWrapper,
                   P2VISC const & p2ViscosityWrapper,
                   FLASH const & flashWrapper,
                   arrayView1d< real64 const > componentMolarWeight,
                   bool const useMass,
                   PhaseProp::ViewType phaseFraction,
                   PhaseProp::ViewType phaseDensity,
                   PhaseProp::ViewType phaseMassDensity,
                   PhaseProp::ViewType phaseViscosity,
                   PhaseComp::ViewType phaseCompFraction,
                   FluidProp::ViewType totalDensity );

    /// Index of the liquid phase
    integer m_p1Index;

    /// Index of the gas phase
    integer m_p2Index;

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
using CO2BrinePhillipsFluid = MultiPhaseMultiComponentFluid< PVTProps::PhillipsBrineDensity,
                                                             PVTProps::PhillipsBrineViscosity,
                                                             PVTProps::SpanWagnerCO2Density,
                                                             PVTProps::FenghourCO2Viscosity,
                                                             PVTProps::CO2Solubility >;

using CO2BrineEzrokhiFluid = MultiPhaseMultiComponentFluid< PVTProps::EzrokhiBrineDensity,
                                                            PVTProps::EzrokhiBrineViscosity,
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
  integer constexpr numComp = 2;
  integer constexpr numPhase = 2;
  integer const ip1 = m_p1Index;
  integer const ip2 = m_p2Index;

  // 1. Convert input mass fractions to mole fractions

  stackArray1d< real64, numComp > compMoleFrac( numComp );
  if( m_useMass )
  {
    // convert mass fractions to mole fractions
    convertToMoleFractions< numComp >( composition,
                                       compMoleFrac );
  }
  else
  {
    for( integer ic = 0; ic < numComp; ++ic )
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

  if( m_useMass )
  {
    // compute the phase molecular weights (ultimately, get that from the PVT function)
    real64 phaseMolecularWeight[numPhase]{};
    real64 phaseMolarDens{};
    m_p1Density.compute( pressure,
                         temperatureInCelsius,
                         phaseCompFraction[ip1].toSliceConst(),
                         phaseMolarDens,
                         false );
    phaseMolecularWeight[ip1] = phaseDensity[ip1] / phaseMolarDens;

    m_p2Density.compute( pressure,
                         temperatureInCelsius,
                         phaseCompFraction[ip2].toSliceConst(),
                         phaseMolarDens,
                         false );
    phaseMolecularWeight[ip2] = phaseDensity[ip2] / phaseMolarDens;

    // copy the densities into the mass densities
    for( integer ip = 0; ip < numPhase; ++ip )
    {
      phaseMassDensity[ip] = phaseDensity[ip];
    }

    // convert mole fractions to mass fractions
    convertToMassFractions< numComp >( phaseMolecularWeight,
                                       phaseFraction,
                                       phaseCompFraction );
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

  // 5. Compute total fluid mass/molar density and derivatives

  computeTotalDensity< numComp, numPhase >( phaseFraction,
                                            phaseDensity,
                                            totalDensity );

}

template< typename P1DENS, typename P1VISC, typename P2DENS, typename P2VISC, typename FLASH >
GEOSX_HOST_DEVICE
inline void
MultiPhaseMultiComponentFluid< P1DENS, P1VISC, P2DENS, P2VISC, FLASH >::KernelWrapper::
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
  integer constexpr numComp = 2;
  integer constexpr numPhase = 2;
  integer const ip1 = m_p1Index;
  integer const ip2 = m_p2Index;

  // 1. Convert input mass fractions to mole fractions and keep derivatives

  stackArray1d< real64, numComp > compMoleFrac( numComp );
  real64 dCompMoleFrac_dCompMassFrac[numComp][numComp]{};

  if( m_useMass )
  {
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

  // 2. Compute phase fractions and phase component fractions

  real64 const temperatureInCelsius = temperature - 273.15;
  m_flash.compute( pressure,
                   temperatureInCelsius,
                   compMoleFrac.toSliceConst(),
                   phaseFraction.value, phaseFraction.dPres, phaseFraction.dTemp, phaseFraction.dComp,
                   phaseCompFraction.value, phaseCompFraction.dPres, phaseCompFraction.dTemp, phaseCompFraction.dComp );

  // 3. Compute phase densities and phase viscosities

  // TODO: these compute functions should take PhaseProp::SliceType in input ...

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

  // 4. Depending on the m_useMass flag, convert to mass variables or simply compute mass density

  // TODO: for now the treatment of molar/mass density requires too many interpolations in the tables, it needs to be fixed
  //       we should modify the PVT functions so that they can return phaseMassDens, phaseDens, and phaseMW in one call

  if( m_useMass )
  {

    // 4.1 Compute the phase molecular weights (ultimately, get that from the PVT function)
    real64 phaseMolecularWeight[numPhase]{};
    real64 dPhaseMolecularWeight_dPres[numPhase]{};
    real64 dPhaseMolecularWeight_dTemp[numPhase]{};
    real64 dPhaseMolecularWeight_dComp[numPhase][numComp]{};

    real64 phaseMolarDens{};
    real64 dPhaseMolarDens_dPres{};
    real64 dPhaseMolarDens_dTemp{};
    stackArray1d< real64, numComp > dPhaseMolarDens_dComp( 2 );
    m_p2Density.compute( pressure,
                         temperatureInCelsius,
                         phaseCompFraction.value[ip2].toSliceConst(), phaseCompFraction.dPres[ip2].toSliceConst(),
                         phaseCompFraction.dTemp[ip2].toSliceConst(), phaseCompFraction.dComp[ip2].toSliceConst(),
                         phaseMolarDens, dPhaseMolarDens_dPres,
                         dPhaseMolarDens_dTemp, dPhaseMolarDens_dComp.toSlice(),
                         false );
    phaseMolecularWeight[ip2] = phaseDensity.value[ip2] / phaseMolarDens;
    dPhaseMolecularWeight_dPres[ip2] = phaseDensity.dPres[ip2] / phaseMolarDens - phaseMolecularWeight[ip2] * dPhaseMolarDens_dPres / phaseMolarDens;
    dPhaseMolecularWeight_dTemp[ip2] = phaseDensity.dTemp[ip2] / phaseMolarDens - phaseMolecularWeight[ip2] * dPhaseMolarDens_dTemp / phaseMolarDens;
    for( integer ic = 0; ic < numComp; ++ic )
    {
      dPhaseMolecularWeight_dComp[ip2][ic] = phaseDensity.dComp[ip2][ic] / phaseMolarDens - phaseMolecularWeight[ip2] * dPhaseMolarDens_dComp[ic] / phaseMolarDens;
    }

    m_p1Density.compute( pressure,
                         temperatureInCelsius,
                         phaseCompFraction.value[ip1].toSliceConst(), phaseCompFraction.dPres[ip1].toSliceConst(),
                         phaseCompFraction.dTemp[ip1].toSliceConst(), phaseCompFraction.dComp[ip1].toSliceConst(),
                         phaseMolarDens, dPhaseMolarDens_dPres,
                         dPhaseMolarDens_dTemp, dPhaseMolarDens_dComp.toSlice(),
                         false );
    phaseMolecularWeight[ip1] = phaseDensity.value[ip1] / phaseMolarDens;
    dPhaseMolecularWeight_dPres[ip1] = phaseDensity.dPres[ip1] / phaseMolarDens - phaseMolecularWeight[ip1] * dPhaseMolarDens_dPres / phaseMolarDens;
    dPhaseMolecularWeight_dTemp[ip1] = phaseDensity.dTemp[ip1] / phaseMolarDens - phaseMolecularWeight[ip1] * dPhaseMolarDens_dTemp / phaseMolarDens;
    for( integer ic = 0; ic < numComp; ++ic )
    {
      dPhaseMolecularWeight_dComp[ip1][ic] = phaseDensity.dComp[ip1][ic] / phaseMolarDens - phaseMolecularWeight[ip1] * dPhaseMolarDens_dComp[ic] / phaseMolarDens;
    }

    // 4.2 Convert the mole fractions to mass fractions
    convertToMassFractions( dCompMoleFrac_dCompMassFrac,
                            phaseMolecularWeight,
                            dPhaseMolecularWeight_dPres,
                            dPhaseMolecularWeight_dTemp,
                            dPhaseMolecularWeight_dComp,
                            phaseFraction,
                            phaseCompFraction,
                            phaseDensity.dComp,
                            phaseViscosity.dComp );


    // 4.3 Copy the phase densities into the phase mass densities
    for( integer ip = 0; ip < numPhase; ++ip )
    {
      phaseMassDensity.value[ip] = phaseDensity.value[ip];
      phaseMassDensity.dPres[ip] = phaseDensity.dPres[ip];
      phaseMassDensity.dTemp[ip] = phaseDensity.dTemp[ip];
      for( integer ic = 0; ic < numComp; ++ic )
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

  // 5. Compute total fluid mass/molar density and derivatives

  computeTotalDensity( phaseFraction,
                       phaseDensity,
                       totalDensity );

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
           m_phaseFraction( k, q ),
           m_phaseDensity( k, q ),
           m_phaseMassDensity( k, q ),
           m_phaseViscosity( k, q ),
           m_phaseCompFraction( k, q ),
           m_totalDensity( k, q ) );
}

} //namespace constitutive

} //namespace geosx

#endif //GEOSX_CONSTITUTIVE_FLUID_MULTIPHASEMULTICOMPONENTFLUID_HPP_
