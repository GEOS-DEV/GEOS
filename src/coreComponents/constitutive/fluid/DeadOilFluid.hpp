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
 * @file DeadOilFluid.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_DEADOILFLUID_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_DEADOILFLUID_HPP_

#include "constitutive/fluid/BlackOilFluidBase.hpp"
#include "functions/TableFunction.hpp"

namespace geosx
{

namespace constitutive
{

/**
 * @brief Kernel wrapper class for DeadOilFluid
 *        This kernel can be called on the GPU
 */
class DeadOilFluidUpdate final : public MultiFluidBaseUpdate
{
public:

  DeadOilFluidUpdate( arrayView1d< integer const > const & phaseTypes,
                      arrayView1d< integer const > const & phaseOrder,
                      arrayView1d< integer const > const & hydrocarbonPhaseOrder,
                      arrayView1d< real64 const > const & surfacePhaseMassDensity,
                      arrayView1d< TableFunction::KernelWrapper const > const & formationVolFactorTables,
                      arrayView1d< TableFunction::KernelWrapper const > const & viscosityTables,
                      real64 const waterRefPressure,
                      real64 const waterFormationVolFactor,
                      real64 const waterCompressibility,
                      real64 const waterViscosity,
                      arrayView1d< real64 const > const & componentMolarWeight,
                      bool useMass,
                      arrayView3d< real64, multifluid::USD_PHASE > const & phaseFraction,
                      arrayView3d< real64, multifluid::USD_PHASE > const & dPhaseFraction_dPressure,
                      arrayView3d< real64, multifluid::USD_PHASE > const & dPhaseFraction_dTemperature,
                      arrayView4d< real64, multifluid::USD_PHASE_DC > const & dPhaseFraction_dGlobalCompFraction,
                      arrayView3d< real64, multifluid::USD_PHASE > const & phaseDensity,
                      arrayView3d< real64, multifluid::USD_PHASE > const & dPhaseDensity_dPressure,
                      arrayView3d< real64, multifluid::USD_PHASE > const & dPhaseDensity_dTemperature,
                      arrayView4d< real64, multifluid::USD_PHASE_DC > const & dPhaseDensity_dGlobalCompFraction,
                      arrayView3d< real64, multifluid::USD_PHASE > const & phaseMassDensity,
                      arrayView3d< real64, multifluid::USD_PHASE > const & dPhaseMassDensity_dPressure,
                      arrayView3d< real64, multifluid::USD_PHASE > const & dPhaseMassDensity_dTemperature,
                      arrayView4d< real64, multifluid::USD_PHASE_DC > const & dPhaseMassDensity_dGlobalCompFraction,
                      arrayView3d< real64, multifluid::USD_PHASE > const & phaseViscosity,
                      arrayView3d< real64, multifluid::USD_PHASE > const & dPhaseViscosity_dPressure,
                      arrayView3d< real64, multifluid::USD_PHASE > const & dPhaseViscosity_dTemperature,
                      arrayView4d< real64, multifluid::USD_PHASE_DC > const & dPhaseViscosity_dGlobalCompFraction,
                      arrayView4d< real64, multifluid::USD_PHASE_COMP > const & phaseCompFraction,
                      arrayView4d< real64, multifluid::USD_PHASE_COMP > const & dPhaseCompFraction_dPressure,
                      arrayView4d< real64, multifluid::USD_PHASE_COMP > const & dPhaseCompFraction_dTemperature,
                      arrayView5d< real64, multifluid::USD_PHASE_COMP_DC > const & dPhaseCompFraction_dGlobalCompFraction,
                      arrayView2d< real64, multifluid::USD_FLUID > const & totalDensity,
                      arrayView2d< real64, multifluid::USD_FLUID > const & dTotalDensity_dPressure,
                      arrayView2d< real64, multifluid::USD_FLUID > const & dTotalDensity_dTemperature,
                      arrayView3d< real64, multifluid::USD_FLUID_DC > const & dTotalDensity_dGlobalCompFraction )
    : MultiFluidBaseUpdate( componentMolarWeight,
                            useMass,
                            phaseFraction,
                            dPhaseFraction_dPressure,
                            dPhaseFraction_dTemperature,
                            dPhaseFraction_dGlobalCompFraction,
                            phaseDensity,
                            dPhaseDensity_dPressure,
                            dPhaseDensity_dTemperature,
                            dPhaseDensity_dGlobalCompFraction,
                            phaseMassDensity,
                            dPhaseMassDensity_dPressure,
                            dPhaseMassDensity_dTemperature,
                            dPhaseMassDensity_dGlobalCompFraction,
                            phaseViscosity,
                            dPhaseViscosity_dPressure,
                            dPhaseViscosity_dTemperature,
                            dPhaseViscosity_dGlobalCompFraction,
                            phaseCompFraction,
                            dPhaseCompFraction_dPressure,
                            dPhaseCompFraction_dTemperature,
                            dPhaseCompFraction_dGlobalCompFraction,
                            totalDensity,
                            dTotalDensity_dPressure,
                            dTotalDensity_dTemperature,
                            dTotalDensity_dGlobalCompFraction ),
    m_phaseTypes( phaseTypes ),
    m_phaseOrder( phaseOrder ),
    m_hydrocarbonPhaseOrder( hydrocarbonPhaseOrder ),
    m_surfacePhaseMassDensity( surfacePhaseMassDensity ),
    m_formationVolFactorTables( formationVolFactorTables ),
    m_viscosityTables( viscosityTables ),
    m_waterRefPressure( waterRefPressure ),
    m_waterFormationVolFactor( waterFormationVolFactor ),
    m_waterCompressibility( waterCompressibility ),
    m_waterViscosity( waterViscosity )
  {}

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
  virtual void update( localIndex const k,
                       localIndex const q,
                       real64 const pressure,
                       real64 const temperature,
                       arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition ) const override
  {
    compute( pressure,
             temperature,
             composition,
             m_phaseFraction[k][q],
             m_dPhaseFraction_dPressure[k][q],
             m_dPhaseFraction_dTemperature[k][q],
             m_dPhaseFraction_dGlobalCompFraction[k][q],
             m_phaseDensity[k][q],
             m_dPhaseDensity_dPressure[k][q],
             m_dPhaseDensity_dTemperature[k][q],
             m_dPhaseDensity_dGlobalCompFraction[k][q],
             m_phaseMassDensity[k][q],
             m_dPhaseMassDensity_dPressure[k][q],
             m_dPhaseMassDensity_dTemperature[k][q],
             m_dPhaseMassDensity_dGlobalCompFraction[k][q],
             m_phaseViscosity[k][q],
             m_dPhaseViscosity_dPressure[k][q],
             m_dPhaseViscosity_dTemperature[k][q],
             m_dPhaseViscosity_dGlobalCompFraction[k][q],
             m_phaseCompFraction[k][q],
             m_dPhaseCompFraction_dPressure[k][q],
             m_dPhaseCompFraction_dTemperature[k][q],
             m_dPhaseCompFraction_dGlobalCompFraction[k][q],
             m_totalDensity[k][q],
             m_dTotalDensity_dPressure[k][q],
             m_dTotalDensity_dTemperature[k][q],
             m_dTotalDensity_dGlobalCompFraction[k][q] );
  }

private:

  GEOSX_HOST_DEVICE
  void computeDensities( real64 pressure,
                         arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseMassDens ) const;

  GEOSX_HOST_DEVICE
  void computeDensities( real64 pressure,
                         arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseMassDens,
                         arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & dPhaseMassDens_dPres,
                         arraySlice2d< real64, multifluid::USD_PHASE_DC - 2 > const & dPhaseMassDens_dGlobalCompFrac ) const;

  GEOSX_HOST_DEVICE
  void computeViscosities( real64 pressure,
                           arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseVisc ) const;

  GEOSX_HOST_DEVICE
  void computeViscosities( real64 pressure,
                           arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseVisc,
                           arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & dPhaseVisc_dPressure,
                           arraySlice2d< real64, multifluid::USD_PHASE_DC - 2 > const & dPhaseVisc_dGlobalCompFraction ) const;

  /// Phase ordering info
  arrayView1d< integer const > m_phaseTypes;
  arrayView1d< integer const > m_phaseOrder;
  arrayView1d< integer const > m_hydrocarbonPhaseOrder;

  /// Surface mass density for each phase
  arrayView1d< real64 const > m_surfacePhaseMassDensity;

  /// Table kernel wrappers to interpolate in the oil and gas (B vs p) tables
  arrayView1d< TableFunction::KernelWrapper const > m_formationVolFactorTables;

  /// Table kernel wrappers to interpolate in the oil and gas (\mu vs p) tables
  arrayView1d< TableFunction::KernelWrapper const > m_viscosityTables;

  /// Water reference pressure
  real64 m_waterRefPressure;

  /// Water formation volume factor
  real64 m_waterFormationVolFactor;

  /// Water compressibility
  real64 m_waterCompressibility;

  /// Water viscosity
  real64 m_waterViscosity;

};

class DeadOilFluid : public BlackOilFluidBase
{
public:

  DeadOilFluid( string const & name, Group * const parent );

  virtual ~DeadOilFluid() override = default;

  virtual std::unique_ptr< ConstitutiveBase >
  deliverClone( string const & name,
                Group * const parent ) const override;

  static string catalogName() { return "DeadOilFluid"; }

  virtual string getCatalogName() const override { return catalogName(); }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = DeadOilFluidUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper()
  {
    return KernelWrapper( m_phaseTypes,
                          m_phaseOrder,
                          m_hydrocarbonPhaseOrder,
                          m_surfacePhaseMassDensity,
                          m_formationVolFactorTables,
                          m_viscosityTables,
                          m_waterRefPressure,
                          m_waterFormationVolFactor,
                          m_waterCompressibility,
                          m_waterViscosity,
                          m_componentMolarWeight,
                          m_useMass,
                          m_phaseFraction,
                          m_dPhaseFraction_dPressure,
                          m_dPhaseFraction_dTemperature,
                          m_dPhaseFraction_dGlobalCompFraction,
                          m_phaseDensity,
                          m_dPhaseDensity_dPressure,
                          m_dPhaseDensity_dTemperature,
                          m_dPhaseDensity_dGlobalCompFraction,
                          m_phaseMassDensity,
                          m_dPhaseMassDensity_dPressure,
                          m_dPhaseMassDensity_dTemperature,
                          m_dPhaseMassDensity_dGlobalCompFraction,
                          m_phaseViscosity,
                          m_dPhaseViscosity_dPressure,
                          m_dPhaseViscosity_dTemperature,
                          m_dPhaseViscosity_dGlobalCompFraction,
                          m_phaseCompFraction,
                          m_dPhaseCompFraction_dPressure,
                          m_dPhaseCompFraction_dTemperature,
                          m_dPhaseCompFraction_dGlobalCompFraction,
                          m_totalDensity,
                          m_dTotalDensity_dPressure,
                          m_dTotalDensity_dTemperature,
                          m_dTotalDensity_dGlobalCompFraction );
  }

private:

  /**
   * @brief Use the TableFunctions provided by the user to get the PVT data
   */
  virtual void useProvidedTableFunctions() override;

  /**
   * @brief Read all the PVT table provided by the user in Eclipse format
   */
  virtual void readInputDataFromPVTFiles() override;

};

GEOSX_HOST_DEVICE
inline void
DeadOilFluidUpdate::computeDensities( real64 const pressure,
                                      arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseMassDens ) const
{
  real64 fvf = 0.0;
  real64 derivative = 0.0;

  // 1. Hydrocarbon phases: look up in the formation vol factor tables

  for( localIndex iph = 0; iph < m_hydrocarbonPhaseOrder.size(); ++iph )
  {
    // get the phase index
    localIndex const ip = m_hydrocarbonPhaseOrder[iph];
    // interpolate in the table to get the phase formation vol factor, and discard the derivative
    m_formationVolFactorTables[iph].compute( &pressure, fvf, &derivative );

    // we are ready to update the densities
    phaseMassDens[ip] = m_surfacePhaseMassDensity[ip] / fvf;
  }

  // 2. Water phase: use the constant formation volume factor and compressibility provided by the user

  using PT = DeadOilFluid::PhaseType;
  integer const ipWater = m_phaseOrder[PT::WATER];

  // if water is present
  if( ipWater >= 0 )
  {
    // note: double check, but I think std::exp is allowed in kernel
    real64 const denom = m_waterFormationVolFactor * std::exp( -m_waterCompressibility * ( pressure - m_waterRefPressure ) );
    phaseMassDens[ipWater] = m_surfacePhaseMassDensity[ipWater] / denom;
  }
}

GEOSX_HOST_DEVICE
inline void
DeadOilFluidUpdate::computeDensities( real64 pressure,
                                      arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseMassDens,
                                      arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & dPhaseMassDens_dPres,
                                      arraySlice2d< real64, multifluid::USD_PHASE_DC - 2 > const & dPhaseMassDens_dGlobalCompFrac ) const
{
  real64 fvf = 0.0;
  real64 derivative = 0.0;

  LvArray::forValuesInSlice( dPhaseMassDens_dGlobalCompFrac, []( real64 & val ){ val = 0.0; } );

  // 1. Hydrocarbon phases: look up in the formation vol factor tables

  for( localIndex iph = 0; iph < m_hydrocarbonPhaseOrder.size(); ++iph )
  {
    // get the phase index
    localIndex const ip = m_hydrocarbonPhaseOrder[iph];
    // interpolate in the table to get the phase formation vol factor and its derivative wrt pressure
    m_formationVolFactorTables[iph].compute( &pressure, fvf, &derivative );

    // we are ready to update the densities
    real64 const fvfInv = 1.0 / fvf;
    phaseMassDens[ip] = m_surfacePhaseMassDensity[ip] * fvfInv;
    dPhaseMassDens_dPres[ip] = -derivative * phaseMassDens[ip] * fvfInv;
  }

  // 2. Water phase: use the constant formation volume factor and compressibility provided by the user

  using PT = DeadOilFluid::PhaseType;
  integer const ipWater = m_phaseOrder[PT::WATER];

  // if water is present
  if( ipWater >= 0 )
  {
    // note: double check, but I think std::exp is allowed in kernel
    real64 const expCompDeltaPres = std::exp( -m_waterCompressibility * ( pressure - m_waterRefPressure ) );
    real64 const dExpCompDeltaPres_dPres = -m_waterCompressibility * expCompDeltaPres;
    real64 const denom = m_waterFormationVolFactor * expCompDeltaPres;
    real64 const dDenom_dPres = m_waterFormationVolFactor * dExpCompDeltaPres_dPres;
    real64 const denomInv = 1.0 / denom;
    phaseMassDens[ipWater] = m_surfacePhaseMassDensity[ipWater] * denomInv;
    dPhaseMassDens_dPres[ipWater] = -dDenom_dPres * phaseMassDens[ipWater] * denomInv;
  }
}

GEOSX_HOST_DEVICE
inline void
DeadOilFluidUpdate::computeViscosities( real64 pressure,
                                        arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseVisc ) const
{
  // this derivative will be discarded
  real64 derivative = 0.0;

  // 1. Hydrocarbon phases: look up in the viscosity tables

  for( localIndex iph = 0; iph < m_hydrocarbonPhaseOrder.size(); ++iph )
  {
    // get the phase index
    localIndex const ip = m_hydrocarbonPhaseOrder[iph];
    // interpolate in the table to get the phase viscosity, and discard the derivative
    m_viscosityTables[iph].compute( &pressure, phaseVisc[ip], &derivative );
  }

  // 2. Water phase: use the constant viscosity provided by the user

  using PT = DeadOilFluid::PhaseType;
  integer const ipWater = m_phaseOrder[PT::WATER];

  // if water is present
  if( ipWater >= 0 )
  {
    // just assign the viscosity value
    phaseVisc[ipWater] = m_waterViscosity;
  }
}

GEOSX_HOST_DEVICE
inline void
DeadOilFluidUpdate::computeViscosities( real64 pressure,
                                        arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseVisc,
                                        arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & dPhaseVisc_dPres,
                                        arraySlice2d< real64, multifluid::USD_PHASE_DC - 2 > const & dPhaseVisc_dGlobalCompFrac ) const
{
  LvArray::forValuesInSlice( dPhaseVisc_dGlobalCompFrac, []( real64 & val ){ val = 0.0; } );

  // 1. Hydrocarbon phases: look up in the viscosity tables

  for( localIndex iph = 0; iph < m_hydrocarbonPhaseOrder.size(); ++iph )
  {
    // get the phase index
    localIndex const ip = m_hydrocarbonPhaseOrder[iph];
    // interpolate in the table to get the phase viscosity and derivatives
    m_viscosityTables[iph].compute( &pressure, phaseVisc[ip], &(dPhaseVisc_dPres)[ip] );
  }

  // 2. Water phase: use the constant viscosity provided by the user

  using PT = DeadOilFluid::PhaseType;
  integer const ipWater = m_phaseOrder[PT::WATER];

  // if water is present
  if( ipWater >= 0 )
  {
    // just assign the viscosity value
    phaseVisc[ipWater] = m_waterViscosity;
    dPhaseVisc_dPres[ipWater] = 0.0;
  }
}

GEOSX_HOST_DEVICE
inline void
DeadOilFluidUpdate::compute( real64 pressure,
                             real64 temperature,
                             arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition,
                             arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseFraction,
                             arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseDensity,
                             arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseMassDensity,
                             arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseViscosity,
                             arraySlice2d< real64, multifluid::USD_PHASE_COMP-2 > const & phaseCompFraction,
                             real64 & totalDens ) const
{
  GEOSX_UNUSED_VAR( temperature );
  localIndex const nComps = m_componentMolarWeight.size();
  localIndex const nPhases = m_phaseTypes.size();

  // 1. Read viscosities and formation volume factors from tables, update mass densities
  computeViscosities( pressure, phaseViscosity );
  computeDensities( pressure, phaseMassDensity );

  // 2. Update phaseDens (mass density if useMass == 1, molar density otherwise)
  for( localIndex ip = 0; ip < nPhases; ++ip )
  {
    real64 const mult = m_useMass ? 1.0 : 1.0 / m_componentMolarWeight[ip];
    phaseDensity[ip] = phaseMassDensity[ip] * mult;
  }

  // 3. Update remaining variables: phaseFrac, phaseCompFrac using Dead-Oil assumptions
  for( localIndex ip = 0; ip < nPhases; ++ip )
  {
    phaseFraction[ip] = composition[ip];
    for( localIndex ic = 0; ic < nComps; ++ic )
    {
      phaseCompFraction[ip][ic] = (ip == ic) ? 1.0 : 0.0;
    }
  }

  // 4. Compute total fluid mass/molar density
  totalDens = 0.0;

  // 4.1. Sum mass/molar fraction/density ratio over all phases to get the inverse of density
  for( localIndex ip = 0; ip < nPhases; ++ip )
  {
    totalDens += phaseFraction[ip] / phaseDensity[ip];
  }
  // 4.2. Invert the previous quantity to get actual density
  totalDens = 1.0 / totalDens;
}

GEOSX_HOST_DEVICE
inline void
DeadOilFluidUpdate::compute( real64 pressure,
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
  GEOSX_UNUSED_VAR( temperature );
  GEOSX_UNUSED_VAR( dPhaseFraction_dTemperature );
  GEOSX_UNUSED_VAR( dPhaseDensity_dTemperature );
  GEOSX_UNUSED_VAR( dPhaseMassDensity_dTemperature );
  GEOSX_UNUSED_VAR( dPhaseViscosity_dTemperature );
  GEOSX_UNUSED_VAR( dPhaseCompFraction_dTemperature );
  GEOSX_UNUSED_VAR( dTotalDensity_dTemperature );

  localIndex const nComps = numComponents();
  localIndex const nPhases = numPhases();

  // 1. Read viscosities and formation volume factors from tables, update mass densities
  computeViscosities( pressure,
                      phaseViscosity,
                      dPhaseViscosity_dPressure,
                      dPhaseViscosity_dGlobalCompFraction );
  computeDensities( pressure,
                    phaseMassDensity,
                    dPhaseMassDensity_dPressure,
                    dPhaseMassDensity_dGlobalCompFraction );

  // 2. Update phaseDens (mass density if useMass == 1, molar density otherwise)
  for( localIndex ip = 0; ip < nPhases; ++ip )
  {
    real64 const mult = m_useMass ? 1.0 : 1.0 / m_componentMolarWeight[ip];
    phaseDensity[ip] = phaseMassDensity[ip] * mult;
    dPhaseDensity_dPressure[ip] = dPhaseMassDensity_dPressure[ip] * mult;
    for( localIndex ic = 0; ic < nComps; ++ic )
    {
      dPhaseDensity_dGlobalCompFraction[ip][ic] = 0.0;
    }
  }

  // 3. Update remaining variables: phaseFrac, phaseCompFrac using Dead-Oil assumptions
  for( localIndex ip = 0; ip < nPhases; ++ip )
  {
    phaseFraction[ip] = composition[ip];
    dPhaseFraction_dPressure[ip] = 0.0;
    for( localIndex ic = 0; ic < nComps; ++ic )
    {
      dPhaseFraction_dGlobalCompFraction[ip][ic] = (ip == ic) ? 1.0 : 0.0;

      phaseCompFraction[ip][ic] = (ip == ic) ? 1.0 : 0.0;
      dPhaseCompFraction_dPressure[ip][ic] = 0.0;
      for( localIndex jc = 0; jc < nComps; ++jc )
      {
        dPhaseCompFraction_dGlobalCompFraction[ip][ic][jc] = 0.0;
      }
    }
  }

  // 4. Compute total fluid mass/molar density and derivatives
  totalDensity = 0.0;
  dTotalDensity_dPressure = 0.0;
  LvArray::forValuesInSlice( dTotalDensity_dGlobalCompFraction, []( real64 & val ){ val = 0.0; } );

  // 4.1. Sum mass/molar fraction/density ratio over all phases to get the inverse of density
  for( localIndex ip = 0; ip < nPhases; ++ip )
  {
    real64 const densInv = 1.0 / phaseDensity[ip];
    real64 const value = phaseFraction[ip] * densInv;

    totalDensity += value;
    dTotalDensity_dPressure += ( dPhaseFraction_dPressure[ip] - value * dPhaseDensity_dPressure[ip] ) * densInv;
    for( localIndex ic = 0; ic < nComps; ++ic )
    {
      dTotalDensity_dGlobalCompFraction[ic] += ( dPhaseFraction_dGlobalCompFraction[ip][ic]
                                                 - value * dPhaseDensity_dGlobalCompFraction[ip][ic] ) * densInv;
    }
  }

  // 4.2. Invert the previous quantity to get actual density
  totalDensity = 1.0 / totalDensity;
  real64 const minusDens2 = -totalDensity * totalDensity;
  dTotalDensity_dPressure *= minusDens2;
  for( localIndex ic = 0; ic < nComps; ++ic )
  {
    dTotalDensity_dGlobalCompFraction[ic] *= minusDens2;
  }
}

} //namespace constitutive

} //namespace geosx

#endif //GEOSX_CONSTITUTIVE_FLUID_DEADOILFLUID_HPP_
