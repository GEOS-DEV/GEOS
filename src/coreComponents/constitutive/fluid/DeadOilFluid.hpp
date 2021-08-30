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

namespace geosx
{

namespace constitutive
{

class DeadOilFluid : public BlackOilFluidBase
{
public:

  DeadOilFluid( string const & name, Group * const parent );

  static string catalogName() { return "DeadOilFluid"; }

  virtual string getCatalogName() const override { return catalogName(); }

  /**
   * @brief Kernel wrapper class for DeadOilFluid
   *        This kernel can be called on the GPU
   */
  class KernelWrapper final : public BlackOilFluidBase::KernelWrapper
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

    friend class DeadOilFluid;

    KernelWrapper( arrayView1d< integer const > const & phaseTypes,
                   arrayView1d< integer const > const & phaseOrder,
                   arrayView1d< integer const > const & hydrocarbonPhaseOrder,
                   arrayView1d< real64 const > const & surfacePhaseMassDensity,
                   arrayView1d< TableFunction::KernelWrapper const > const & formationVolFactorTables,
                   arrayView1d< TableFunction::KernelWrapper const > const & viscosityTables,
                   BlackOilFluidBase::WaterParams waterParams,
                   arrayView1d< real64 const > const & componentMolarWeight,
                   bool useMass,
                   PhasePropViews const & phaseFraction,
                   PhasePropViews const & phaseDensity,
                   PhasePropViews const & phaseMassDensity,
                   PhasePropViews const & phaseViscosity,
                   PhaseCompViews const & phaseCompFraction,
                   FluidPropViews const & totalDensity );

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

  };

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper();

private:

  /**
   * @brief Use the TableFunctions provided by the user to get the PVT data
   */
  virtual void readInputDataFromTableFunctions() override;

  /**
   * @brief Read all the PVT table provided by the user in Eclipse format
   */
  virtual void readInputDataFromPVTFiles() override;

};

GEOSX_HOST_DEVICE
inline void
DeadOilFluid::KernelWrapper::
  computeDensities( real64 const pressure,
                    arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseMassDens ) const
{
  // 1. Hydrocarbon phases: look up in the formation vol factor tables

  for( integer iph = 0; iph < m_hydrocarbonPhaseOrder.size(); ++iph )
  {
    // get the phase index
    integer const ip = m_hydrocarbonPhaseOrder[iph];
    // interpolate in the table to get the phase formation vol factor, and discard the derivative
    real64 const fvf = m_formationVolFactorTables[iph].compute( &pressure );

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
    real64 const denom = m_waterParams.formationVolFactor * std::exp( -m_waterParams.compressibility * ( pressure - m_waterParams.referencePressure ) );
    phaseMassDens[ipWater] = m_surfacePhaseMassDensity[ipWater] / denom;
  }
}

GEOSX_HOST_DEVICE
inline void
DeadOilFluid::KernelWrapper::
  computeDensities( real64 pressure,
                    arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseMassDens,
                    arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & dPhaseMassDens_dPres,
                    arraySlice2d< real64, multifluid::USD_PHASE_DC - 2 > const & dPhaseMassDens_dGlobalCompFrac ) const
{
  LvArray::forValuesInSlice( dPhaseMassDens_dGlobalCompFrac, []( real64 & val ){ val = 0.0; } );

  // 1. Hydrocarbon phases: look up in the formation vol factor tables

  for( integer iph = 0; iph < m_hydrocarbonPhaseOrder.size(); ++iph )
  {
    // get the phase index
    integer const ip = m_hydrocarbonPhaseOrder[iph];
    // interpolate in the table to get the phase formation vol factor and its derivative wrt pressure
    real64 derivative;
    real64 const fvf = m_formationVolFactorTables[iph].compute( &pressure, &derivative );

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
    real64 const expCompDeltaPres = std::exp( -m_waterParams.compressibility * ( pressure - m_waterParams.referencePressure ) );
    real64 const dExpCompDeltaPres_dPres = -m_waterParams.compressibility * expCompDeltaPres;
    real64 const denom = m_waterParams.formationVolFactor * expCompDeltaPres;
    real64 const dDenom_dPres = m_waterParams.formationVolFactor * dExpCompDeltaPres_dPres;
    real64 const denomInv = 1.0 / denom;
    phaseMassDens[ipWater] = m_surfacePhaseMassDensity[ipWater] * denomInv;
    dPhaseMassDens_dPres[ipWater] = -dDenom_dPres * phaseMassDens[ipWater] * denomInv;
  }
}

GEOSX_HOST_DEVICE
inline void
DeadOilFluid::KernelWrapper::
  computeViscosities( real64 const pressure,
                      arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseVisc ) const
{
  // 1. Hydrocarbon phases: look up in the viscosity tables

  for( integer iph = 0; iph < m_hydrocarbonPhaseOrder.size(); ++iph )
  {
    // get the phase index
    integer const ip = m_hydrocarbonPhaseOrder[iph];
    // interpolate in the table to get the phase viscosity, and discard the derivative
    phaseVisc[ip] = m_viscosityTables[iph].compute( &pressure );
  }

  // 2. Water phase: use the constant viscosity provided by the user

  using PT = DeadOilFluid::PhaseType;
  integer const ipWater = m_phaseOrder[PT::WATER];

  // if water is present
  if( ipWater >= 0 )
  {
    // just assign the viscosity value
    phaseVisc[ipWater] = m_waterParams.viscosity;
  }
}

GEOSX_HOST_DEVICE
inline void
DeadOilFluid::KernelWrapper::
  computeViscosities( real64 const pressure,
                      arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseVisc,
                      arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & dPhaseVisc_dPres,
                      arraySlice2d< real64, multifluid::USD_PHASE_DC - 2 > const & dPhaseVisc_dGlobalCompFrac ) const
{
  LvArray::forValuesInSlice( dPhaseVisc_dGlobalCompFrac, []( real64 & val ){ val = 0.0; } );

  // 1. Hydrocarbon phases: look up in the viscosity tables

  for( integer iph = 0; iph < m_hydrocarbonPhaseOrder.size(); ++iph )
  {
    // get the phase index
    integer const ip = m_hydrocarbonPhaseOrder[iph];
    // interpolate in the table to get the phase viscosity and derivatives
    phaseVisc[ip] = m_viscosityTables[iph].compute( &pressure, &(dPhaseVisc_dPres)[ip] );
  }

  // 2. Water phase: use the constant viscosity provided by the user

  using PT = DeadOilFluid::PhaseType;
  integer const ipWater = m_phaseOrder[PT::WATER];

  // if water is present
  if( ipWater >= 0 )
  {
    // just assign the viscosity value
    phaseVisc[ipWater] = m_waterParams.viscosity;
    dPhaseVisc_dPres[ipWater] = 0.0;
  }
}

GEOSX_HOST_DEVICE
inline void
DeadOilFluid::KernelWrapper::
  compute( real64 const pressure,
           real64 const temperature,
           arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition,
           arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseFraction,
           arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseDensity,
           arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseMassDensity,
           arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseViscosity,
           arraySlice2d< real64, multifluid::USD_PHASE_COMP - 2 > const & phaseCompFraction,
           real64 & totalDens ) const
{
  GEOSX_UNUSED_VAR( temperature );
  integer const nComps = numComponents();
  integer const nPhases = numPhases();

  // 1. Read viscosities and formation volume factors from tables, update mass densities
  computeViscosities( pressure, phaseViscosity );
  computeDensities( pressure, phaseMassDensity );

  // 2. Update phaseDens (mass density if useMass == 1, molar density otherwise)
  for( integer ip = 0; ip < nPhases; ++ip )
  {
    real64 const mult = m_useMass ? 1.0 : 1.0 / m_componentMolarWeight[ip];
    phaseDensity[ip] = phaseMassDensity[ip] * mult;
  }

  // 3. Update remaining variables: phaseFrac, phaseCompFrac using Dead-Oil assumptions
  for( integer ip = 0; ip < nPhases; ++ip )
  {
    phaseFraction[ip] = composition[ip];
    for( integer ic = 0; ic < nComps; ++ic )
    {
      phaseCompFraction[ip][ic] = (ip == ic) ? 1.0 : 0.0;
    }
  }

  // 4. Compute total fluid mass/molar density
  totalDens = 0.0;

  // 4.1. Sum mass/molar fraction/density ratio over all phases to get the inverse of density
  for( integer ip = 0; ip < nPhases; ++ip )
  {
    totalDens += phaseFraction[ip] / phaseDensity[ip];
  }
  // 4.2. Invert the previous quantity to get actual density
  totalDens = 1.0 / totalDens;
}

GEOSX_HOST_DEVICE
inline void
DeadOilFluid::KernelWrapper::
  compute( real64 const pressure,
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
           arraySlice2d< real64, multifluid::USD_PHASE_COMP - 2 > const & phaseCompFraction,
           arraySlice2d< real64, multifluid::USD_PHASE_COMP - 2 > const & dPhaseCompFraction_dPressure,
           arraySlice2d< real64, multifluid::USD_PHASE_COMP - 2 > const & dPhaseCompFraction_dTemperature,
           arraySlice3d< real64, multifluid::USD_PHASE_COMP_DC - 2 > const & dPhaseCompFraction_dGlobalCompFraction,
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

  integer const nComps = numComponents();
  integer const nPhases = numPhases();

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
  for( integer ip = 0; ip < nPhases; ++ip )
  {
    real64 const mult = m_useMass ? 1.0 : 1.0 / m_componentMolarWeight[ip];
    phaseDensity[ip] = phaseMassDensity[ip] * mult;
    dPhaseDensity_dPressure[ip] = dPhaseMassDensity_dPressure[ip] * mult;
    for( integer ic = 0; ic < nComps; ++ic )
    {
      dPhaseDensity_dGlobalCompFraction[ip][ic] = 0.0;
    }
  }

  // 3. Update remaining variables: phaseFrac, phaseCompFrac using Dead-Oil assumptions
  for( integer ip = 0; ip < nPhases; ++ip )
  {
    phaseFraction[ip] = composition[ip];
    dPhaseFraction_dPressure[ip] = 0.0;
    for( integer ic = 0; ic < nComps; ++ic )
    {
      dPhaseFraction_dGlobalCompFraction[ip][ic] = (ip == ic) ? 1.0 : 0.0;

      phaseCompFraction[ip][ic] = (ip == ic) ? 1.0 : 0.0;
      dPhaseCompFraction_dPressure[ip][ic] = 0.0;
      for( integer jc = 0; jc < nComps; ++jc )
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
  for( integer ip = 0; ip < nPhases; ++ip )
  {
    real64 const densInv = 1.0 / phaseDensity[ip];
    real64 const value = phaseFraction[ip] * densInv;

    totalDensity += value;
    dTotalDensity_dPressure += ( dPhaseFraction_dPressure[ip] - value * dPhaseDensity_dPressure[ip] ) * densInv;
    for( integer ic = 0; ic < nComps; ++ic )
    {
      dTotalDensity_dGlobalCompFraction[ic] += ( dPhaseFraction_dGlobalCompFraction[ip][ic]
                                                 - value * dPhaseDensity_dGlobalCompFraction[ip][ic] ) * densInv;
    }
  }

  // 4.2. Invert the previous quantity to get actual density
  totalDensity = 1.0 / totalDensity;
  real64 const minusDens2 = -totalDensity * totalDensity;
  dTotalDensity_dPressure *= minusDens2;
  for( integer ic = 0; ic < nComps; ++ic )
  {
    dTotalDensity_dGlobalCompFraction[ic] *= minusDens2;
  }
}

GEOSX_HOST_DEVICE
inline void
DeadOilFluid::KernelWrapper::
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

#endif //GEOSX_CONSTITUTIVE_FLUID_DEADOILFLUID_HPP_
