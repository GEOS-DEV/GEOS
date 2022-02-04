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

    friend class DeadOilFluid;

    KernelWrapper( arrayView1d< integer const > phaseTypes,
                   arrayView1d< integer const > phaseOrder,
                   arrayView1d< integer const > hydrocarbonPhaseOrder,
                   arrayView1d< real64 const > surfacePhaseMassDensity,
                   arrayView1d< TableFunction::KernelWrapper const > formationVolFactorTables,
                   arrayView1d< TableFunction::KernelWrapper const > viscosityTables,
                   BlackOilFluidBase::WaterParams const waterParams,
                   arrayView1d< real64 const > componentMolarWeight,
                   bool useMass,
                   PhaseProp::ViewType phaseFraction,
                   PhaseProp::ViewType phaseDensity,
                   PhaseProp::ViewType phaseMassDensity,
                   PhaseProp::ViewType phaseViscosity,
                   PhaseComp::ViewType phaseCompFraction,
                   FluidProp::ViewType totalDensity );

    GEOSX_HOST_DEVICE
    void computeDensities( real64 const pressure,
                           arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseMassDens ) const;

    GEOSX_HOST_DEVICE
    void computeDensities( real64 const pressure,
                           PhaseProp::SliceType const & phaseMassDens ) const;

    GEOSX_HOST_DEVICE
    void computeViscosities( real64 const pressure,
                             arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseVisc ) const;

    GEOSX_HOST_DEVICE
    void computeViscosities( real64 const pressure,
                             PhaseProp::SliceType const & phaseVisc ) const;

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
  computeDensities( real64 const pressure,
                    PhaseProp::SliceType const & phaseMassDens ) const
{
  LvArray::forValuesInSlice( phaseMassDens.dComp, []( real64 & val ){ val = 0.0; } );

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
    phaseMassDens.value[ip] = m_surfacePhaseMassDensity[ip] * fvfInv;
    phaseMassDens.dPres[ip] = -derivative * phaseMassDens.value[ip] * fvfInv;
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
    phaseMassDens.value[ipWater] = m_surfacePhaseMassDensity[ipWater] * denomInv;
    phaseMassDens.dPres[ipWater] = -dDenom_dPres * phaseMassDens.value[ipWater] * denomInv;
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
                      PhaseProp::SliceType const & phaseVisc ) const
{
  LvArray::forValuesInSlice( phaseVisc.dComp, []( real64 & val ){ val = 0.0; } );

  // 1. Hydrocarbon phases: look up in the viscosity tables

  for( integer iph = 0; iph < m_hydrocarbonPhaseOrder.size(); ++iph )
  {
    // get the phase index
    integer const ip = m_hydrocarbonPhaseOrder[iph];
    // interpolate in the table to get the phase viscosity and derivatives
    phaseVisc.value[ip] = m_viscosityTables[iph].compute( &pressure, &(phaseVisc.dPres)[ip] );
  }

  // 2. Water phase: use the constant viscosity provided by the user

  using PT = DeadOilFluid::PhaseType;
  integer const ipWater = m_phaseOrder[PT::WATER];

  // if water is present
  if( ipWater >= 0 )
  {
    // just assign the viscosity value
    phaseVisc.value[ipWater] = m_waterParams.viscosity;
    phaseVisc.dPres[ipWater] = 0.0;
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
  computeTotalDensity( phaseFraction,
                       phaseDensity,
                       totalDens );

}

GEOSX_HOST_DEVICE
inline void
DeadOilFluid::KernelWrapper::
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
  GEOSX_UNUSED_VAR( temperature );

  integer const nComps = numComponents();
  integer const nPhases = numPhases();

  // 1. Read viscosities and formation volume factors from tables, update mass densities
  computeViscosities( pressure,
                      phaseViscosity );
  computeDensities( pressure,
                    phaseMassDensity );

  // 2. Update phaseDens (mass density if useMass == 1, molar density otherwise)
  for( integer ip = 0; ip < nPhases; ++ip )
  {
    real64 const mult = m_useMass ? 1.0 : 1.0 / m_componentMolarWeight[ip];
    phaseDensity.value[ip] = phaseMassDensity.value[ip] * mult;
    phaseDensity.dPres[ip] = phaseMassDensity.dPres[ip] * mult;
    for( integer ic = 0; ic < nComps; ++ic )
    {
      phaseDensity.dComp[ip][ic] = 0.0;
    }
  }

  // 3. Update remaining variables: phaseFrac, phaseCompFrac using Dead-Oil assumptions
  for( integer ip = 0; ip < nPhases; ++ip )
  {
    phaseFraction.value[ip] = composition[ip];
    phaseFraction.dPres[ip] = 0.0;
    for( integer ic = 0; ic < nComps; ++ic )
    {
      phaseFraction.dComp[ip][ic] = (ip == ic) ? 1.0 : 0.0;

      phaseCompFraction.value[ip][ic] = (ip == ic) ? 1.0 : 0.0;
      phaseCompFraction.dPres[ip][ic] = 0.0;
      for( integer jc = 0; jc < nComps; ++jc )
      {
        phaseCompFraction.dComp[ip][ic][jc] = 0.0;
      }
    }
  }

  // 4. Compute total fluid mass/molar density and derivatives
  computeTotalDensity( phaseFraction,
                       phaseDensity,
                       totalDensity );

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
           m_phaseFraction( k, q ),
           m_phaseDensity( k, q ),
           m_phaseMassDensity( k, q ),
           m_phaseViscosity( k, q ),
           m_phaseCompFraction( k, q ),
           m_totalDensity( k, q ) );
}

} //namespace constitutive

} //namespace geosx

#endif //GEOSX_CONSTITUTIVE_FLUID_DEADOILFLUID_HPP_
