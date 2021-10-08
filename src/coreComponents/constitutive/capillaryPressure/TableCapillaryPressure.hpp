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
 * @file TableCapillaryPressure.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_TABLECAPILLARYPRESSURE_HPP
#define GEOSX_CONSTITUTIVE_TABLECAPILLARYPRESSURE_HPP

#include "constitutive/capillaryPressure/CapillaryPressureBase.hpp"

#include "functions/TableFunction.hpp"

namespace geosx
{
namespace constitutive
{

class TableCapillaryPressure : public CapillaryPressureBase
{
public:

  TableCapillaryPressure( std::string const & name, dataRepository::Group * const parent );

  static std::string catalogName() { return "TableCapillaryPressure"; }

  virtual string getCatalogName() const override { return catalogName(); }

  /// Type of kernel wrapper for in-kernel update
  class KernelWrapper final : public CapillaryPressureBaseUpdate
  {
public:

    KernelWrapper( arrayView1d< TableFunction::KernelWrapper const > const & capPresKernelWrappers,
                   arrayView1d< integer const > const & phaseTypes,
                   arrayView1d< integer const > const & phaseOrder,
                   arrayView3d< real64, cappres::USD_CAPPRES > const & phaseCapPres,
                   arrayView4d< real64, cappres::USD_CAPPRES_DS > const & dPhaseCapPres_dPhaseVolFrac );

    GEOSX_HOST_DEVICE
    virtual void compute( arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFraction,
                          arraySlice1d< real64, cappres::USD_CAPPRES - 2 > const & phaseCapPres,
                          arraySlice2d< real64, cappres::USD_CAPPRES_DS - 2 > const & dPhaseCapPres_dPhaseVolFrac ) const override;

    GEOSX_HOST_DEVICE
    virtual void update( localIndex const k,
                         localIndex const q,
                         arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFraction ) const override;

private:

    /// Array of kernel wrappers for the capillary pressures
    /// Is of size 1 for two-phase flow, and of size 2 for three-phase flow
    arrayView1d< TableFunction::KernelWrapper const > const m_capPresKernelWrappers;

  };

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper();

  struct viewKeyStruct : CapillaryPressureBase::viewKeyStruct
  {
    static constexpr char const * phaseMinVolumeFractionString() { return "phaseMinVolumeFraction"; }
    static constexpr char const * wettingNonWettingCapPresTableNameString() { return "wettingNonWettingCapPressureTableName"; }
    static constexpr char const * wettingIntermediateCapPresTableNameString() { return "wettingIntermediateCapPressureTableName"; }
    static constexpr char const * nonWettingIntermediateCapPresTableNameString() { return "nonWettingIntermediateCapPressureTableName"; }
  };

private:

  virtual void postProcessInput() override;

  virtual void initializePreSubGroups() override;

  /**
   * @brief Create all the table kernel wrappers needed for the simulation (for all the phases present)
   */
  void createAllTableKernelWrappers();

  /**
   * @brief Validate the capillary pressure table provided in input (increasing phase vol frac and cap pressure, etc)
   * @param[in] capPresTable the capillary pressure table (pc vs s) for a given phase)
   * @param[in] capPresMustBeIncreasing flag saying that we expect an increasing cap pressure (otherwise, we expect a decreasing cap
   * pressure)
   */
  void validateCapillaryPressureTable( TableFunction const & capPresTable,
                                       bool const capPresMustBeIncreasing ) const;

  /// Relative permeability table names (one for each phase in the wetting-non-wetting pair)
  string m_wettingNonWettingCapPresTableName;

  /// Relative permeability table names (one for each phase in the wetting-intermediate pair)
  string m_wettingIntermediateCapPresTableName;

  /// Relative permeability table names (one for each phase in the non-wetting-intermediate pair)
  string m_nonWettingIntermediateCapPresTableName;

  /// Relative permeability kernel wrapper for the first pair (wetting-intermediate if NP=3, wetting-non-wetting otherwise)
  array1d< TableFunction::KernelWrapper > m_capPresKernelWrappers;

};

GEOSX_HOST_DEVICE
inline void
TableCapillaryPressure::KernelWrapper::
  compute( arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFraction,
           arraySlice1d< real64, cappres::USD_CAPPRES - 2 > const & phaseCapPres,
           arraySlice2d< real64, cappres::USD_CAPPRES_DS - 2 > const & dPhaseCapPres_dPhaseVolFrac ) const
{
  LvArray::forValuesInSlice( dPhaseCapPres_dPhaseVolFrac, []( real64 & val ){ val = 0.0; } );

  using PT = CapillaryPressureBase::PhaseType;
  integer const ipWater = m_phaseOrder[PT::WATER];
  integer const ipOil   = m_phaseOrder[PT::OIL];
  integer const ipGas   = m_phaseOrder[PT::GAS];

  if( ipWater >= 0 && ipOil >= 0 && ipGas >= 0 )
  {
    // water-oil capillary pressure
    phaseCapPres[ipWater] =
      m_capPresKernelWrappers[0].compute( &(phaseVolFraction)[ipWater],
                                          &(dPhaseCapPres_dPhaseVolFrac)[ipWater][ipWater] );

    // gas-oil capillary pressure
    phaseCapPres[ipGas] =
      m_capPresKernelWrappers[1].compute( &(phaseVolFraction)[ipGas],
                                          &(dPhaseCapPres_dPhaseVolFrac)[ipGas][ipGas] );

    // when pc is on the gas phase, we need to must multiply user input by -1
    // because CompositionalMultiphaseFVM does: pres_gas = pres_oil - pc_og, so we need a negative pc_og
    phaseCapPres[ipGas] *= -1;
    dPhaseCapPres_dPhaseVolFrac[ipGas][ipGas] *= -1;
  }
  else if( ipWater < 0 )
  {
    // put capillary pressure on the non-wetting phase
    phaseCapPres[ipGas] =
      m_capPresKernelWrappers[0].compute( &(phaseVolFraction)[ipGas],
                                          &(dPhaseCapPres_dPhaseVolFrac)[ipGas][ipGas] );

    // when pc is on the gas phase, we need to must multiply user input by -1
    // because CompositionalMultiphaseFVM does: pres_gas = pres_oil - pc_og, so we need a negative pc_og
    phaseCapPres[ipGas] *= -1;
    dPhaseCapPres_dPhaseVolFrac[ipGas][ipGas] *= -1;
  }
  else if( ipOil < 0 )
  {
    // put capillary pressure on the wetting phase
    phaseCapPres[ipWater] =
      m_capPresKernelWrappers[0].compute( &(phaseVolFraction)[ipWater],
                                          &(dPhaseCapPres_dPhaseVolFrac)[ipWater][ipWater] );
  }
  else if( ipGas < 0 )
  {
    // put capillary pressure on the wetting phase
    phaseCapPres[ipWater] =
      m_capPresKernelWrappers[0].compute( &(phaseVolFraction)[ipWater],
                                          &(dPhaseCapPres_dPhaseVolFrac)[ipWater][ipWater] );
  }
}

GEOSX_HOST_DEVICE
inline void
TableCapillaryPressure::KernelWrapper::
  update( localIndex const k,
          localIndex const q,
          arraySlice1d< geosx::real64 const, compflow::USD_PHASE - 1 > const & phaseVolFraction ) const
{
  compute( phaseVolFraction,
           m_phaseCapPressure[k][q],
           m_dPhaseCapPressure_dPhaseVolFrac[k][q] );
}

} // namespace constitutive

} // namespace geosx

#endif //GEOSX_CONSTITUTIVE_TABLECAPILLARYPRESSURE_HPP
