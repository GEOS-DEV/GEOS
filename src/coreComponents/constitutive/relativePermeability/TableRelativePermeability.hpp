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
 * @file TableRelativePermeability.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_TABLERELATIVEPERMEABILITY_HPP
#define GEOSX_CONSTITUTIVE_TABLERELATIVEPERMEABILITY_HPP

#include "constitutive/relativePermeability/RelativePermeabilityBase.hpp"
#include "constitutive/relativePermeability/RelativePermeabilityInterpolators.hpp"
#include "functions/TableFunction.hpp"

namespace geosx
{
namespace constitutive
{

class TableRelativePermeabilityUpdate final : public RelativePermeabilityBaseUpdate
{
public:

  TableRelativePermeabilityUpdate( arrayView1d< TableFunction::KernelWrapper const > const & waterOilRelPermTableKernelWrappers,
                                   arrayView1d< TableFunction::KernelWrapper const > const & gasOilRelPermTableKernelWrappers,
                                   arrayView1d< real64 const > const & phaseMinVolumeFraction,
                                   arrayView1d< integer const > const & phaseTypes,
                                   arrayView1d< integer const > const & phaseOrder,
                                   arrayView3d< real64 > const & phaseRelPerm,
                                   arrayView4d< real64 > const & dPhaseRelPerm_dPhaseVolFrac )
    : RelativePermeabilityBaseUpdate( phaseTypes,
                                      phaseOrder,
                                      phaseRelPerm,
                                      dPhaseRelPerm_dPhaseVolFrac ),
    m_waterOilRelPermTableKernelWrappers( waterOilRelPermTableKernelWrappers ),
    m_gasOilRelPermTableKernelWrappers( gasOilRelPermTableKernelWrappers ),
    m_phaseMinVolumeFraction( phaseMinVolumeFraction )
  {}

  /// Default copy constructor
  TableRelativePermeabilityUpdate( TableRelativePermeabilityUpdate const & ) = default;

  /// Default move constructor
  TableRelativePermeabilityUpdate( TableRelativePermeabilityUpdate && ) = default;

  /// Deleted copy assignment operator
  TableRelativePermeabilityUpdate & operator=( TableRelativePermeabilityUpdate const & ) = delete;

  /// Deleted move assignment operator
  TableRelativePermeabilityUpdate & operator=( TableRelativePermeabilityUpdate && ) = delete;

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  virtual void compute( arraySlice1d< real64 const > const & phaseVolFraction,
                        arraySlice1d< real64 > const & phaseRelPerm,
                        arraySlice2d< real64 > const & dPhaseRelPerm_dPhaseVolFrac ) const override;

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  virtual void update( localIndex const k,
                       localIndex const q,
                       arraySlice1d< real64 const > const & phaseVolFraction ) const override
  {
    compute( phaseVolFraction,
             m_phaseRelPerm[k][q],
             m_dPhaseRelPerm_dPhaseVolFrac[k][q] );
  }

private:

  /// Kernel wrappers for the water-oil relative permeabilities
  arrayView1d< TableFunction::KernelWrapper const > m_waterOilRelPermTableKernelWrappers;

  /// Kernel wrappers for the gas-oil relative permeabilities
  arrayView1d< TableFunction::KernelWrapper const > m_gasOilRelPermTableKernelWrappers;

  /// Minimum volume fraction for each phase (deduced from the table)
  arrayView1d< real64 const > m_phaseMinVolumeFraction;

};

class TableRelativePermeability : public RelativePermeabilityBase
{
public:

  TableRelativePermeability( std::string const & name, dataRepository::Group * const parent );

  virtual ~TableRelativePermeability() override;

  std::unique_ptr< ConstitutiveBase > deliverClone( string const & name, Group * const parent ) const override;

  static std::string catalogName() { return "TableRelativePermeability"; }

  virtual string getCatalogName() const override { return catalogName(); }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = TableRelativePermeabilityUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper();

  struct viewKeyStruct : RelativePermeabilityBase::viewKeyStruct
  {
    static constexpr char const * phaseMinVolumeFractionString() { return "phaseMinVolumeFraction"; }
    static constexpr char const * waterOilRelPermTableNamesString() { return "waterOilRelPermTableNames"; }
    static constexpr char const * gasOilRelPermTableNamesString() { return "gasOilRelPermTableNames"; }

    using ViewKey = dataRepository::ViewKey;
    ViewKey phaseMinVolumeFraction = { phaseMinVolumeFractionString() };
    ViewKey waterOilRelPermTableNames = { waterOilRelPermTableNamesString() };
    ViewKey gasOilRelPermTableNames = { gasOilRelPermTableNamesString() };

  } vieKeysTableRelativePermeability;

private:

  virtual void postProcessInput() override;

  virtual void initializePreSubGroups() override;

  /**
   * @brief Create all the table kernel wrappers needed for the simulation (for all the phases present)
   */
  void createAllTableKernelWrappers();

  /**
   * @brief Validate the relative permeability table provided in input (increasing phase vol frac and rel perm, etc)
   * @param[in] relPermTable the relative permeability table (kr vs s) for a given phase)
   */
  real64 validateRelativePermeabilityTable( TableFunction const & relPermTable ) const;

  /// Relative permeability table names (one for each phase in the oil-water pair)
  array1d< string > m_waterOilRelPermTableNames;

  /// Relative permeability table names (one for each phase in the oil-gas pair)
  array1d< string > m_gasOilRelPermTableNames;

  /// Relative permeability kernel wrapper for two-phase oil water data
  array1d< TableFunction::KernelWrapper > m_waterOilRelPermTableKernelWrappers;

  /// Relative permeability kernel wrapper for two-phase oil gas data
  array1d< TableFunction::KernelWrapper > m_gasOilRelPermTableKernelWrappers;

  /// Min phase volume fractions (deduced from the tables). With Baker, only the water phase entry is used
  array1d< real64 > m_phaseMinVolumeFraction;

};

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
TableRelativePermeabilityUpdate::
  compute( arraySlice1d< real64 const > const & phaseVolFraction,
           arraySlice1d< real64 > const & phaseRelPerm,
           arraySlice2d< real64 > const & dPhaseRelPerm_dPhaseVolFrac ) const
{
  for( real64 & val : dPhaseRelPerm_dPhaseVolFrac )
  {
    val = 0.0;
  }

  using PT = RelativePermeabilityBase::PhaseType;
  integer const ipWater = m_phaseOrder[PT::WATER];
  integer const ipOil   = m_phaseOrder[PT::OIL];
  integer const ipGas   = m_phaseOrder[PT::GAS];

  real64 oilRelPerm_wo = 0; // oil rel perm using two-phase gas-oil data
  real64 dOilRelPerm_wo_dOilVolFrac = 0; // derivative w.r.t to So
  real64 oilRelPerm_go = 0; // oil rel perm using two-phase gas-oil data
  real64 dOilRelPerm_go_dOilVolFrac = 0; // derivative w.r.t to So

  // this function assumes that the oil phase can always be present (i.e., ipOil > 0)

  // 1) Water and oil phase relative permeabilities using water-oil data
  if( ipWater >= 0 )
  {
    using WOPT = RelativePermeabilityBase::WaterOilPairPhaseType;

    // water rel perm
    m_waterOilRelPermTableKernelWrappers[WOPT::WATER].compute( &(phaseVolFraction)[ipWater],
                                                               phaseRelPerm[ipWater],
                                                               &(dPhaseRelPerm_dPhaseVolFrac)[ipWater][ipWater] );

    // oil rel perm
    m_waterOilRelPermTableKernelWrappers[WOPT::OIL].compute( &(phaseVolFraction)[ipOil],
                                                             oilRelPerm_wo,
                                                             &dOilRelPerm_wo_dOilVolFrac );
  }


  // 2) Gas and oil phase relative permeabilities using gas-oil data
  if( ipGas >= 0 )
  {
    using GOPT = RelativePermeabilityBase::GasOilPairPhaseType;

    // gas rel perm
    m_gasOilRelPermTableKernelWrappers[GOPT::GAS].compute( &(phaseVolFraction)[ipGas],
                                                           phaseRelPerm[ipGas],
                                                           &(dPhaseRelPerm_dPhaseVolFrac)[ipGas][ipGas] );

    // oil rel perm
    m_gasOilRelPermTableKernelWrappers[GOPT::OIL].compute( &(phaseVolFraction)[ipOil],
                                                           oilRelPerm_go,
                                                           &dOilRelPerm_go_dOilVolFrac );
  }


  // 3) Compute the "three-phase" oil relperm

  // if no gas, use water-oil data
  if( ipGas < 0 )
  {
    phaseRelPerm[ipOil] = oilRelPerm_wo;
    dPhaseRelPerm_dPhaseVolFrac[ipOil][ipOil] = dOilRelPerm_wo_dOilVolFrac;
  }
  // if no water, use gas-oil data
  else if( ipWater < 0 )
  {
    phaseRelPerm[ipOil] = oilRelPerm_go;
    dPhaseRelPerm_dPhaseVolFrac[ipOil][ipOil] = dOilRelPerm_go_dOilVolFrac;
  }
  // if water and oil and gas can be present, use saturation-weighted interpolation
  else
  {
    real64 const shiftedWaterVolFrac = (phaseVolFraction[ipWater] - m_phaseMinVolumeFraction[ipWater]);

    // TODO: add template to choose the interpolator from the XML file
    relpermInterpolators::Baker::compute( shiftedWaterVolFrac,
                                          phaseVolFraction[ipGas],
                                          m_phaseOrder,
                                          oilRelPerm_wo,
                                          dOilRelPerm_wo_dOilVolFrac,
                                          oilRelPerm_go,
                                          dOilRelPerm_go_dOilVolFrac,
                                          phaseRelPerm[ipOil],
                                          dPhaseRelPerm_dPhaseVolFrac[ipOil] );
  }
}

} // namespace constitutive

} // namespace geosx

#endif //GEOSX_CONSTITUTIVE_TABLERELATIVEPERMEABILITY_HPP
