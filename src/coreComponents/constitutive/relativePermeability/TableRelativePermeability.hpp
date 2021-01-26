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
#include "managers/Functions/TableFunction.hpp"

namespace geosx
{
namespace constitutive
{

class TableRelativePermeabilityUpdate final : public RelativePermeabilityBaseUpdate
{
public:

  TableRelativePermeabilityUpdate( arrayView1d< TableFunctionKernelWrapper const > const & waterOilRelPermTableKernelWrappers,
                                   arrayView1d< TableFunctionKernelWrapper const > const & gasOilRelPermTableKernelWrappers,
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
  virtual void Compute( arraySlice1d< real64 const > const & phaseVolFraction,
                        arraySlice1d< real64 > const & phaseRelPerm,
                        arraySlice2d< real64 > const & dPhaseRelPerm_dPhaseVolFrac ) const override;

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  virtual void Update( localIndex const k,
                       localIndex const q,
                       arraySlice1d< real64 const > const & phaseVolFraction ) const override
  {
    Compute( phaseVolFraction,
             m_phaseRelPerm[k][q],
             m_dPhaseRelPerm_dPhaseVolFrac[k][q] );
  }

private:

  /// Kernel wrappers for the water-oil relative permeabilities
  arrayView1d< TableFunctionKernelWrapper const > m_waterOilRelPermTableKernelWrappers;

  /// Kernel wrappers for the gas-oil relative permeabilities
  arrayView1d< TableFunctionKernelWrapper const > m_gasOilRelPermTableKernelWrappers;

  /// Minimum volume fraction for each phase (deduced from the table)
  arrayView1d< real64 const > m_phaseMinVolumeFraction;

};

class TableRelativePermeability : public RelativePermeabilityBase
{
public:

  TableRelativePermeability( std::string const & name, dataRepository::Group * const parent );

  virtual ~TableRelativePermeability() override;

  std::unique_ptr< ConstitutiveBase > deliverClone( string const & name, Group * const parent ) const override;

  static std::string CatalogName() { return "TableRelativePermeability"; }

  virtual string getCatalogName() const override { return CatalogName(); }

  /**
   * @brief Create all the table kernel wrappers needed for the simulation (for all the phases present)
   */
  void CreateAllTableKernelWrappers();

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = TableRelativePermeabilityUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper();

  struct viewKeyStruct : RelativePermeabilityBase::viewKeyStruct
  {
    static constexpr auto phaseMinVolumeFractionString = "phaseMinVolumeFraction";
    static constexpr auto waterOilRelPermTableNamesString = "waterOilRelPermTableNames";
    static constexpr auto gasOilRelPermTableNamesString = "gasOilRelPermTableNames";

    using ViewKey = dataRepository::ViewKey;
    ViewKey phaseMinVolumeFraction = { phaseMinVolumeFractionString };
    ViewKey waterOilRelPermTableNames = { waterOilRelPermTableNamesString };
    ViewKey gasOilRelPermTableNames = { gasOilRelPermTableNamesString };

  } vieKeysTableRelativePermeability;

protected:

  virtual void PostProcessInput() override;

  virtual void InitializePreSubGroups( Group * const ) override;

private:

  /**
   * @brief Validate the relative permeability table provided in input (increasing phase vol frac and rel perm, etc)
   * @param[in] relPermTable the relative permeability table (kr vs s) for a given phase)
   */
  real64 ValidateRelativePermeabilityTable( TableFunction const & relPermTable ) const;

  /// Relative permeability table names (one for each phase in the oil-water pair)
  array1d< string > m_waterOilRelPermTableNames;

  /// Relative permeability table names (one for each phase in the oil-gas pair)
  array1d< string > m_gasOilRelPermTableNames;

  /// Relative permeability kernel wrapper for two-phase oil water data
  array1d< TableFunctionKernelWrapper > m_waterOilRelPermTableKernelWrappers;

  /// Relative permeability kernel wrapper for two-phase oil gas data
  array1d< TableFunctionKernelWrapper > m_gasOilRelPermTableKernelWrappers;

  /// Min phase volume fractions (deduced from the tables). With Baker, only the water phase entry is used
  array1d< real64 > m_phaseMinVolumeFraction;

};

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
TableRelativePermeabilityUpdate::
  Compute( arraySlice1d< real64 const > const & phaseVolFraction,
           arraySlice1d< real64 > const & phaseRelPerm,
           arraySlice2d< real64 > const & dPhaseRelPerm_dPhaseVolFrac ) const
{
  real64 relPermDerivative[ 1 ] = { 0.0 };

  localIndex const NP = numPhases();

  for( localIndex ip = 0; ip < NP; ++ip )
  {
    for( localIndex jp = 0; jp < NP; ++jp )
    {
      dPhaseRelPerm_dPhaseVolFrac[ip][jp] = 0.0;
    }
  }

  integer const ip_water = m_phaseOrder[RelativePermeabilityBase::PhaseType::WATER];
  integer const ip_oil   = m_phaseOrder[RelativePermeabilityBase::PhaseType::OIL];
  integer const ip_gas   = m_phaseOrder[RelativePermeabilityBase::PhaseType::GAS];

  real64 oilRelPerm_wo = 0; // oil rel perm using two-phase gas-oil data
  real64 dOilRelPerm_wo_dOilVolFrac = 0; // derivative w.r.t to So
  real64 oilRelPerm_go = 0; // oil rel perm using two-phase gas-oil data
  real64 dOilRelPerm_go_dOilVolFrac = 0; // derivative w.r.t to So

  // this function assumes that the oil phase can always be present (i.e., ip_oil > 0)

  // 1) Water and oil phase relative permeabilities using water-oil data
  if( ip_water >= 0 )
  {
    // water rel perm
    m_waterOilRelPermTableKernelWrappers[RelativePermeabilityBase::WaterOilPairPhaseType::WATER].Compute( &(phaseVolFraction)[ip_water],
                                                                                                          phaseRelPerm[ip_water],
                                                                                                          relPermDerivative );
    dPhaseRelPerm_dPhaseVolFrac[ip_water][ip_water] = relPermDerivative[0];

    // oil rel perm
    m_waterOilRelPermTableKernelWrappers[RelativePermeabilityBase::WaterOilPairPhaseType::OIL].Compute( &(phaseVolFraction)[ip_oil],
                                                                                                        phaseRelPerm[ip_oil],
                                                                                                        relPermDerivative );
    dPhaseRelPerm_dPhaseVolFrac[ip_oil][ip_oil] = relPermDerivative[0];
  }


  // 2) Gas and oil phase relative permeabilities using gas-oil data
  if( ip_gas >= 0 )
  {
    // gas rel perm
    m_gasOilRelPermTableKernelWrappers[RelativePermeabilityBase::GasOilPairPhaseType::GAS].Compute( &(phaseVolFraction)[ip_gas],
                                                                                                    phaseRelPerm[ip_gas],
                                                                                                    relPermDerivative );
    dPhaseRelPerm_dPhaseVolFrac[ip_gas][ip_gas] = relPermDerivative[0];

    // oil rel perm
    m_gasOilRelPermTableKernelWrappers[RelativePermeabilityBase::GasOilPairPhaseType::OIL].Compute( &(phaseVolFraction)[ip_oil],
                                                                                                    phaseRelPerm[ip_oil],
                                                                                                    relPermDerivative );
    dPhaseRelPerm_dPhaseVolFrac[ip_oil][ip_oil] = relPermDerivative[0];
  }


  // 3) Compute the "three-phase" oil relperm

  // if no gas, use water-oil data
  if( ip_gas < 0 )
  {
    phaseRelPerm[ip_oil] = oilRelPerm_wo;
    dPhaseRelPerm_dPhaseVolFrac[ip_oil][ip_oil] = dOilRelPerm_wo_dOilVolFrac;
  }
  // if no water, use gas-oil data
  else if( ip_water < 0 )
  {
    phaseRelPerm[ip_oil] = oilRelPerm_go;
    dPhaseRelPerm_dPhaseVolFrac[ip_oil][ip_oil] = dOilRelPerm_go_dOilVolFrac;
  }
  // if water and oil and gas can be present, use saturation-weighted interpolation
  else
  {
    real64 const shiftedWaterVolFrac = (phaseVolFraction[ip_water] - m_phaseMinVolumeFraction[ip_water]);

    // TODO: add template to choose the interpolator from the XML file
    relpermInterpolators::Baker::Compute( shiftedWaterVolFrac,
                                          phaseVolFraction[ip_gas],
                                          m_phaseOrder,
                                          oilRelPerm_wo,
                                          dOilRelPerm_wo_dOilVolFrac,
                                          oilRelPerm_go,
                                          dOilRelPerm_go_dOilVolFrac,
                                          phaseRelPerm[ip_oil],
                                          dPhaseRelPerm_dPhaseVolFrac[ip_oil] );
  }

}

} // namespace constitutive

} // namespace geosx

#endif //GEOSX_CONSTITUTIVE_TABLERELATIVEPERMEABILITY_HPP
