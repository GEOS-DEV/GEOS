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
 * @file RelativePermeabilityBase.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_RELATIVEPERMEABILITYBASE_HPP
#define GEOSX_CONSTITUTIVE_RELATIVEPERMEABILITYBASE_HPP

#include "common/DataLayouts.hpp"
#include "constitutive/ConstitutiveBase.hpp"
#include "constitutive/relativePermeability/layouts.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

namespace geosx
{

namespace constitutive
{

class RelativePermeabilityBaseUpdate
{
public:

  /**
   * @brief Get number of elements in this wrapper.
   * @return number of elements
   */
  GEOSX_HOST_DEVICE
  localIndex numElems() const { return m_phaseRelPerm.size( 0 ); }

  /**
   * @brief Get number of gauss points per element.
   * @return number of gauss points per element
   */
  GEOSX_HOST_DEVICE
  localIndex numGauss() const { return m_phaseRelPerm.size( 1 ); }

  /**
   * @brief Get number of fluid phases.
   * @return number of phases
   */
  GEOSX_HOST_DEVICE
  integer numPhases() const { return LvArray::integerConversion< integer >( m_phaseTypes.size() ); }

protected:

  RelativePermeabilityBaseUpdate( arrayView1d< integer const > const & phaseTypes,
                                  arrayView1d< integer const > const & phaseOrder,
                                  arrayView1d< real64 const > const & phaseMinVolFrac,
                                  arrayView3d< real64, relperm::USD_RELPERM > const & phaseRelPerm,
                                  arrayView4d< real64, relperm::USD_RELPERM_DS > const & dPhaseRelPerm_dPhaseVolFrac,
                                  arrayView3d< real64, relperm::USD_RELPERM > const & phaseTrappedVolFrac )
    : m_phaseTypes( phaseTypes ),
    m_phaseOrder( phaseOrder ),
    m_phaseMinVolumeFraction( phaseMinVolFrac ),
    m_phaseRelPerm( phaseRelPerm ),
    m_dPhaseRelPerm_dPhaseVolFrac( dPhaseRelPerm_dPhaseVolFrac ),
    m_phaseTrappedVolFrac( phaseTrappedVolFrac ) {}

  arrayView1d< integer const > m_phaseTypes;
  arrayView1d< integer const > m_phaseOrder;

  arrayView1d< real64 const > m_phaseMinVolumeFraction;

  arrayView3d< real64, relperm::USD_RELPERM > m_phaseRelPerm;
  arrayView4d< real64, relperm::USD_RELPERM_DS > m_dPhaseRelPerm_dPhaseVolFrac;

  arrayView3d< real64, relperm::USD_RELPERM > m_phaseTrappedVolFrac;


private:

  GEOSX_HOST_DEVICE
  virtual void update( localIndex const k,
                       localIndex const q,
                       arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFraction ) const = 0;
};

class RelativePermeabilityBase : public ConstitutiveBase
{
public:

  static constexpr integer MAX_NUM_PHASES = 3;

  struct PhaseType
  {
    enum : integer
    {
      OIL = 0,
      GAS = 1,
      WATER = 2,
    };
  };

  /// order of the phase properties in the water-oil data
  struct WaterOilPairPhaseType
  {
    enum : integer
    {
      WATER = 0, ///< first water phase property
      OIL = 1,   ///< second oil phase property
    };
  };

  /// order of the phase properties in the gas-oil data
  struct GasOilPairPhaseType
  {
    enum : integer
    {
      GAS = 0, ///< first gas phase property
      OIL = 1, ///< second oil phase property
    };
  };

  RelativePermeabilityBase( string const & name, dataRepository::Group * const parent );

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  integer numFluidPhases() const { return LvArray::integerConversion< integer >( m_phaseNames.size() ); }

  arrayView1d< string const > phaseNames() const { return m_phaseNames; }

  arrayView1d< real64 const > phaseMinVolumeFraction() const { return m_phaseMinVolumeFraction; }
  arrayView3d< real64 const, relperm::USD_RELPERM > phaseTrapped() const { return m_phaseTrappedVolFrac; }

  arrayView3d< real64 const, relperm::USD_RELPERM > phaseRelPerm() const { return m_phaseRelPerm; }
  arrayView4d< real64 const, relperm::USD_RELPERM_DS > dPhaseRelPerm_dPhaseVolFraction() const { return m_dPhaseRelPerm_dPhaseVolFrac; }

  /**
   * @brief Save converged phase volume fraction at the end of a time step (needed for hysteresis)
   * @param[in] phaseVolFraction an array containing the phase volume fractions at the end of a converged time step
   */
  virtual void saveConvergedPhaseVolFractionState( arrayView2d< real64 const, compflow::USD_PHASE > const & phaseVolFraction ) const
  { GEOSX_UNUSED_VAR( phaseVolFraction ); }

  /**
   * @brief init and save trapped volfrac at the end of time step in each element
   * @param[in] phaseVolFraction an array containing the phase volume fractions at the end of a converged time step
   */
  virtual void updateTrappedPhaseVolFraction( arrayView2d< real64 const, compflow::USD_PHASE > const & phaseVolFraction ) const;

  struct viewKeyStruct : ConstitutiveBase::viewKeyStruct
  {
    static constexpr char const * phaseNamesString() { return "phaseNames"; }
    static constexpr char const * phaseTypesString() { return "phaseTypes"; }
    static constexpr char const * phaseOrderString() { return "phaseOrder"; }
    static constexpr char const * phaseMinVolumeFractionString() { return "phaseMinVolumeFraction"; }
  };

private:

  /**
   * @brief Called internally to set array dim labels.
   */
  void setLabels();

protected:

  /**
   * @brief Function called internally to resize member arrays
   * @param size primary dimension (e.g. number of cells)
   * @param numPts secondary dimension (e.g. number of gauss points per cell)
   */
  virtual void resizeFields( localIndex const size, localIndex const numPts );

  virtual void postProcessInput() override;

  // phase names read from input
  string_array m_phaseNames;

  // phase ordering info
  array1d< integer > m_phaseTypes;
  array1d< integer > m_phaseOrder;
  array1d< real64 > m_phaseMinVolumeFraction;

  // output quantities
  array3d< real64, relperm::LAYOUT_RELPERM > m_phaseRelPerm;
  array4d< real64, relperm::LAYOUT_RELPERM_DS > m_dPhaseRelPerm_dPhaseVolFrac;
  array3d< real64, relperm::LAYOUT_RELPERM > m_phaseTrappedVolFrac;

};

} // namespace constitutive

} // namespace geosx


#endif //GEOSX_CONSTITUTIVE_RELATIVEPERMEABILITYBASE_HPP
