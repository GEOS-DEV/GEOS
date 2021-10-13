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
                                  arrayView3d< real64, relperm::USD_RELPERM > const & phaseRelPerm,
                                  arrayView4d< real64, relperm::USD_RELPERM_DS > const & dPhaseRelPerm_dPhaseVolFrac )
    : m_phaseTypes( phaseTypes ),
    m_phaseOrder( phaseOrder ),
    m_phaseRelPerm( phaseRelPerm ),
    m_dPhaseRelPerm_dPhaseVolFrac( dPhaseRelPerm_dPhaseVolFrac )
  {}

  arrayView1d< integer const > m_phaseTypes;
  arrayView1d< integer const > m_phaseOrder;

  arrayView3d< real64, relperm::USD_RELPERM > m_phaseRelPerm;
  arrayView4d< real64, relperm::USD_RELPERM_DS > m_dPhaseRelPerm_dPhaseVolFrac;

private:

  GEOSX_HOST_DEVICE
  virtual void compute( arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFraction,
                        arraySlice1d< real64, relperm::USD_RELPERM - 2 > const & phaseRelPerm,
                        arraySlice2d< real64, relperm::USD_RELPERM_DS - 2 > const & dPhaseRelPerm_dPhaseVolFrac ) const = 0;

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

  arrayView3d< real64 const, relperm::USD_RELPERM > phaseRelPerm() const { return m_phaseRelPerm; }
  arrayView4d< real64 const, relperm::USD_RELPERM_DS > dPhaseRelPerm_dPhaseVolFraction() const { return m_dPhaseRelPerm_dPhaseVolFrac; }

  struct viewKeyStruct : ConstitutiveBase::viewKeyStruct
  {
    static constexpr char const * phaseNamesString() { return "phaseNames"; }
    static constexpr char const * phaseTypesString() { return "phaseTypes"; }
    static constexpr char const * phaseOrderString() { return "phaseOrder"; }

    static constexpr char const * phaseRelPermString() { return "phaseRelPerm"; }                                       // Kr
    static constexpr char const * dPhaseRelPerm_dPhaseVolFractionString() { return "dPhaseRelPerm_dPhaseVolFraction"; } // dKr_p/dS_p
  };

private:

  /**
   * @brief Function called internally to resize member arrays
   * @param size primary dimension (e.g. number of cells)
   * @param numPts secondary dimension (e.g. number of gauss points per cell)
   */
  void resizeFields( localIndex const size, localIndex const numPts );

  /**
   * @brief Called internally to set array dim labels.
   */
  void setLabels();

protected:

  virtual void postProcessInput() override;

  // phase names read from input
  string_array m_phaseNames;

  // phase ordering info
  array1d< integer > m_phaseTypes;
  array1d< integer > m_phaseOrder;

  // output quantities
  array3d< real64, relperm::LAYOUT_RELPERM >  m_phaseRelPerm;
  array4d< real64, relperm::LAYOUT_RELPERM_DS >  m_dPhaseRelPerm_dPhaseVolFrac;
};

} // namespace constitutive

} // namespace geosx


#endif //GEOSX_CONSTITUTIVE_RELATIVEPERMEABILITYBASE_HPP
