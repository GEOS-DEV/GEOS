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
 * @file CompositionalMultiphaseFluid.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_COMPOSITIONALMULTIPHASEFLUID_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_COMPOSITIONALMULTIPHASEFLUID_HPP_

#include "constitutive/fluid/multifluid/compositional/CompositionalMultiphaseFluidUpdates.hpp"
#include "constitutive/fluid/multifluid/compositional/CompositionalMultiphaseFluidParameters.hpp"
#include "constitutive/fluid/multifluid/compositional/models/ConstantViscosity.hpp"
#include "constitutive/fluid/multifluid/compositional/models/CubicEOSDensity.hpp"
#include "constitutive/fluid/multifluid/compositional/models/NegativeTwoPhaseFlashModel.hpp"
#include "constitutive/fluid/multifluid/compositional/models/NullModel.hpp"
#include "constitutive/fluid/multifluid/compositional/models/PhaseModel.hpp"

#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "constitutive/fluid/multifluid/MultiFluidUtils.hpp"

namespace geos
{
namespace constitutive
{

/**
 * @brief A general compositional fluid model.
 * @tparam FLASH Class describing the phase equilibrium model
 * @tparam PHASE1 Class describing the phase property models for the first phase.
 * @tparam PHASE2 Class describing the phase property models for the second phase.
 * @tparam PHASE3 Class describing the phase property models for the possible third phase.
 */
template< typename FLASH, typename PHASE1, typename PHASE2, typename PHASE3 = compositional::NullPhaseModel >
class CompositionalMultiphaseFluid : public MultiFluidBase
{
public:
  using FlashModel = FLASH;
  using Phase1Model = PHASE1;
  using Phase2Model = PHASE2;
  using Phase3Model = PHASE3;

  // Get the number of phases
  static constexpr integer NUM_PHASES = FlashModel::getNumberOfPhases();
  // Currently restrict to 2 or 3 phases
  static_assert( NUM_PHASES == 2 || NUM_PHASES == 3 );

  using exec_policy = parallelDevicePolicy<>;

public:
  CompositionalMultiphaseFluid( string const & name, Group * const parent );

  virtual std::unique_ptr< ConstitutiveBase >
  deliverClone( string const & name,
                Group * const parent ) const override;

  static string catalogName() { return "CompositionalMultiphaseFluidGeos"; }

  virtual string getCatalogName() const override { return catalogName(); }

  virtual integer getWaterPhaseIndex() const override final;

  struct viewKeyStruct : MultiFluidBase::viewKeyStruct
  {
    static constexpr char const * equationsOfStateString() { return "equationsOfState"; }
    static constexpr char const * componentCriticalPressureString() { return "componentCriticalPressure"; }
    static constexpr char const * componentCriticalTemperatureString() { return "componentCriticalTemperature"; }
    static constexpr char const * componentAcentricFactorString() { return "componentAcentricFactor"; }
    static constexpr char const * componentVolumeShiftString() { return "componentVolumeShift"; }
    static constexpr char const * componentBinaryCoeffString() { return "componentBinaryCoeff"; }
  };

public:
  using KernelWrapper = CompositionalMultiphaseFluidUpdates< FLASH, PHASE1, PHASE2, PHASE3 >;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper();

protected:

  virtual void postProcessInput() override;

  virtual void initializePostSubGroups() override;

private:
  PhaseType getPhaseType( string const & name ) const;

  // Create the fluid models
  void createModels();

  // Flash model
  std::unique_ptr< FLASH > m_flash{};

  // Phase models
  std::unique_ptr< PHASE1 > m_phase1;
  std::unique_ptr< PHASE2 > m_phase2;
  std::unique_ptr< PHASE3 > m_phase3;

  std::unique_ptr< compositional::ComponentProperties > m_componentProperties{};

  // names of equations of state to use for each phase
  string_array m_equationsOfState;

  // standard EOS component input
  array1d< real64 > m_componentCriticalPressure;
  array1d< real64 > m_componentCriticalTemperature;
  array1d< real64 > m_componentAcentricFactor;
  array1d< real64 > m_componentVolumeShift;
  array2d< real64 > m_componentBinaryCoeff;
};

using CompositionalTwoPhaseFluidPengRobinson =  CompositionalMultiphaseFluid<
  compositional::NegativeTwoPhaseFlashPRPR,
  compositional::PhaseModel< compositional::CubicEOSDensityPR, compositional::ConstantViscosity, compositional::NullModel >,
  compositional::PhaseModel< compositional::CubicEOSDensityPR, compositional::ConstantViscosity, compositional::NullModel > >;

} /* namespace constitutive */

} /* namespace geos */

#endif //GEOS_CONSTITUTIVE_FLUID_COMPOSITIONALMULTIPHASEFLUID_HPP_
