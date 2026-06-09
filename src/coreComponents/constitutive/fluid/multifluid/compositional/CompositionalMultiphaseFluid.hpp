/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
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
#include "constitutive/fluid/multifluid/compositional/models/ConstantViscosity.hpp"
#include "constitutive/fluid/multifluid/compositional/models/CompositionalDensity.hpp"
#include "constitutive/fluid/multifluid/compositional/models/CompositionalEnthalpy.hpp"
#include "constitutive/fluid/multifluid/compositional/models/ImmiscibleWaterDensity.hpp"
#include "constitutive/fluid/multifluid/compositional/models/ImmiscibleWaterFlashModel.hpp"
#include "constitutive/fluid/multifluid/compositional/models/ImmiscibleWaterViscosity.hpp"
#include "constitutive/fluid/multifluid/compositional/models/KValueFlashModel.hpp"
#include "constitutive/fluid/multifluid/compositional/models/LohrenzBrayClarkViscosity.hpp"
#include "constitutive/fluid/multifluid/compositional/models/NegativeTwoPhaseFlashModel.hpp"
#include "constitutive/fluid/multifluid/compositional/models/NullModel.hpp"
#include "constitutive/fluid/multifluid/compositional/models/PhaseModel.hpp"
#include "constitutive/fluid/multifluid/compositional/models/PhillipsBrineDensity.hpp"
#include "constitutive/fluid/multifluid/compositional/models/PhillipsBrineViscosity.hpp"
#include "constitutive/fluid/multifluid/compositional/parameters/ModelParameters.hpp"

namespace geos
{
namespace constitutive
{

/**
 * @brief A general compositional fluid model.
 * @tparam FLASH Class describing the phase equilibrium model
 * @tparam PHASES Classes describing the phase property models for each phase.
 */
template< typename FLASH, typename ... PHASES >
class CompositionalMultiphaseFluid : public MultiFluidBase
{
public:
  using FlashModel = FLASH;

  // Get the number of phases
  static constexpr integer NUM_PHASES = static_cast< integer >(sizeof...(PHASES));

  // Ensure that the number of phases matches the flash object
  static_assert( NUM_PHASES == FlashModel::KernelWrapper::getNumberOfPhases(),
                 "Number of phases should match the flash" );

  // Check if all phase models have a valid enthalpy model
  static constexpr bool isThermalType()
  {
    return (!std::is_same_v< typename PHASES::Enthalpy, compositional::NullModel > && ...);
  }

  // Check if all phase models do not provide an enthalpy model
  static constexpr bool isIsoThermalType()
  {
    return (std::is_same_v< typename PHASES::Enthalpy, compositional::NullModel > && ...);
  }

  // Either all phases are thermal or all are not
  static_assert( ( isThermalType() || isIsoThermalType() ),
                 "All phase models must either use NullModel for Enthalpy, or none should." );

public:
  CompositionalMultiphaseFluid( string const & name, dataRepository::Group * const parent );

  virtual std::unique_ptr< ConstitutiveBase >
  deliverClone( string const & name,
                dataRepository::Group * const parent ) const override;

  static string catalogName();

  virtual string getCatalogName() const override { return catalogName(); }

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numPts ) override;

  // TODO: This method should be implemented if an incorrect extrapolation of the pressure and temperature is encountered in the kernel
  /**
   * @copydoc MultiFluidBase::checkTablesParameters( real64 pressure, real64 temperature )
   */
  virtual void checkTablesParameters( real64 pressure, real64 temperature ) const override final
  {
    GEOS_UNUSED_VAR( pressure, temperature );
  }

  virtual void initializeState() const override;

  virtual integer getWaterPhaseIndex() const override final;

  struct viewKeyStruct : MultiFluidBase::viewKeyStruct
  {
    static constexpr char const * componentCriticalPressureString() { return "componentCriticalPressure"; }
    static constexpr char const * componentCriticalTemperatureString() { return "componentCriticalTemperature"; }
    static constexpr char const * componentAcentricFactorString() { return "componentAcentricFactor"; }
    static constexpr char const * componentVolumeShiftString() { return "componentVolumeShift"; }
    static constexpr char const * componentBinaryCoeffString() { return "componentBinaryCoeff"; }
  };

public:
  using KernelWrapper = CompositionalMultiphaseFluidUpdates< FLASH, PHASES... >;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper();

protected:

  virtual void postInputInitialization() override;

  virtual void initializePostSubGroups() override;

private:
  // Create the fluid models
  void createModels();

  // Create each of the phase models
  template< std::size_t... Is >
  void createPhaseModels( std::index_sequence< Is... > );

  // Helper to create the kernel wrappers
  template< std::size_t... Is >
  KernelWrapper createKernelWrapper( std::index_sequence< Is... > );

  array1d< integer > getPhaseTypes() const;

  static std::unique_ptr< compositional::ModelParameters > createModelParameters();

  // Flash model
  std::unique_ptr< FLASH > m_flash{};

  // Phase ordering
  array1d< integer > m_phaseOrder{};
  array1d< integer > m_phaseType{};

  // Phase models
  camp::tuple< std::unique_ptr< PHASES >... > m_phases{};

  // Standard EOS component input
  std::unique_ptr< compositional::ComponentProperties > m_componentProperties{};

  // Extra parameters specific to this model
  std::unique_ptr< compositional::ModelParameters > m_parameters{};

  // backup data
  PhaseComp::ValueType m_kValues;
};

using CompositionalTwoPhaseConstantViscosity = CompositionalMultiphaseFluid<
  compositional::NegativeTwoPhaseFlashModel,
  compositional::PhaseModel< compositional::CompositionalDensity, compositional::ConstantViscosity >,
  compositional::PhaseModel< compositional::CompositionalDensity, compositional::ConstantViscosity > >;
using CompositionalTwoPhaseLohrenzBrayClarkViscosity = CompositionalMultiphaseFluid<
  compositional::NegativeTwoPhaseFlashModel,
  compositional::PhaseModel< compositional::CompositionalDensity, compositional::LohrenzBrayClarkViscosity >,
  compositional::PhaseModel< compositional::CompositionalDensity, compositional::LohrenzBrayClarkViscosity > >;
using CompositionalTwoPhasePhillipsBrine = CompositionalMultiphaseFluid<
  compositional::NegativeTwoPhaseFlashModel,
  compositional::PhaseModel< compositional::PhillipsBrineDensity, compositional::PhillipsBrineViscosity >,
  compositional::PhaseModel< compositional::CompositionalDensity, compositional::LohrenzBrayClarkViscosity > >;
using CompositionalThreePhaseLohrenzBrayClarkViscosity = CompositionalMultiphaseFluid<
  compositional::ImmiscibleWaterFlashModel,
  compositional::PhaseModel< compositional::CompositionalDensity, compositional::LohrenzBrayClarkViscosity >,
  compositional::PhaseModel< compositional::CompositionalDensity, compositional::LohrenzBrayClarkViscosity >,
  compositional::PhaseModel< compositional::ImmiscibleWaterDensity, compositional::ImmiscibleWaterViscosity > >;
using CompositionalKValueLohrenzBrayClarkViscosity = CompositionalMultiphaseFluid<
  compositional::KValueFlashModel< 2 >,
  compositional::PhaseModel< compositional::CompositionalDensity, compositional::LohrenzBrayClarkViscosity >,
  compositional::PhaseModel< compositional::CompositionalDensity, compositional::LohrenzBrayClarkViscosity > >;
using CompositionalKValuePhillipsBrine = CompositionalMultiphaseFluid<
  compositional::KValueFlashModel< 2 >,
  compositional::PhaseModel< compositional::PhillipsBrineDensity, compositional::PhillipsBrineViscosity >,
  compositional::PhaseModel< compositional::CompositionalDensity, compositional::LohrenzBrayClarkViscosity > >;
using CompositionalThermalTwoPhaseLohrenzBrayClarkViscosity = CompositionalMultiphaseFluid<
  compositional::NegativeTwoPhaseFlashModel,
  compositional::PhaseModel< compositional::CompositionalDensity, compositional::LohrenzBrayClarkViscosity, compositional::CompositionalEnthalpy >,
  compositional::PhaseModel< compositional::CompositionalDensity, compositional::LohrenzBrayClarkViscosity, compositional::CompositionalEnthalpy > >;
using CompositionalThermalTwoPhasePhillipsBrine = CompositionalMultiphaseFluid<
  compositional::NegativeTwoPhaseFlashModel,
  compositional::PhaseModel< compositional::PhillipsBrineDensity, compositional::PhillipsBrineViscosity, compositional::CompositionalEnthalpy >,
  compositional::PhaseModel< compositional::CompositionalDensity, compositional::LohrenzBrayClarkViscosity, compositional::CompositionalEnthalpy > >;
using CompositionalThermalKValuePhillipsBrine = CompositionalMultiphaseFluid<
  compositional::KValueFlashModel< 2 >,
  compositional::PhaseModel< compositional::PhillipsBrineDensity, compositional::PhillipsBrineViscosity, compositional::CompositionalEnthalpy >,
  compositional::PhaseModel< compositional::CompositionalDensity, compositional::LohrenzBrayClarkViscosity, compositional::CompositionalEnthalpy > >;

} /* namespace constitutive */

} /* namespace geos */

#endif //GEOS_CONSTITUTIVE_FLUID_COMPOSITIONALMULTIPHASEFLUID_HPP_
