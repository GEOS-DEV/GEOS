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
 * @file ReactiveCompressibleSinglePhaseFluid.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_REACTIVECOMPRESSIBLESINGLEPHASEFLUID_HPP_
#define GEOS_CONSTITUTIVE_FLUID_REACTIVECOMPRESSIBLESINGLEPHASEFLUID_HPP_

#include "constitutive/fluid/singlefluid/reactive/ReactiveSingleFluid.hpp"
#include "constitutive/ExponentialRelation.hpp"

#include <cmath>

namespace geos
{

namespace constitutive
{

/**
 * @brief Update class for the model suitable for lambda capture.
 * @tparam DENS_EAT type of density exponent approximation
 * @tparam VISC_EAT type of viscosity exponent approximation
 */
template< ExponentApproximationType DENS_EAT, ExponentApproximationType VISC_EAT >
class ReactiveCompressibleSinglePhaseUpdate : public ReactiveSingleFluidUpdate
{
public:

  using DensRelationType  = ExponentialRelation< real64, DENS_EAT >;
  using ViscRelationType  = ExponentialRelation< real64, VISC_EAT >;
  using DerivOffset = constitutive::singlefluid::DerivativeOffsetC< 1 >;

  ReactiveCompressibleSinglePhaseUpdate( DensRelationType const & densRelation,
                                         ViscRelationType const & viscRelation,
                                         arrayView2d< real64, constitutive::singlefluid::USD_FLUID > const & density,
                                         arrayView3d< real64, constitutive::singlefluid::USD_FLUID_DER > const & dDensity,
                                         arrayView2d< real64, constitutive::singlefluid::USD_FLUID > const & viscosity,
                                         arrayView3d< real64, constitutive::singlefluid::USD_FLUID_DER > const & dViscosity,
                                         integer const numPrimarySpecies,
                                         //  chemicalReactions::EquilibriumReactions const & equilibriumReactions,
                                         //  chemicalReactions::KineticReactions const & kineticReactions,
                                         arrayView2d< real64, compflow::USD_COMP > const & primarySpeciesConcentration,
                                         arrayView2d< real64, compflow::USD_COMP > const & secondarySpeciesConcentration,
                                         arrayView2d< real64, compflow::USD_COMP > const & primarySpeciesAggregateConcentration,
                                         arrayView3d< real64, compflow::USD_COMP_DC > const & dPrimarySpeciesAggregateConcentration_dLogPrimaryConc,
                                         arrayView2d< real64, compflow::USD_COMP > const & kineticReactionRates )
    : ReactiveSingleFluidUpdate( density, dDensity, viscosity, dViscosity, numPrimarySpecies,
                                 //  equilibriumReactions, kineticReactions,
                                 primarySpeciesConcentration, secondarySpeciesConcentration, primarySpeciesAggregateConcentration,
                                 dPrimarySpeciesAggregateConcentration_dLogPrimaryConc, kineticReactionRates ),
    m_densRelation( densRelation ),
    m_viscRelation( viscRelation )
  {}

  /// Default copy constructor
  ReactiveCompressibleSinglePhaseUpdate( ReactiveCompressibleSinglePhaseUpdate const & ) = default;

  /// Default move constructor
  ReactiveCompressibleSinglePhaseUpdate( ReactiveCompressibleSinglePhaseUpdate && ) = default;

  /// Deleted copy assignment operator
  ReactiveCompressibleSinglePhaseUpdate & operator=( ReactiveCompressibleSinglePhaseUpdate const & ) = delete;

  /// Deleted move assignment operator
  ReactiveCompressibleSinglePhaseUpdate & operator=( ReactiveCompressibleSinglePhaseUpdate && ) = delete;

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  virtual void compute( real64 const pressure,
                        real64 & density,
                        real64 & dDensity_dPressure,
                        real64 & viscosity,
                        real64 & dViscosity_dPressure ) const override
  {
    m_densRelation.compute( pressure, density, dDensity_dPressure );
    m_viscRelation.compute( pressure, viscosity, dViscosity_dPressure );
  }

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  virtual void compute( real64 const pressure,
                        real64 const GEOS_UNUSED_PARAM( temperature ),
                        real64 & density,
                        real64 & dDensity_dPressure,
                        real64 & GEOS_UNUSED_PARAM( dDensity_dTemperature ),
                        real64 & viscosity,
                        real64 & dViscosity_dPressure,
                        real64 & GEOS_UNUSED_PARAM( dViscosity_dTemperature ),
                        real64 & GEOS_UNUSED_PARAM( internalEnergy ),
                        real64 & GEOS_UNUSED_PARAM( dInternalEnergy_dPressure ),
                        real64 & GEOS_UNUSED_PARAM( dInternalEnergy_dTemperature ),
                        real64 & GEOS_UNUSED_PARAM( enthalpy ),
                        real64 & GEOS_UNUSED_PARAM( dEnthalpy_dPressure ),
                        real64 & GEOS_UNUSED_PARAM( dEnthalpy_dTemperature ) ) const override
  {
    m_densRelation.compute( pressure, density, dDensity_dPressure );
    m_viscRelation.compute( pressure, viscosity, dViscosity_dPressure );
  }


  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  virtual void update( localIndex const k,
                       localIndex const q,
                       real64 const pressure ) const override
  {
    compute( pressure,
             m_density[k][q],
             m_dDensity[k][q][DerivOffset::dP],
             m_viscosity[k][q],
             m_dViscosity[k][q][DerivOffset::dP] );
  }

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  virtual void update( localIndex const k,
                       localIndex const q,
                       real64 const pressure,
                       real64 const GEOS_UNUSED_PARAM( temperature ) ) const override
  {
    compute( pressure,
             m_density[k][q],
             m_dDensity[k][q][DerivOffset::dP],
             m_viscosity[k][q],
             m_dViscosity[k][q][DerivOffset::dP] );
  }

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  virtual void update( localIndex const k,
                       localIndex const q,
                       real64 const pressure,
                       real64 const GEOS_UNUSED_PARAM( temperature ),
                       arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & GEOS_UNUSED_PARAM( composition ) ) const override
  {
    compute( pressure,
             m_density[k][q],
             m_dDensity[k][q][DerivOffset::dP],
             m_viscosity[k][q],
             m_dViscosity[k][q][DerivOffset::dP] );
  }

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  virtual void updateChemistry( localIndex const k,
                                localIndex const GEOS_UNUSED_PARAM( q ),
                                real64 const pressure,
                                real64 const temperature,
                                arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition ) const override
  {
    for( int i=0; i < m_numPrimarySpecies; i++ )
    {
      m_primarySpeciesAggregateConcentration[k][i] = composition[i];
    }

    computeChemistry( pressure,
                      temperature,
                      m_primarySpeciesAggregateConcentration[k],
                      m_primarySpeciesConcentration[k],
                      m_secondarySpeciesConcentration[k],
                      m_kineticReactionRates[k] );
  }

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  virtual void updateChemistryLogConc( localIndex const k,
                                       localIndex const GEOS_UNUSED_PARAM( q ),
                                       real64 const GEOS_UNUSED_PARAM( pressure ),
                                       real64 const GEOS_UNUSED_PARAM( temperature ),
                                       arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & logPrimaryConc ) const override
  {
    for( int i=0; i < m_numPrimarySpecies; i++ )
    {
      m_primarySpeciesConcentration[k][i] = std::exp( logPrimaryConc[i] );

      m_primarySpeciesAggregateConcentration[k][i] = m_primarySpeciesConcentration[k][i];

      m_dPrimarySpeciesAggregateConcentration_dLogPrimaryConc[k][i][i] = m_primarySpeciesConcentration[k][i];
    }

    // computeChemistry( pressure,
    //                   temperature,
    //                   m_primarySpeciesAggregateConcentration[k],
    //                   m_primarySpeciesConcentration[k],
    //                   m_secondarySpeciesConcentration[k],
    //                   m_kineticReactionRates[k] );
  }

private:

  /// Relationship between the fluid density and pressure
  DensRelationType m_densRelation;

  /// Relationship between the fluid viscosity and pressure
  ViscRelationType m_viscRelation;

};

class ReactiveCompressibleSinglePhase : public ReactiveSingleFluid
{
public:
  using DerivOffset = singlefluid::DerivativeOffset;
  ReactiveCompressibleSinglePhase( string const & name, Group * const parent );

  virtual ~ReactiveCompressibleSinglePhase() override;

  static string catalogName() { return "ReactiveCompressibleSinglePhase"; }

  virtual string getCatalogName() const override { return catalogName(); }

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  /// Type of kernel wrapper for in-kernel update (TODO: support multiple EAT, not just linear)
  using KernelWrapper = ReactiveCompressibleSinglePhaseUpdate< ExponentApproximationType::Linear, ExponentApproximationType::Linear >;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper();

  struct viewKeyStruct
  {
    static constexpr char const * defaultDensityString() { return "defaultDensity"; }
    static constexpr char const * defaultViscosityString() { return "defaultViscosity"; }
    static constexpr char const * compressibilityString() { return "compressibility"; }
    static constexpr char const * viscosibilityString() { return "viscosibility"; }
    static constexpr char const * referencePressureString() { return "referencePressure"; }
    static constexpr char const * referenceDensityString() { return "referenceDensity"; }
    static constexpr char const * referenceViscosityString() { return "referenceViscosity"; }
    static constexpr char const * densityModelTypeString() { return "densityModelType"; }
    static constexpr char const * viscosityModelTypeString() { return "viscosityModelType"; }
  };

  real64 defaultDensity() const override final { return m_defaultDensity; }
  real64 defaultViscosity() const override final { return m_defaultViscosity; }

protected:

  virtual void postInputInitialization() override;

  /// default density value
  real64 m_defaultDensity;

  /// default viscosity value
  real64 m_defaultViscosity;

  /// scalar fluid bulk modulus parameter
  real64 m_compressibility;

  /// scalar fluid viscosity exponential coefficient
  real64 m_viscosibility;

  /// reference pressure parameter
  real64 m_referencePressure;

  /// reference density parameter
  real64 m_referenceDensity;

  /// reference viscosity parameter
  real64 m_referenceViscosity;

  /// type of density model in terms of pressure
  ExponentApproximationType m_densityModelType;

  /// type of viscosity model
  ExponentApproximationType m_viscosityModelType;

};

} /* namespace constitutive */

} /* namespace geos */

#endif /* GEOS_CONSTITUTIVE_FLUID_REACTIVECOMPRESSIBLESINGLEPHASEFLUID_HPP_ */
