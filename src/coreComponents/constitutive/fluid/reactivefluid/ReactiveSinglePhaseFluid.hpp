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
 * @file ReactiveSinglePhaseFluid.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_REACTIVEFLUID_REACTIVESINGLEPHASEFLUID_HPP_
#define GEOS_CONSTITUTIVE_FLUID_REACTIVEFLUID_REACTIVESINGLEPHASEFLUID_HPP_

#include "common/format/EnumStrings.hpp"

#include "constitutive/ConstitutiveBase.hpp"
#include "constitutive/fluid/reactivefluid/ReactiveFluidLayouts.hpp"
#include "constitutive/fluid/singlefluid/CompressibleSinglePhaseFluid.hpp"
#include "constitutive/fluid/singlefluid/ThermalCompressibleSinglePhaseFluid.hpp"

#include "constitutive/fluid/reactivefluid/ReactiveFluidSystemSelector.hpp"
#include "constitutive/HPCReact/src/reactions/reactionsSystems/EquilibriumReactions.hpp"
#include "constitutive/HPCReact/src/reactions/reactionsSystems/MixedEquilibriumKineticReactions.hpp"
#include <memory>
#include <optional>

namespace geos
{

namespace constitutive
{

namespace reactivefluid
{

using namespace hpcReact::reactionsSystems;

template< typename BASE >
class ReactiveSinglePhaseFluid : public BASE
{
public:

  ReactiveSinglePhaseFluid( string const & name,
                            dataRepository::Group * const parent );

  virtual std::unique_ptr< ConstitutiveBase >
  deliverClone( string const & name,
                dataRepository::Group * const parent ) const override;

  static string catalogName() { return string( "Reactive" ) + BASE::catalogName(); }
  virtual string getCatalogName() const override { return catalogName(); }

  virtual void saveConvergedState() const override;

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  static constexpr integer MAX_NUM_SPECIES = 20;
  static constexpr integer MAX_NUM_KINETIC_REACTIONS = 10;

  arrayView3d< real64 const, reactivefluid::USD_SPECIES > primarySpeciesAggregateConcentration() const
  { return m_primarySpeciesAggregateConcentration; }

  arrayView3d< real64 const, reactivefluid::USD_SPECIES > primarySpeciesAggregateConcentration_n() const
  { return m_primarySpeciesAggregateConcentration_n; }

  arrayView3d< real64 const, reactivefluid::USD_SPECIES > primarySpeciesMobileAggregateConcentration() const
  { return m_primarySpeciesMobileAggregateConcentration; }

  arrayView4d< real64 const, reactivefluid::USD_SPECIES_DC > dPrimarySpeciesAggregateConcentration_dLogPrimarySpeciesConcentrations() const
  { return m_dPrimarySpeciesAggregateConcentration_dLogPrimarySpeciesConcentrations; }

  arrayView4d< real64 const, reactivefluid::USD_SPECIES_DC > dPrimarySpeciesMobileAggregateConcentration_dLogPrimarySpeciesConcentrations() const
  { return m_dPrimarySpeciesMobileAggregateConcentration_dLogPrimarySpeciesConcentrations; }

  arrayView3d< real64 const, reactivefluid::USD_SPECIES > secondarySpeciesConcentration() const
  { return m_secondarySpeciesConcentration; }

  arrayView3d< real64 const, reactivefluid::USD_SPECIES > aggregateSpeciesRates() const
  { return m_aggregateSpeciesRates; }

  arrayView3d< real64 const, reactivefluid::USD_SPECIES > kineticReactionRates() const
  { return m_kineticReactionRates; }

  arrayView4d< real64 const, reactivefluid::USD_SPECIES_DC > dAggregateSpeciesRates_dLogPrimarySpeciesConcentrations() const
  { return m_dAggregateSpeciesRates_dLogPrimarySpeciesConcentrations; }

  integer numPrimarySpecies() const { return m_numPrimarySpecies; }

  integer numSecondarySpecies() const { return m_numSecondarySpecies; }

  integer numKineticReactions() const { return m_numKineticReactions; }

  /**
   * @brief Mass of solvent per unit volume of solution [kg/m^3].
   *
   * Converts species molality [mol/kg solvent] to molarity [mol/m^3 solution]. HPCReact is a
   * molality-based library: concentrations, equilibrium constants and mass-action quotients are all
   * on the molal scale.
   */
  real64 solventMassPerSolutionVolume() const { return m_solventMassPerSolutionVolume; }

  /**
   * @brief Kernel wrapper class for ReactiveSinglePhaseFluid.
   */
  template< typename REACTION_PARAMS_TYPE, typename ACTIVITY_MODEL >
  class ReactionKernelWrapper
  {

public:

    ReactionKernelWrapper( arrayView3d< real64, reactivefluid::USD_SPECIES > const & primarySpeciesAggregateConcentration,
                           arrayView3d< real64, reactivefluid::USD_SPECIES > const & primarySpeciesMobileAggregateConcentration,
                           arrayView4d< real64, reactivefluid::USD_SPECIES_DC > const & dPrimarySpeciesAggregateConcentration_dLogPrimarySpeciesConcentrations,
                           arrayView4d< real64, reactivefluid::USD_SPECIES_DC > const & dPrimarySpeciesMobileAggregateConcentration_dLogPrimarySpeciesConcentrations,
                           arrayView3d< real64 const, reactivefluid::USD_SPECIES > const & initialPrimarySpeciesConcentration,
                           arrayView3d< real64, reactivefluid::USD_SPECIES > const & secondarySpeciesConcentration,
                           arrayView3d< real64, reactivefluid::USD_SPECIES > const & kineticReactionRates,
                           arrayView3d< real64, reactivefluid::USD_SPECIES > const & aggregateSpeciesRates,
                           arrayView4d< real64, reactivefluid::USD_SPECIES_DC > const & dAggregateSpeciesRates_dLogPrimarySpeciesConcentrations,
                           integer const numPrimarySpecies,
                           integer const numSecondarySpecies,
                           integer const numKineticReactions,
                           REACTION_PARAMS_TYPE params,
                           typename ACTIVITY_MODEL::Params activityParams ):
      m_numPrimarySpecies( numPrimarySpecies ),
      m_numSecondarySpecies( numSecondarySpecies ),
      m_numKineticReactions( numKineticReactions ),
      m_primarySpeciesAggregateConcentration( primarySpeciesAggregateConcentration ),
      m_primarySpeciesMobileAggregateConcentration( primarySpeciesMobileAggregateConcentration ),
      m_dPrimarySpeciesAggregateConcentration_dLogPrimarySpeciesConcentrations( dPrimarySpeciesAggregateConcentration_dLogPrimarySpeciesConcentrations ),
      m_dPrimarySpeciesMobileAggregateConcentration_dLogPrimarySpeciesConcentrations( dPrimarySpeciesMobileAggregateConcentration_dLogPrimarySpeciesConcentrations ),
      m_initialPrimarySpeciesConcentration( initialPrimarySpeciesConcentration ),
      m_secondarySpeciesConcentration( secondarySpeciesConcentration ),
      m_kineticReactionRates( kineticReactionRates ),
      m_aggregateSpeciesRates( aggregateSpeciesRates ),
      m_dAggregateSpeciesRates_dLogPrimarySpeciesConcentrations( dAggregateSpeciesRates_dLogPrimarySpeciesConcentrations ),
      m_params( params ),
      m_activityParams( activityParams )
    {}

    using EquilibriumReactionsType = hpcReact::reactionsSystems::EquilibriumReactions< real64, integer, localIndex, ACTIVITY_MODEL >;

    /**
     * @brief Get number of elements in this wrapper.
     * @return number of elements
     */
    GEOS_HOST_DEVICE
    localIndex numElems() const { return m_secondarySpeciesConcentration.size( 0 ); }

    /**
     * @brief Speciate cell @p k at equilibrium.
     * @return whether the equilibrium solve converged
     */
    GEOS_HOST_DEVICE
    bool updateEquilibriumReaction( localIndex const k,
                                    real64 const pressure,
                                    real64 const temperature,
                                    arraySlice1d< real64, compflow::USD_COMP - 1 > const & logPrimarySpeciesConcentration ) const;

    /**
     * @brief Solve for the primary and secondary concentrations at the target aggregates.
     * @return whether the solve converged
     */
    GEOS_HOST_DEVICE
    bool enforceEquilibrium( real64 const pressure,
                             real64 const temperature,
                             arraySlice1d< real64 const, reactivefluid::USD_SPECIES - 2 > const & targetPrimarySpeciesAggregateConcentration,
                             arraySlice1d< real64 const, reactivefluid::USD_SPECIES - 2 > const & initialPrimarySpeciesConcentration,
                             arraySlice1d< real64, compflow::USD_COMP - 1 > const & logPrimarySpeciesConcentration,
                             arraySlice1d< real64 > const & logSecondarySpeciesConcentration ) const;

    GEOS_HOST_DEVICE
    void updateMixedReactionSystem( localIndex const k,
                                    real64 const pressure,
                                    real64 const temperature,
                                    arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & logPrimarySpeciesConcentration,
                                    arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & surfaceArea ) const;

    GEOS_HOST_DEVICE
    void computeAggregateConcentrationsAndRates( real64 const pressure,
                                                 real64 const temperature,
                                                 arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & logPrimarySpeciesConcentration,
                                                 arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & surfaceArea,
                                                 arraySlice1d< real64 > const & logSecondarySpeciesConcentration,
                                                 arraySlice1d< real64, reactivefluid::USD_SPECIES - 2 > const & primarySpeciesAggregateConcentration,
                                                 arraySlice1d< real64, reactivefluid::USD_SPECIES - 2 > const & primarySpeciesMobileAggregateConcentration,
                                                 arraySlice2d< real64, reactivefluid::USD_SPECIES_DC - 2 > const & dPrimarySpeciesAggregateConcentration_dLogPrimarySpeciesConcentrations,
                                                 arraySlice2d< real64, reactivefluid::USD_SPECIES_DC - 2 > const & dPrimarySpeciesMobileAggregateConcentration_dLogPrimarySpeciesConcentrations,
                                                 arraySlice1d< real64, reactivefluid::USD_SPECIES - 2 > const & reactionRates,
                                                 arraySlice2d< real64 > const & dReactionRates_dLogPrimarySpeciesConcentrations,
                                                 arraySlice1d< real64, reactivefluid::USD_SPECIES - 2 > const & aggregateSpeciesRates,
                                                 arraySlice2d< real64, reactivefluid::USD_SPECIES_DC - 2 > const & dAggregateSpeciesRates_dLogPrimarySpeciesConcentrations ) const;

protected:

    integer m_numPrimarySpecies;

    integer m_numSecondarySpecies;

    integer m_numKineticReactions;

    arrayView3d< real64, reactivefluid::USD_SPECIES >  m_primarySpeciesAggregateConcentration;

    arrayView3d< real64, reactivefluid::USD_SPECIES >  m_primarySpeciesMobileAggregateConcentration;

    arrayView4d< real64, reactivefluid::USD_SPECIES_DC >  m_dPrimarySpeciesAggregateConcentration_dLogPrimarySpeciesConcentrations;

    arrayView4d< real64, reactivefluid::USD_SPECIES_DC >  m_dPrimarySpeciesMobileAggregateConcentration_dLogPrimarySpeciesConcentrations;

    arrayView3d< real64 const, reactivefluid::USD_SPECIES > const m_initialPrimarySpeciesConcentration;

    arrayView3d< real64, reactivefluid::USD_SPECIES >  m_secondarySpeciesConcentration;

    arrayView3d< real64, reactivefluid::USD_SPECIES >  m_kineticReactionRates;

    arrayView3d< real64, reactivefluid::USD_SPECIES >  m_aggregateSpeciesRates;

    arrayView4d< real64, reactivefluid::USD_SPECIES_DC >  m_dAggregateSpeciesRates_dLogPrimarySpeciesConcentrations;

    REACTION_PARAMS_TYPE m_params;

    typename ACTIVITY_MODEL::Params m_activityParams;
  };

  /// The kernel wrapper for one system of reactivefluid::ReactionSystemList.
  template< typename SYSTEM >
  using WrapperFor = ReactionKernelWrapper< typename SYSTEM::ReactionParamsType, typename SYSTEM::ActivityType >;

  /// @cond DO_NOT_DOCUMENT
  template< typename LIST >
  struct WrapperVariantHelper;

  template< typename ... SYSTEMS >
  struct WrapperVariantHelper< std::variant< SYSTEMS... > >
  {
    using type = std::variant< WrapperFor< SYSTEMS > ... >;
  };
  /// @endcond

  /// One alternative per system of reactivefluid::ReactionSystemList.
  using ReactionKernelWrapperVariant = typename WrapperVariantHelper< reactivefluid::ReactionSystemList >::type;

  /**
   * @brief Build the kernel wrapper for the chemical system and activity model this fluid was given.
   * @return the wrapper, as the alternative of ReactionKernelWrapperVariant matching that pairing
   *
   * postInputInitialization has already rejected a pairing ReactionSystemList does not hold.
   */
  ReactionKernelWrapperVariant createReactionKernelWrapper() const
  {
    std::optional< ReactionKernelWrapperVariant > wrapper;

    reactivefluid::forEachReactionSystem( [&]( auto system )
    {
      using System = decltype( system );
      if( System::chemicalSystem == m_chemicalSystemType && System::activityModel == m_activityModelType )
      {
        wrapper.emplace( makeReactionKernelWrapper< WrapperFor< System > >( System::reactionParams(),
                                                                            System::activityParams() ) );
      }
    } );

    return std::move( wrapper.value() );
  }

  struct viewKeyStruct : ConstitutiveBase::viewKeyStruct
  {
    static constexpr char const * chemicalSystemNameString() { return "chemicalSystemType"; }
    static constexpr char const * activityModelNameString() { return "activityModelType"; }
    static constexpr char const * solventMassPerSolutionVolumeString() { return "solventMassPerSolutionVolume"; }
  };

protected:

  virtual void postInputInitialization() override;

  virtual void resizeFields( localIndex const size, localIndex const numPts );

  /**
   * @brief Build one kernel wrapper for the given reaction system and activity model.
   */
  template< typename WRAPPER_TYPE, typename REACTION_PARAMS_TYPE, typename ACTIVITY_PARAMS_TYPE >
  WRAPPER_TYPE makeReactionKernelWrapper( REACTION_PARAMS_TYPE const & params,
                                          ACTIVITY_PARAMS_TYPE const & activityParams ) const
  {
    return WRAPPER_TYPE( m_primarySpeciesAggregateConcentration,
                         m_primarySpeciesMobileAggregateConcentration,
                         m_dPrimarySpeciesAggregateConcentration_dLogPrimarySpeciesConcentrations,
                         m_dPrimarySpeciesMobileAggregateConcentration_dLogPrimarySpeciesConcentrations,
                         m_initialPrimarySpeciesConcentration,
                         m_secondarySpeciesConcentration,
                         m_kineticReactionRates,
                         m_aggregateSpeciesRates,
                         m_dAggregateSpeciesRates_dLogPrimarySpeciesConcentrations,
                         m_numPrimarySpecies,
                         m_numSecondarySpecies,
                         m_numKineticReactions,
                         params,
                         activityParams );
  }

  integer m_numPrimarySpecies;

  integer m_numSecondarySpecies;

  integer m_numKineticReactions;

  array3d< real64, constitutive::reactivefluid::LAYOUT_SPECIES >  m_initialPrimarySpeciesConcentration;

  array3d< real64, constitutive::reactivefluid::LAYOUT_SPECIES >  m_secondarySpeciesConcentration;

  array3d< real64, constitutive::reactivefluid::LAYOUT_SPECIES >  m_primarySpeciesAggregateConcentration;

  array3d< real64, constitutive::reactivefluid::LAYOUT_SPECIES >  m_primarySpeciesAggregateConcentration_n;

  array3d< real64, constitutive::reactivefluid::LAYOUT_SPECIES >  m_primarySpeciesMobileAggregateConcentration;

  array4d< real64, constitutive::reactivefluid::LAYOUT_SPECIES_DC >  m_dPrimarySpeciesAggregateConcentration_dLogPrimarySpeciesConcentrations;

  array4d< real64, constitutive::reactivefluid::LAYOUT_SPECIES_DC >  m_dPrimarySpeciesMobileAggregateConcentration_dLogPrimarySpeciesConcentrations;

  array3d< real64, constitutive::reactivefluid::LAYOUT_SPECIES >  m_kineticReactionRates;

  array3d< real64, constitutive::reactivefluid::LAYOUT_SPECIES >  m_aggregateSpeciesRates;

  array4d< real64, constitutive::reactivefluid::LAYOUT_SPECIES_DC >  m_dAggregateSpeciesRates_dLogPrimarySpeciesConcentrations;

  ChemicalSystemType m_chemicalSystemType;

  ActivityModelType m_activityModelType;

  /// TODO: prescribed as a constant for now. The exact factor is
  ///
  ///         rho_s = rho * w
  ///
  ///       where rho_s is this quantity [kg/m^3], rho the solution density [kg/m^3] and w the
  ///       solvent mass fraction [-]. For the carbonate brine EQ3/6 gives 1070.9 * 0.898 = 961.6,
  ///       not the 1000 defaulted here. Ideally rho is a function of pressure, temperature and
  ///       species concentration, and w a function of concentration. The update methods and where
  ///       they should be launched are TBD.
  real64 m_solventMassPerSolutionVolume;
};

// these aliases are useful in constitutive dispatch
using ReactiveCompressibleSinglePhaseFluid = ReactiveSinglePhaseFluid< CompressibleSinglePhaseFluid >;

using ReactiveThermalCompressibleSinglePhaseFluid = ReactiveSinglePhaseFluid< ThermalCompressibleSinglePhaseFluid >;

template< typename BASE >
template< typename REACTION_PARAMS_TYPE, typename ACTIVITY_MODEL >
GEOS_HOST_DEVICE
inline bool
ReactiveSinglePhaseFluid< BASE >::ReactionKernelWrapper< REACTION_PARAMS_TYPE, ACTIVITY_MODEL >::
updateEquilibriumReaction( localIndex const k,
                           real64 const pressure,
                           real64 const temperature,
                           arraySlice1d< real64, compflow::USD_COMP - 1 > const & logPrimarySpeciesConcentration ) const
{
  constexpr integer numSecondarySpecies = REACTION_PARAMS_TYPE::numSecondarySpecies();

  if constexpr ( numSecondarySpecies > 0 )
  {
    stackArray1d< real64, numSecondarySpecies > logSecondarySpeciesConcentration( numSecondarySpecies );

    bool const converged = enforceEquilibrium( pressure, temperature, m_primarySpeciesAggregateConcentration[k][0],
                                               m_initialPrimarySpeciesConcentration[k][0], logPrimarySpeciesConcentration,
                                               logSecondarySpeciesConcentration.toSlice() );

    for( integer i=0; i < numSecondarySpecies; ++i )
    {
      m_secondarySpeciesConcentration[k][0][i] =  LvArray::math::exp( logSecondarySpeciesConcentration[i] );
    }

    return converged;
  }
  else
  {
    GEOS_UNUSED_VAR( k, pressure, temperature, logPrimarySpeciesConcentration );
    return true;
  }

}

template< typename BASE >
template< typename REACTION_PARAMS_TYPE, typename ACTIVITY_MODEL >
GEOS_HOST_DEVICE
inline bool
ReactiveSinglePhaseFluid< BASE >::ReactionKernelWrapper< REACTION_PARAMS_TYPE, ACTIVITY_MODEL >::
enforceEquilibrium( real64 const pressure,
                    real64 const temperature,
                    arraySlice1d< real64 const, reactivefluid::USD_SPECIES - 2 > const & targetPrimarySpeciesAggregateConcentration,
                    arraySlice1d< real64 const, reactivefluid::USD_SPECIES - 2 > const & initialPrimarySpeciesConcentration,
                    arraySlice1d< real64, compflow::USD_COMP - 1 > const & logPrimarySpeciesConcentration,
                    arraySlice1d< real64 > const & logSecondarySpeciesConcentration ) const
{
  GEOS_UNUSED_VAR( pressure );

  integer const numPrimarySpecies = m_numPrimarySpecies;

  stackArray1d< real64, MAX_NUM_SPECIES > logPrimarySpeciesConcentration0( numPrimarySpecies );
  stackArray1d< real64, MAX_NUM_SPECIES > targetPrimarySpeciesAggregateConc( numPrimarySpecies );

  for( integer i=0; i < numPrimarySpecies; ++i )
  {
    targetPrimarySpeciesAggregateConc[i] = targetPrimarySpeciesAggregateConcentration[i];
    logPrimarySpeciesConcentration0[i] = LvArray::math::log( initialPrimarySpeciesConcentration[i] );
  }

  // Solve for the primary and secondary concentrations with equilibrium enforced at the target aggregates
  return EquilibriumReactionsType::enforceEquilibrium_Aggregate( temperature,
                                                                 m_params,
                                                                 m_activityParams,
                                                                 targetPrimarySpeciesAggregateConc,
                                                                 logPrimarySpeciesConcentration0,
                                                                 logPrimarySpeciesConcentration,
                                                                 logSecondarySpeciesConcentration );
}

template< typename BASE >
template< typename REACTION_PARAMS_TYPE, typename ACTIVITY_MODEL >
GEOS_HOST_DEVICE
inline void
ReactiveSinglePhaseFluid< BASE >::ReactionKernelWrapper< REACTION_PARAMS_TYPE, ACTIVITY_MODEL >::
updateMixedReactionSystem( localIndex const k,
                           real64 const pressure,
                           real64 const temperature,
                           arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & logPrimarySpeciesConcentration,
                           arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & surfaceArea ) const
{
  integer const numPrimarySpecies = m_numPrimarySpecies;
  integer const numSecondarySpecies = m_numSecondarySpecies;
  integer const numKineticReactions = m_numKineticReactions;

  stackArray1d< real64, MAX_NUM_SPECIES > logSecondarySpeciesConcentration( numSecondarySpecies );
  stackArray2d< real64, MAX_NUM_KINETIC_REACTIONS * MAX_NUM_SPECIES > dReactionRates_dLogPrimarySpeciesConcentrations( numKineticReactions, numPrimarySpecies );

  computeAggregateConcentrationsAndRates( pressure,
                                          temperature,
                                          logPrimarySpeciesConcentration,
                                          surfaceArea,
                                          logSecondarySpeciesConcentration.toSlice(),
                                          m_primarySpeciesAggregateConcentration[k][0],
                                          m_primarySpeciesMobileAggregateConcentration[k][0],
                                          m_dPrimarySpeciesAggregateConcentration_dLogPrimarySpeciesConcentrations[k][0],
                                          m_dPrimarySpeciesMobileAggregateConcentration_dLogPrimarySpeciesConcentrations[k][0],
                                          m_kineticReactionRates[k][0],
                                          dReactionRates_dLogPrimarySpeciesConcentrations.toSlice(),
                                          m_aggregateSpeciesRates[k][0],
                                          m_dAggregateSpeciesRates_dLogPrimarySpeciesConcentrations[k][0] );

  for( integer i=0; i < numSecondarySpecies; ++i )
  {
    m_secondarySpeciesConcentration[k][0][i] =  LvArray::math::exp( logSecondarySpeciesConcentration[i] );
  }
}

template< typename BASE >
template< typename REACTION_PARAMS_TYPE, typename ACTIVITY_MODEL >
GEOS_HOST_DEVICE
inline void
ReactiveSinglePhaseFluid< BASE >::ReactionKernelWrapper< REACTION_PARAMS_TYPE, ACTIVITY_MODEL >::
computeAggregateConcentrationsAndRates( real64 const pressure,
                                        real64 const temperature,
                                        arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & logPrimarySpeciesConcentration,
                                        arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & surfaceArea,
                                        arraySlice1d< real64 > const & logSecondarySpeciesConcentration,
                                        arraySlice1d< real64, reactivefluid::USD_SPECIES - 2 > const & primarySpeciesAggregateConcentration,
                                        arraySlice1d< real64, reactivefluid::USD_SPECIES - 2 > const & primarySpeciesMobileAggregateConcentration,
                                        arraySlice2d< real64, reactivefluid::USD_SPECIES_DC - 2 > const & dPrimarySpeciesAggregateConcentration_dLogPrimarySpeciesConcentrations,
                                        arraySlice2d< real64, reactivefluid::USD_SPECIES_DC - 2 > const & dPrimarySpeciesMobileAggregateConcentration_dLogPrimarySpeciesConcentrations,
                                        arraySlice1d< real64, reactivefluid::USD_SPECIES - 2 > const & reactionRates,
                                        arraySlice2d< real64 > const & dReactionRates_dLogPrimarySpeciesConcentrations,
                                        arraySlice1d< real64, reactivefluid::USD_SPECIES - 2 > const & aggregateSpeciesRates,
                                        arraySlice2d< real64, reactivefluid::USD_SPECIES_DC - 2 > const & dAggregateSpeciesRates_dLogPrimarySpeciesConcentrations ) const
{
  GEOS_UNUSED_VAR( pressure );

  MixedEquilibriumKineticReactions< real64, localIndex, localIndex, ACTIVITY_MODEL, true >::
  updateMixedSystem( temperature,
                     m_params,
                     m_activityParams,
                     logPrimarySpeciesConcentration,
                     surfaceArea,
                     logSecondarySpeciesConcentration,
                     primarySpeciesAggregateConcentration,
                     primarySpeciesMobileAggregateConcentration,
                     dPrimarySpeciesAggregateConcentration_dLogPrimarySpeciesConcentrations,
                     dPrimarySpeciesMobileAggregateConcentration_dLogPrimarySpeciesConcentrations,
                     reactionRates,
                     dReactionRates_dLogPrimarySpeciesConcentrations,
                     aggregateSpeciesRates,
                     dAggregateSpeciesRates_dLogPrimarySpeciesConcentrations );
}

} // namespace reactivefluid

} // namespace constitutive

} // namespace geos

#endif // GEOS_CONSTITUTIVE_FLUID_REACTIVEFLUID_REACTIVESINGLEPHASEFLUID_HPP_
