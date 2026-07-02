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

#include "constitutive/HPCReact/src/reactions/geochemistry/GeochemicalSystems.hpp"
#include "constitutive/HPCReact/src/reactions/exampleSystems/BulkGeneric.hpp"
#include "constitutive/HPCReact/src/reactions/exampleSystems/ChainGeneric.hpp"
#include "constitutive/HPCReact/src/reactions/exampleSystems/MoMasBenchmark.hpp"
#include "constitutive/HPCReact/src/reactions/reactionsSystems/EquilibriumReactions.hpp"
#include "constitutive/HPCReact/src/reactions/reactionsSystems/MixedEquilibriumKineticReactions.hpp"
#include "constitutive/HPCReact/src/reactions/massActions/MassActions.hpp"
#include <memory>

namespace geos
{

namespace constitutive
{

namespace reactivefluid
{

using namespace hpcReact::reactionsSystems;

enum class ChemicalSystemType : integer
{
  carbonate,
  carbonateAllEquilibrium,
  ultramafic,
  serpentinization,
  kineticCarbonate,
  momasEasy,
  momasMedium,
  chainSerialAllKinetic
};

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

  real64 solventDensity() const { return m_solventDensity; }

  /**
   * @brief Kernel wrapper class for ReactiveSinglePhaseFluid.
   */
  template< typename REACTION_PARAMS_TYPE >
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
                           REACTION_PARAMS_TYPE params ):
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
      m_params( params )
    {}

    using EquilibriumReactionsType = hpcReact::reactionsSystems::EquilibriumReactions< real64, integer, localIndex >;

    /**
     * @brief Get number of elements in this wrapper.
     * @return number of elements
     */
    GEOS_HOST_DEVICE
    localIndex numElems() const { return m_secondarySpeciesConcentration.size( 0 ); }

    GEOS_HOST_DEVICE
    void updateEquilibriumReaction( localIndex const k,
                                    real64 const pressure,
                                    real64 const temperature,
                                    arraySlice1d< real64, compflow::USD_COMP - 1 > const & logPrimarySpeciesConcentration ) const;

    GEOS_HOST_DEVICE
    void enforceEquilibrium( real64 const pressure,
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
  };

  std::variant<
    typename ReactiveSinglePhaseFluid< BASE >::template ReactionKernelWrapper< hpcReact::geochemistry::ultramaficSystemType >,
    typename ReactiveSinglePhaseFluid< BASE >::template ReactionKernelWrapper< hpcReact::geochemistry::carbonateSystemType >,
    typename ReactiveSinglePhaseFluid< BASE >::template ReactionKernelWrapper< hpcReact::geochemistry::carbonateSystemAllEquilibriumType >,
    typename ReactiveSinglePhaseFluid< BASE >::template ReactionKernelWrapper< hpcReact::ChainGeneric::serialAllKineticType >,
    typename ReactiveSinglePhaseFluid< BASE >::template ReactionKernelWrapper< hpcReact::MoMasBenchmark::mediumCaseType >,
    typename ReactiveSinglePhaseFluid< BASE >::template ReactionKernelWrapper< hpcReact::MoMasBenchmark::easyCaseType >,
    typename ReactiveSinglePhaseFluid< BASE >::template ReactionKernelWrapper< hpcReact::geochemistry::kineticCarbonateSystemType > >
  createReactionKernelWrapper() const
  {
    using namespace hpcReact::geochemistry;
    using namespace hpcReact::MoMasBenchmark;
    using namespace hpcReact::bulkGeneric;
    using namespace hpcReact::ChainGeneric;
    switch( m_chemicalSystemType )
    {
      case ChemicalSystemType::ultramafic:
        return ReactionKernelWrapper< ultramaficSystemType >( m_primarySpeciesAggregateConcentration,
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
                                                              ultramaficSystem );

      case ChemicalSystemType::carbonate:
        return ReactionKernelWrapper< carbonateSystemType >( m_primarySpeciesAggregateConcentration,
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
                                                             carbonateSystem );
      case ChemicalSystemType::carbonateAllEquilibrium:
        return ReactionKernelWrapper< carbonateSystemAllEquilibriumType >( m_primarySpeciesAggregateConcentration,
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
                                                                           carbonateSystemAllEquilibrium );
      case ChemicalSystemType::serpentinization:
        return ReactionKernelWrapper< serpentinizationSystemType >( m_primarySpeciesAggregateConcentration,
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
                                                                    serpentinizationSystem );
      case ChemicalSystemType::kineticCarbonate:
        return ReactionKernelWrapper< kineticCarbonateSystemType >( m_primarySpeciesAggregateConcentration,
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
                                                                     kineticCarbonateSystem );
      case ChemicalSystemType::chainSerialAllKinetic:
        return ReactionKernelWrapper< serialAllKineticType >( m_primarySpeciesAggregateConcentration,
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
                                                              serialAllKineticParams );
      case ChemicalSystemType::momasMedium:
        return ReactionKernelWrapper< mediumCaseType >( m_primarySpeciesAggregateConcentration,
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
                                                        mediumCaseParams );
      default:
        return ReactionKernelWrapper< easyCaseType >( m_primarySpeciesAggregateConcentration,
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
                                                      easyCaseParams );
    }
  }

  struct viewKeyStruct : ConstitutiveBase::viewKeyStruct
  {
    static constexpr char const * chemicalSystemNameString() { return "chemicalSystemType"; }
  };

protected:

  virtual void postInputInitialization() override;

  virtual void resizeFields( localIndex const size, localIndex const numPts );

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

  real64 m_solventDensity;
};

// these aliases are useful in constitutive dispatch
using ReactiveCompressibleSinglePhaseFluid = ReactiveSinglePhaseFluid< CompressibleSinglePhaseFluid >;

using ReactiveThermalCompressibleSinglePhaseFluid = ReactiveSinglePhaseFluid< ThermalCompressibleSinglePhaseFluid >;

template< typename BASE >
template< typename REACTION_PARAMS_TYPE >
GEOS_HOST_DEVICE
inline void
ReactiveSinglePhaseFluid< BASE >::ReactionKernelWrapper< REACTION_PARAMS_TYPE >::
updateEquilibriumReaction( localIndex const k,
                           real64 const pressure,
                           real64 const temperature,
                           arraySlice1d< real64, compflow::USD_COMP - 1 > const & logPrimarySpeciesConcentration ) const
{
  integer const numSecondarySpecies = m_numSecondarySpecies;

  if( numSecondarySpecies > 0 )
  {
    stackArray1d< real64, MAX_NUM_SPECIES > logSecondarySpeciesConcentration( numSecondarySpecies );

    enforceEquilibrium( pressure, temperature, m_primarySpeciesAggregateConcentration[k][0], m_initialPrimarySpeciesConcentration[k][0], logPrimarySpeciesConcentration,
                        logSecondarySpeciesConcentration.toSlice() );

    for( integer i=0; i < numSecondarySpecies; ++i )
    {
      m_secondarySpeciesConcentration[k][0][i] =  LvArray::math::exp( logSecondarySpeciesConcentration[i] );
    }
  }
  else
  {
    GEOS_UNUSED_VAR( k, pressure, temperature, logPrimarySpeciesConcentration );
  }

}

template< typename BASE >
template< typename REACTION_PARAMS_TYPE >
GEOS_HOST_DEVICE
inline void
ReactiveSinglePhaseFluid< BASE >::ReactionKernelWrapper< REACTION_PARAMS_TYPE >::
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

  // 1. We enforce equilibrium
  EquilibriumReactionsType::enforceEquilibrium_Aggregate( temperature, m_params, targetPrimarySpeciesAggregateConc, logPrimarySpeciesConcentration0, logPrimarySpeciesConcentration );

  // 2. We calculate the secondary species concentration
  hpcReact::massActions::calculateLogSecondarySpeciesConcentration< real64,
                                                                    localIndex,
                                                                    localIndex >( m_params, logPrimarySpeciesConcentration, logSecondarySpeciesConcentration );
}

template< typename BASE >
template< typename REACTION_PARAMS_TYPE >
GEOS_HOST_DEVICE
inline void
ReactiveSinglePhaseFluid< BASE >::ReactionKernelWrapper< REACTION_PARAMS_TYPE >::
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
template< typename REACTION_PARAMS_TYPE >
GEOS_HOST_DEVICE
inline void
ReactiveSinglePhaseFluid< BASE >::ReactionKernelWrapper< REACTION_PARAMS_TYPE >::
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

  MixedEquilibriumKineticReactions< real64, localIndex, localIndex, true >::
  updateMixedSystem( temperature,
                     m_params,
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

ENUM_STRINGS( ChemicalSystemType,
              "carbonate",
              "carbonateAllEquilibrium",
              "ultramafic",
              "serpentinization",
              "kineticCarbonate",
              "momasEasy",
              "momasMedium",
              "chainSerialAllKinetic" );

} // namespace reactivefluid

} // namespace constitutive

} // namespace geos

#endif // GEOS_CONSTITUTIVE_FLUID_REACTIVEFLUID_REACTIVESINGLEPHASEFLUID_HPP_
