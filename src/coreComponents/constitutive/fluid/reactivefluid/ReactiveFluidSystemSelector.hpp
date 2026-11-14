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
 * @file ReactiveFluidSystemSelector.hpp
 * @brief The chemical systems ReactiveSinglePhaseFluid can instantiate, and the activity model each
 *        may be paired with.
 *
 * Adding a system, or an activity model for one, is one entry below plus its name in
 * ReactionSystemList. The kernel wrapper variant, the dispatch, and the check that rejects an
 * unsupported pairing all follow from that list.
 */
#ifndef GEOS_CONSTITUTIVE_FLUID_REACTIVEFLUID_REACTIVEFLUIDSYSTEMSELECTOR_HPP_
#define GEOS_CONSTITUTIVE_FLUID_REACTIVEFLUID_REACTIVEFLUIDSYSTEMSELECTOR_HPP_

#include "common/DataTypes.hpp"
#include "common/format/EnumStrings.hpp"

#include "constitutive/HPCReact/src/reactions/geochemistry/GeochemicalSystems.hpp"
#include "constitutive/HPCReact/src/reactions/exampleSystems/BulkGeneric.hpp"
#include "constitutive/HPCReact/src/reactions/exampleSystems/ChainGeneric.hpp"
#include "constitutive/HPCReact/src/reactions/exampleSystems/MoMasBenchmark.hpp"

#include <utility>
#include <variant>

namespace geos
{

namespace constitutive
{

namespace reactivefluid
{

/**
 * @brief The chemical system a ReactiveSinglePhaseFluid solves.
 */
enum class ChemicalSystemType : integer
{
  carbonate,
  carbonateAllEquilibrium,
  ultramafic,
  momasEasy,
  momasMedium,
  chainSerialAllKinetic
};

/**
 * @brief The activity model applied to it.
 * @details ``identity`` is the identity activity model, gamma = 1 and a_w = 1. ``bdot`` needs ion
 *          size and b-dot parameters, which only the geochemical systems carry.
 */
enum class ActivityModelType : integer
{
  identity,
  bdot
};

ENUM_STRINGS( ChemicalSystemType,
              "carbonate",
              "carbonateAllEquilibrium",
              "ultramafic",
              "momasEasy",
              "momasMedium",
              "chainSerialAllKinetic" );

ENUM_STRINGS( ActivityModelType,
              "identity",
              "bdot" );

/**
 * @brief What every entry below shares: the two HPCReact types, and the two enum values the XML
 *        input selects it by.
 * @tparam SYSTEM the ChemicalSystemType it answers to
 * @tparam MODEL the ActivityModelType it answers to
 * @tparam REACTION_PARAMS the HPCReact reaction parameters type
 * @tparam ACTIVITY_MODEL the HPCReact activity model type
 *
 * Each entry adds the two parameter objects, as accessors returning by value. They cannot be
 * template arguments: the HPCReact objects are namespace-scope constexpr and so have internal
 * linkage, which would make the entry a different type in every translation unit. Returning a
 * copy also matches how the kernel wrapper stores them.
 */
template< ChemicalSystemType SYSTEM,
          ActivityModelType MODEL,
          typename REACTION_PARAMS,
          typename ACTIVITY_MODEL >
struct ReactionSystemEntry
{
  /// the HPCReact reaction parameters type
  using ReactionParamsType = REACTION_PARAMS;

  /// the HPCReact activity model type
  using ActivityType = ACTIVITY_MODEL;

  /// the chemical system this answers to
  static constexpr ChemicalSystemType chemicalSystem = SYSTEM;

  /// the activity model this answers to
  static constexpr ActivityModelType activityModel = MODEL;
};

/// Carbonate, 16 species with calcite a pure phase, identity activity model.
struct CarbonateIdentityEntry : ReactionSystemEntry< ChemicalSystemType::carbonate,
                                                     ActivityModelType::identity,
                                                     hpcReact::geochemistry::carbonateSystemType,
                                                     hpcReact::geochemistry::carbonateNosolidIdentityActivityType >
{
  static ReactionParamsType reactionParams() { return hpcReact::geochemistry::carbonateSystem; }
  static ActivityType::Params activityParams() { return hpcReact::geochemistry::carbonateNosolidIdentityActivityParams; }
};

/// Carbonate, 16 species with calcite a pure phase, B-dot activity model on EQ3/6 parameters.
struct CarbonateBdotEntry : ReactionSystemEntry< ChemicalSystemType::carbonate,
                                                 ActivityModelType::bdot,
                                                 hpcReact::geochemistry::carbonateSystemType,
                                                 hpcReact::geochemistry::carbonateNosolidActivityType >
{
  static ReactionParamsType reactionParams() { return hpcReact::geochemistry::carbonateSystem; }
  static ActivityType::Params activityParams() { return hpcReact::geochemistry::carbonateNosolidActivityParamsEQ36; }
};

/// Carbonate, 17 species with calcite among them, identity activity model.
struct CarbonateAllEquilibriumIdentityEntry : ReactionSystemEntry< ChemicalSystemType::carbonateAllEquilibrium,
                                                                   ActivityModelType::identity,
                                                                   hpcReact::geochemistry::carbonateSystemAllEquilibriumType,
                                                                   hpcReact::geochemistry::carbonateIdentityActivityType >
{
  static ReactionParamsType reactionParams() { return hpcReact::geochemistry::carbonateSystemAllEquilibrium; }
  static ActivityType::Params activityParams() { return hpcReact::geochemistry::carbonateIdentityActivityParams; }
};

/// Carbonate, 17 species with calcite among them, B-dot activity model on EQ3/6 parameters.
struct CarbonateAllEquilibriumBdotEntry : ReactionSystemEntry< ChemicalSystemType::carbonateAllEquilibrium,
                                                               ActivityModelType::bdot,
                                                               hpcReact::geochemistry::carbonateSystemAllEquilibriumType,
                                                               hpcReact::geochemistry::carbonateActivityType >
{
  static ReactionParamsType reactionParams() { return hpcReact::geochemistry::carbonateSystemAllEquilibrium; }
  static ActivityType::Params activityParams() { return hpcReact::geochemistry::carbonateActivityParamsEQ36; }
};

/// Ultramafic, 20 species with the five minerals pure phases, identity activity model.
struct UltramaficIdentityEntry : ReactionSystemEntry< ChemicalSystemType::ultramafic,
                                                      ActivityModelType::identity,
                                                      hpcReact::geochemistry::ultramaficSystemType,
                                                      hpcReact::geochemistry::ultramaficIdentityActivityType >
{
  static ReactionParamsType reactionParams() { return hpcReact::geochemistry::ultramaficSystem; }
  static ActivityType::Params activityParams() { return hpcReact::geochemistry::ultramaficIdentityActivityParams; }
};

/// Ultramafic, 20 species with the five minerals pure phases, B-dot activity model on EQ3/6 parameters.
struct UltramaficBdotEntry : ReactionSystemEntry< ChemicalSystemType::ultramafic,
                                                  ActivityModelType::bdot,
                                                  hpcReact::geochemistry::ultramaficSystemType,
                                                  hpcReact::geochemistry::ultramaficActivityType >
{
  static ReactionParamsType reactionParams() { return hpcReact::geochemistry::ultramaficSystem; }
  static ActivityType::Params activityParams() { return hpcReact::geochemistry::ultramaficActivityParamsEQ36; }
};

/// MoMaS easy benchmark, 12 abstract species, identity activity model.
struct MomasEasyIdentityEntry : ReactionSystemEntry< ChemicalSystemType::momasEasy,
                                                     ActivityModelType::identity,
                                                     hpcReact::MoMasBenchmark::easyCaseType,
                                                     hpcReact::MoMasBenchmark::easyCaseIdentityActivityType >
{
  static ReactionParamsType reactionParams() { return hpcReact::MoMasBenchmark::easyCaseParams; }
  static ActivityType::Params activityParams() { return hpcReact::MoMasBenchmark::easyCaseIdentityActivityParams; }
};

/// MoMaS medium benchmark, 14 abstract species, identity activity model.
struct MomasMediumIdentityEntry : ReactionSystemEntry< ChemicalSystemType::momasMedium,
                                                       ActivityModelType::identity,
                                                       hpcReact::MoMasBenchmark::mediumCaseType,
                                                       hpcReact::MoMasBenchmark::mediumCaseIdentityActivityType >
{
  static ReactionParamsType reactionParams() { return hpcReact::MoMasBenchmark::mediumCaseParams; }
  static ActivityType::Params activityParams() { return hpcReact::MoMasBenchmark::mediumCaseIdentityActivityParams; }
};

/// Generic serial chain, 3 species all kinetic, identity activity model.
struct ChainSerialIdentityEntry : ReactionSystemEntry< ChemicalSystemType::chainSerialAllKinetic,
                                                       ActivityModelType::identity,
                                                       hpcReact::ChainGeneric::serialAllKineticType,
                                                       hpcReact::ChainGeneric::serialAllKineticIdentityActivityType >
{
  static ReactionParamsType reactionParams() { return hpcReact::ChainGeneric::serialAllKineticParams; }
  static ActivityType::Params activityParams() { return hpcReact::ChainGeneric::serialAllKineticIdentityActivityParams; }
};

/**
 * @brief Every pairing the fluid model can build. B-dot appears only for the geochemical systems,
 *        the ones carrying ion size and b-dot parameters.
 */
using ReactionSystemList = std::variant< CarbonateIdentityEntry,
                                         CarbonateBdotEntry,
                                         CarbonateAllEquilibriumIdentityEntry,
                                         CarbonateAllEquilibriumBdotEntry,
                                         UltramaficIdentityEntry,
                                         UltramaficBdotEntry,
                                         MomasEasyIdentityEntry,
                                         MomasMediumIdentityEntry,
                                         ChainSerialIdentityEntry >;

/**
 * @brief Apply @p func to one instance of every system in ReactionSystemList.
 * @tparam FUNC the callable type
 * @tparam INDICES the alternatives of ReactionSystemList
 * @param func a generic callable taking one entry
 */
template< typename FUNC, std::size_t ... INDICES >
void forEachReactionSystem( FUNC && func, std::index_sequence< INDICES... > )
{
  ( func( std::variant_alternative_t< INDICES, ReactionSystemList >{} ), ... );
}

/**
 * @copydoc forEachReactionSystem
 */
template< typename FUNC >
void forEachReactionSystem( FUNC && func )
{
  forEachReactionSystem( std::forward< FUNC >( func ),
                         std::make_index_sequence< std::variant_size_v< ReactionSystemList > >{} );
}

/**
 * @brief Whether @p chemicalSystem may be paired with @p activityModel.
 * @param chemicalSystem the chemical system
 * @param activityModel the activity model
 * @return true if ReactionSystemList holds that pairing
 */
inline bool isSupportedReactionSystem( ChemicalSystemType const chemicalSystem,
                                       ActivityModelType const activityModel )
{
  bool supported = false;
  forEachReactionSystem( [&]( auto system )
  {
    using System = decltype( system );
    supported = supported || ( System::chemicalSystem == chemicalSystem &&
                               System::activityModel == activityModel );
  } );
  return supported;
}

} // namespace reactivefluid

} // namespace constitutive

} // namespace geos

#endif // GEOS_CONSTITUTIVE_FLUID_REACTIVEFLUID_REACTIVEFLUIDSYSTEMSELECTOR_HPP_
