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

#include "constitutive/HPCReact/src/reactions/bulkGeneric/ParametersPredefined.hpp"
#include "constitutive/HPCReact/src/reactions/bulkGeneric/EquilibriumReactions.hpp"
#include "constitutive/HPCReact/src/reactions/bulkGeneric/SpeciesUtilities.hpp"
#include <memory>

namespace geos
{

namespace constitutive
{

namespace reactivefluid
{

using namespace hpcReact::bulkGeneric;

enum class ChemicalSystemType : integer
{
  carbonate,
  ultramafic,
  simple
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

  // virtual bool isThermal() const override;

  arrayView2d< real64 const, reactivefluid::USD_COMP > primarySpeciesConcentration() const
  { return m_primarySpeciesConcentration; }

  arrayView2d< real64 const, reactivefluid::USD_COMP > primarySpeciesAggregateConcentration() const
  { return m_primarySpeciesAggregateConcentration; }

  arrayView2d< real64 const, reactivefluid::USD_COMP > primarySpeciesAggregateConcentration_n() const
  { return m_primarySpeciesAggregateConcentration_n; }

  arrayView3d< real64 const, reactivefluid::USD_COMP_DC > dPrimarySpeciesAggregateConcentration_dLogPrimaryConc() const
  { return m_dPrimarySpeciesAggregateConcentration_dLogPrimaryConc; }

  arrayView2d< real64 const, reactivefluid::USD_COMP > secondarySpeciesConcentration() const
  { return m_secondarySpeciesConcentration; }

  integer numPrimarySpecies() const { return m_numPrimarySpecies; }

  integer numSecondarySpecies() const { return m_numSecondarySpecies; }

  /**
   * @brief Kernel wrapper class for ReactiveSinglePhaseFluid.
   */
  template< typename REACTION_PARAMS_TYPE >
  class ReactionKernelWrapper
  {

public:

    ReactionKernelWrapper( arrayView2d< real64, reactivefluid::USD_COMP > const & primarySpeciesConcentration,
                           arrayView2d< real64, reactivefluid::USD_COMP > const & primarySpeciesAggregateConcentration,
                           arrayView2d< real64, reactivefluid::USD_COMP > const & secondarySpeciesConcentration,
                           integer const numPrimarySpecies,
                           integer const numSecondarySpecies,
                           REACTION_PARAMS_TYPE params ):
      m_numPrimarySpecies( numPrimarySpecies ),
      m_numSecondarySpecies( numSecondarySpecies ),
      m_primarySpeciesConcentration( primarySpeciesConcentration ),
      m_primarySpeciesAggregateConcentration( primarySpeciesAggregateConcentration ),
      m_secondarySpeciesConcentration( secondarySpeciesConcentration ),
      m_params( params )
    {}

    using EquilibriumReactionsType = hpcReact::bulkGeneric::EquilibriumReactions< real64, integer, localIndex >;

    /**
     * @brief Get number of elements in this wrapper.
     * @return number of elements
     */
    GEOS_HOST_DEVICE
    localIndex numElems() const { return m_secondarySpeciesConcentration.size( 0 ); }

    void updateEquilibriumReaction( localIndex const k,
                                    real64 const pressure,
                                    real64 const temperature,
                                    arraySlice1d< real64 const, reactivefluid::USD_COMP - 1 > const & primarySpeciesTotalConcentration,
                                    arraySlice1d< real64, reactivefluid::USD_COMP - 1 > const & logPrimarySpeciesConcentration ) const;

    void enforceEquilibrium( real64 const pressure,
                             real64 const temperature,
                             arraySlice1d< real64 const, reactivefluid::USD_COMP - 1 > const & primarySpeciesTotalConcentration,
                             arraySlice1d< real64, reactivefluid::USD_COMP - 1 > const & primarySpeciesConcentration,
                             arraySlice1d< real64, reactivefluid::USD_COMP - 1 > const & secondarySpeciesConcentration ) const;

    void computeReactionRates( real64 const pressure,
                               real64 const temperature,
                               arraySlice1d< real64 const, reactivefluid::USD_COMP - 1 > const & primarySpeciesTotalConcentration,
                               arraySlice1d< real64, reactivefluid::USD_COMP - 1 > const & primarySpeciesConcentration,
                               arraySlice1d< real64, reactivefluid::USD_COMP - 1 > const & secondarySpeciesConcentration ) const;

protected:

    integer m_numPrimarySpecies;

    integer m_numSecondarySpecies;

    arrayView2d< real64, reactivefluid::USD_COMP >  m_primarySpeciesConcentration;

    arrayView2d< real64, reactivefluid::USD_COMP >  m_primarySpeciesAggregateConcentration;

    arrayView2d< real64, reactivefluid::USD_COMP >  m_secondarySpeciesConcentration;

    arrayView2d< real64, reactivefluid::USD_COMP >  m_kineticReactionRates;

    REACTION_PARAMS_TYPE m_params;
  };

  std::variant<
    typename ReactiveSinglePhaseFluid< BASE >::template ReactionKernelWrapper< hpcReact::bulkGeneric::carbonateSystemType >,
    typename ReactiveSinglePhaseFluid< BASE >::template ReactionKernelWrapper< hpcReact::bulkGeneric::simpleTestType > >
  createReactionKernelWrapper() const
  {
    switch( m_chemicalSystemType )
    {
      case ChemicalSystemType::carbonate:
        return ReactionKernelWrapper< hpcReact::bulkGeneric::carbonateSystemType >( m_primarySpeciesConcentration,
                                                                                    m_primarySpeciesAggregateConcentration,
                                                                                    m_secondarySpeciesConcentration,
                                                                                    m_numPrimarySpecies,
                                                                                    m_numSecondarySpecies,
                                                                                    carbonateSystem );
      default:
        return ReactionKernelWrapper< hpcReact::bulkGeneric::simpleTestType >( m_primarySpeciesConcentration,
                                                                               m_primarySpeciesAggregateConcentration,
                                                                               m_secondarySpeciesConcentration,
                                                                               m_numPrimarySpecies,
                                                                               m_numSecondarySpecies,
                                                                               simpleTestRateParams );
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

  array2d< real64, constitutive::reactivefluid::LAYOUT_FLUID >  m_primarySpeciesConcentration;

  array2d< real64, constitutive::reactivefluid::LAYOUT_FLUID >  m_secondarySpeciesConcentration;

  array2d< real64, constitutive::reactivefluid::LAYOUT_FLUID >  m_primarySpeciesAggregateConcentration;

  array2d< real64, constitutive::reactivefluid::LAYOUT_FLUID >  m_primarySpeciesAggregateConcentration_n;

  array3d< real64, constitutive::reactivefluid::LAYOUT_FLUID_DC >  m_dPrimarySpeciesAggregateConcentration_dLogPrimaryConc;

  array2d< real64, constitutive::reactivefluid::LAYOUT_FLUID >  m_kineticReactionRates;

  ChemicalSystemType m_chemicalSystemType;
};

// these aliases are useful in constitutive dispatch
using ReactiveSinglePhaseFluidBase = ReactiveSinglePhaseFluid< SingleFluidBase >; 

using ReactiveCompressibleSinglePhaseFluid = ReactiveSinglePhaseFluid< CompressibleSinglePhaseFluid >; 

using ReactiveThermalCompressibleSinglePhaseFluid = ReactiveSinglePhaseFluid< ThermalCompressibleSinglePhaseFluid >; 

template< typename BASE >
template< typename REACTION_PARAMS_TYPE >
inline void
ReactiveSinglePhaseFluid< BASE >::ReactionKernelWrapper< REACTION_PARAMS_TYPE >::
updateEquilibriumReaction( localIndex const k,
                           real64 const pressure,
                           real64 const temperature,
                           arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & primarySpeciesAggregateConcentration,
                           arraySlice1d< real64, compflow::USD_COMP - 1 > const & logPrimarySpeciesConcentration ) const
{
  enforceEquilibrium( pressure, temperature, primarySpeciesAggregateConcentration, m_primarySpeciesConcentration[k], m_secondarySpeciesConcentration[k] );

  for( int i=0; i < m_numPrimarySpecies; i++ )
  {
    real64 const primarySpeciesConc_i = m_primarySpeciesConcentration[k][i];

    logPrimarySpeciesConcentration[i] = LvArray::math::log( primarySpeciesConc_i );
  }
}

template< typename BASE >
template< typename REACTION_PARAMS_TYPE >
inline void
ReactiveSinglePhaseFluid< BASE >::ReactionKernelWrapper< REACTION_PARAMS_TYPE >::
enforceEquilibrium( real64 const pressure,
                    real64 const temperature,
                    arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & primarySpeciesAggregateConcentration,
                    arraySlice1d< real64, compflow::USD_COMP - 1 > const & primarySpeciesConcentration,
                    arraySlice1d< real64, compflow::USD_COMP - 1 > const & secondarySpeciesConcentration ) const
{
  GEOS_UNUSED_VAR( pressure, temperature, secondarySpeciesConcentration );
  // 1. We enforce equilibrium
  EquilibriumReactionsType::enforceEquilibrium_Extents( 298.15, m_params, primarySpeciesAggregateConcentration, primarySpeciesConcentration );
  // // 2. We calculate the secondary species concentration
  // hpcReact::bulkGeneric::calculateLogSecondarySpeciesConcentration< real64,
  //                                                                   localIndex,
  //                                                                   localIndex >( m_params, primarySpeciesConcentration, secondarySpeciesConcentration );
}

template< typename BASE >
template< typename REACTION_PARAMS_TYPE >
inline void
ReactiveSinglePhaseFluid< BASE >::ReactionKernelWrapper< REACTION_PARAMS_TYPE >::
computeReactionRates( real64 const pressure,
                      real64 const temperature,
                      arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & primarySpeciesAggregateConcentration,
                      arraySlice1d< real64, compflow::USD_COMP - 1 > const & primarySpeciesConcentration,
                      arraySlice1d< real64, compflow::USD_COMP - 1 > const & secondarySpeciesConcentration ) const
{
  GEOS_UNUSED_VAR( pressure );
  // 1. We enforce equilibrium
  EquilibriumReactionsType::enforceEquilibrium_Extents( temperature, m_params, primarySpeciesAggregateConcentration, primarySpeciesConcentration );
  // 2. We calculate the secondary species concentration
  hpcReact::bulkGeneric::utilities_impl::calculateLogSecondarySpeciesConcentration( m_params, primarySpeciesConcentration, secondarySpeciesConcentration );
}

ENUM_STRINGS( ChemicalSystemType,
              "carbonate",
              "ultramafic",
              "simple" );

} // namespace reactivefluid

} // namespace constitutive

} // namespace geos

#endif // GEOS_CONSTITUTIVE_FLUID_REACTIVEFLUID_REACTIVESINGLEPHASEFLUID_HPP_
