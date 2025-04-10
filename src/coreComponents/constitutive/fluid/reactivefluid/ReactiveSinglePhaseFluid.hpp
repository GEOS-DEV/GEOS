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

#include "constitutive/HPCReact/src/reactions/bulkGeneric/ParametersPredefined.hpp"
#include "constitutive/HPCReact/src/reactions/bulkGeneric/EquilibriumReactions.hpp"
#include <memory>

namespace geos
{

namespace constitutive
{
  
namespace reactivefluid
{

enum class ChemicalSystemType : integer
{
  carbonate,
  ultramafic
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

  virtual bool isThermal() const override;

  arrayView2d< real64 const, compflow::USD_COMP > primarySpeciesConcentration() const
  { return m_primarySpeciesConcentration; }

  arrayView2d< real64 const, compflow::USD_COMP > secondarySpeciesConcentration() const
  { return m_secondarySpeciesConcentration; }

  integer numPrimarySpecies() const { return m_numPrimarySpecies; }

  integer numSecondarySpecies() const { return m_numSecondarySpecies; }

  auto createKernelWraper() const 
  {

  
  std::variant< KernelWrapper< hpcReact::bubulkGeneric::carbonateSystemType >,
                KernelWrapper< hpcReact::bubulkGeneric::simpleTestType > >
  createKernelWrapper() const;
  
  /**
   * @brief Kernel wrapper class for ReactiveSinglePhaseFluid.
   */
  template< typename REACTION_PARAMS_TYPE >
  class KernelWrapper : public BASE::KernelWrapper
  {

public:

    KernelWrapper( arrayView2d< real64, reactivefluid::USD_COMP > const & secondarySpeciesConcentration,
                   integer const numPrimarySpecies,
                   integer const numSecondarySpecies,
                   REACTION_PARAMS_TYPE params ) :
      m_secondarySpeciesConcentration( secondarySpeciesConcentration ),
      m_numPrimarySpecies( numPrimarySpecies ),
      m_numSecondarySpecies( numSecondarySpecies ),
      m_params( params )
    {}

    using EquilibriumReactionsType = hpcReact::bulkGeneric::EquilibriumReactions< real64, integer, localIndex >;

    void computeChemistry( real64 const pressure,
                           real64 const temperature,
                           arraySlice1d< real64 const, reactivefluid::USD_COMP - 1 > const & primarySpeciesTotalConcentration,
                           arraySlice1d< real64, reactivefluid::USD_COMP - 1 > const & primarySpeciesConcentration,
                           arraySlice1d< real64, reactivefluid::USD_COMP - 1 > const & secondarySpeciesConcentration ) const;

protected:

    integer m_numPrimarySpecies;

    integer m_numSecondarySpecies;

    arrayView2d< real64, reactivefluid::USD_COMP >  m_secondarySpeciesConcentration;

    arrayView2d< real64, reactivefluid::USD_COMP >  m_kineticReactionRates;

    REACTION_PARAMS_TYPE m_params;
  };

  struct viewKeyStruct : ConstitutiveBase::viewKeyStruct
  {
    static constexpr char const * chemicalSystemNameString() { return "chemicalSystemType"; }
  };

protected:

  virtual void postInputInitialization() override;

  virtual void resizeFields( localIndex const size, localIndex const numPts ) override;

  integer m_numPrimarySpecies;
    
  integer m_numSecondarySpecies;

  array2d< real64, constitutive::reactivefluid::LAYOUT_FLUID >  m_secondarySpeciesConcentration;

  array2d< real64, constitutive::reactivefluid::LAYOUT_FLUID >  m_kineticReactionRates;

  ChemicalSystemType m_chemicalSystemType;
};

template< typename BASE >
inline void
ReactiveSinglePhaseFluid< BASE >::KernelWrapper::
  enforceEquilibrium( real64 const pressure,
                      real64 const temperature,
                      arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & primarySpeciesAggregateConcentration,
                      arraySlice1d< real64, compflow::USD_COMP - 1 > const & primarySpeciesConcentration,
                      arraySlice1d< real64, compflow::USD_COMP - 1 > const & secondarySpeciesConcentration ) const
{
  GEOS_UNUSED_VAR( pressure );
  // 1. We enforce equilibrium 
  EquilibriumReactionsType::enforceEquilibrium_Extents( temperature, m_params, primarySpeciesAggregateConcentration, primarySpeciesConcentration );
  // 2. We calculate the secondary species concentration 
  speciesUtilities::calculateLogSecondarySpeciesConcentration( m_params, primarySpeciesConcentration, secondarySpeciesConcentration );
}

template< typename BASE >
inline void
ReactiveSinglePhaseFluid< BASE >::KernelWrapper::
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
  speciesUtilities::calculateLogSecondarySpeciesConcentration( m_params, primarySpeciesConcentration, secondarySpeciesConcentration );
}

ENUM_STRINGS( ChemicalSystemType,
              "carbonate",
              "ultramafic" );

} // namespace reactivefluid

} // namespace constitutive

} // namespace geos

#endif // GEOS_CONSTITUTIVE_FLUID_REACTIVEFLUID_REACTIVESINGLEPHASEFLUID_HPP_
