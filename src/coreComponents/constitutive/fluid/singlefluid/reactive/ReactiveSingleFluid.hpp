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
 * @file ReactiveSingleFluid.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_SINGLEFLUID_REACTIVE_REACTIVESINGLEFLUID_HPP_
#define GEOS_CONSTITUTIVE_FLUID_SINGLEFLUID_REACTIVE_REACTIVESINGLEFLUID_HPP_


#include "common/format/EnumStrings.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidBase.hpp"
#include "constitutive/fluid/multifluid/reactive/chemicalReactions/EquilibriumReactions.hpp"
#include "constitutive/fluid/multifluid/reactive/chemicalReactions/KineticReactions.hpp"

#include <memory>

namespace geos
{

namespace constitutive
{

class ReactiveSingleFluidUpdate : public SingleFluidBaseUpdate
{
public:

  GEOS_HOST_DEVICE
  void computeChemistry( real64 const pressure,
                         real64 const temperature,
                         arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & primarySpeciesAggregateConcentration,
                         arraySlice1d< real64, compflow::USD_COMP - 1 > const & primarySpeciesConcentration,
                         arraySlice1d< real64, compflow::USD_COMP - 1 > const & secondarySpeciesConcentration,
                         arraySlice1d< real64, compflow::USD_COMP - 1 > const & kineticReactionRates ) const;

  GEOS_HOST_DEVICE
  virtual void updateChemistry( localIndex const k,
                                localIndex const q,
                                real64 const pressure,
                                real64 const temperature,
                                arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition ) const = 0;

  GEOS_HOST_DEVICE
  virtual void updateChemistryLogConc( localIndex const k,
                                       localIndex const q,
                                       real64 const pressure,
                                       real64 const temperature,
                                       arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & logPrimaryConc ) const = 0;

protected:

  ReactiveSingleFluidUpdate( arrayView2d< real64, constitutive::singlefluid::USD_FLUID >  const & density,
                             arrayView3d< real64, constitutive::singlefluid::USD_FLUID_DER >  const & dDensity,
                             arrayView2d< real64, constitutive::singlefluid::USD_FLUID > const & viscosity,
                             arrayView3d< real64, constitutive::singlefluid::USD_FLUID_DER >  const & dViscosity,
                             integer const numPrimarySpecies,
                             //  chemicalReactions::EquilibriumReactions const & equilibriumReactions,
                             //  chemicalReactions::KineticReactions const & kineticReactions,
                             arrayView2d< real64, compflow::USD_COMP > const & primarySpeciesConcentration,
                             arrayView2d< real64, compflow::USD_COMP > const & secondarySpeciesConcentration,
                             arrayView2d< real64, compflow::USD_COMP > const & primarySpeciesAggregateConcentration,
                             arrayView3d< real64, compflow::USD_COMP_DC > const & dPrimarySpeciesAggregateConcentration_dLogPrimaryConc,
                             arrayView2d< real64, compflow::USD_COMP > const & kineticReactionRates )
    : SingleFluidBaseUpdate( density,
                             dDensity,
                             viscosity,
                             dViscosity ),
    m_numPrimarySpecies( numPrimarySpecies ),
    // m_equilibriumReactions( equilibriumReactions.createKernelWrapper() ),
    // m_kineticReactions( kineticReactions.createKernelWrapper() ),
    m_primarySpeciesConcentration( primarySpeciesConcentration ),
    m_secondarySpeciesConcentration( secondarySpeciesConcentration ),
    m_primarySpeciesAggregateConcentration( primarySpeciesAggregateConcentration ),
    m_dPrimarySpeciesAggregateConcentration_dLogPrimaryConc( dPrimarySpeciesAggregateConcentration_dLogPrimaryConc ),
    m_kineticReactionRates( kineticReactionRates )
  {}

  /**
   * @brief Copy constructor.
   */
  ReactiveSingleFluidUpdate( ReactiveSingleFluidUpdate const & ) = default;

  /**
   * @brief Move constructor.
   */
  ReactiveSingleFluidUpdate( ReactiveSingleFluidUpdate && ) = default;

  /**
   * @brief Deleted copy assignment operator
   * @return reference to this object
   */
  ReactiveSingleFluidUpdate & operator=( ReactiveSingleFluidUpdate const & ) = delete;

  /**
   * @brief Deleted move assignment operator
   * @return reference to this object
   */
  ReactiveSingleFluidUpdate & operator=( ReactiveSingleFluidUpdate && ) = delete;

  /// Reaction related terms
  integer m_numPrimarySpecies;

  // chemicalReactions::EquilibriumReactions::KernelWrapper m_equilibriumReactions;

  // chemicalReactions::KineticReactions::KernelWrapper m_kineticReactions;

  arrayView2d< real64, compflow::USD_COMP >  m_primarySpeciesConcentration;

  arrayView2d< real64, compflow::USD_COMP >  m_secondarySpeciesConcentration;

  arrayView2d< real64, compflow::USD_COMP >  m_primarySpeciesAggregateConcentration;

  arrayView3d< real64, compflow::USD_COMP_DC > m_dPrimarySpeciesAggregateConcentration_dLogPrimaryConc;

  arrayView2d< real64, compflow::USD_COMP >  m_kineticReactionRates;
};

class ReactiveSingleFluid : public SingleFluidBase
{
public:

  using exec_policy = serialPolicy;

  ReactiveSingleFluid( string const & name,
                       Group * const parent );

  virtual void saveConvergedState() const override;

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  arrayView2d< real64 const, compflow::USD_COMP > primarySpeciesConcentration() const
  { return m_primarySpeciesConcentration; }

  arrayView2d< real64 const, compflow::USD_COMP > primarySpeciesAggregateConcentration() const
  { return m_primarySpeciesAggregateConcentration; }

  arrayView2d< real64 const, compflow::USD_COMP > primarySpeciesAggregateConcentration_n() const
  { return m_primarySpeciesAggregateConcentration_n; }

  arrayView3d< real64 const, compflow::USD_COMP_DC > dPrimarySpeciesAggregateConcentration_dLogPrimaryConc() const
  { return m_dPrimarySpeciesAggregateConcentration_dLogPrimaryConc; }

  arrayView2d< real64 const, compflow::USD_COMP > secondarySpeciesConcentration() const
  { return m_secondarySpeciesConcentration; }

  arrayView2d< real64 const, compflow::USD_COMP > kineticReactionRates() const
  { return m_kineticReactionRates; }

  integer numPrimarySpecies() const { return m_numPrimarySpecies; }

  integer numSecondarySpecies() const { return m_numSecondarySpecies; }

  integer numKineticReactions() const { return m_numKineticReactions; }


  struct viewKeyStruct : ConstitutiveBase::viewKeyStruct
  {
    static constexpr char const * primarySpeciesNamesString() { return "primarySpeciesNames"; }
  };

protected:

  virtual void postInputInitialization() override;

  void createChemicalReactions();

  virtual void resizeFields( localIndex const size, localIndex const numPts );

  /// Reaction related terms
  array1d< string > m_primarySpeciesNames;

  integer m_numPrimarySpecies;

  integer m_numSecondarySpecies;

  integer m_numKineticReactions;

  // std::unique_ptr< chemicalReactions::EquilibriumReactions > m_equilibriumReactions;

  // std::unique_ptr< chemicalReactions::KineticReactions > m_kineticReactions;

  array2d< real64, constitutive::multifluid::LAYOUT_FLUID >  m_primarySpeciesConcentration;

  array2d< real64, constitutive::multifluid::LAYOUT_FLUID >  m_secondarySpeciesConcentration;

  array2d< real64, constitutive::multifluid::LAYOUT_FLUID >  m_primarySpeciesAggregateConcentration;

  array2d< real64, constitutive::multifluid::LAYOUT_FLUID >  m_primarySpeciesAggregateConcentration_n;

  array3d< real64, constitutive::multifluid::LAYOUT_FLUID_DC >  m_dPrimarySpeciesAggregateConcentration_dLogPrimaryConc;

  array2d< real64, constitutive::multifluid::LAYOUT_FLUID >  m_kineticReactionRates;
};

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void ReactiveSingleFluidUpdate::
  computeChemistry( real64 const pressure,
                    real64 const temperature,
                    arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & primarySpeciesAggregateConcentration,
                    arraySlice1d< real64, compflow::USD_COMP - 1 > const & primarySpeciesConcentration,
                    arraySlice1d< real64, compflow::USD_COMP - 1 > const & secondarySpeciesConcentration,
                    arraySlice1d< real64, compflow::USD_COMP - 1 > const & kineticReactionRates ) const
{
  GEOS_UNUSED_VAR( pressure, temperature, primarySpeciesAggregateConcentration, primarySpeciesConcentration, secondarySpeciesConcentration, kineticReactionRates );

  // // 2. solve for equilibrium
  // m_equilibriumReactions.updateConcentrations( temperature,
  //                                              primarySpeciesAggregateConcentration,
  //                                              primarySpeciesConcentration,
  //                                              secondarySpeciesConcentration );

  // // 3. compute kinetic reaction rates
  // m_kineticReactions.computeReactionRates( temperature,
  //                                          primarySpeciesConcentration,
  //                                          secondarySpeciesConcentration,
  //                                          kineticReactionRates );
}

} // namespace constitutive

} // namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_REACTIVEMULTIFLUID_HPP
