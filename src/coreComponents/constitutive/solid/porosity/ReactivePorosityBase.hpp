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
 * @file ReactivePorosityBase.hpp
 */

#ifndef GEOS_CONSTITUTIVE_POROSITY_REACTIVEPOROSITYBASE_HPP_
#define GEOS_CONSTITUTIVE_POROSITY_REACTIVEPOROSITYBASE_HPP_

#include "PorosityBase.hpp"

#include "constitutive/fluid/reactivefluid/ReactiveFluidLayouts.hpp"

namespace geos
{
namespace constitutive
{

class ReactivePorosityBaseUpdates : public PorosityBaseUpdates
{
public:

  ReactivePorosityBaseUpdates( arrayView2d< real64 > const & newPorosity,
                               arrayView2d< real64 const > const & porosity_n,
                               arrayView2d< real64 > const & dPorosity_dPressure,
                               arrayView2d< real64 > const & dPorosity_dTemperature,
                               arrayView2d< real64 const > const & initialPorosity,
                               arrayView1d< real64 const > const & referencePorosity,
                               arrayView3d< real64, reactivefluid::USD_SPECIES > const & volumeFractions,
                               arrayView3d< real64 const, reactivefluid::USD_SPECIES > const & initialVolumeFractions,
                               arrayView3d< real64 const, reactivefluid::USD_SPECIES > const & volumeFractions_n,
                               integer const numKineticReactions,
                               arrayView1d< real64 const > const & molarWeights,
                               arrayView1d< real64 const > const & mineralDensities,
                               arrayView1d< real64 > const & bulkModulus,
                               arrayView2d< real64 > const & meanEffectiveStressIncrement_k,
                               integer const fixedPorosity,
                               arrayView2d< real64 > const & intactPorosity,
                               arrayView2d< real64 const > const & intactPorosity_n,
                               real64 const crackPorosity ):
    PorosityBaseUpdates( newPorosity,
                         porosity_n,
                         dPorosity_dPressure,
                         dPorosity_dTemperature,
                         initialPorosity,
                         referencePorosity ),
    m_volumeFractions( volumeFractions ),
    m_initialVolumeFractions( initialVolumeFractions ),
    m_volumeFractions_n( volumeFractions_n ),
    m_numKineticReactions( numKineticReactions ),
    m_molarWeights( molarWeights ),
    m_mineralDensities( mineralDensities ),
    m_bulkModulus( bulkModulus ),
    m_meanEffectiveStressIncrement_k( meanEffectiveStressIncrement_k ),
    m_fixedPorosity( fixedPorosity ),
    m_intactPorosity( intactPorosity ),
    m_intactPorosity_n( intactPorosity_n ),
    m_crackPorosity( crackPorosity )
  {}

  GEOS_HOST_DEVICE
  void computePorosityFromReaction( arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & kineticReactionMolarIncrements,
                                    arraySlice1d< real64, reactivefluid::USD_SPECIES - 2 > const & volumeFractions,
                                    arraySlice1d< real64 const, reactivefluid::USD_SPECIES - 2 > const & volumeFractions_n,
                                    real64 & reactionPorosityIncrement,
                                    integer const & numKineticReactions,
                                    arrayView1d< real64 const > const & molarWeights,
                                    arrayView1d< real64 const > const & mineralDensities ) const
  {
    for( integer r=0; r < numKineticReactions; ++r )
    {
      real64 const volumeFractionIncrement = -kineticReactionMolarIncrements[r] * molarWeights[r]/mineralDensities[r];
      volumeFractions[r] = volumeFractions_n[r] + volumeFractionIncrement;

      reactionPorosityIncrement -= volumeFractionIncrement;
    }
  }

  GEOS_HOST_DEVICE
  void updateFromReactions( localIndex const k,
                            localIndex const q,
                            arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & kineticReactionMolarIncrements ) const
  {
    real64 reactionPorosityIncrement = 0.0;

    computePorosityFromReaction( kineticReactionMolarIncrements,
                                 m_volumeFractions[k][q],
                                 m_volumeFractions_n[k][q],
                                 reactionPorosityIncrement,
                                 m_numKineticReactions,
                                 m_molarWeights,
                                 m_mineralDensities );

    if( !m_fixedPorosity )
    {
      m_intactPorosity[k][q] = LvArray::math::min( LvArray::math::max( m_intactPorosity_n[k][q] + reactionPorosityIncrement, 0.0 ), 1.0 );
    }
    m_newPorosity[k][q] = m_intactPorosity[k][q];
  }

  // this function is used in flow solver
  // it uses average stress increment (element-based)
  GEOS_HOST_DEVICE
  virtual void updateFixedStress( localIndex const k,
                                  localIndex const q,
                                  real64 const & pressure,
                                  real64 const & pressure_k,
                                  real64 const & pressure_n,
                                  real64 const & temperature,
                                  real64 const & temperature_k,
                                  real64 const & temperature_n,
                                  arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & kineticReactionMolarIncrements ) const
  {
    // For now, we ignore the pressure or temperature dependence but just follow Evans et al. for this eigenstrain approach
    GEOS_UNUSED_VAR( pressure, pressure_k, pressure_n, temperature, temperature_k, temperature_n );

    // 1. Update the porosity due to reactions
    real64 reactionPorosityIncrement = 0.0;

    computePorosityFromReaction( kineticReactionMolarIncrements,
                                 m_volumeFractions[k][q],
                                 m_volumeFractions_n[k][q],
                                 reactionPorosityIncrement,
                                 m_numKineticReactions,
                                 m_molarWeights,
                                 m_mineralDensities );

    // 2. Update the porosity due to solid deformation
    if( !m_fixedPorosity )
    {
      // + m_meanEffectiveStressIncrement_k[k][q]/m_bulkModulus[k];
      m_intactPorosity[k][q] = LvArray::math::min( LvArray::math::max( m_intactPorosity_n[k][q] + reactionPorosityIncrement, 0.0 ), 1.0 );
    }
    m_newPorosity[k][q] = m_intactPorosity[k][q];
  }

  /**
   * @brief Add the damage contribution to the porosity: phi = max( phi_f, phi_f + ( phi_c - phi_f ) d ), with phi_f the
   *        damage-free (reactive/poromechanical) porosity. The result is the porosity of the material, used for the fluid
   *        storage, the permeability and the diffusive flux; phi_f is kept separately so that the damage contribution
   *        does not accumulate in the incremental porosity update. Must be called after the porosity update.
   *        No-op if the crack porosity phi_c is not set (phi_c <= 0).
   * @param[in] k the element index
   * @param[in] q the quadrature point index
   * @param[in] damage the damage in [0, 1]
   */
  GEOS_HOST_DEVICE
  void updateFromDamage( localIndex const k,
                         localIndex const q,
                         real64 const damage ) const
  {
    if( m_crackPorosity <= 0.0 )
    {
      return;
    }
    real64 const intactPorosity = m_intactPorosity[k][q];
    real64 const damagePorosity = intactPorosity + ( m_crackPorosity - intactPorosity ) * damage;
    // dPorosity_dPressure keeps d phi_f / dp (exact derivative ( 1 - d ) d phi_f / dp): only the Jacobian is approximate
    if( damagePorosity > intactPorosity )
    {
      m_newPorosity[k][q] = damagePorosity;
    }
  }

  /// @return the damage-free (reactive/poromechanical) porosity
  GEOS_HOST_DEVICE
  real64 getIntactPorosity( localIndex const k,
                            localIndex const q ) const
  {
    return m_intactPorosity[k][q];
  }

  GEOS_HOST_DEVICE
  void updateSolidBulkModulus( localIndex const k,
                               real64 const bulkModulus ) const
  {
    m_bulkModulus[k] = bulkModulus;
  }

  GEOS_HOST_DEVICE
  void updateMeanEffectiveStressIncrement( localIndex const k,
                                           localIndex const q,
                                           real64 const & meanEffectiveStressIncrement ) const
  {
    m_meanEffectiveStressIncrement_k[k][q] = meanEffectiveStressIncrement;
  }

  GEOS_HOST_DEVICE
  inline
  real64 getVolumeFractionForMineral( localIndex const k,
                                      localIndex const q,
                                      localIndex const r ) const
  {
    return m_volumeFractions[k][q][r];
  }

  GEOS_HOST_DEVICE
  inline
  real64 getInitialVolumeFractionForMineral( localIndex const k,
                                             localIndex const q,
                                             localIndex const r ) const
  {
    return m_initialVolumeFractions[k][q][r];
  }

  GEOS_HOST_DEVICE
  inline
  real64 getMolarWeights( localIndex const r ) const
  {
    return m_molarWeights[r];
  }

  GEOS_HOST_DEVICE
  inline
  real64 getMineralDensities( localIndex const r ) const
  {
    return m_mineralDensities[r];
  }

protected:

  arrayView3d< real64, reactivefluid::USD_SPECIES > m_volumeFractions;
  arrayView3d< real64 const, reactivefluid::USD_SPECIES > m_initialVolumeFractions;
  arrayView3d< real64 const, reactivefluid::USD_SPECIES > const m_volumeFractions_n;

  integer const m_numKineticReactions;
  arrayView1d< real64 const > const m_molarWeights;
  arrayView1d< real64 const > const m_mineralDensities;

  arrayView1d< real64 > const m_bulkModulus;
  arrayView2d< real64 > const m_meanEffectiveStressIncrement_k;

  integer const m_fixedPorosity;

  /// damage-free (reactive/poromechanical) porosity and its value at the previous converged step
  arrayView2d< real64 > const m_intactPorosity;
  arrayView2d< real64 const > const m_intactPorosity_n;

  /// porosity phi_c of the fully damaged material (0 = no damage contribution)
  real64 const m_crackPorosity;
};


class ReactivePorosityBase : public PorosityBase
{
public:
  ReactivePorosityBase( string const & name, Group * const parent );

  virtual std::unique_ptr< ConstitutiveBase >
  deliverClone( string const & name,
                dataRepository::Group * const parent ) const override;

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  static string catalogName() { return "ReactivePorosity"; }

  virtual string getCatalogName() const override { return catalogName(); }

  virtual void saveConvergedState() const override;
  virtual void ignoreConvergedState() const override;

  integer numKineticReactions() const { return m_numKineticReactions; }

  virtual void initializeState() const override;

  struct viewKeyStruct : public PorosityBase::viewKeyStruct
  {
    static constexpr char const * defaultInitialVolumeFractionsString() { return "defaultInitialVolumeFractions"; }
    static constexpr char const * initialVolumeFractionsString() { return "initialVolumeFractions"; }
    static constexpr char const * volumeFractionsString() { return "volumeFractions"; }
    static constexpr char const * volumeFractions_nString() { return "volumeFractions_n"; }
    static constexpr char const * molarWeightsString() { return "molarWeights"; }
    static constexpr char const * mineralDensitiesString() { return "mineralDensities"; }
    static constexpr char const * solidBulkModulusString() { return "solidBulkModulus"; }
    static constexpr char const * meanEffectiveStressIncrement_kString() { return "meanEffectiveStressIncrement_k"; }
    static constexpr char const * fixedPorosityString() { return "fixedPorosity"; }
    static constexpr char const * intactPorosityString() { return "intactPorosity"; }
    static constexpr char const * intactPorosity_nString() { return "intactPorosity_n"; }
    static constexpr char const * crackPorosityString() { return "crackPorosity"; }
  } viewKeys;


  using KernelWrapper = ReactivePorosityBaseUpdates;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelUpdates() const
  {
    return KernelWrapper( m_newPorosity,
                          m_porosity_n,
                          m_dPorosity_dPressure,
                          m_dPorosity_dTemperature,
                          m_initialPorosity,
                          m_referencePorosity,
                          m_volumeFractions,
                          m_initialVolumeFractions,
                          m_volumeFractions_n,
                          m_numKineticReactions,
                          m_molarWeights,
                          m_mineralDensities,
                          m_bulkModulus,
                          m_meanEffectiveStressIncrement_k,
                          m_fixedPorosity,
                          m_intactPorosity,
                          m_intactPorosity_n,
                          m_crackPorosity );
  }


protected:
  virtual void postInputInitialization() override;

  virtual void resizeFields( localIndex const size, localIndex const numPts );

  array1d< real64 > m_defaultInitialVolumeFractions;

  array3d< real64, constitutive::reactivefluid::LAYOUT_SPECIES > m_volumeFractions;
  array3d< real64, constitutive::reactivefluid::LAYOUT_SPECIES > m_initialVolumeFractions;
  array3d< real64, constitutive::reactivefluid::LAYOUT_SPECIES > m_volumeFractions_n;

  integer m_numKineticReactions;
  array1d< real64 > m_molarWeights;
  array1d< real64 > m_mineralDensities;

  array1d< real64 > m_bulkModulus;
  array2d< real64 > m_meanEffectiveStressIncrement_k;

  integer m_fixedPorosity;

  /// damage-free (reactive/poromechanical) porosity, equal to the porosity without damage
  array2d< real64 > m_intactPorosity;
  array2d< real64 > m_intactPorosity_n;

  /// porosity phi_c of the fully damaged material (0 = no damage contribution)
  real64 m_crackPorosity;
};


}/* namespace constitutive */

} /* namespace geos */


#endif //GEOS_CONSTITUTIVE_POROSITY_REACTIVEPOROSITYBASE_HPP_
