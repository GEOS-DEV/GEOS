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
 * @file BiotReactivePorosity.hpp
 */

#ifndef GEOS_CONSTITUTIVE_POROSITY_BIOTREACTIVEPOROSITY_HPP_
#define GEOS_CONSTITUTIVE_POROSITY_BIOTREACTIVEPOROSITY_HPP_

#include "ReactivePorosityBase.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geos
{
namespace constitutive
{

class BiotReactivePorosityUpdates : public ReactivePorosityBaseUpdates
{
public:

  BiotReactivePorosityUpdates( arrayView2d< real64 > const & newPorosity,
                               arrayView2d< real64 > const & porosity_n,
                               arrayView2d< real64 > const & dPorosity_dPressure,
                               arrayView2d< real64 > const & dPorosity_dTemperature,
                               arrayView2d< real64 > const & initialPorosity,
                               arrayView1d< real64 > const & referencePorosity,
                               arrayView3d< real64, reactivefluid::USD_SPECIES > const & volumeFractions,
                               arrayView3d< real64 const, reactivefluid::USD_SPECIES > const & initialVolumeFractions,
                               arrayView3d< real64 const, reactivefluid::USD_SPECIES > const & volumeFractions_n,
                               integer const numKineticReactions,
                               arrayView1d< real64 const > const & molarWeights,
                               arrayView1d< real64 const > const & mineralDensities,
                               arrayView1d< real64 > const & biotCoefficient,
                               arrayView1d< real64 > const & bulkModulus,
                               arrayView1d< real64 > const & grainBulkModulus,
                               arrayView1d< real64 > const & mineralBulkModulus,
                               arrayView1d< real64 > const & mineralPressure,
                               arrayView1d< real64 > const & mineralPressure_n,
                               real64 const fluidCompressibility,
                               arrayView2d< real64 > const & meanEffectiveStressIncrement_k,
                               integer const fixedPorosity ): ReactivePorosityBaseUpdates( newPorosity,
                                                                                            porosity_n,
                                                                                            dPorosity_dPressure,
                                                                                            dPorosity_dTemperature,
                                                                                            initialPorosity,
                                                                                            referencePorosity,
                                                                                            volumeFractions,
                                                                                            initialVolumeFractions,
                                                                                            volumeFractions_n,
                                                                                            numKineticReactions,
                                                                                            molarWeights,
                                                                                            mineralDensities,
                                                                                            bulkModulus,
                                                                                            meanEffectiveStressIncrement_k,
                                                                                            fixedPorosity ),
    m_grainBulkModulus( grainBulkModulus ),
    m_biotCoefficient( biotCoefficient ),
    m_mineralBulkModulus( mineralBulkModulus ),
    m_mineralPressure( mineralPressure ),
    m_mineralPressure_n( mineralPressure_n ),
    m_fluidCompressibility( fluidCompressibility )
  {}

  GEOS_HOST_DEVICE
  real64 getBiotCoefficient( localIndex const k ) const { return m_biotCoefficient[k]; }

  GEOS_HOST_DEVICE
  real64 getGrainBulkModulus( localIndex const k ) const { return m_grainBulkModulus[k]; }

  /// @return the fluid bulk modulus K_f = 1/compressibility [Pa] (0 if incompressible)
  GEOS_HOST_DEVICE
  real64 getFluidBulkModulus( localIndex const k ) const
  {
    GEOS_UNUSED_VAR( k );
    return m_fluidCompressibility > 0.0 ? 1.0 / m_fluidCompressibility : 0.0;
  }

  GEOS_HOST_DEVICE
  real64 getPoreMineralPressure( localIndex const k ) const { return m_mineralPressure[k]; }

  GEOS_HOST_DEVICE
  void computePorosityFixedStress( real64 const & pressure,
                                   real64 const & pressure_k,
                                   real64 const & pressure_n,
                                   real64 const & porosity_n,
                                   real64 const & referencePorosity,
                                   real64 & porosity,
                                   real64 & dPorosity_dPressure,
                                   real64 const & biotCoefficient,
                                   real64 const & meanEffectiveStressIncrement_k,
                                   real64 const & bulkModulus,
                                   real64 const & grainBulkModulus,
                                   real64 const & mineralBulkModulus,
                                   real64 const & fluidBulkModulus,
                                   real64 const & reactionPorosityIncrement ) const
  {
    GEOS_UNUSED_VAR( pressure_k );

    real64 const biotSkeletonModulusInverse = (biotCoefficient - referencePorosity) / grainBulkModulus;
    real64 const porosityMultiplierInverse = 1 / ( 1 + biotSkeletonModulusInverse*mineralBulkModulus/referencePorosity );
    real64 const poreFluidMineralBulkModRatio = ( fluidBulkModulus + mineralBulkModulus ) / fluidBulkModulus;

    real64 const anelasticStrainIncrement = -reactionPorosityIncrement;

    porosity = porosity_n
               // change due to stress increment
               + biotCoefficient * meanEffectiveStressIncrement_k / bulkModulus * porosityMultiplierInverse
               // change due to pressure increment
               + biotSkeletonModulusInverse * poreFluidMineralBulkModRatio * ( pressure - pressure_n ) * porosityMultiplierInverse
               // change due to mineral pressure increment
               + biotSkeletonModulusInverse * 3 * anelasticStrainIncrement * mineralBulkModulus * porosityMultiplierInverse
               // change due to mineral volume fraction increment
               + reactionPorosityIncrement;  

    dPorosity_dPressure = biotSkeletonModulusInverse * poreFluidMineralBulkModRatio * porosityMultiplierInverse;
  }

  GEOS_HOST_DEVICE
  void computePoreMineralPressure( real64 const & mineralPressure_n,
                                   real64 & mineralPressure,
                                   real64 & dMineralPres_dMeanEffStressIncre,
                                   real64 const & referencePorosity,
                                   real64 const & biotCoefficient,
                                   real64 const & deltaPressureFromLastStep,
                                   real64 const & meanEffectiveStressIncrement,
                                   real64 const & bulkModulus,
                                   real64 const & grainBulkModulus,
                                   real64 const & mineralBulkModulus,
                                   real64 const & fluidBulkModulus,
                                   real64 const & anelasticStrainIncrement ) const
  {
    // GEOS_UNUSED_VAR( meanEffectiveStressIncrement, bulkModulus );
    real64 const biotSkeletonModulusInverse = (biotCoefficient - referencePorosity) / grainBulkModulus;
    real64 const mineralPressureMultiplier = mineralBulkModulus / ( referencePorosity + biotSkeletonModulusInverse*mineralBulkModulus );

    mineralPressure = mineralPressure_n
                      // change due to inelastic strain increment
                      + mineralPressureMultiplier * referencePorosity * anelasticStrainIncrement * 3
                      // change due to stress increment
                      - mineralPressureMultiplier * biotCoefficient * meanEffectiveStressIncrement / bulkModulus
                      // change due to pressure increment
                      - mineralPressureMultiplier * biotSkeletonModulusInverse * deltaPressureFromLastStep
                      + mineralPressureMultiplier * referencePorosity / fluidBulkModulus * deltaPressureFromLastStep;

    dMineralPres_dMeanEffStressIncre = -mineralPressureMultiplier * biotCoefficient / bulkModulus;
    // dMineralPres_dMeanEffStressIncre = 0.0;
  }

  GEOS_HOST_DEVICE
  void updatePoreMineralPressure( localIndex const k,
                                  localIndex const q,
                                  real64 const & deltaPressureFromLastStep,
                                  real64 const & meanEffectiveStressIncrement,
                                  real64 const & anelasticStrainIncrement,
                                  real64 & dMineralPres_dMeanEffStressIncre ) const
  {
    GEOS_UNUSED_VAR( q );

    computePoreMineralPressure( m_mineralPressure_n[k],
                                m_mineralPressure[k],
                                dMineralPres_dMeanEffStressIncre,
                                m_referencePorosity[k],
                                m_biotCoefficient[k],
                                deltaPressureFromLastStep,
                                meanEffectiveStressIncrement,
                                m_bulkModulus[k],
                                m_grainBulkModulus[k],
                                m_mineralBulkModulus[k],
                                getFluidBulkModulus( k ),
                                anelasticStrainIncrement );
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
                                  arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & kineticReactionMolarIncrements ) const override final
  {
    // 1. Update the porosity due to reactions
    real64 reactionPorosityIncrement = 0.0;

    ReactivePorosityBaseUpdates::computePorosityFromReaction( kineticReactionMolarIncrements,
                                                              m_volumeFractions[k][q],
                                                              m_volumeFractions_n[k][q],
                                                              reactionPorosityIncrement,
                                                              m_numKineticReactions,
                                                              m_molarWeights,
                                                              m_mineralDensities );

    // 2. Update the porosity due to solid, pore mineral, and pore fluid pressure
    // Currently ignore thermal effects
    GEOS_UNUSED_VAR( temperature, temperature_k, temperature_n );
    if( !m_fixedPorosity )
    {
      computePorosityFixedStress( pressure, pressure_k, pressure_n,
                                  m_porosity_n[k][q],
                                  m_referencePorosity[k],
                                  m_newPorosity[k][q],
                                  m_dPorosity_dPressure[k][q],
                                  m_biotCoefficient[k],
                                  m_meanEffectiveStressIncrement_k[k][q],
                                  m_bulkModulus[k],
                                  m_grainBulkModulus[k],
                                  m_mineralBulkModulus[k],
                                  getFluidBulkModulus( k ),
                                  reactionPorosityIncrement );
    }
  }

  GEOS_HOST_DEVICE
  void updateBiotCoefficientAndAssignModuli( localIndex const k,
                                             real64 const bulkModulus ) const
  {
    m_bulkModulus[k] = bulkModulus;

    m_biotCoefficient[k] =  1.0 - bulkModulus / m_grainBulkModulus[k];
  }



protected:

  /// View on the grain bulk modulus (read from XML)
  arrayView1d< real64 > const m_grainBulkModulus;

  /// View on the Biot coefficient (updated by PorousSolid)
  arrayView1d< real64 > const m_biotCoefficient;

  /// View on the mineral bulk modulus (read from XML)
  arrayView1d< real64 > const m_mineralBulkModulus;

  /// View on the mineral pressure
  arrayView1d< real64 > const m_mineralPressure;

  /// View on the mineral pressure at the previous timestep
  arrayView1d< real64 > const m_mineralPressure_n;

  /// Fluid compressibility [1/Pa] (constant, pushed in from the flow solver's fluid model)
  real64 const m_fluidCompressibility;
};

class BiotReactivePorosity : public ReactivePorosityBase
{
public:
  BiotReactivePorosity( string const & name, dataRepository::Group * const parent );

  static string catalogName() { return "BiotReactivePorosity"; }

  virtual string getCatalogName() const override { return catalogName(); }

  struct viewKeyStruct : public ReactivePorosityBase::viewKeyStruct
  {
    static constexpr char const *defaultMineralBulkModulusString() { return "defaultMineralBulkModulus"; }

    static constexpr char const *mineralBulkModulusString() { return "mineralBulkModulus"; }

    static constexpr char const *mineralPressure_nString() { return "mineralPressure_n"; }

    static constexpr char const *mineralPressureString() { return "mineralPressure"; }

    static constexpr char const *defaultGrainBulkModulusString() { return "defaultGrainBulkModulus"; }
  };

  virtual void initializeState() const override final;

  virtual void saveConvergedState() const override final;

  using KernelWrapper = BiotReactivePorosityUpdates;

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
                          m_biotCoefficient,
                          m_bulkModulus,
                          m_grainBulkModulus,
                          m_mineralBulkModulus,
                          m_mineralPressure,
                          m_mineralPressure_n,
                          m_fluidCompressibility,
                          m_meanEffectiveStressIncrement_k,
                          m_fixedPorosity );
  }

  /**
   * @brief Store the (constant) fluid compressibility read from the flow solver's fluid model.
   * @param fluidCompressibility the fluid compressibility [1/Pa]
   */
  void setFluidCompressibility( real64 const fluidCompressibility )
  { m_fluidCompressibility = fluidCompressibility; }

protected:
  virtual void postInputInitialization() override;

  /// Biot coefficients (update in the update class, not read in input)
  array1d< real64 > m_biotCoefficient;

  /// Grain bulk modulus (read from XML)
  real64 m_defaultGrainBulkModulus;

  /// Grain bulk modulus (can be specified in XML)
  array1d< real64 > m_grainBulkModulus;

  /// Mineral bulk modulus (read from XML)
  real64 m_defaultMineralBulkModulus;

  /// Mineral bulk modulus (can be specified in XML)
  array1d< real64 > m_mineralBulkModulus;

  /// Mineral pressure
  array1d< real64 > m_mineralPressure;

  /// Mineral pressure at the previous timestep
  array1d< real64 > m_mineralPressure_n;

  /// Fluid compressibility [1/Pa], pushed in from the flow solver's fluid model (0 until set)
  real64 m_fluidCompressibility = 0.0;
};

}   /* namespace constitutive */

} /* namespace geos */

#endif //GEOS_CONSTITUTIVE_POROSITY_BIOTREACTIVEPOROSITY_HPP_
