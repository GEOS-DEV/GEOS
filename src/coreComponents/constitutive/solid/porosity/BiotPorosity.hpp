/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file BiotPorosity.hpp
 */

#ifndef GEOS_CONSTITUTIVE_POROSITY_BIOTPOROSITY_HPP_
#define GEOS_CONSTITUTIVE_POROSITY_BIOTPOROSITY_HPP_

#include "PorosityBase.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geos
{
namespace constitutive
{

class BiotPorosityUpdates : public PorosityBaseUpdates
{
public:
  /**
   * @brief Get number of elements in this wrapper.
   * @return number of elements
   */
  GEOS_HOST_DEVICE
  localIndex numElems() const { return m_newPorosity.size( 0 ); }

  /**
   * @brief Get number of gauss points per element.
   * @return number of gauss points per element
   */
  GEOS_HOST_DEVICE
  localIndex numGauss() const { return m_newPorosity.size( 1 ); }

  BiotPorosityUpdates( arrayView2d< real64 > const & newPorosity,
                       arrayView2d< real64 > const & porosity_n,
                       arrayView2d< real64 > const & dPorosity_dPressure,
                       arrayView2d< real64 > const & dPorosity_dTemperature,
                       arrayView2d< real64 > const & initialPorosity,
                       arrayView1d< real64 > const & referencePorosity,
                       arrayView1d< real64 > const & biotCoefficient,
                       arrayView1d< real64 > const & thermalExpansionCoefficient,
                       arrayView2d< real64 > const & meanEffectiveStressIncrement_k,
                       arrayView1d< real64 const > const & averageMeanEffectiveStressIncrement_k,
                       arrayView1d< real64 > const & bulkModulus,
                       real64 const & grainBulkModulus ): PorosityBaseUpdates( newPorosity,
                                                                               porosity_n,
                                                                               dPorosity_dPressure,
                                                                               dPorosity_dTemperature,
                                                                               initialPorosity,
                                                                               referencePorosity ),
    m_grainBulkModulus( grainBulkModulus ),
    m_thermalExpansionCoefficient( thermalExpansionCoefficient ),
    m_biotCoefficient( biotCoefficient ),
    m_bulkModulus( bulkModulus ),
    m_meanEffectiveStressIncrement_k( meanEffectiveStressIncrement_k ),
    m_averageMeanEffectiveStressIncrement_k( averageMeanEffectiveStressIncrement_k )
  {}

  GEOS_HOST_DEVICE
  real64 getBiotCoefficient( localIndex const k ) const { return m_biotCoefficient[k]; }

  GEOS_HOST_DEVICE
  real64 getGrainBulkModulus() const { return m_grainBulkModulus; }

  GEOS_HOST_DEVICE
  real64 dGrainDensity_dPressure() const { return 1.0 / m_grainBulkModulus; }

  GEOS_HOST_DEVICE
  void updateFromPressureTemperatureAndStrain( localIndex const k,
                                               localIndex const q,
                                               real64 const & deltaPressure,
                                               real64 const & deltaTemperature,
                                               real64 const (&strainIncrement)[6],
                                               real64 & dPorosity_dVolStrain,
                                               real64 & dPorosity_dPressure,
                                               real64 & dPorosity_dTemperature ) const
  {
    real64 const biotSkeletonModulusInverse = (m_biotCoefficient[k] - m_referencePorosity[k]) / m_grainBulkModulus;
    real64 const porosityThermalExpansion = 3 * m_thermalExpansionCoefficient[k] * ( m_biotCoefficient[k] - m_referencePorosity[k] );

    real64 const porosity = m_porosity_n[k][q]
                            + m_biotCoefficient[k] * LvArray::tensorOps::symTrace< 3 >( strainIncrement )
                            + biotSkeletonModulusInverse * deltaPressure
                            - porosityThermalExpansion * deltaTemperature;

    dPorosity_dVolStrain = m_biotCoefficient[k];
    dPorosity_dPressure = biotSkeletonModulusInverse;
    dPorosity_dTemperature = -porosityThermalExpansion;

    savePorosity( k, q, porosity, biotSkeletonModulusInverse );
  }

  GEOS_HOST_DEVICE
  void computePorosityFixedStress( real64 const & deltaPressureFromBeginningOfTimeStep,
                                   real64 const & deltaPressureFromLastIteration,
                                   real64 const & deltaTemperatureFromBeginningOfTimeStep,
                                   real64 const & deltaTemperatureFromLastIteration,
                                   real64 const & porosity_n,
                                   real64 const & referencePorosity,
                                   real64 & porosity,
                                   real64 & dPorosity_dPressure,
                                   real64 & dPorosity_dTemperature,
                                   real64 & dPorosity_dVolStrain,
                                   real64 const & biotCoefficient,
                                   real64 const & thermalExpansionCoefficient,
                                   real64 const & meanEffectiveStressIncrement_k,
                                   real64 const & bulkModulus ) const
  {
    real64 const biotSkeletonModulusInverse = (biotCoefficient - referencePorosity) / m_grainBulkModulus;
    real64 const porosityThermalExpansion = 3 * thermalExpansionCoefficient * ( biotCoefficient - referencePorosity );
    real64 const fixedStressPressureCoefficient = biotCoefficient * biotCoefficient / bulkModulus;
    real64 const fixedStressTemperatureCoefficient = 3 * biotCoefficient * thermalExpansionCoefficient;

    porosity = porosity_n
               + biotCoefficient * meanEffectiveStressIncrement_k / bulkModulus // change due to stress increment (at the previous
                                                                                // sequential iteration)
               + biotSkeletonModulusInverse * deltaPressureFromBeginningOfTimeStep // change due to pressure increment
               - porosityThermalExpansion * deltaTemperatureFromBeginningOfTimeStep; // change due to temperature increment
    dPorosity_dPressure = biotSkeletonModulusInverse;
    dPorosity_dTemperature = -porosityThermalExpansion;
    dPorosity_dVolStrain = biotCoefficient / bulkModulus; // ???

    // Fixed-stress part
    porosity += fixedStressPressureCoefficient * deltaPressureFromLastIteration   // fixed-stress pressure term
                + fixedStressTemperatureCoefficient * deltaTemperatureFromLastIteration;   // fixed-stress temperature term
    dPorosity_dPressure += fixedStressPressureCoefficient;
    dPorosity_dTemperature += fixedStressTemperatureCoefficient;
  }

  GEOS_HOST_DEVICE
  virtual void updateFixedStress( localIndex const k,
                                  localIndex const q,
                                  real64 const & pressure,                // current
                                  real64 const & pressure_k,                // last iteration (for sequential)
                                  real64 const & pressure_n,                // last time step
                                  real64 const & temperature,
                                  real64 const & temperature_k,
                                  real64 const & temperature_n,
                                  real64 & dPorosity_dVolStrain ) const
  {
    real64 const deltaPressureFromBeginningOfTimeStep = pressure - pressure_n;
    real64 const deltaPressureFromLastIteration = pressure - pressure_k;
    real64 const deltaTemperatureFromBeginningOfTimeStep = temperature - temperature_n;
    real64 const deltaTemperatureFromLastIteration = temperature - temperature_k;

    computePorosityFixedStress( deltaPressureFromBeginningOfTimeStep,
                                deltaPressureFromLastIteration,
                                deltaTemperatureFromBeginningOfTimeStep,
                                deltaTemperatureFromLastIteration,
                                m_porosity_n[k][q],
                                m_referencePorosity[k],
                                m_newPorosity[k][q],
                                m_dPorosity_dPressure[k][q],
                                m_dPorosity_dTemperature[k][q],
                                dPorosity_dVolStrain,
                                m_biotCoefficient[k],
                                m_thermalExpansionCoefficient[k],
                                m_averageMeanEffectiveStressIncrement_k[k],
                                m_bulkModulus[k] );
  }

  GEOS_HOST_DEVICE
  void updateBiotCoefficientAndAssignBulkModulus( localIndex const k,
                                                  real64 const bulkModulus ) const
  {
    m_bulkModulus[k] = bulkModulus;
    m_biotCoefficient[k] = 1 - bulkModulus / m_grainBulkModulus;
  }

  GEOS_HOST_DEVICE
  void updateMeanEffectiveStressIncrement( localIndex const k,
                                           localIndex const q,
                                           real64 const & meanEffectiveStressIncrement ) const
  {
    m_meanEffectiveStressIncrement_k[k][q] = meanEffectiveStressIncrement;
  }

protected:

  /// Grain bulk modulus (read from XML)
  real64 const m_grainBulkModulus;

  /// View on the thermal expansion coefficients (read from XML)
  arrayView1d< real64 const > const m_thermalExpansionCoefficient;

  /// View on the Biot coefficient (updated by PorousSolid)
  arrayView1d< real64 > const m_biotCoefficient;

  /// View on the bulk modulus (updated by PorousSolid)
  arrayView1d< real64 > const m_bulkModulus;

  /// View on the mean stress increment at quadrature points (updated by PorousSolid)
  arrayView2d< real64 > const m_meanEffectiveStressIncrement_k;

  /// View on the average mean stress increment
  arrayView1d< real64 const > const m_averageMeanEffectiveStressIncrement_k;

};

class BiotPorosity : public PorosityBase
{
public:
  BiotPorosity( string const & name, Group * const parent );

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  static string catalogName() { return "BiotPorosity"; }

  virtual string getCatalogName() const override { return catalogName(); }

  struct viewKeyStruct : public PorosityBase::viewKeyStruct
  {
    static constexpr char const *grainBulkModulusString() { return "grainBulkModulus"; }

    static constexpr char const *meanEffectiveStressIncrementString() { return "meanEffectiveStressIncrement"; }

    static constexpr char const *averageMeanEffectiveStressIncrementString() { return "averageMeanEffectiveStressIncrement"; }

    static constexpr char const *solidBulkModulusString() { return "solidBulkModulus"; }

    static constexpr char const *defaultThermalExpansionCoefficientString() { return "defaultThermalExpansionCoefficient"; }
  } viewKeys;

  virtual void initializeState() const override final;

  virtual void saveConvergedState() const override final;

  virtual void ignoreConvergedState() const override final;

  virtual arrayView1d< real64 const > const getBiotCoefficient() const override final
  {
    return m_biotCoefficient.toViewConst();
  }

  virtual arrayView1d< real64 > const getAverageMeanEffectiveStressIncrement_k() override final
  {
    return m_averageMeanEffectiveStressIncrement_k.toView();
  }

  virtual arrayView2d< real64 const > const getMeanEffectiveStressIncrement_k() const override final
  {
    return m_meanEffectiveStressIncrement_k.toViewConst();
  }

  using KernelWrapper = BiotPorosityUpdates;

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
                          m_biotCoefficient,
                          m_thermalExpansionCoefficient,
                          m_meanEffectiveStressIncrement_k,
                          m_averageMeanEffectiveStressIncrement_k,
                          m_bulkModulus,
                          m_grainBulkModulus );
  }

protected:
  virtual void postProcessInput() override;


  /// Default thermal expansion coefficients (read from XML)
  real64 m_defaultThermalExpansionCoefficient;

  /// Thermal expansion coefficients (read from XML)
  array1d< real64 > m_thermalExpansionCoefficient;

  /// Biot coefficients (update in the update class, not read in input)
  array1d< real64 > m_biotCoefficient;

  /// Bulk modulus (updated in the update class, not read in input)
  array1d< real64 > m_bulkModulus;

  /// Mean stress increment (updated in the update class, not read in input)
  array2d< real64 > m_meanEffectiveStressIncrement_k;

  /// Average mean stress increment (not read in input)
  array1d< real64 > m_averageMeanEffectiveStressIncrement_k;

  /// Grain bulk modulus (read from XML)
  real64 m_grainBulkModulus;
};

}   /* namespace constitutive */

} /* namespace geos */

#endif //GEOS_CONSTITUTIVE_POROSITY_BIOTPOROSITY_HPP_
