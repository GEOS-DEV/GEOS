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

#ifndef GEOSX_CONSTITUTIVE_POROSITY_BIOTPOROSITY_HPP_
#define GEOSX_CONSTITUTIVE_POROSITY_BIOTPOROSITY_HPP_

#include "PorosityBase.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geosx
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
  GEOSX_HOST_DEVICE
  localIndex numElems() const { return m_newPorosity.size( 0 ); }

  /**
   * @brief Get number of gauss points per element.
   * @return number of gauss points per element
   */
  GEOSX_HOST_DEVICE
  localIndex numGauss() const { return m_newPorosity.size( 1 ); }

  BiotPorosityUpdates( arrayView2d< real64 > const & newPorosity,
                       arrayView2d< real64 > const & porosity_n,
                       arrayView2d< real64 > const & dPorosity_dPressure,
                       arrayView2d< real64 > const & dPorosity_dTemperature,
                       arrayView2d< real64 > const & initialPorosity,
                       arrayView1d< real64 > const & referencePorosity,
                       arrayView1d< real64 > const & biotCoefficient,
                       arrayView1d< real64 > const & thermalExpansionCoefficient,
                       arrayView2d< real64 > const & meanStressIncrement,
                       arrayView1d< real64 > const & bulkModulus,
                       real64 const & grainBulkModulus ): PorosityBaseUpdates( newPorosity,
                                                                               porosity_n,
                                                                               dPorosity_dPressure,
                                                                               dPorosity_dTemperature,
                                                                               initialPorosity,
                                                                               referencePorosity ),
    m_biotCoefficient( biotCoefficient ),
    m_thermalExpansionCoefficient( thermalExpansionCoefficient ),
    m_meanStressIncrement( meanStressIncrement ),
    m_bulkModulus( bulkModulus ),
    m_grainBulkModulus( grainBulkModulus )
  {}

  GEOSX_HOST_DEVICE
  real64 getBiotCoefficient( localIndex const k ) const { return m_biotCoefficient[k]; }

  GEOSX_HOST_DEVICE
  real64 getGrainBulkModulus() const { return m_grainBulkModulus; }

  GEOSX_HOST_DEVICE
  real64 dGrainDensity_dPressure() const { return 1.0 / m_grainBulkModulus; }

  GEOSX_HOST_DEVICE
  void updateFromPressureTemperatureAndStrain( localIndex const k,
                                               localIndex const q,
                                               real64 const & deltaPressure,
                                               real64 const & deltaTemperature,
                                               real64 const (&strainIncrement)[6],
                                               real64 const & thermalExpansionCoefficient,
                                               real64 & dPorosity_dVolStrain,
                                               real64 & dPorosity_dPressure,
                                               real64 & dPorosity_dTemperature ) const
  {
    real64 const biotSkeletonModulusInverse = (m_biotCoefficient[k] - m_referencePorosity[k]) / m_grainBulkModulus;
    real64 const porosityThermalExpansion = 3 * thermalExpansionCoefficient * m_biotCoefficient[k];

    real64 porosity = m_porosity_n[k][q]
                            + m_biotCoefficient[k] * LvArray::tensorOps::symTrace< 3 >( strainIncrement )
                            + biotSkeletonModulusInverse * deltaPressure
                            - porosityThermalExpansion * deltaTemperature;

    dPorosity_dVolStrain = m_biotCoefficient[k];
    dPorosity_dPressure = biotSkeletonModulusInverse;
    dPorosity_dTemperature = -porosityThermalExpansion;

    real64 constexpr eps = 1e-5;
    if( porosity < eps )
    {
      porosity = eps;
      dPorosity_dVolStrain   = 0.0;
      dPorosity_dPressure    = 0.0;
      dPorosity_dTemperature = 0.0;
    }

    savePorosity( k, q, porosity, biotSkeletonModulusInverse );
  }

  GEOSX_HOST_DEVICE
  void computePorosity( real64 const & pressure,
                        real64 const & temperature,
                        real64 & porosity,
                        real64 & dPorosity_dPressure,
                        real64 & dPorosity_dTemperature,
                        real64 const & biotCoefficient,
                        real64 const & thermalExpansionCoefficient,
                        real64 const & meanStressIncrement,
                        real64 const & bulkModulus,
                        real64 const & porosity_n ) const
  {
    real64 const biotSkeletonModulusInverse = (biotCoefficient - porosity_n) / m_grainBulkModulus;
    real64 const porosityThermalExpansion = 3 * thermalExpansionCoefficient * biotCoefficient;

    porosity = porosity_n + biotSkeletonModulusInverse * pressure + biotCoefficient * biotCoefficient / bulkModulus * pressure
               + porosityThermalExpansion * temperature
               + biotCoefficient * meanStressIncrement / bulkModulus;

    dPorosity_dPressure = biotSkeletonModulusInverse + biotCoefficient * biotCoefficient / bulkModulus;
    dPorosity_dTemperature = porosityThermalExpansion;
  }

  GEOSX_HOST_DEVICE
  virtual void updateFromPressureAndTemperature( localIndex const k,
                                                 localIndex const q,
                                                 real64 const & pressure,
                                                 real64 const & pressure_n,
                                                 real64 const & temperature,
                                                 real64 const & temperature_n ) const override final
  {
    real64 const deltaPressure    = pressure - pressure_n;
    real64 const deltaTemperature = temperature - temperature_n;

    computePorosity( deltaPressure,
                     deltaTemperature,
                     m_newPorosity[k][q],
                     m_dPorosity_dPressure[k][q],
                     m_dPorosity_dTemperature[k][q],
                     m_biotCoefficient[k],
                     m_thermalExpansionCoefficient[k],
                     m_meanStressIncrement[k][q],
                     m_bulkModulus[k],
                     m_porosity_n[k][q] );
  }

  GEOSX_HOST_DEVICE
  void updateBiotCoefficient( localIndex const k,
                              real64 const bulkModulus ) const
  {
    m_bulkModulus[k] = bulkModulus;

    m_biotCoefficient[k] = 1 - bulkModulus / m_grainBulkModulus;
  }

  GEOSX_HOST_DEVICE
  void updateThermalExpansionCoefficient( localIndex const k,
                                          real64 const thermalExpansionCoefficient ) const
  {
    m_thermalExpansionCoefficient[k] = thermalExpansionCoefficient;
  }

  GEOSX_HOST_DEVICE
  void updateTotalMeanStressIncrement( localIndex const k,
                                       localIndex const q,
                                       real64 const & totalMeanStressIncrement ) const
  {
    m_meanStressIncrement[k][q] = totalMeanStressIncrement;
  }

protected:
  arrayView1d< real64 > m_biotCoefficient;

  arrayView1d< real64 > m_thermalExpansionCoefficient;

  arrayView2d< real64 > m_meanStressIncrement;

  arrayView1d< real64 > m_bulkModulus;

  real64 m_grainBulkModulus;
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

    static constexpr char const *thermalExpansionCoefficientString() { return "thermalExpansionCoefficient"; }

    static constexpr char const *meanStressIncrementString() { return "meanStressIncrement"; }

    static constexpr char const *solidBulkModulusString() { return "solidBulkModulus"; }
  } viewKeys;

  virtual void initializeState() const override final;

  virtual arrayView1d< real64 const > const getBiotCoefficient() const override final
  {
    return m_biotCoefficient.toViewConst();
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
                          m_meanStressIncrement,
                          m_bulkModulus,
                          m_grainBulkModulus );
  }

protected:
  virtual void postProcessInput() override;

  array1d< real64 > m_biotCoefficient;

  array1d< real64 > m_thermalExpansionCoefficient;

  array2d< real64 > m_meanStressIncrement;

  array1d< real64 > m_bulkModulus;

  real64 m_grainBulkModulus;
};

}   /* namespace constitutive */

} /* namespace geosx */

#endif //GEOSX_CONSTITUTIVE_POROSITY_BIOTPOROSITY_HPP_
