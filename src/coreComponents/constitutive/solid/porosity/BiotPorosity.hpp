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
                       arrayView2d< real64 > const & meanStressIncrement,
                       arrayView1d< real64 > const & bulkModulus,
                       arrayView1d< real64 > const & pressure,
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
    m_grainBulkModulus( grainBulkModulus ),
    m_pressure(pressure)
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
    real64 const porosityThermalExpansion = 3 * m_thermalExpansionCoefficient[k] * m_biotCoefficient[k];

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
  void computePorosity( real64 const & pressure,
                        real64 const & temperature,
                        real64 const & porosity_n,
                        real64 const & referencePorosity,
                        real64 & porosity,
                        real64 & dPorosity_dPressure,
                        real64 & dPorosity_dTemperature,
                        real64 & biotCoefficient,
                        real64 const & thermalExpansionCoefficient,
                        real64 const & meanStressIncrement,
                        real64 const & bulkModulus ) const
  {

    biotCoefficient= 1.0;

    real64 const biotSkeletonModulusInverse = (biotCoefficient - referencePorosity) / m_grainBulkModulus;
    real64 const porosityThermalExpansion = 3 * thermalExpansionCoefficient * biotCoefficient;

    GEOS_UNUSED_VAR(biotSkeletonModulusInverse, porosityThermalExpansion, pressure, temperature, meanStressIncrement, bulkModulus, porosity_n);

    porosity = porosity_n + biotSkeletonModulusInverse * pressure
               - porosityThermalExpansion * temperature
               + biotCoefficient * (meanStressIncrement + biotCoefficient*pressure + 3*thermalExpansionCoefficient*bulkModulus*temperature) / bulkModulus;


    //std::cout << porosity << "\t" << porosity_n << "\t" << pressure << "\t" << meanStressIncrement << std::endl;

    //porosity = porosity_n + (biotCoefficient - porosity_n) * pressure / bulkModulus + 3.0 * thermalExpansionCoefficient * temperature + biotCoefficient * meanStressIncrement / bulkModulus;

    dPorosity_dPressure = biotSkeletonModulusInverse + biotCoefficient * biotCoefficient / bulkModulus;
    //dPorosity_dPressure = (biotCoefficient - porosity_n) / bulkModulus;
    
    
    dPorosity_dTemperature = -porosityThermalExpansion;
  }

  GEOS_HOST_DEVICE
  virtual void updateFromPressureAndTemperature( localIndex const k,
                                                 localIndex const q,
                                                 real64 const & pressure,
                                                 real64 const & pressure_n,
                                                 real64 const & temperature,
                                                 real64 const & temperature_n ) const override final
  {
    real64 const deltaPressure    = pressure - m_pressure[k];
    real64 const deltaTemperature = temperature - temperature_n;

    GEOS_UNUSED_VAR(pressure_n);

    //std::cout << k << std::endl;

    computePorosity( deltaPressure,
                     deltaTemperature,
                     m_porosity_n[k][q],
                     m_referencePorosity[k],
                     m_newPorosity[k][q],
                     m_dPorosity_dPressure[k][q],
                     m_dPorosity_dTemperature[k][q],
                     m_biotCoefficient[k],
                     m_thermalExpansionCoefficient[k],
                     m_meanStressIncrement[k][q],
                     m_bulkModulus[k] );
  }

  GEOS_HOST_DEVICE
  void updateBiotCoefficient( localIndex const k,
                              real64 const bulkModulus ) const
  {
    m_bulkModulus[k] = bulkModulus;

    m_biotCoefficient[k] = 1 - bulkModulus / m_grainBulkModulus;
  }

  GEOS_HOST_DEVICE
  void updateTotalMeanStressIncrement( localIndex const k,
                                       localIndex const q,
                                       real64 const & totalMeanStressIncrement ) const
  {
    m_meanStressIncrement[k][q] = totalMeanStressIncrement;
    //

    //if (q == 0) {m_meanStressIncrement[k][0] = 0.0;}
    //m_meanStressIncrement[k][0] += totalMeanStressIncrement / 8.0;
    //GEOS_UNUSED_VAR(q);

    //std::cout << "k = " << k << "\t q = " << q << "\t totalMeanStressIncrement = " << totalMeanStressIncrement << std::endl;
  }

  GEOS_HOST_DEVICE
  void updateSequentialPressure(real64 const & pressure, localIndex const k) const {
    m_pressure[k] = pressure;
  }

protected:
  arrayView1d< real64 > m_biotCoefficient;

  arrayView1d< real64 > m_thermalExpansionCoefficient;

  arrayView2d< real64 > m_meanStressIncrement;

  arrayView1d< real64 > m_bulkModulus;

  real64 m_grainBulkModulus;

  arrayView1d <real64> m_pressure;
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

    static constexpr char const *meanStressIncrementString() { return "meanStressIncrement"; }

    static constexpr char const *solidBulkModulusString() { return "solidBulkModulus"; }

    static constexpr char const *pressureString() { return "sequentialPressure"; }

    static constexpr char const *defaultThermalExpansionCoefficientString() { return "defaultThermalExpansionCoefficient"; }
  } viewKeys;

  virtual void initializeState() const override final;

  virtual void saveConvergedState() const override final                                                                                                                                 
  {                                                                                                                                                                                      
       PorosityBase::saveConvergedState();                                                                                                                                                  
          m_meanStressIncrement.zero();                                                                                                                                                        
  }  

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
                          m_pressure,
                          m_grainBulkModulus );
  }

protected:
  virtual void postProcessInput() override;

  array1d< real64 > m_biotCoefficient;

  array1d< real64 > m_thermalExpansionCoefficient;

  array1d< real64 > m_bulkModulus;

  array2d< real64 > m_meanStressIncrement;

  real64 m_grainBulkModulus;

  real64 m_defaultThermalExpansionCoefficient;

  array1d< real64 > m_pressure;   
};

}   /* namespace constitutive */

} /* namespace geos */

#endif //GEOS_CONSTITUTIVE_POROSITY_BIOTPOROSITY_HPP_
