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
                       arrayView2d< real64 > const & initialPorosity,
                       arrayView1d< real64 > const & referencePorosity,
                       arrayView1d< real64 > const & biotCoefficient,
                       real64 const & grainBulkModulus ): PorosityBaseUpdates( newPorosity,
                                                                               porosity_n,
                                                                               dPorosity_dPressure,
                                                                               initialPorosity,
                                                                               referencePorosity ),
    m_biotCoefficient( biotCoefficient ),
    m_grainBulkModulus( grainBulkModulus )
  {}

  GEOSX_HOST_DEVICE
  real64 getBiotCoefficient( localIndex const k ) const { return m_biotCoefficient[k]; }

  GEOSX_HOST_DEVICE
  real64 getGrainBulkModulus() const { return m_grainBulkModulus; }

  GEOSX_HOST_DEVICE
  real64 dGrainDensity_dPressure() const { return 1.0 / m_grainBulkModulus; }

  GEOSX_HOST_DEVICE
  void updateFromPressureAndStrain( localIndex const k,
                                    localIndex const q,
                                    real64 const & deltaPressure,
                                    real64 const (&strainIncrement)[6],
                                    real64 & dPorosity_dPressure,
                                    real64 & dPorosity_dVolStrain ) const
  {
    real64 const biotSkeletonModulusInverse = (m_biotCoefficient[k] - m_referencePorosity[k]) / m_grainBulkModulus;

    real64 const porosity = m_porosity_n[k][q] +
                            +m_biotCoefficient[k] * LvArray::tensorOps::symTrace< 3 >( strainIncrement ) + biotSkeletonModulusInverse * deltaPressure;

    dPorosity_dPressure = biotSkeletonModulusInverse;

    dPorosity_dVolStrain = m_biotCoefficient[k];

    savePorosity( k, q, porosity, biotSkeletonModulusInverse );
  }

  GEOSX_HOST_DEVICE
  void updateBiotCoefficient( localIndex const k,
                              real64 const bulkModulus ) const
  {
    m_biotCoefficient[k] = 1 - bulkModulus / m_grainBulkModulus;
  }

protected:
  arrayView1d< real64 > m_biotCoefficient;

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
  } viewKeys;

  virtual void initializeState() const override final;

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
                          m_initialPorosity,
                          m_referencePorosity,
                          m_biotCoefficient,
                          m_grainBulkModulus );
  }

protected:
  virtual void postProcessInput() override;

  array1d< real64 > m_biotCoefficient;

  real64 m_grainBulkModulus;
};

}   /* namespace constitutive */

} /* namespace geosx */

#endif //GEOSX_CONSTITUTIVE_POROSITY_BIOTPOROSITY_HPP_
