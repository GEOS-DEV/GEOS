/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
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
#include "constitutive/solid/PropertyConversions.hpp"

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
                       arrayView2d< real64 > const & oldPorosity,
                       arrayView2d< real64 > const & dPorosity_dPressure,
                       arrayView1d< real64 > const & referencePorosity,
                       arrayView1d< real64 > const & biotCoefficient,
                       real64 const & grainBulkModulus,
                       real64 const & grainShearModulus ):
    PorosityBaseUpdates( newPorosity,
                         oldPorosity,
                         dPorosity_dPressure,
                         referencePorosity ),
    m_biotCoefficient( biotCoefficient ),
    m_grainBulkModulus( grainBulkModulus ),
    m_grainShearModulus( grainShearModulus )
  {}


  GEOSX_HOST_DEVICE
  void updatePorosity( localIndex const k,
                       localIndex const q,
                       real64 const & deltaPressure,
                       real64 const ( &strainIncrement )[6],
                       real64 & dPorosity_dPressure,
                       real64 & dPorosity_dVolStrain,
                       real64 & dTotalStress_dPressure ) const
  {
    real64 const biotSkeletonModulusInverse = ( m_biotCoefficient[k] - m_referencePorosity[k] ) / m_grainBulkModulus;

    real64 const porosity = m_oldPorosity[k][q] +
                            +m_biotCoefficient[k] * LvArray::tensorOps::symTrace< 3 >( strainIncrement )
                            + biotSkeletonModulusInverse * deltaPressure;

    dPorosity_dPressure = biotSkeletonModulusInverse;

    dPorosity_dVolStrain = m_biotCoefficient[k];

    savePorosity( k, q, porosity, biotSkeletonModulusInverse );

    dTotalStress_dPressure = m_biotCoefficient[k];
  }


  GEOSX_HOST_DEVICE
  void updateBiotCoefficient( localIndex const k,
                              real64 const bulkModulus ) const
  {
    m_biotCoefficient[k] = 1.0 - bulkModulus / m_grainBulkModulus;
  }

  /**
   * @brief update transversely isotropic Biot's tensor.
   * @param c11 (1,1) component of the stiffness matrix
   * @param c12 (1,2) component of the stiffness matrix
   * @param c13 (1,3) component of the stiffness matrix
   * @param c33 (3,3) component of the stiffness matrix
   */
  GEOSX_HOST_DEVICE
  void updateBiotCoefficient( localIndex const k,
                              real64 const c11,
                              real64 const c12,
                              real64 const c13,
                              real64 const c33 ) const
  {
    // Components of grain compliance matrix (isotropic grain is assumed)
    real64 const grainPoissonRatio = conversions::BulkModAndShearMod::
                                       toPoissonRatio( m_grainBulkModulus, m_grainShearModulus );

    real64 const grainYoungModulus = conversions::BulkModAndShearMod::
                                       toYoungsMod( m_grainBulkModulus, m_grainShearModulus );

    real64 const s11Grain = 1.0 / grainYoungModulus;
    real64 const s12Grain = - grainPoissonRatio / grainYoungModulus;
    real64 const s13Grain = s12Grain;
    real64 const s33Grain = s11Grain;

    // Components of Biot's matrix
    real64 const b11 = 1.0 - ( s11Grain + s12Grain ) * ( c11 + c12 )
                       - s13Grain * ( c11 + c12 + 2.0 * c13 )
                       - s33Grain * c13;

    real64 const b33 = 1.0 - 2.0 * s11Grain * c13
                       - 2.0 * s12Grain * c13
                       - 2.0 * s13Grain * ( c13 + c33 )
                       - s33Grain * c33;

    real64 const b22 = b11;

    // Biot's coefficient is approximated by the average of the Biot's tensor diagonal.
    m_biotCoefficient[k] = ( b11 + b22 + b33 ) / 3.0;
  }

protected:

  arrayView1d< real64 > m_biotCoefficient;

  real64 m_grainBulkModulus;

  real64 m_grainShearModulus;
};


class BiotPorosity : public PorosityBase
{
public:
  BiotPorosity( string const & name, Group * const parent );

  virtual ~BiotPorosity() override;

  std::unique_ptr< ConstitutiveBase > deliverClone( string const & name,
                                                    Group * const parent ) const override;

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  static string catalogName() { return "BiotPorosity"; }

  virtual string getCatalogName() const override { return catalogName(); }

  struct viewKeyStruct : public PorosityBase::viewKeyStruct
  {
    static constexpr char const * biotCoefficientString() { return "biotCoefficient"; }
    static constexpr char const * grainBulkModulusString() { return "grainBulkModulus"; }
    static constexpr char const * grainShearModulusString() { return "grainShearModulus"; }
  } viewKeys;

  using KernelWrapper = BiotPorosityUpdates;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelUpdates() const
  {
    return KernelWrapper( m_newPorosity,
                          m_oldPorosity,
                          m_dPorosity_dPressure,
                          m_referencePorosity,
                          m_biotCoefficient,
                          m_grainBulkModulus,
                          m_grainShearModulus );
  }


protected:
  virtual void postProcessInput() override;

  array1d< real64 > m_biotCoefficient;

  real64 m_grainBulkModulus; ///< Bulk modulus of the solid grain

  real64 m_grainShearModulus; ///< Shear modulus of the solid grain
};

}/* namespace constitutive */

} /* namespace geosx */


#endif //GEOSX_CONSTITUTIVE_POROSITY_BIOTPOROSITY_HPP_
