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
 * @file RockBase.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_POROSITY_ROCKBASE_HPP_
#define GEOSX_CONSTITUTIVE_POROSITY_ROCKBASE_HPP_

#include "constitutive/ConstitutiveBase.hpp"

namespace geosx
{
namespace constitutive
{

class RockBaseUpdates
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

  RockBaseUpdates( arrayView2d< real64 > const & newPorosity,
                   arrayView2d< real64 > const & oldPorosity,
                   arrayView2d< real64 > const & dPorosity_dPressure,
                   real64 const & compressibility,
                   real64 const & grainBulkModulus,
                   real64 const & grainDensity ):
    m_newPorosity( newPorosity ),
    m_oldPorosity( oldPorosity ),
    m_dPorosity_dPressure( dPorosity_dPressure ),
    m_compressibility ( compressibility ),
    m_grainBulkModulus( grainBulkModulus ),
    m_grainDensity( grainDensity )
  {}

  /// Default copy constructor
  RockBaseUpdates( RockBaseUpdates const & ) = default;

  /// Default move constructor
  RockBaseUpdates( RockBaseUpdates && ) = default;

  /// Deleted copy assignment operator
  RockBaseUpdates & operator=( RockBaseUpdates const & ) = delete;

  /// Deleted move assignment operator
  RockBaseUpdates & operator=( RockBaseUpdates && ) = delete;


protected:
  arrayView2d< real64 > m_newPorosity;

  arrayView2d< real64 > m_oldPorosity;

  arrayView2d< real64 > m_dPorosity_dPressure;

  real64 m_compressibility;

  real64 m_grainBulkModulus;

  real64 m_grainDensity;
};


class RockBase : public ConstitutiveBase
{
public:
  RockBase( string const & name, Group * const parent );

  virtual ~RockBase() override;

  std::unique_ptr< ConstitutiveBase > deliverClone( string const & name,
                                                    Group * const parent ) const override;

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  static string catalogName() { return "RockBase"; }

  virtual string getCatalogName() const override { return catalogName(); }

  struct viewKeyStruct : public ConstitutiveBase::viewKeyStruct
  {
    static constexpr char const * newPorosityString() { return "porosity"; }
    static constexpr char const * oldPorosityString() { return "oldPorosity"; }
    static constexpr char const * dPorosity_dPressureString() { return "dPorosity_dPressure"; }
    static constexpr char const * compressibilityString() { return "compressibility"; }
    static constexpr char const * grainBulkModulusString() { return "grainBulkModulus"; }
    static constexpr char const * grainDensityString() { return "grainDensity"; }
  } viewKeys;

  arrayView2d< real64 const > const  getPorosity() const { return m_newPorosity; }

  arrayView2d< real64 const > const  getOldPorosity() const { return m_oldPorosity; }

  arrayView2d< real64 > const getOldPorosity() { return m_oldPorosity; }

  arrayView2d< real64 const > const  dPorosity_dPressure() const { return m_dPorosity_dPressure; }

protected:
  virtual void postProcessInput() override;

  array2d< real64 > m_newPorosity;

  array2d< real64 > m_oldPorosity;

  array2d< real64 > m_dPorosity_dPressure;

  real64 m_compressibility;

  real64 m_grainBulkModulus;

  real64 m_grainDensity;
};

}/* namespace constitutive */

} /* namespace geosx */


#endif //GEOSX_CONSTITUTIVE_POROSITY_ROCKBASE_HPP_
