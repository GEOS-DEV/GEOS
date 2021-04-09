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
 * @file CompressibleRock.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_POROSITY_COMPRESSIBLEROCK_HPP_
#define GEOSX_CONSTITUTIVE_POROSITY_COMPRESSIBLEROCK_HPP_

#include "constitutive/ExponentialRelation.hpp"
#include "RockBase.hpp"

namespace geosx
{
namespace constitutive
{

class CompressibleRockUpdates : RockBaseUpdates
{
public:

  CompressibleRockUpdates( arrayView2d< real64 > const & newPorosity,
                           arrayView2d< real64 > const & oldPorosity,
                           arrayView2d< real64 > const & dPorosity_dPressure,
                           real64 const & compressibility,
                           real64 const & grainBulkModulus,
                           real64 const & grainDensity,
                           arrayView1d< real64 > const & referencePorosity,
                           real64 const & referencePressure ):
   RockBaseUpdates( newPorosity,
                    oldPorosity,
                    dPorosity_dPressure,
                    compressibility,
                    grainBulkModulus,
                    grainDensity ),
    m_referencePorosity( referencePorosity ),
    m_referencePressure( referencePressure )
  {}

  /// Default copy constructor
  CompressibleRockUpdates( CompressibleRockUpdates const & ) = default;

  /// Default move constructor
  CompressibleRockUpdates( CompressibleRockUpdates && ) = default;

  /// Deleted copy assignment operator
  CompressibleRockUpdates & operator=( CompressibleRockUpdates const & ) = delete;

  /// Deleted move assignment operator
  CompressibleRockUpdates & operator=( CompressibleRockUpdates && ) = delete;

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void compute( real64 const & pressure,
                real64 & porosity,
                real64 & dPorosity_dPressure,
                real64 const & referencePorosity ) const
  {

    porosity            =  referencePorosity * exp( m_compressibility * (pressure - m_referencePressure) );
    dPorosity_dPressure =  m_compressibility * porosity;
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void update( localIndex const k,
               localIndex const q,
               real64 const & pressure ) const
  {
    compute( pressure,
             m_newPorosity[k][q],
             m_dPorosity_dPressure[k][q],
             m_referencePorosity[k] );
  }

private:

  arrayView1d< real64 > m_referencePorosity;

  real64 m_referencePressure;

};


class CompressibleRock : public RockBase
{
public:
  CompressibleRock( string const & name, Group * const parent );

  virtual ~CompressibleRock() override;

  std::unique_ptr< ConstitutiveBase > deliverClone( string const & name,
                                                    Group * const parent ) const override;

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  static string catalogName() { return "CompressibleRock"; }

  virtual string getCatalogName() const override { return catalogName(); }

  struct viewKeyStruct : public RockBase::viewKeyStruct
  {
    static constexpr char const * referencePressureString() { return "referencePressure"; }
    static constexpr char const * referencePorosityString() { return "referencePorosity"; }
    static constexpr char const * defaultRefererencePorosityString() { return "defaultReferencePorosity"; }
  } viewKeys;

  using KernelWrapper = CompressibleRockUpdates;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper()
  {
    return KernelWrapper( m_newPorosity,
                          m_oldPorosity,
                          m_dPorosity_dPressure,
                          m_compressibility,
                          m_grainBulkModulus,
                          m_grainDensity,
                          m_referencePorosity,
                          m_referencePressure );
  }


protected:
  virtual void postProcessInput() override;

  real64 m_referencePressure;
  array1d< real64 > m_referencePorosity;
  real64 m_defaultReferencePorosity;

};


}/* namespace constitutive */

} /* namespace geosx */


#endif //GEOSX_CONSTITUTIVE_POROSITY_COMPRESSIBLEROCK_HPP_
