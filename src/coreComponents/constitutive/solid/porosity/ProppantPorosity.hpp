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
 * @file PressurePorosity.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_POROSITY_PROPPANTPOROSITY_HPP_
#define GEOSX_CONSTITUTIVE_POROSITY_PROPPANTPOROSITY_HPP_

#include "PorosityBase.hpp"

namespace geosx
{
namespace constitutive
{

class ProppantPorosityUpdates : public PorosityBaseUpdates
{
public:

  ProppantPorosityUpdates( arrayView2d< real64 > const & newPorosity,
                           arrayView2d< real64 > const & oldPorosity,
                           arrayView2d< real64 > const & dPorosity_dPressure,
                           arrayView1d< real64 > const & referencePorosity,
                           real64 const & maxProppantConcentration ):
    PorosityBaseUpdates( newPorosity,
                         oldPorosity,
                         dPorosity_dPressure,
                         referencePorosity ),
    m_maxProppantConcentration( maxProppantConcentration )
  {}

  GEOSX_HOST_DEVICE
  void computePorosity( real64 const & proppantPackVolumeFraction,
                        real64 & porosity ) const
  {
    porosity = 1.0 - m_maxProppantConcentration * proppantPackVolumeFraction;
  }

  GEOSX_HOST_DEVICE
  void updateFromProppantVolumeFraction( localIndex const k,
                                         localIndex const q,
                                         real64 const & proppantPackVolumeFraction ) const
  {
    computePorosity( proppantPackVolumeFraction,
                     m_newPorosity[k][q] );
  }

private:

  real64 m_maxProppantConcentration;

};


class ProppantPorosity : public PorosityBase
{
public:
  ProppantPorosity( string const & name, Group * const parent );

  virtual ~ProppantPorosity() override;

  std::unique_ptr< ConstitutiveBase > deliverClone( string const & name,
                                                    Group * const parent ) const override;

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  static string catalogName() { return "ProppantPorosity"; }

  virtual string getCatalogName() const override { return catalogName(); }

  struct viewKeyStruct : public PorosityBase::viewKeyStruct
  {
    static constexpr char const * maxProppantConcentrationString() { return "maxProppantConcentration"; }
  } viewKeys;

  using KernelWrapper = ProppantPorosityUpdates;

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
                          m_maxProppantConcentration );
  }


private:
  virtual void postProcessInput() override;

  real64 m_maxProppantConcentration;


};


}/* namespace constitutive */

} /* namespace geosx */


#endif //GEOSX_CONSTITUTIVE_POROSITY_PROPPANTPOROSITY_HPP_
