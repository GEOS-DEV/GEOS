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
 * @file PressurePorosityModel.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_POROSITY_PRESSUREPOROSITYMODEL_HPP_
#define GEOSX_CONSTITUTIVE_POROSITY_PRESSUREPOROSITYMODEL_HPP_

#include "constitutive/ExponentialRelation.hpp"
#include "PorosityBase.hpp"

namespace geosx
{
namespace constitutive
{

class PressurePorosityModelUpdates : public PorosityBaseUpdates
{
public:

  PressurePorosityModelUpdates( arrayView2d< real64 > const & newPorosity,
                           arrayView2d< real64 > const & oldPorosity,
                           arrayView2d< real64 > const & dPorosity_dPressure,
                           arrayView1d< real64 > const & referencePorosity,
                           real64 const & referencePressure,
                           real64 const & compressibility ):
    PorosityBaseUpdates( newPorosity,
                         oldPorosity,
                         dPorosity_dPressure,
                         referencePorosity ),
    m_referencePressure( referencePressure ),
    m_compressibility( compressibility )
  {}

  /// Default copy constructor
  PressurePorosityModelUpdates( PressurePorosityModelUpdates const & ) = default;

  /// Default move constructor
  PressurePorosityModelUpdates( PressurePorosityModelUpdates && ) = default;

  /// Deleted copy assignment operator
  PressurePorosityModelUpdates & operator=( PressurePorosityModelUpdates const & ) = delete;

  /// Deleted move assignment operator
  PressurePorosityModelUpdates & operator=( PressurePorosityModelUpdates && ) = delete;

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void computePorosity( real64 const & pressure,
                        real64 & porosity,
                        real64 & dPorosity_dPressure,
                        real64 const & referencePorosity ) const
  {

    porosity            =  referencePorosity * exp( m_compressibility * (pressure - m_referencePressure) );
    dPorosity_dPressure =  m_compressibility * porosity;
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  virtual void updatePorosity( localIndex const k,
                               localIndex const q,
                               real64 const & pressure ) const override
  {
    computePorosity( pressure,
                     m_newPorosity[k][q],
                     m_dPorosity_dPressure[k][q],
                     m_referencePorosity[k] );
  }

private:

  real64 m_referencePressure;

  real64 m_compressibility;

};


class PressurePorosityModel : public PorosityBase
{
public:
  PressurePorosityModel( string const & name, Group * const parent );

  virtual ~PressurePorosityModel() override;

  std::unique_ptr< ConstitutiveBase > deliverClone( string const & name,
                                                    Group * const parent ) const override;

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  static string catalogName() { return "PressurePorosityModel"; }

  virtual string getCatalogName() const override { return catalogName(); }

  struct viewKeyStruct : public RockBase::viewKeyStruct
  {
    static constexpr char const * referencePressureString() { return "referencePressure"; }
    static constexpr char const * compressibilityString() { return "compressibility"; }
  } viewKeys;

  using KernelWrapper = PressurePorosityModelUpdates;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelUpdates()
  {
    return KernelWrapper( m_newPorosity,
                          m_oldPorosity,
                          m_dPorosity_dPressure,
                          m_referencePorosity,
                          m_referencePressure,
                          m_compressibility );
  }


protected:
  virtual void postProcessInput() override;

  real64 m_referencePressure;

  real64 m_compressibility;

};


}/* namespace constitutive */

} /* namespace geosx */


#endif //GEOSX_CONSTITUTIVE_POROSITY_PRESSUREPOROSITYMODEL_HPP_
