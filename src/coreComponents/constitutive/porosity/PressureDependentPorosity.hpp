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
 * @file PressureDependentPorosity.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_POROSITY_PRESSUREDEPENDENTPOROSITY_HPP_
#define GEOSX_CONSTITUTIVE_POROSITY_PRESSUREDEPENDENTPOROSITY_HPP_

#include "PorosityBase.hpp"

#include "constitutive/ExponentialRelation.hpp"

namespace geosx
{
namespace constitutive
{

class PressureDependentPorosityUpdate
{
public:

  /**
   * @brief Get number of elements in this wrapper.
   * @return number of elements
   */
  GEOSX_HOST_DEVICE
  localIndex numElems() const { return m_porosity.size( 0 ); }

  /**
   * @brief Get number of gauss points per element.
   * @return number of gauss points per element
   */
  GEOSX_HOST_DEVICE
  localIndex numGauss() const { return m_porosity.size( 1 ); }

  PressureDependentPorosityUpdate( arrayView2d< real64 > const & porosity,
                                   arrayView2d< real64 > const & dPorosity_dPressure,
                                   arrayView1d< real64 > const & referencePorosity,
                                   real64 const & referencePressure,
                                   real64 const & compressibility )
    : m_porosity( porosity ),
    m_dPorosity_dPressure( dPorosity_dPressure ),
    m_referencePorosity( referencePorosity ),
    m_referencePressure( referencePressure ),
    m_compressibility ( compressibility )
  {}

  /// Default copy constructor
  PressureDependentPorosityUpdate( PressureDependentPorosityUpdate const & ) = default;

  /// Default move constructor
  PressureDependentPorosityUpdate( PressureDependentPorosityUpdate && ) = default;

  /// Deleted copy assignment operator
  PressureDependentPorosityUpdate & operator=( PressureDependentPorosityUpdate const & ) = delete;

  /// Deleted move assignment operator
  PressureDependentPorosityUpdate & operator=( PressureDependentPorosityUpdate && ) = delete;

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
             m_porosity[k][q],
             m_dPorosity_dPressure[k][q],
             m_referencePorosity[k] );
  }

private:
  arrayView2d< real64 > m_porosity;

  arrayView2d< real64 > m_porosityOld;

  arrayView2d< real64 > m_dPorosity_dPressure;

  arrayView1d< real64 > m_referencePorosity;

  real64 m_referencePressure;

  real64 m_compressibility;
};


class PressureDependentPorosity : public PorosityBase
{
public:
  PressureDependentPorosity( string const & name, Group * const parent );

  virtual ~PressureDependentPorosity() override;

  std::unique_ptr< ConstitutiveBase > deliverClone( string const & name,
                                                    Group * const parent ) const override;

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  static string catalogName() { return "PressureDependentPorosity"; }

  virtual string getCatalogName() const override { return catalogName(); }

  struct viewKeyStruct : public PorosityBase::viewKeyStruct
  {
    static constexpr char const * compressibilityString() { return "compressibility"; }
    static constexpr char const * referencePressureString() { return "referencePressure"; }
    static constexpr char const * referencePorosityString() { return "referencePorosity"; }
    static constexpr char const * defaultRefererencePorosityString() { return "defaultReferencePorosity"; }
  } viewKeys;

  using KernelWrapper = PressureDependentPorosityUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper()
  {
    return KernelWrapper( m_porosity,
                          m_dPorosity_dPressure,
                          m_referencePorosity,
                          m_referencePressure,
                          m_compressibility );
  }


protected:
  virtual void postProcessInput() override;

  real64 m_referencePressure;
  real64 m_compressibility;
  array1d< real64 > m_referencePorosity;
  real64 m_defaultReferencePorosity;

};


}/* namespace constitutive */

} /* namespace geosx */


#endif //GEOSX_CONSTITUTIVE_POROSITY_PRESSUREDEPENDENTPOROSITY_HPP_
