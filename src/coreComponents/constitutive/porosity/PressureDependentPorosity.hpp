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
 * @file PorosityBase.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_PERMEABILITY_POROSITYBASE_HPP_
#define GEOSX_CONSTITUTIVE_PERMEABILITY_POROSITYBASE_HPP_

#include "PorosityBase.hpp"

#include "constitutive/ExponentialRelation.hpp"

namespace geosx
{
namespace constitutive
{

template< ExponentApproximationType POR_EAT >
class PressureDependentPorosityUpdate
{
public:

  using PorRelationType = ExponentialRelation< real64, POR_EAT >;

  /**
   * @brief Get number of elements in this wrapper.
   * @return number of elements
   */
  GEOSX_HOST_DEVICE
  localIndex numElems() const { return m_permeability.size( 0 ); }

  /**
   * @brief Get number of gauss points per element.
   * @return number of gauss points per element
   */
  GEOSX_HOST_DEVICE
  localIndex numGauss() const { return m_permeability.size( 1 ); }

protected:

  PressureDependentPorosityUpdate( PorRelationType const & porRelation,
                                   arrayView3d< real64 > const & porosity,
                                   arrayView3d< real64 > const & dPorosity_dPressure )
    : m_porosity( porosity ),
      m_dPorosity_dPressure( dPorosity_dPressure ),
      m_porRelation( porRelation )
  {}

  /// Default copy constructor
  PressureDependentPorosityUpdate( PermeabilityBaseUpdate const & ) = default;

  /// Default move constructor
  PressureDependentPorosityUpdate( PermeabilityBaseUpdate && ) = default;

  /// Deleted copy assignment operator
  PressureDependentPorosityUpdate & operator=( PermeabilityBaseUpdate const & ) = delete;

  /// Deleted move assignment operator
  PressureDependentPorosityUpdate & operator=( PermeabilityBaseUpdate && ) = delete;

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void compute( real64 const & pressure,
                arraySlice1d< real64 > const & porosity,
                arraySlice1d< real64 > const & dPorosity_dPressure ) const;

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void update( localIndex const k,
               localIndex const q,
               real64 const & pressure ) const
  {
    compute( pressure,
             m_porosity[k][q],
             m_dPorosity_dPressure[k][q] );
  }

private:
  arrayView1d< real64 > m_porosity;

  arrayView1d< real64 > m_dPorosity_dPressure;

  PorRelationType m_porRelation;
};


class PressureDependentPorosity : public PorosityBase
{
public:
  PressureDependentPorosity( string const & name, Group * const parent );

  virtual ~PressureDependentPorosity() override;

  std::unique_ptr< ConstitutiveBase > deliverClone( string const & name,
                                                    Group * const parent ) const override;

  virtual void allocateConstitutiveData( dataRepository::Group * const parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  static string catalogName() { return "PressureDependentPorosity"; }

  virtual string getCatalogName() const override { return catalogName(); }

  struct viewKeyStruct : public ConstitutiveBase::viewKeyStruct
  {
  } viewKeys;

  using KernelWrapper = PressureDependentPorosityUpdate< ExponentApproximationType::Linear >;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper()
  {
      return KernelWrapper( KernelWrapper::PorRelationType( m_referencePressure,
                                                            m_referencePorosity,
                                                            m_compressibility ),
                            m_porosity,
                            m_dPorosity_dPressure );
  }


protected:
  virtual void postProcessInput() override;

  real64 m_referencePressure;
  real64 m_referencePorosity;
  real64 m_compressibility;


};


GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void PressureDependentPorosityUpdate::compute( real64 const & pressure,
                                         arraySlice1d< real64 > const & porosity,
                                         arraySlice1d< real64 > const & dPorosity_dPressure ) const
{
  m_porRelation.compute( pressure, porosity, dPorosity_dPressure );
}




}/* namespace constitutive */

} /* namespace geosx */


#endif //GEOSX_CONSTITUTIVE_PERMEABILITY_POROSITYBASE_HPP_
