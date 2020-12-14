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
 * @file Damage.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_DAMAGE_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_DAMAGE_HPP_

#include "constitutive/solid/SolidBase.hpp"

namespace geosx
{
namespace constitutive
{

// DAMAGE MODEL UPDATES
//
// NOTE: This model uses the m_newStress array to represent the stress in an
//       elastic, "undamaged" configuration.  We then scale the results
//       by the damage factor whenever the true stress is requested through a getter.
//
// TODO: This approach is probably error prone---e.g. when stress data is
//       accessed directly through an array view.  We should refactor so
//       the true state is saved, similar to the plasticity models

template< typename UPDATE_BASE >
class DamageUpdates : public UPDATE_BASE
{
public:
  template< typename ... PARAMS >
  DamageUpdates( arrayView2d< real64 > const & inputDamage,
                 arrayView2d< real64 > const & inputStrainEnergyDensity,
                 PARAMS && ... baseParams ):
    UPDATE_BASE( std::forward< PARAMS >( baseParams )... ),
    m_damage( inputDamage ),
    m_strainEnergyDensity( inputStrainEnergyDensity )
  {}

  using DiscretizationOps = typename UPDATE_BASE::DiscretizationOps;

  using UPDATE_BASE::smallStrainNoStateUpdate;
  using UPDATE_BASE::smallStrainUpdate;
  using UPDATE_BASE::smallStrainNoStateUpdate_StressOnly;
  using UPDATE_BASE::smallStrainUpdate_StressOnly;

  /**
   * Compute current damage scaling
   * @param k Element index
   * @param q Quadrature point index
   * @return Current damage scaling
   */
  GEOSX_HOST_DEVICE
  real64 damageFactor( localIndex const k,
                       localIndex const q ) const
  {
    return ( 1.0 - m_damage[k][q] )*( 1.0 - m_damage[k][q] );
  }


  GEOSX_HOST_DEVICE
  virtual void smallStrainUpdate( localIndex const k,
                                  localIndex const q,
                                  real64 const ( &strainIncrement )[6],
                                  real64 ( & stress )[6],
                                  DiscretizationOps & stiffness ) const final
  {
    UPDATE_BASE::smallStrainUpdate( k, q, strainIncrement, stress, stiffness );
    real64 factor = damageFactor( k, q );
    LvArray::tensorOps::scale< 6 >( stress, factor );
    stiffness.scaleParams( factor );
  }

  // TODO: The code below assumes the strain energy density will never be
  //       evaluated in a non-converged / garbage configuration.
  GEOSX_HOST_DEVICE
  virtual real64 getStrainEnergyDensity( localIndex const k,
                                         localIndex const q ) const override
  {
    real64 sed = SolidBaseUpdates::getStrainEnergyDensity( k, q );
    if( sed > m_strainEnergyDensity[k][q] )
    {
      m_strainEnergyDensity[k][q] = sed;
    }
    return m_strainEnergyDensity[k][q];
  }

  arrayView2d< real64 > const m_damage;
  arrayView2d< real64 > const m_strainEnergyDensity;
};



class DamageBase : public SolidBase
{};

template< typename BASE >
class Damage : public BASE
{
public:

  /// @typedef Alias for LinearElasticIsotropicUpdates
  using KernelWrapper = DamageUpdates< typename BASE::KernelWrapper >;

  Damage( string const & name, dataRepository::Group * const parent );
  virtual ~Damage() override;

  static std::string CatalogName() { return string( "Damage" ) + BASE::m_catalogNameString; }
  virtual string getCatalogName() const override { return CatalogName(); }

  virtual void PostProcessInput() override;

  virtual void allocateConstitutiveData( dataRepository::Group * const parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;


  KernelWrapper createKernelUpdates()
  {
    return BASE::template createDerivedKernelUpdates< KernelWrapper >( m_damage.toView(),
                                                                       m_strainEnergyDensity.toView() );
  }

  struct viewKeyStruct : public BASE::viewKeyStruct
  {
    static constexpr auto damageString =  "damage";
    static constexpr auto strainEnergyDensityString =  "strainEnergyDensity";

  };


protected:
  array2d< real64 > m_damage;
  array2d< real64 > m_strainEnergyDensity;
};

}
} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_SOLID_DAMAGE_HPP_ */
