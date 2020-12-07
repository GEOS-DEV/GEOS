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

  using UPDATE_BASE::getStress;
  using UPDATE_BASE::getElasticStiffness;
  using UPDATE_BASE::smallStrainNoStateUpdate;
  using UPDATE_BASE::smallStrainUpdate;
  using UPDATE_BASE::hypoUpdate;
  using UPDATE_BASE::hyperUpdate;

  using UPDATE_BASE::setDiscretizationOps;
  using UPDATE_BASE::GetStiffness;
  using UPDATE_BASE::SmallStrainNoState;
  using UPDATE_BASE::SmallStrain;
  using UPDATE_BASE::HypoElastic;
  using UPDATE_BASE::HyperElastic;

/*
  GEOSX_HOST_DEVICE inline
  virtual void GetStiffness( localIndex const k,
                             localIndex const q,
                             real64 (& c)[6][6] ) const override final
  {
    UPDATE_BASE::GetStiffness( k, q, c );
    real64 const damageFactor = ( 1.0 - m_damage( k, q ) )*( 1.0 - m_damage( k, q ) );
    for( localIndex i=0; i<6; ++i )
    {
      for( localIndex j=0; j<6; ++j )
      {
        c[i][j] *= damageFactor;
      }
    }
  }
*/

  GEOSX_FORCE_INLINE
  GEOSX_HOST_DEVICE
  void setDiscretizationOps( localIndex const k,
                             localIndex const q,
                             typename UPDATE_BASE::DiscretizationOps & discOps ) const
  {
    UPDATE_BASE::setDiscretizationOps( k, q, discOps );
    real64 const damageFactor = ( 1.0 - m_damage( k, q ) )*( 1.0 - m_damage( k, q ) );
    discOps.scaleParams( damageFactor );
  }


  GEOSX_HOST_DEVICE
  virtual real64 calculateStrainEnergyDensity( localIndex const k,
                                               localIndex const q ) const override final
  {
    real64 const sed = UPDATE_BASE::calculateStrainEnergyDensity( k, q );
    if( sed > m_strainEnergyDensity( k, q ) )
    {
      m_strainEnergyDensity( k, q ) = sed;
    }
    return m_strainEnergyDensity( k, q );
  }

  GEOSX_HOST_DEVICE
  virtual void getStress( localIndex const k,
                          localIndex const q,
                          real64 (& stress)[6] ) const override
  {
    real64 const damageFactor = ( 1.0 - m_damage( k, q ) )*( 1.0 - m_damage( k, q ) );

    stress[0] = this->m_newStress( k, q, 0 ) * damageFactor;
    stress[1] = this->m_newStress( k, q, 1 ) * damageFactor;
    stress[2] = this->m_newStress( k, q, 2 ) * damageFactor;
    stress[3] = this->m_newStress( k, q, 3 ) * damageFactor;
    stress[4] = this->m_newStress( k, q, 4 ) * damageFactor;
    stress[5] = this->m_newStress( k, q, 5 ) * damageFactor;
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
