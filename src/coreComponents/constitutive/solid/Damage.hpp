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
#define LORENTZ 0
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
                 real64 const & inputLengthScale,
                 real64 const & inputCriticalFractureEnergy,
                 real64 const & inputcriticalStrainEnergy,
                 PARAMS && ... baseParams ):
    UPDATE_BASE( std::forward< PARAMS >( baseParams )... ),
    m_damage( inputDamage ),
    m_strainEnergyDensity( inputStrainEnergyDensity ),
    m_lengthScale( inputLengthScale),
    m_criticalFractureEnergy( inputCriticalFractureEnergy),
    m_criticalStrainEnergy( inputcriticalStrainEnergy)
  {}


  using UPDATE_BASE::GetStiffness;
  using UPDATE_BASE::SmallStrainNoState;
  using UPDATE_BASE::SmallStrain;
  using UPDATE_BASE::HypoElastic;
  using UPDATE_BASE::HyperElastic;

  //Quasi-Quadratic Lorentz Degradation Function
  #if LORENTZ

  GEOSX_HOST_DEVICE
  virtual real64 GetDegradationValue( localIndex const k,
                              localIndex const q) const {
     real64 m = 3*m_criticalFractureEnergy/(8*m_lengthScale*m_criticalStrainEnergy);
     real64 p = 1;
     return pow(1 - m_damage( k,q ),2) /( pow(1 - m_damage( k,q ),2) + m * m_damage( k,q ) * (1 + p*m_damage( k,q )) );
  }

  GEOSX_HOST_DEVICE
  virtual real64 GetDegradationDerivative( real64 const d) const {
     real64 m = 3*m_criticalFractureEnergy/(8*m_lengthScale*m_criticalStrainEnergy);
     real64 p = 1;
     return -m*(1 - d)*(1 + (2*p + 1)*d) / pow( pow(1-d,2) + m*d*(1+p*d), 2);
  }

  GEOSX_HOST_DEVICE
  virtual real64 GetDegradationSecondDerivative( real64 const d) const {
     real64 m = 3*m_criticalFractureEnergy/(8*m_lengthScale*m_criticalStrainEnergy);
     real64 p = 1;
     return -2*m*( pow(d,3)*(2*m*p*p + m*p + 2*p + 1) + pow(d,2)*(-3*m*p*p -3*p) + d*(-3*m*p - 3) + (-m+p+2) )/pow( pow(1-d,2) + m*d*(1+p*d), 3);
  }

  #else
  //Standard Quadratic Degradation Function

  GEOSX_HOST_DEVICE
  virtual real64 GetDegradationValue( localIndex const k,
                              localIndex const q) const {
     return (1 - m_damage( k,q ))*(1 - m_damage( k,q ));
  }

  GEOSX_HOST_DEVICE
  virtual real64 GetDegradationDerivative( real64 const d) const {
     return -2*(1 - d);
  }

  GEOSX_HOST_DEVICE
  virtual real64 GetDegradationSecondDerivative( real64 const d) const {
     return 2 * (d - d + 1);
  }
  #endif

  GEOSX_HOST_DEVICE inline
  virtual void GetStiffness( localIndex const k,
                             localIndex const q,
                             real64 (& c)[6][6] ) const override
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

  GEOSX_HOST_DEVICE
  virtual real64 calculateActiveStrainEnergyDensity( localIndex const k,
                                               localIndex const q ) const 
  {
    real64 const sed = UPDATE_BASE::calculateStrainEnergyDensity(k,q);

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
    //no tension-compression asymmetry
    
    real64 const damageFactor = GetDegradationValue( k,q );
  
    stress[0] = this->m_stress(k,q,0) * damageFactor;
    stress[1] = this->m_stress(k,q,1) * damageFactor;
    stress[2] = this->m_stress(k,q,2) * damageFactor;
    stress[3] = this->m_stress(k,q,3) * damageFactor;
    stress[4] = this->m_stress(k,q,4) * damageFactor;
    stress[5] = this->m_stress(k,q,5) * damageFactor;

  }

  GEOSX_HOST_DEVICE
  real64 getRegularizationLength() const {
    return m_lengthScale;
  }

  GEOSX_HOST_DEVICE
  real64 getCriticalFractureEnergy() const {
    return m_criticalFractureEnergy;
  }

  GEOSX_HOST_DEVICE
  real64 getEnergyThreshold() const {
    #if LORENTZ
    return m_criticalStrainEnergy;
    #else
    return 3*m_criticalFractureEnergy/(16 * m_lengthScale);
    #endif
  }

  arrayView2d< real64 > const m_damage;
  arrayView2d< real64 > const m_strainEnergyDensity;
  real64 const m_lengthScale;
  real64 const m_criticalFractureEnergy;
  real64 const m_criticalStrainEnergy;

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
                                                                       m_strainEnergyDensity.toView(),
                                                                       m_lengthScale,
                                                                       m_criticalFractureEnergy,
                                                                       m_criticalStrainEnergy);
  }

  struct viewKeyStruct : public BASE::viewKeyStruct
  {
    static constexpr auto damageString =  "damage";
    static constexpr auto strainEnergyDensityString =  "strainEnergyDensity";
    /// string/key for regularization length
    static constexpr auto lengthScaleString  = "lengthScale";
    /// string/key for Gc
    static constexpr auto criticalFractureEnergyString = "criticalFractureEnergy";
    /// string/key for sigma_c
    static constexpr auto criticalStrainEnergyString = "criticalStrainEnergy";

  };


protected:
  array2d< real64 > m_damage;
  array2d< real64 > m_strainEnergyDensity;
  real64 m_lengthScale;
  real64 m_criticalFractureEnergy;
  real64 m_criticalStrainEnergy;
};

}
} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_SOLID_DAMAGE_HPP_ */
