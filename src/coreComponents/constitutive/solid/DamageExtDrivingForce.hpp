/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */


/**
 * @file DamageExtDrivingForce.hpp
 * @brief This class overrides the SSLE constitutive updates to account for a Damage field
 * @brief This class also introduces a new external driving force in the Damage field evolution equation 

 * Reference: 
 *
 * Kumar, A., Bourdin, B., Francfort, G. A., & Lopez-Pamies, O. (2020). Revisiting nucleation in the phase-field 
 * approach to brittle fracture. Journal of the Mechanics and Physics of Solids, 142, 104027.
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_DAMAGEEXTDRIVFORCE_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_DAMAGEEXTDRIVFORCE_HPP_

#include "Damage.hpp"
#include "InvariantDecompositions.hpp"
#include "SolidBase.hpp"

namespace geosx
{
namespace constitutive
{

template< typename UPDATE_BASE >
class DamageExtDrivForceUpdates : public DamageUpdates< UPDATE_BASE >
{
public:
  template< typename ... PARAMS >
  DamageExtDrivForceUpdates( arrayView2d< real64 > const & inputDamage,
                             arrayView2d< real64 > const & inputStrainEnergyDensity,
                             arrayView2d< real64 > const & inputExtDrivingForce, 
                             real64 const & inputLengthScale,
                             real64 const & inputCriticalFractureEnergy,
                             real64 const & inputcriticalStrainEnergy,
                             real64 const & inputTensileStrength, 
                             real64 const & inputCompressStrength,
                             real64 const & inputDeltaCoefficient,
                             PARAMS && ... baseParams ):
    DamageUpdates< UPDATE_BASE >( inputDamage, inputStrainEnergyDensity, inputExtDrivingForce, inputLengthScale,
                                  inputCriticalFractureEnergy, inputcriticalStrainEnergy,
                                  std::forward< PARAMS >( baseParams )... ), 
    m_tensileStrength( inputTensileStrength ), 
    m_compressStrength( inputCompressStrength ),
    m_deltaCoefficient( inputDeltaCoefficient )
  {}

  using DiscretizationOps = typename DamageUpdates< UPDATE_BASE >::DiscretizationOps;

  using DamageUpdates< UPDATE_BASE >::smallStrainUpdate;
  using DamageUpdates< UPDATE_BASE >::saveConvergedState;

  using DamageUpdates< UPDATE_BASE >::getDegradationValue;
  using DamageUpdates< UPDATE_BASE >::getDegradationDerivative;
  using DamageUpdates< UPDATE_BASE >::getDegradationSecondDerivative;
  using DamageUpdates< UPDATE_BASE >::getEnergyThreshold;

  using DamageUpdates< UPDATE_BASE >::m_strainEnergyDensity;
  using DamageUpdates< UPDATE_BASE >::m_criticalStrainEnergy;
  using DamageUpdates< UPDATE_BASE >::m_extDrivingForce; 
  using DamageUpdates< UPDATE_BASE >::m_criticalFractureEnergy;
  using DamageUpdates< UPDATE_BASE >::m_lengthScale;

  using UPDATE_BASE::m_bulkModulus; 
  using UPDATE_BASE::m_shearModulus;

  GEOSX_HOST_DEVICE
  virtual void smallStrainUpdate( localIndex const k,
                                  localIndex const q,
                                  real64 const ( &strainIncrement )[6],
                                  real64 ( & stress )[6],
                                  DiscretizationOps & stiffness ) const final
  {
    UPDATE_BASE::smallStrainUpdate( k, q, strainIncrement, stress, stiffness );
    real64 factor = getDegradationValue( k, q );

    real64 stressP;
    real64 stressQ;
    real64 deviator[6];

    twoInvariant::stressDecomposition( stress,
                                       stressP,
                                       stressQ,
                                       deviator );

    real64 const mu    = stiffness.m_shearModulus; 
    real64 const kappa = stiffness.m_bulkModulus; 

    // compute volumetric and deviatoric strain invariants
    real64 strain[6];
    UPDATE_BASE::getElasticStrain( k, q, strain );

    // compute invariants of degraded stress 
    real64 I1 = factor * stressP * 3.; 
    real64 sqrt_J2 = factor * stressQ / sqrt(3.); 

    // Calculate the external driving force according to Kumar et al. 
    real64 beta0 = m_deltaCoefficient * 0.375 * m_criticalFractureEnergy / m_lengthScale; 
    
    real64 beta1 = - 0.375 * m_criticalFractureEnergy / m_lengthScale * ((1 + m_deltaCoefficient)*(m_compressStrength - m_tensileStrength)/2./m_compressStrength/m_tensileStrength)
                   - (8*mu + 24*kappa - 27*m_tensileStrength) * (m_compressStrength - m_tensileStrength) / 144. / mu / kappa
                   - m_lengthScale / m_criticalFractureEnergy * ((mu + 3*kappa)*(pow(m_compressStrength, 3) - pow(m_tensileStrength, 3))*m_tensileStrength/18/(mu*mu)/(kappa*kappa)); 
    
    real64 beta2 = - 0.375 * m_criticalFractureEnergy / m_lengthScale * (sqrt(3.)*(1 + m_deltaCoefficient)*(m_compressStrength + m_tensileStrength)/2./m_compressStrength/m_tensileStrength)
                   - (8*mu + 24*kappa - 27*m_tensileStrength)*(m_compressStrength + m_tensileStrength) / 48. / sqrt(3.) / mu / kappa
                   - m_lengthScale / m_criticalFractureEnergy * ((mu + 3*kappa)*(pow(m_compressStrength,3) + pow(m_tensileStrength,3))*m_tensileStrength/6./sqrt(3.)/(mu*mu)/(kappa*kappa)); 

    real64 beta3 = m_lengthScale * (m_tensileStrength/mu/kappa) / m_criticalFractureEnergy; 

    m_extDrivingForce( k, q ) = 1. / (1 + beta3*I1*I1) * (beta2 * sqrt_J2 + beta1*I1 + beta0); 

    LvArray::tensorOps::scale< 6 >( stress, factor );

    stiffness.scaleParams( factor );
  }

  GEOSX_HOST_DEVICE
  virtual real64 getExtDrivingForce( localIndex const k, 
                                     localIndex const q ) const override final
  {
    return m_extDrivingForce( k, q );  
  }

  real64 const m_tensileStrength;
  real64 const m_compressStrength;
  real64 const m_deltaCoefficient; 
};

template< typename BASE >
class DamageExtDrivingForce : public Damage< BASE >
{
public:

  /// @typedef Alias for LinearElasticIsotropicUpdates
  using KernelWrapper = DamageExtDrivForceUpdates< typename BASE::KernelWrapper >;

  using Damage< BASE >::m_damage;
  using Damage< BASE >::m_strainEnergyDensity;
  using Damage< BASE >::m_extDrivingForce; 
  using Damage< BASE >::m_criticalFractureEnergy;
  using Damage< BASE >::m_lengthScale;
  using Damage< BASE >::m_criticalStrainEnergy;

  DamageExtDrivingForce( string const & name, dataRepository::Group * const parent );
  virtual ~DamageExtDrivingForce() override = default;

  static string catalogName() { return string( "DamageExtDrivingForce" ) + BASE::m_catalogNameString; }
  virtual string getCatalogName() const override { return catalogName(); }

  struct viewKeyStruct : public BASE::viewKeyStruct
  {
    /// string/key for the uniaxial tensile strength 
    static constexpr char const * tensileStrengthString() { return "tensileStrength"; }
    /// string/key for the uniaxial compressive strength 
    static constexpr char const * compressStrengthString() { return "compressiveStrength"; }
    /// string/key for a delta coefficient in computing the external driving force  
    static constexpr char const * deltaCoefficientString() { return "deltaCoefficient"; }
  };

  KernelWrapper createKernelUpdates() const
  {
    return BASE::template createDerivedKernelUpdates< KernelWrapper >( m_damage.toView(),
                                                                       m_strainEnergyDensity.toView(),
                                                                       m_extDrivingForce.toView(), 
                                                                       m_lengthScale,
                                                                       m_criticalFractureEnergy,
                                                                       m_criticalStrainEnergy,  
                                                                       m_tensileStrength, 
                                                                       m_compressStrength,
                                                                       m_deltaCoefficient);
  }


protected:
  real64 m_tensileStrength; 
  real64 m_compressStrength; 
  real64 m_deltaCoefficient; 
};

}
} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_SOLID_DAMAGE_HPP_ */
