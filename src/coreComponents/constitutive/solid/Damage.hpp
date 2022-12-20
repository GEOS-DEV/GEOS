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
 * @file Damage.hpp
 * @brief This class overrides the SSLE constitutive updates to account for a Damage field
 *
 * In a phase-field for fracture model, the damage variable affects the Elasticity equation
 * with the degradation of the stresses. Instead of sigma = C : epsilon, we have sigma = g(d)*C:epsilon,
 * where g(d) is the degradation function. This degradation function can either be a quadratic one
 * (set LORENTZ 0) or a quasi-quadratic one (set LORENTZ 1). In general, the quadratic one will give you
 * brittle fracture behaviour. The quasi-quadratic one, combined with linear dissipation, will give you
 * cohesive fracture behaviour, with a user-defined critical stress. If you use quadratic dissipation in
 * your damage solver, set QUADRATIC_DISSIPATION to 1.
 *
 * References:
 *
 * Miehe, Christian; Hofacker, Martina; Welschinger, Fabian. A phase field model for rate-independent crack
 * propagation: Robust algorithmic implementation based on operator splits.
 * Computer Methods in Applied Mechianics and Engineering, v. 199, n. 45-48, p. 2765-2778, 2010
 *
 * Borden, Micheal J., et al. A phase-field description of dynamic brittle fracture.
 * Computer Methods in Applied Mechanics and Engineering, v. 217, p. 77-95, 2012
 *
 * Bourdin, Blaise; Francfort, Gille A.; Marigo, Jean-Jacques. The variational approach to fracture.
 * Journal of Elasticity, v. 91, n. 1-3, p. 5-148, 2008.
 *
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_DAMAGE_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_DAMAGE_HPP_

#include "constitutive/solid/SolidBase.hpp"
#include "InvariantDecompositions.hpp"

namespace geosx
{
namespace constitutive
{

// DAMAGE MODEL UPDATES
//
// NOTE: This model uses the m_newStress array to represent the stress in an
//       elastic, "undamaged" configuration.  We then scale the results
//       by the damage factor whenever the true stress is requested through an update
//       function.  The developer should be very cautious if accessing the stress
//       directly through an arrayView, as it does not represent the true stress.
//
// NOTE: This model is designed to work with phase field implementation where m_damage
//       is updated externally to the material model routine.  A modified implementation
//       would be required to use this directly within a SolidMechanics-only solver, to
//       internally update the damage variable and consistently linearize the system.

template< typename UPDATE_BASE >
class DamageUpdates : public UPDATE_BASE
{
public:
  template< typename ... PARAMS >
  DamageUpdates( arrayView2d< real64 > const & inputDamage,
                 arrayView2d< real64 > const & inputStrainEnergyDensity,
                 arrayView2d< real64 > const & inputExtDrivingForce,
                 real64 const & inputLengthScale,
                 real64 const & inputCriticalFractureEnergy,
                 real64 const & inputcriticalStrainEnergy,
                 real64 const & inputDegradationLowerLimit,
                 integer const & inputExtDrivingForceFlag,
                 real64 const & inputTensileStrength,
                 real64 const & inputCompressStrength,
                 real64 const & inputDeltaCoefficient,
                 PARAMS && ... baseParams ):
    UPDATE_BASE( std::forward< PARAMS >( baseParams )... ),
    m_damage( inputDamage ),
    m_strainEnergyDensity( inputStrainEnergyDensity ),
    m_extDrivingForce ( inputExtDrivingForce ),
    m_lengthScale( inputLengthScale ),
    m_criticalFractureEnergy( inputCriticalFractureEnergy ),
    m_criticalStrainEnergy( inputcriticalStrainEnergy ),
    m_degradationLowerLimit( inputDegradationLowerLimit ),
    m_extDrivingForceFlag( inputExtDrivingForceFlag ),
    m_tensileStrength( inputTensileStrength ),
    m_compressStrength( inputCompressStrength ),
    m_deltaCoefficient( inputDeltaCoefficient )
  {}

  using DiscretizationOps = typename UPDATE_BASE::DiscretizationOps;

  using UPDATE_BASE::smallStrainNoStateUpdate;
  using UPDATE_BASE::smallStrainUpdate;
  using UPDATE_BASE::smallStrainNoStateUpdate_StressOnly;
  using UPDATE_BASE::smallStrainUpdate_StressOnly;
  using UPDATE_BASE::hypoUpdate;
  using UPDATE_BASE::hypoUpdate_StressOnly;
  using UPDATE_BASE::hyperUpdate;
  using UPDATE_BASE::hyperUpdate_StressOnly;
  using UPDATE_BASE::saveConvergedState;

  using UPDATE_BASE::m_disableInelasticity;

  //Standard quadratic degradation functions

  inline
  GEOSX_HOST_DEVICE
  virtual real64 getDegradationValue( localIndex const k,
                                      localIndex const q ) const
  {
    real64 pf;

    if( m_extDrivingForceFlag )
    {
      pf = fmax( fmin( 1.0, m_damage( k, q )), 0.0 );
    }
    else
    {
      pf = m_damage( k, q );
    }

    // Set a lower bound tolerance for the degradation
    real64 const eps = m_degradationLowerLimit;

    return ((1 - eps)*(1 - pf)*(1 - pf) + eps);
  }


  inline
  GEOSX_HOST_DEVICE
  virtual real64 getDegradationDerivative( real64 const d ) const
  {
    return -2*(1 - d);
  }


  inline
  GEOSX_HOST_DEVICE
  virtual real64 getDegradationSecondDerivative( real64 const d ) const
  {
    GEOSX_UNUSED_VAR( d );

    return 2.0;
  }

  GEOSX_HOST_DEVICE
  virtual void smallStrainUpdate( localIndex const k,
                                  localIndex const q,
                                  real64 const ( &strainIncrement )[6],
                                  real64 ( & stress )[6],
                                  DiscretizationOps & stiffness ) const override
  {
    UPDATE_BASE::smallStrainUpdate( k, q, strainIncrement, stress, stiffness );

    if( m_disableInelasticity )
    {
      return;
    }

    real64 factor = getDegradationValue( k, q );

    if( m_extDrivingForceFlag )
    {
      real64 stressP;
      real64 stressQ;
      real64 deviator[6];

      twoInvariant::stressDecomposition( stress,
                                         stressP,
                                         stressQ,
                                         deviator );

      real64 const mu    = stiffness.m_shearModulus;
      real64 const kappa = stiffness.m_bulkModulus;

      // compute invariants of degraded stress
      real64 I1 = factor * stressP * 3.;
      real64 sqrt_J2 = factor * stressQ / sqrt( 3. );

      // Calculate the external driving force according to Kumar et al.
      real64 beta0 = m_deltaCoefficient * 0.375 * m_criticalFractureEnergy / m_lengthScale;

      real64 beta1 = -0.375 * m_criticalFractureEnergy / m_lengthScale * ((1 + m_deltaCoefficient)*(m_compressStrength - m_tensileStrength)/2./m_compressStrength/m_tensileStrength)
                     - (8*mu + 24*kappa - 27*m_tensileStrength) * (m_compressStrength - m_tensileStrength) / 144. / mu / kappa
                     - m_lengthScale / m_criticalFractureEnergy * ((mu + 3*kappa)*(pow( m_compressStrength, 3 ) - pow( m_tensileStrength, 3 ))*m_tensileStrength/18/(mu*mu)/(kappa*kappa));

      real64 beta2 = -0.375 * m_criticalFractureEnergy / m_lengthScale * (sqrt( 3. )*(1 + m_deltaCoefficient)*(m_compressStrength + m_tensileStrength)/2./m_compressStrength/m_tensileStrength)
                     + (8*mu + 24*kappa - 27*m_tensileStrength)*(m_compressStrength + m_tensileStrength) / 48. / sqrt( 3. ) / mu / kappa
                     + m_lengthScale / m_criticalFractureEnergy * ((mu + 3*kappa)*(pow( m_compressStrength, 3 ) + pow( m_tensileStrength, 3 ))*m_tensileStrength/6./sqrt( 3. )/(mu*mu)/(kappa*kappa));

      real64 beta3 = m_lengthScale * (m_tensileStrength/mu/kappa) / m_criticalFractureEnergy;

      m_extDrivingForce( k, q ) = 1. / (1 + beta3*I1*I1) * (beta2 * sqrt_J2 + beta1*I1 + beta0);
    }

    LvArray::tensorOps::scale< 6 >( stress, factor );

    stiffness.scaleParams( factor );
  }


  // TODO: The code below assumes the strain energy density will never be
  //       evaluated in a non-converged / garbage configuration.

  GEOSX_HOST_DEVICE
  virtual real64 getStrainEnergyDensity( localIndex const k,
                                         localIndex const q ) const override
  {
    real64 const sed = SolidBaseUpdates::getStrainEnergyDensity( k, q );

    if( sed > m_strainEnergyDensity( k, q ) )
    {
      m_strainEnergyDensity( k, q ) = sed;
    }

    return m_strainEnergyDensity( k, q );
  }

  GEOSX_HOST_DEVICE
  real64 getRegularizationLength() const
  {
    return m_lengthScale;
  }

  GEOSX_HOST_DEVICE
  real64 getCriticalFractureEnergy() const
  {
    return m_criticalFractureEnergy;
  }

  GEOSX_HOST_DEVICE
  virtual real64 getEnergyThreshold( localIndex const k,
                                     localIndex const q ) const
  {
    #if LORENTZ
    return m_criticalStrainEnergy;
    #else
    if( m_extDrivingForceFlag )
      return 3*m_criticalFractureEnergy/(16 * m_lengthScale) + 0.5 * m_extDrivingForce( k, q );
    else
      return 3*m_criticalFractureEnergy/(16 * m_lengthScale);

    #endif


  }

  arrayView2d< real64 > const m_damage;
  arrayView2d< real64 > const m_strainEnergyDensity;
  arrayView2d< real64 > const m_extDrivingForce;
  real64 const m_lengthScale;
  real64 const m_criticalFractureEnergy;
  real64 const m_criticalStrainEnergy;
  real64 const m_degradationLowerLimit;
  integer const m_extDrivingForceFlag;
  real64 const m_tensileStrength;
  real64 const m_compressStrength;
  real64 const m_deltaCoefficient;
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
  virtual ~Damage() override = default;

  static string catalogName() { return string( "Damage" ) + BASE::m_catalogNameString; }
  virtual string getCatalogName() const override { return catalogName(); }

  virtual void postProcessInput() override;

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  /// *** The interface to get member variables
  arrayView2d< real64 const > getDamage() const { return m_damage; }

  arrayView2d< real64 const > getExtDrivingForce() const { return m_extDrivingForce; }


  KernelWrapper createKernelUpdates() const
  {
    return BASE::template createDerivedKernelUpdates< KernelWrapper >( m_damage.toView(),
                                                                       m_strainEnergyDensity.toView(),
                                                                       m_extDrivingForce.toView(),
                                                                       m_lengthScale,
                                                                       m_criticalFractureEnergy,
                                                                       m_criticalStrainEnergy,
                                                                       m_degradationLowerLimit,
                                                                       m_extDrivingForceFlag,
                                                                       m_tensileStrength,
                                                                       m_compressStrength,
                                                                       m_deltaCoefficient );
  }

  struct viewKeyStruct : public BASE::viewKeyStruct
  {
    static constexpr char const * damageString() { return "damage"; }
    static constexpr char const * strainEnergyDensityString() { return "strainEnergyDensity"; }
    static constexpr char const * extDrivingForceString() { return "extDrivingForce"; }
    /// string/key for regularization length
    static constexpr char const * lengthScaleString() { return "lengthScale"; }
    /// string/key for Gc
    static constexpr char const * criticalFractureEnergyString() { return "criticalFractureEnergy"; }
    /// string/key for sigma_c
    static constexpr char const * criticalStrainEnergyString() { return "criticalStrainEnergy"; }
    /// string/key for degradation lower limit
    static constexpr char const * degradationLowerLimitString() { return "degradationLowerLimit"; }
    // string/key for c_e switch
    static constexpr char const * extDrivingForceFlagString() { return "extDrivingForceFlag"; }
    /// string/key for the uniaxial tensile strength
    static constexpr char const * tensileStrengthString() { return "tensileStrength"; }
    /// string/key for the uniaxial compressive strength
    static constexpr char const * compressStrengthString() { return "compressiveStrength"; }
    /// string/key for a delta coefficient in computing the external driving force
    static constexpr char const * deltaCoefficientString() { return "deltaCoefficient"; }
  };


protected:
  array2d< real64 > m_damage;
  array2d< real64 > m_strainEnergyDensity;
  array2d< real64 > m_extDrivingForce;
  real64 m_lengthScale;
  real64 m_criticalFractureEnergy;
  real64 m_criticalStrainEnergy;
  real64 m_degradationLowerLimit;
  integer m_extDrivingForceFlag;
  real64 m_tensileStrength;
  real64 m_compressStrength;
  real64 m_deltaCoefficient;
};

}
} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_SOLID_DAMAGE_HPP_ */
