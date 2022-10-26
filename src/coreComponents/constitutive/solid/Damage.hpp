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
#include "SolidModelDiscretizationOpsFullyAnisotroipic.hpp"
#include "constitutive/solid/damage/DegradationFunction.hpp"
#include "constitutive/solid/damage/Decomposition.hpp"
#include "constitutive/solid/damage/PressureFunction.hpp"
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
// NOTE: This model is designed to work with phase field implementation where m_newDamage
//       is updated externally to the material model routine.  A modified implementation
//       would be required to use this directly within a SolidMechanics-only solver, to
//       internally update the damage variable and consistently linearize the system.

template< typename UPDATE_BASE >
class DamageUpdates : public UPDATE_BASE
{
public:
  template< typename ... PARAMS >
  DamageUpdates( arrayView2d< real64 > const & inputNewDamage,
                 arrayView3d< real64 > const & inputDamageGrad,
                 arrayView2d< real64 > const & inputStrainEnergyDensity,
                 arrayView2d< real64 > const & inputVolumetricStrain,
                 arrayView2d< real64 > const & inputExtDrivingForce,
                 string const & inputDegradationFunction,
                 string const & inputDecomposition,
                 string const & inputPressureIndicatorFunction,
                 real64 const & inputLengthScale,
                 real64 const & inputCriticalFractureEnergy,
                 real64 const & inputcriticalStrainEnergy,
                 real64 const & inputDegradationLowerLimit,
                 integer const & inputExtDrivingForceFlag,
                 real64 const & inputTensileStrength,
                 real64 const & inputCompressStrength,
                 real64 const & inputDeltaCoefficient,
                 arrayView1d< real64 > const & inputBiotCoefficient,
                 PARAMS && ... baseParams ):
    UPDATE_BASE( std::forward< PARAMS >( baseParams )... ),
    m_newDamage( inputNewDamage ),
    m_damageGrad( inputDamageGrad ),
    m_strainEnergyDensity( inputStrainEnergyDensity ),
    m_volStrain( inputVolumetricStrain ),
    m_extDrivingForce ( inputExtDrivingForce ),
    m_validDegradations({"Quadratic","QuasiQuadratic"}),
    m_validDecompositions({"None","VolDev","Spectral"}),
    m_validPressureIndicators({"Linear","Cosine"}),
    m_degradationOption(-1),
    m_decompositionOption(-1),
    m_pressureIndicatorOption(-1),
    m_degradationFunction( inputDegradationFunction ),
    m_decomposition( inputDecomposition ),
    m_pressureIndicatorFunction( inputPressureIndicatorFunction ),
    m_lengthScale( inputLengthScale ),
    m_criticalFractureEnergy( inputCriticalFractureEnergy ),
    m_criticalStrainEnergy( inputcriticalStrainEnergy ),
    m_degradationLowerLimit( inputDegradationLowerLimit ),
    m_extDrivingForceFlag( inputExtDrivingForceFlag ),
    m_tensileStrength( inputTensileStrength ),
    m_compressStrength( inputCompressStrength ),
    m_deltaCoefficient( inputDeltaCoefficient ),
    m_biotCoefficient( inputBiotCoefficient )
  {
    //set degradation function option
    for (size_t i = 0; i < m_validDegradations.size(); i++) 
    {
      if(m_validDegradations[i]==m_degradationFunction)
      {
        m_degradationOption = i;
      }
    }
    if(m_degradationOption==-1)
    {
      GEOSX_ERROR("In Damage.hpp: "<<m_degradationFunction<<" is not a valid degradation function type.");
    }
    //set decomposition option
    for (size_t i = 0; i < m_validDecompositions.size(); i++) 
    {
      if(m_validDecompositions[i]==m_decomposition)
      {
        m_decompositionOption= i;
      }
    }
    if(m_decompositionOption==-1)
    {
      GEOSX_ERROR("In Damage.hpp: "<<m_decomposition<<" is not a valid strain decomposition.");
    }
    //set pressure indicator option
    for (size_t i = 0; i < m_validPressureIndicators.size(); i++) 
    {
      if(m_validPressureIndicators[i]==m_pressureIndicatorFunction)
      {
        m_pressureIndicatorOption = i;
      }
    }
    if(m_pressureIndicatorOption==-1)
    {
      GEOSX_ERROR("In Damage.hpp: "<<m_pressureIndicatorFunction<<" is not a valid pressure indicator function.");
    }
    
  }

  //using DiscretizationOps = typename UPDATE_BASE::DiscretizationOps;
  using DiscretizationOps = SolidModelDiscretizationOpsFullyAnisotroipic; // could maybe optimize, but general for now


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

  using UPDATE_BASE::m_bulkModulus;  // TODO: model below strongly assumes iso elasticity, templating not so useful
  using UPDATE_BASE::m_shearModulus;

  //Degradation function and its derivatives

  GEOSX_FORCE_INLINE
  GEOSX_HOST_DEVICE
  virtual real64 getDegradationValue( localIndex const k,
                                      localIndex const q ) const
  {
    //this part enforces damage bounds for the nucleation model
    //if everything is correct, it can be applied without the nucleation model as well
    //however, it could hide some bugs that would be detected by seeing negative damage
    real64 pf;
    if( m_extDrivingForceFlag )
    {
      pf = fmax( fmin( 1.0, m_newDamage( k, q )), 0.0 );
    }
    else
    {
      pf = m_newDamage( k,q );
    }

    real64 eps = 0.0;//residual stiffness
    return DegradationFunction::getValue( pf, eps, m_lengthScale, m_criticalFractureEnergy, m_criticalStrainEnergy, m_degradationOption );
  }

  GEOSX_FORCE_INLINE
  GEOSX_HOST_DEVICE
  virtual real64 getDegradationDerivative( real64 const d ) const
  {

    real64 eps = 0.0;//residual stiffness
    return DegradationFunction::getDerivative( d, eps, m_lengthScale, m_criticalFractureEnergy, m_criticalStrainEnergy, m_degradationOption );
  }

  GEOSX_FORCE_INLINE
  GEOSX_HOST_DEVICE
  virtual real64 getDegradationSecondDerivative( real64 const d ) const
  //virtual real64 getDegradationSecondDerivative( real64 const GEOSX_UNUSED_PARAM(d) ) const
  {
    real64 eps = 0.0;//residual stiffness
    return DegradationFunction::getSecondDerivative( d, eps, m_lengthScale, m_criticalFractureEnergy, m_criticalStrainEnergy, m_degradationOption );
  }

  //End degradation functions

  //Pressure indicator function and its derivatives
  GEOSX_FORCE_INLINE
  GEOSX_HOST_DEVICE
  virtual real64 pressureDamageFunction( localIndex const k,
                                         localIndex const q ) const
  {
    return PressureFunction::getValue( m_newDamage( k, q ), m_pressureIndicatorOption );
  }

  GEOSX_FORCE_INLINE
  GEOSX_HOST_DEVICE
  virtual real64 pressureDamageFunctionDerivative( real64 const d ) const
  {
    return PressureFunction::getDerivative( d, m_pressureIndicatorOption );
  }

  GEOSX_FORCE_INLINE
  GEOSX_HOST_DEVICE
  virtual real64 pressureDamageFunctionSecondDerivative( real64 const d ) const
  {
    return PressureFunction::getSecondDerivative( d, m_pressureIndicatorOption );
  }
  //End pressure indicator function 

  //stress and stiffness update

  GEOSX_HOST_DEVICE
  virtual void smallStrainUpdate( localIndex const k,
                                  localIndex const q,
                                  real64 const ( &strainIncrement )[6],
                                  real64 ( & stress )[6],
                                  DiscretizationOps & stiffness ) const final

  {

    //perform undamaged update
    UPDATE_BASE::smallStrainUpdate( k, q, strainIncrement, stress, stiffness.m_c );
    if( m_disableInelasticity )
    {
      return;
    }
    real64 strain[6] = {0};
    UPDATE_BASE::getElasticStrain( k, q, strain );
    //GEOSX_LOG_RANK_0(k);
    //GEOSX_LOG_RANK_0(q);
    real64 traceOfStrain = strain[0] + strain[1] + strain[2];
    //GEOSX_LOG_RANK_0( "m_volStrain: "<<m_volStrain( k, q ));
    //m_volStrain( k, q ) = traceOfStrain;
    //m_volStrain( k, q ) =  strain[0] + strain[1] + strain[2];
    //GEOSX_LOG_RANK_0("After m_volStrain");
    //get degradation value
    real64 factor = getDegradationValue( k, q );
    real64 const mu = m_shearModulus[k];
    real64 const lambda = conversions::bulkModAndShearMod::toFirstLame( m_bulkModulus[k], mu );
    real64 const kappa = m_bulkModulus[k];        

    //compute c_e for nucleation model using intact stress and stiffness
    if( m_extDrivingForceFlag )
    {
      real64 stressP;
      real64 stressQ;
      real64 deviator[6];

      twoInvariant::stressDecomposition( stress,
                                         stressP,
                                         stressQ,
                                         deviator );

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

    real64 sed; //placehold for strain energy density
    //use template class to get decomposition
    //this will compute the damage stresses and stiffness at the current k,q, and also update sed

    Decomposition::stressCalculator(strain, factor, stress, stiffness.m_c, mu, lambda, sed, m_decompositionOption);

    //this enforces the history approach
    if( sed > m_strainEnergyDensity( k, q ) )
    {
      m_strainEnergyDensity( k, q ) = sed;
    }

  }

  //several accessors

  GEOSX_FORCE_INLINE
  GEOSX_HOST_DEVICE
  virtual real64 getDamage( localIndex const k,
                            localIndex const q ) const
  {
    return m_newDamage( k, q );
  }

  GEOSX_FORCE_INLINE
  GEOSX_HOST_DEVICE
  virtual void getDamageGrad( localIndex const k,
                              localIndex const q,
                              real64 ( & damageGrad )[3] ) const
  {
    for( int dim=0; dim < 3; ++dim )
    {
      damageGrad[dim] = m_damageGrad[k][q][dim];
    }
  }

  GEOSX_FORCE_INLINE
  GEOSX_HOST_DEVICE
  virtual void updateBiotCoefficient( localIndex const k,
                                      real64 const biotCoefficient ) const
  {
    m_biotCoefficient[k] = biotCoefficient;
  }

  // TODO: The code below assumes the strain energy density will never be
  //       evaluated in a non-converged / garbage configuration.

  GEOSX_HOST_DEVICE
  virtual real64 getStrainEnergyDensity( localIndex const k,
                                         localIndex const q ) const override
  {
    // real64 const sed = SolidBaseUpdates::getStrainEnergyDensity( k, q );

    // if( sed > m_strainEnergyDensity( k, q ) )
    // {
    //   m_strainEnergyDensity( k, q ) = sed;
    // }

    return m_strainEnergyDensity( k, q );
  }

  GEOSX_HOST_DEVICE
  virtual real64 getVolStrain( localIndex const k,
                               localIndex const q ) const
  {
    return m_volStrain( k, q );
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
    if(m_degradationFunction == "QuasiQuadratic")
    {    
      return m_criticalStrainEnergy;
    }
    if( m_extDrivingForceFlag )
      return 3*m_criticalFractureEnergy/(16 * m_lengthScale) + 0.5 * m_extDrivingForce( k, q );
    else
      return 3*m_criticalFractureEnergy/(16 * m_lengthScale);

  }

  GEOSX_HOST_DEVICE
  virtual real64 getBiotCoefficient( localIndex const k ) const
  {
    return m_biotCoefficient[k];
  }

  GEOSX_HOST_DEVICE
  real64 getParamByName(string const name) const
  {
    if(name=="lengthScale")
    {
        return m_lengthScale;
    }
    if(name=="criticalFractureEnergy")
    {
        return m_criticalFractureEnergy;
    }
    if(name=="criticalStrainEnergy")
    {
        return m_criticalStrainEnergy;
    }
    GEOSX_ERROR("no parameter with name "<< name << " was found in DamageUpdates.");
  }

  arrayView2d< real64 > const m_newDamage;
  arrayView3d< real64 > const m_damageGrad;
  arrayView2d< real64 > const m_strainEnergyDensity;
  arrayView2d< real64 > const m_volStrain;
  arrayView2d< real64 > const m_extDrivingForce;
  std::vector<std::string> m_validDegradations;
  std::vector<std::string> m_validDecompositions;
  std::vector<std::string> m_validPressureIndicators;
  localIndex m_degradationOption;
  localIndex m_decompositionOption;
  localIndex m_pressureIndicatorOption;
  string const m_degradationFunction;
  string const m_decomposition;
  string const m_pressureIndicatorFunction;
  real64 const m_lengthScale;
  real64 const m_criticalFractureEnergy;
  real64 const m_criticalStrainEnergy;
  real64 const m_degradationLowerLimit;
  integer const m_extDrivingForceFlag;
  real64 const m_tensileStrength;
  real64 const m_compressStrength;
  real64 const m_deltaCoefficient;
  arrayView1d< real64 > const m_biotCoefficient;
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
  arrayView2d< real64 const > getDamage() const { return m_newDamage; }

  arrayView2d< real64 const > getExtDrivingForce() const { return m_extDrivingForce; }


  KernelWrapper createKernelUpdates() const
  {
    return BASE::template createDerivedKernelUpdates< KernelWrapper >( m_newDamage.toView(),
                                                                       m_damageGrad.toView(),
                                                                       m_strainEnergyDensity.toView(),
                                                                       m_volStrain.toView(),
                                                                       m_extDrivingForce.toView(),
                                                                       m_degradationFunction,
                                                                       m_decomposition,
                                                                       m_pressureIndicatorFunction,
                                                                       m_lengthScale,
                                                                       m_criticalFractureEnergy,
                                                                       m_criticalStrainEnergy,
                                                                       m_degradationLowerLimit,
                                                                       m_extDrivingForceFlag,
                                                                       m_tensileStrength,
                                                                       m_compressStrength,
                                                                       m_deltaCoefficient,
                                                                       m_biotCoefficient.toView() );
  }

  struct viewKeyStruct : public BASE::viewKeyStruct
  {
    static constexpr char const * damageString() { return "damage"; }
    static constexpr char const * damageGradString() { return "damageGrad"; }
    static constexpr char const * strainEnergyDensityString() { return "strainEnergyDensity"; }
    static constexpr char const * extDrivingForceString() { return "extDrivingForce"; }
    static constexpr char const * volumetricStrainString() { return "volumetricStrain"; }
    /// string/key for degradation function
    static constexpr char const * degradationFunctionString() { return "degradationFunction";}
    /// string/key for strain decomposition
    static constexpr char const * decompositionString() { return "decomposition";}
    /// string/key for pressure indicator function
    static constexpr char const * pressureIndicatorFunctionString() { return "pressureIndicatorFunction";}
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
    /// string/key for the Biot coefficient
    static constexpr char const * biotCoefficientString() { return "biotCoefficient"; }
  };


protected:
  array2d< real64 > m_newDamage;
  array3d< real64 > m_damageGrad;
  array2d< real64 > m_strainEnergyDensity;
  array2d< real64 > m_volStrain;
  array2d< real64 > m_extDrivingForce;
  string m_degradationFunction;
  string m_decomposition;
  string m_pressureIndicatorFunction;
  real64 m_lengthScale;
  real64 m_criticalFractureEnergy;
  real64 m_criticalStrainEnergy;
  real64 m_degradationLowerLimit;
  integer m_extDrivingForceFlag;
  real64 m_tensileStrength;
  real64 m_compressStrength;
  real64 m_deltaCoefficient;
  array1d< real64 > m_biotCoefficient;
};

}
} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_SOLID_DAMAGE_HPP_ */