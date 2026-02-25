/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file CeramicDamage.hpp
 * @brief Simple damage model for modeling material failure in brittle materials.
 *
 * This damage model is intended for use with damage-field partitioning (DFG) within the
 * MPM solver, but can also be used without DFG by any solver. It is only appropriate for
 * schemes implementing explicit time integration. The model is really a hybrid plasticity/
 * damage model in the sense that we assume damaged material behaves like granular material
 * and hence follows a modified Mohr-Coulomb law. The modifications are that at low pressures,
 * the shape of the yield surface is modified to resemble a maximum principal stress criterion,
 * and at high pressures the shape converges on the von Mises yield surface. The damage
 * parameter results in softening of the deviatoric response i.e. causes the yield surface to
 * contract. Furthermore, damage is used to scale back tensile pressure: p = (1 - d) * pTrial.
 * pTrial is calculatd as pTrial = -k * log(J), where the Jacobian J of the material motion is
 * integrated and tracked by this model.
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_CERAMICDAMAGE_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_CERAMICDAMAGE_HPP_

#include "ElasticIsotropic.hpp"
#include "InvariantDecompositions.hpp"
#include "PropertyConversions.hpp"
#include "SolidModelDiscretizationOpsFullyAnisotroipic.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geos
{

namespace constitutive
{

/**
 * @class CeramicDamageUpdates
 *
 * Class to provide material updates that may be
 * called from a kernel function.
 */
class CeramicDamageUpdates : public ElasticIsotropicUpdates
{
public:

  /**
   * @brief Constructor
   * @param[in] damage The ArrayView holding the damage for each quardrature point.
   * @param[in] jacobian The ArrayView holding the jacobian for each quardrature point.
   * @param[in] lengthScale The ArrayView holding the length scale for each element.
   * @param[in] strengthScale The ArrayView holding the strength scale for each element/particle.
   * @param[in] tensileStrength The unconfined tensile strength.
   * @param[in] compressiveStrength The unconfined compressive strength.
   * @param[in] maximumStrength The theoretical maximum strength.
   * @param[in] crackSpeed The crack speed.
   * @param[in] bulkModulus The ArrayView holding the bulk modulus data for each element.
   * @param[in] shearModulus The ArrayView holding the shear modulus data for each element.
   * @param[in] thermalExpansionCoefficient The ArrayView holding the thermal expansion coefficient data for each element.
   * @param[in] newStress The ArrayView holding the new stress data for each quadrature point.
   * @param[in] oldStress The ArrayView holding the old stress data for each quadrature point. 
   * @param[in] refStrainRate The reference strain rate for rate-dependent strength scaling.
   * @param[in] strainRate The ArrayView holding the strain rate data for each quadrature point (for rate-dependent scaling). //added here by SG
   * @param[in] instTensileStrength Instantaneous Tensile Strength values.
   * @param[in] instCompressiveStrength Instantaneous Compressive Strength values.
   * @param[in] instStrength Instantaneous Strength values.
   * @param[in] instPressure Instantaneous pressure values.
   * @param[in] rateSensitivity The rate sensitivity for rate-dependent strength scaling.
   * @param[in] m2 The second slope for pressure dependent strength scaling.
   */
  CeramicDamageUpdates( arrayView2d< real64 > const & damage,
                        arrayView2d< real64 > const & jacobian,
                        arrayView1d< real64 > const & lengthScale,
                        arrayView1d< real64 > const & strengthScale,
                        arrayView1d< real64 > const & porosity,
                        arrayView1d< real64 > const & referencePorosity,
                        real64 const & tensileStrength,
                        real64 const & compressiveStrength,
                        real64 const & maximumStrength,
                        real64 const & crackSpeed,
                        real64 const & damagedMaterialFrictionSlope,
                        int const & thirdInvariantDependence,
                        arrayView3d< real64 > const & velocityGradient,
                        arrayView3d< real64 > const & plasticStrain,
                        //arrayView1d< real64 const > const & bulkModulus,
                        arrayView1d< real64 > const & bulkModulus,  //added by SG or rather removed real64 const
                        arrayView1d< real64 const > const & shearModulus,
                        arrayView1d< real64 const > const & thermalExpansionCoefficient,
                        arrayView3d< real64, solid::STRESS_USD > const & newStress,
                        arrayView3d< real64, solid::STRESS_USD > const & oldStress,
                        arrayView2d< real64 > const & density,
                        arrayView2d< real64 > const & wavespeed,                       
                        real64  const & refStrainRate,
                        arrayView2d< real64 > const & strainRate, //added by SG
                        arrayView2d< real64 > const & instTensileStrength, //added by SG
                        arrayView2d< real64 > const & instCompressiveStrength, //added by SG
                        arrayView2d< real64 > const & instStrength, //added by SG
                        arrayView2d< real64 > const & instPressure, //added by SG
                        real64  const & rateSensitivity,
                        real64  const & m2,
                        bool const & disableInelasticity ):
    ElasticIsotropicUpdates( bulkModulus,
                             shearModulus,
                             thermalExpansionCoefficient,
                             newStress,
                             oldStress,
                             density,
                             wavespeed,
                             disableInelasticity ),
    m_damage( damage ),
    m_jacobian( jacobian ),
    m_lengthScale( lengthScale ),
    m_strengthScale( strengthScale ),
    m_porosity( porosity ),
    m_referencePorosity( referencePorosity ),
    m_tensileStrength( tensileStrength ),
    m_compressiveStrength( compressiveStrength ),
    m_maximumStrength( maximumStrength ),
    m_crackSpeed( crackSpeed ),
    m_damagedMaterialFrictionSlope( damagedMaterialFrictionSlope ),
    m_thirdInvariantDependence( thirdInvariantDependence ),
    m_velocityGradient( velocityGradient ),
    m_plasticStrain( plasticStrain ),
    m_refStrainRate( refStrainRate ),
    m_strainRate( strainRate ), // <-- Added by SG
    m_instTensileStrength( instTensileStrength ), // <-- Added by SG
    m_instCompressiveStrength( instCompressiveStrength ), // <-- Added by SG
    m_instStrength( instStrength ), // <-- Added by SG
    m_instPressure( instPressure ), // <-- Added by SG
    m_rateSensitivity( rateSensitivity ),  // <-- Added by SG
    m_m2( m2 )  // <-- Added by SG

  {}

  /// Default copy constructor
  CeramicDamageUpdates( CeramicDamageUpdates const & ) = default;

  /// Default move constructor
  CeramicDamageUpdates( CeramicDamageUpdates && ) = default;

  /// Deleted default constructor
  CeramicDamageUpdates() = delete;

  /// Deleted copy assignment operator
  CeramicDamageUpdates & operator=( CeramicDamageUpdates const & ) = delete;

  /// Deleted move assignment operator
  CeramicDamageUpdates & operator=( CeramicDamageUpdates && ) =  delete;

  /// Use the uncompressed version of the stiffness bilinear form
  using DiscretizationOps = SolidModelDiscretizationOpsFullyAnisotroipic; // TODO: typo in anistropic (fix in DiscOps PR)

  // Bring in base implementations to prevent hiding warnings
  using ElasticIsotropicUpdates::smallStrainUpdate;

  GEOS_HOST_DEVICE
  void smallStrainUpdate( localIndex const k,
                          localIndex const q,
                          real64 const & timeIncrement,
                          real64 const ( &strainIncrement )[6],
                          real64 ( &stress )[6],
                          real64 ( &stiffness )[6][6] ) const;

  GEOS_HOST_DEVICE
  virtual void smallStrainUpdate( localIndex const k,
                                  localIndex const q,
                                  real64 const & timeIncrement,
                                  real64 const ( &strainIncrement )[6],
                                  real64 ( &stress )[6],
                                  DiscretizationOps & stiffness ) const;

  GEOS_HOST_DEVICE
  virtual void smallStrainUpdate_StressOnly( localIndex const k,
                                             localIndex const q,
                                             real64 const & timeIncrement,
                                             real64 const ( &strainIncrement )[6],
                                             real64 ( &stress )[6] ) const override;

  GEOS_HOST_DEVICE
  virtual void smallStrainUpdate_StressOnly( localIndex const k,
                                             localIndex const q,
                                             real64 const & timeIncrement,
                                             real64 const ( & beginningRotation )[3][3],
                                             real64 const ( & endRotation )[3][3],
                                             real64 const ( &strainIncrement )[6],
                                             real64 ( &stress )[6] ) const override;

  GEOS_HOST_DEVICE
  void smallStrainUpdateHelper( localIndex const k,
                                localIndex const q,
                                real64 const timeIncrement,
                                real64 const ( & beginningRotation )[3][3],
                                real64 const ( & endRotation )[3][3],
                                real64 const ( &strainIncrement )[6], //added by SG
                                real64 ( &stress )[6] ) const; 

  GEOS_HOST_DEVICE
  real64 getStrength( const real64 damage,      // damage
                      const real64 pressure,    // pressure
                      const real64 J2,          // J2 invariant of stress
                      const real64 J3,          // J3 invariant of stress
                      const real64 mu,          // friction slope
                      const real64 Yc,          // Compressive strength
                      const real64 Yt0,         // Tensile parameter
                      const real64 Ymax ) const; // Max strength

GEOS_HOST_DEVICE
real64 ceramicY10( const real64 pLocal,   // pressure
                                         const real64 dLocal,   // damage,
                                         const real64 muLocal,  // friction slope
                                         const real64 Yt0Local, // tensile strength parameter
                                         const real64 YcLocal ) const;

GEOS_HOST_DEVICE
real64 ceramicdY10dp(const real64 d, // damage,
                                const real64 mu, // friction slope
                                const real64 Yc, // unconfined compressive strength
                                const real64 Yt0 ) const; // unconfined tensile strength before 3rd invariant scaling

GEOS_HOST_DEVICE
real64 ceramicdY20dp(const real64 p, // pressure
                                const real64 d,   // damage,
                                const real64 mu,  // friction slope
                                const real64 Yc,  // unconfined compressive strength
                                const real64 Yt0,  // unconfined tensile strength before 3rd invariant scaling
                                const real64 Ymax, // max shear stress
                                const real64 m2 ) const; 

GEOS_HOST_DEVICE
real64 smoothStep(const real64 x,
                             const real64 xmin,
                             const real64 xmax) const;

GEOS_HOST_DEVICE
real64 thirdInvariantStrengthScaling( const real64 J2,
                                      const real64 J3,
                                      const real64 dfdp ) const;

                                

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  virtual void saveConvergedState( localIndex const k,
                                   localIndex const q ) const override final
  {
    ElasticIsotropicUpdates::saveConvergedState( k, q );
  }

private:
  /// A reference to the ArrayView holding the damage for each quadrature point.
  arrayView2d< real64 > const m_damage;

  /// A reference to the ArrayView holding the jacobian for each quadrature point.
  arrayView2d< real64 > const m_jacobian;

  /// A reference to the ArrayView holding the length scale for each element/particle.
  arrayView1d< real64 > const m_lengthScale;

  /// A reference to the ArrayView holding the strength scale.
  arrayView1d< real64 > const m_strengthScale;

  /// A reference to the ArrayView holding the porosity
  arrayView1d< real64 > const m_porosity;

  ///A reference to the ArrayView holding the reference porosity
  arrayView1d< real64 > const m_referencePorosity;

  /// The tensile strength
  real64 const m_tensileStrength;

  /// The compressive strength
  real64 const m_compressiveStrength;

  /// The maximum theoretical strength
  real64 const m_maximumStrength;

  // The crack speed
  real64 const m_crackSpeed;

  /// The damaged material friction slope
  real64 const m_damagedMaterialFrictionSlope;

  // The third invariant dependence flag
  int const m_thirdInvariantDependence;

  /// State variable: The velocity gradient for each element/particle
  arrayView3d< real64 > const m_velocityGradient;

  /// State variable: The plastic strain values for each quadrature point
  arrayView3d< real64 > const m_plasticStrain;

  real64 m_refStrainRate;       // <-- ADD THIS

  arrayView2d< real64 > const m_strainRate; //Added by SG

  arrayView2d< real64 > const m_instTensileStrength; //Added by SG

  arrayView2d< real64 > const m_instCompressiveStrength; //Added by SG

  arrayView2d< real64 > const m_instStrength; //Added by SG

  arrayView2d< real64 > const m_instPressure; //Added by SG

  real64 m_rateSensitivity;     // <-- Added by SG

  real64 m_m2;     // <-- Added by SG

  

};


GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void CeramicDamageUpdates::smallStrainUpdate( localIndex const k,
                                              localIndex const q,
                                              real64 const & timeIncrement,
                                              real64 const ( &strainIncrement )[6],
                                              real64 ( & stress )[6],
                                              real64 ( & stiffness )[6][6] ) const
{

  m_jacobian[k][q] *= exp( strainIncrement[0] + strainIncrement[1] + strainIncrement[2] ); //changed

  real64 epsilon_v = log(m_jacobian[k][q]);

  real64 S=1.61;  //for concrete
  real64 C0=3.7;
  real64 rho0=2.4;
  real64 K_0= rho0*C0*C0; 

  real64 num=C0*C0*rho0*exp(epsilon_v)*(S+1-S*exp(epsilon_v));
  real64 den=pow((S*exp(epsilon_v)-S+1),3);


  m_bulkModulus[k] = fmax(K_0,num/den); // to be changed by SG
  
  // Elastic trial update (assume strainIncrement is all elastic)
  ElasticIsotropicUpdates::smallStrainUpdate( k, 
                                              q, 
                                              timeIncrement,
                                              strainIncrement, 
                                              stress, 
                                              stiffness );
  

  if( m_disableInelasticity )
  {
    return;
  }

  real64 beginningRotation[3][3] = { { 0 } };
  beginningRotation[0][0] = 1.0;
  beginningRotation[1][1] = 1.0;
  beginningRotation[2][2] = 1.0;

  real64 endRotation[3][3] = { { 0 } }; 
  endRotation[0][0] = 1.0;
  endRotation[1][1] = 1.0;
  endRotation[2][2] = 1.0;

  CeramicDamageUpdates::smallStrainUpdateHelper( k, 
                                                 q, 
                                                 timeIncrement, 
                                                 beginningRotation, 
                                                 endRotation, 
                                                 strainIncrement, //added by SG
                                                 stress ); 

  // It doesn't make sense to modify stiffness with this model

  // save new stress and return
  saveStress( k, q, stress );
  return;
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void CeramicDamageUpdates::smallStrainUpdate( localIndex const k,
                                              localIndex const q,
                                              real64 const & timeIncrement,
                                              real64 const ( &strainIncrement )[6],
                                              real64 ( & stress )[6],
                                              DiscretizationOps & stiffness ) const
{
  smallStrainUpdate( k, 
                     q, 
                     timeIncrement,
                     strainIncrement, 
                     stress, 
                     stiffness.m_c );
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void CeramicDamageUpdates::smallStrainUpdate_StressOnly( localIndex const k,
                                                         localIndex const q,
                                                         real64 const & timeIncrement,
                                                         real64 const ( & strainIncrement )[6],
                                                         real64 ( & stress )[6] ) const
{
  GEOS_UNUSED_VAR( k );
  GEOS_UNUSED_VAR( q );
  GEOS_UNUSED_VAR( timeIncrement );
  GEOS_UNUSED_VAR( strainIncrement );
  GEOS_UNUSED_VAR( stress );
  GEOS_ERROR( "smallStrainUpdateStressOnly overload not implemented for CeramicDamage" );
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void CeramicDamageUpdates::smallStrainUpdate_StressOnly( localIndex const k,
                                                         localIndex const q,
                                                         real64 const & timeIncrement,
                                                         real64 const ( & beginningRotation )[3][3],
                                                         real64 const ( & endRotation )[3][3],
                                                         real64 const ( &strainIncrement )[6],
                                                         real64 ( & stress )[6] ) const
{
  GEOS_UNUSED_VAR( beginningRotation );
  GEOS_UNUSED_VAR( endRotation );

  // elastic predictor (assume strainIncrement is all elastic)
  ElasticIsotropicUpdates::smallStrainUpdate_StressOnly( k, 
                                                         q, 
                                                         timeIncrement,
                                                         strainIncrement, 
                                                         stress );

  m_jacobian[k][q] *= exp( strainIncrement[0] + strainIncrement[1] + strainIncrement[2] );

  if( m_disableInelasticity )
  {
    return;
  }

  // Call the constitutive model
  CeramicDamageUpdates::smallStrainUpdateHelper( k, 
                                                 q, 
                                                 timeIncrement,
                                                 beginningRotation, 
                                                 endRotation, 
                                                 strainIncrement,
                                                 stress );

  // Save new stress and return
  saveStress( k, q, stress );
  return;
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void CeramicDamageUpdates::smallStrainUpdateHelper( localIndex const k,
                                                    localIndex const q,
                                                    real64 const timeIncrement,
                                                    real64 const ( & beginningRotation )[3][3],
                                                    real64 const ( & endRotation )[3][3],
                                                    real64 const ( &strainIncrement )[6], //added by SG
                                                    real64 ( & stress )[6] ) const 
{
  GEOS_UNUSED_VAR( beginningRotation );
  GEOS_UNUSED_VAR( endRotation );


  //static real64 maxStrainRate = 0.0;

  // cohesion slope
  real64 mu = m_damagedMaterialFrictionSlope;

  // Scaled strengths
  real64 gamma = m_compressiveStrength / m_tensileStrength; // TXC/TXE

  real64 Yt = m_strengthScale[k] * m_tensileStrength * ( 1.0 - m_porosity[k] );
  real64 Yc = m_strengthScale[k] * m_compressiveStrength * ( 1.0 - m_porosity[k] );

  real64 Ycmax = m_maximumStrength;
  real64 Ytmax = Ycmax / gamma;

  (void) Ycmax;
  (void) Ytmax;

  //strain rate

  real64 Yt_baseline=Yt;
  real64 Yc_baseline=Yc;

  real64 SR =
  sqrt( (strainIncrement[0]*strainIncrement[0] +
  strainIncrement[1]*strainIncrement[1] +
  strainIncrement[2]*strainIncrement[2] +
  2*strainIncrement[3]*strainIncrement[3] +
  2*strainIncrement[4]*strainIncrement[4] +
  2*strainIncrement[5]*strainIncrement[5] )) / timeIncrement ;

  //  //added by SG
  m_strainRate[k][q] = SR;

  real64 rateScaling_T=std::pow( (SR / (m_refStrainRate)), m_rateSensitivity );
  real64 rateScaling_C=std::pow( (SR / (m_refStrainRate)), m_rateSensitivity );

  Yt *= (rateScaling_T);
  Yc *= (rateScaling_C);

  Yt = std::max(Yt_baseline,Yt); //added for cutoff strength cannot go below initial
  Yc = std::max(Yc_baseline,Yc); //added for cutoff strength cannot go below initial


    // Yt = std::min(Yt, 0.999*Ytmax);
    // Yc = std::min(Yc, 0.999*Ycmax);

  //added by SG

  m_instTensileStrength[k][q]=Yt;
  m_instCompressiveStrength[k][q]=Yc;
  
  // get failure time
  real64 tFail = m_lengthScale[k] / m_crackSpeed;
  //tFail /= rateFactor;  //Faster degradation 

  // get trial pressure
  real64 pressure = -m_bulkModulus[k] * log( m_jacobian[k][q] );

  // Intermediate strength parameter
  real64 Yt0 = m_thirdInvariantDependence == 1 ? fmax( 0.5 * Yt, std::min( 2.0 * Yt, (3.0 * Yc * Yt ) / ( 2.0 * Yc + Yt + 1.0e-16 ) ) ) : Yt;
 
  // real64 Yt0 = fmax( 0.5 * Yt, fmin( 2.0 * Yt, (3.0 * Yc * Yt) / (2.0 * Yc + Yt) ) );
  Yt0 = fmin( Yt0, ( 3.0 * Yc - Yc * mu ) / ( 3.0 + mu ) );

  // Compute the vertex pressure (should be pmin0 < 0) for the undamaged yield surface.
  real64 pmin0 = -( 2.0 * Yc * Yt0 ) / ( 3.0 * ( Yc - Yt0 ) );
  pmin0 = fmin( pmin0, -1.0e-12 );
  real64 pmin = ( 1.0 - m_damage[k][q] ) * pmin0;


  //added by SG
  m_instPressure[k][q]=pressure;

  // Enforce vertex solution
  // if( trialPressure < 0 )
  if( m_jacobian[k][q] > 1.0 )
  {
    pressure *= ( 1.0 - m_damage[k][q] ); // Tensile cutoff pressure (negative value in tension), scaled by damage. Goes to 0
                                                         // as damage -> 1.
  }
  
  if( pressure < pmin ) // TODO: pressure or trial pressure?
  {
    // Increment damage
    m_damage[k][q] = fmin( m_damage[k][q] + timeIncrement / tFail, 1.0 );

    // Pressure is on the vertex
    pressure = ( 1.0 - m_damage[k][q] ) * pmin0;

    // updated stress is isotropic, at the vertex:
    stress[0] = -pressure;
    stress[1] = -pressure;
    stress[2] = -pressure;
    stress[3] = 0.0;
    stress[4] = 0.0;
    stress[5] = 0.0;
  }
  else // Enforce strength solution
  {
    real64 meanStress;    // negative of pressure
    real64 vonMises;      // von Mises stress
    real64 deviator[6] = { 0 };   // direction of stress deviator
    twoInvariant::stressDecomposition( stress,
                                       meanStress,
                                       vonMises,
                                       deviator );


    real64 brittleDuctileTransitionPressure = Ycmax / mu;
    real64 J2 = vonMises * vonMises / 3.0;
    real64 J3 = vonMises * vonMises * vonMises *
                ( deviator[0] * deviator[1] * deviator[2] +
                  2.0 * deviator[3] * deviator[4] * deviator[5] -
                  deviator[0] * deviator[3] * deviator[3] -
                  deviator[1] * deviator[4] * deviator[4] -
                  deviator[2] * deviator[5] * deviator[5] );

    // Find the strength
    real64 strength = CeramicDamageUpdates::getStrength( m_damage[k][q], pressure, J2, J3, mu, Yc, Yt0, Ycmax );

    // Increment damage and get new associated yield surface
    real64 newDeviatorMagnitude = vonMises;
    if( vonMises > strength )
    {
      if( pressure <= brittleDuctileTransitionPressure )
      {
        m_damage[k][q] = fmin( m_damage[k][q] + timeIncrement / tFail, 1.0 );
        strength = CeramicDamageUpdates::getStrength( m_damage[k][q], pressure, J2, J3, mu, Yc, Yt0, Ycmax );
      }
      newDeviatorMagnitude = strength;
    }


    //added by SG
    m_instStrength[k][q]=strength;

    // Radial return
    twoInvariant::stressRecomposition( -pressure,
                                       newDeviatorMagnitude,
                                       deviator,
                                       stress );
  }
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 CeramicDamageUpdates::getStrength( const real64 damage,     // damage
                                          const real64 pressure,   // pressure
                                          const real64 J2,         // J2 invariant of stress
                                          const real64 J3,         // J3 invariant of stress
                                          const real64 mu,         // friction slope
                                          const real64 Yc,
                                          const real64 Yt0,
                                          const real64 Ymax ) const // strength parameter
{
  real64 dfdp = 0.0;
  real64 oneOverGamma = 1.0;

  real64 p1 = Yc / 3.0;
  //real64 p2 = Ymax / mu;

  real64 m2=m_m2; // slope m2 which should be smaller than m1
  //real64 p2=(Yc-p1*m2)/(mu-m2); //added by SG

  real64 p2=(Ymax)/(mu-m2); 
  
  // Determine scaled strength
  if( pressure <= p1 )
  {
    dfdp = ceramicdY10dp( damage, mu, Yc, Yt0 );
    oneOverGamma = m_thirdInvariantDependence == 1 ? thirdInvariantStrengthScaling( J2, J3, dfdp ) : 1.0;
    return oneOverGamma * ceramicY10( pressure, damage, mu, Yt0, Yc );
  }
  
  if( pressure < p2 )
  {
    dfdp = ceramicdY20dp( pressure, damage, mu, Yc, Yt0, Ymax, m2 ); //m2 was added by SG
	  oneOverGamma = m_thirdInvariantDependence == 1 ? thirdInvariantStrengthScaling( J2, J3, dfdp ) : 1.0;

    real64 m1 = oneOverGamma * ceramicdY10dp( damage, mu, Yc, Yt0 );
    real64 y1 = oneOverGamma * ceramicY10( p1, damage, mu, Yt0, Yc);
   // real64 y2 = oneOverGamma * Ymax;
    //real64 y2 = y1+ oneOverGamma * (p2-p1)*m2; //added by SG
    real64 y2 = oneOverGamma * (Ymax + m2 * pressure); //added by SG
    return pow((pressure - p2) / (p1 - p2), m1 * (p1 - p2) / (y1 - y2)) * (y1 - y2) + y2;
  }
  else
  {
    oneOverGamma = m_thirdInvariantDependence == 1 ? thirdInvariantStrengthScaling( J2, J3, m2 ) : 1.0; //changed dfdp to m2
    //return oneOverGamma * Ymax;
    return oneOverGamma * (Ymax+m2 * (pressure)); //added by SG -p2 add this
  }
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 CeramicDamageUpdates::ceramicY10( const real64 pLocal,   // pressure
                                         const real64 dLocal,   // damage,
                                         const real64 muLocal,  // friction slope
                                         const real64 Yt0Local, // tensile strength parameter
                                         const real64 YcLocal ) const 
{
  return (((3.0 + dLocal * (-3.0 + muLocal)) * YcLocal + 
              (-3.0 + dLocal * (3.0 + muLocal)) * Yt0Local) * (pLocal - 
              (2.0 * (dLocal - 1.0) * YcLocal * Yt0Local) / 
              (3.0 * (YcLocal - Yt0Local)))) / 
              (YcLocal + Yt0Local);
};

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 CeramicDamageUpdates::ceramicdY10dp( const real64 d, // damage,
                                            const real64 mu, // friction slope
                                            const real64 Yc, // unconfined compressive strength
                                            const real64 Yt0 ) const // unconfined tensile strength before 3rd invariant scaling
{
  // Linear portion of the Y'(p) yield strength for the ceramic model in the linear region,
  // before application of 3rd invariant scaling.
  return ((3 + d*(-3 + mu))*Yc + (-3 + d*(3 + mu))*Yt0)/(Yc + Yt0);
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 CeramicDamageUpdates::ceramicdY20dp( const real64 p, // pressure
                                            const real64 d,   // damage,
                                            const real64 mu,  // friction slope
                                            const real64 Yc,  // unconfined compressive strength
                                            const real64 Yt0,  // unconfined tensile strength before 3rd invariant scaling
                                            const real64 Ymax, // max shear stress
                                            const real64 m2 ) const // new: slope of high-pressure asymptote
                                             
{
  // This slope is just used to define the third invatiant dependence scaling, rather than use the actual
  // dfdp, which is discontinuous at Ymax/mu for d=1, we use a smooth blending function.

  real64 dfdp1 = ceramicdY10dp( d, mu, Yc, Yt0 );
  real64 p1 = Yc/3;
  //real64 p2 = (Yc-p1*m2)/(mu-m2); 
  //real64 p2 = Ymax/mu;
  real64 p2=(Ymax)/(mu-m2); 
  //return dfdp1*( 1.0 - smoothStep(p,p1,p2) );

  //(void)Ymax; //added by SG

  // Blend from dfdp1 at p1 to m2 at p2
    real64 s = smoothStep(p, p1, p2); //added by SG
    return dfdp1 * (1.0 - s) + m2 * s;//added by SG

  
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 CeramicDamageUpdates::smoothStep( const real64 x,
                                         const real64 xmin,
                                         const real64 xmax ) const
{
  // Smooth blending function from 0 to 1 as
  // x goes from xmin to xmax.
  //
  // will fail if xmax=xmin, so don't do that.

  if(x <= xmin)
  {
    return 0.0;
  }
  else if(x >= xmax)
  {
    return 1.0;
  }
  else
  {
    real64 xi = (x - xmin)/(xmax - xmin);
    return (3.0*xi*xi - 2.0*xi*xi*xi );
  }
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 CeramicDamageUpdates::thirdInvariantStrengthScaling( const real64 J2,         // J2 invariant of stress
                                                            const real64 J3,         // J3 invariant of stress
                                                            const real64 dfdp ) const // strength parameter
{
  real64 oneOverGamma = 1.0; // This is the ratio of strength relative to the unscaled value.

  // pressure dependent scaling based on slop of strength vs. pressure
  // This ignores the friction cutoff for failed material.
  real64 psi = std::min( 2.0, std::max( 0.5, 1.0 / ( 1.0 + dfdp / 3. ) ) );

  // Compute Lode angle
  if( J2 > 1e-12 )
  {
    // Lode angle
    real64 theta = ( 1.0 / 3.0 ) * asin( std::min( 1.0, std::max( -1.0, -0.5 * J3 * std::pow( 3.0 / J2, 1.5 ) ) ) );

    // This is the Willam-Warnke third-invariant scale function as defined in the Kayenta manual.
    real64 cosPi6plusTheta = cos( 0.5235987755982989 + theta );
    real64 num = 2 * ( 1 - psi * psi ) * cosPi6plusTheta + ( 2.0 * psi - 1.0 ) * sqrt( std::max( 0., -4.0 * psi + 5.0 * psi * psi + 4.0 * ( 1.0 - psi * psi ) * cosPi6plusTheta * cosPi6plusTheta ) );
    real64 denom = ( 2 * psi - 1.0 ) * ( 2 * psi - 1.0 ) + 4 * ( 1 - psi * psi ) * cosPi6plusTheta * cosPi6plusTheta;

    if( denom > 1e-12 )
    {
      oneOverGamma = num / denom;
    }
  }

  return oneOverGamma;
}

/**
 * @class CeramicDamage
 *
 * Ceramic damage material model.
 */
class CeramicDamage : public ElasticIsotropic
{
public:

  /// @typedef Alias for CeramicDamageUpdates
  using KernelWrapper = CeramicDamageUpdates;

  /**
   * constructor
   * @param[in] name name of the instance in the catalog
   * @param[in] parent the group which contains this instance
   */
  CeramicDamage( string const & name, Group * const parent );

  /**
   * Default Destructor
   */
  virtual ~CeramicDamage() override;


  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  virtual void saveConvergedState() const override;

  /**
   * @name Static Factory Catalog members and functions
   */
  ///@{

  /// string name to use for this class in the catalog
  static constexpr auto m_catalogNameString = "CeramicDamage";

  /**
   * @return A string that is used to register/lookup this class in the registry
   */
  static string catalogName() { return m_catalogNameString; }

  virtual string getCatalogName() const override { return catalogName(); }

  ///@}

  /**
   * Keys for data specified in this class.
   */
  struct viewKeyStruct : public SolidBase::viewKeyStruct
  {
    /// string/key for strength scale value
    static constexpr char const * strengthScaleString() { return "strengthScale"; }

    /// string/key for porosity value
    static constexpr char const * porosityString() { return "porosity"; }

    /// string/key for reference porosity value
    static constexpr char const * referencePorosityString() { return "referencePorosity"; }

    /// string/key for quadrature point damage value
    static constexpr char const * damageString() { return "damage"; }

    /// string/key for quadrature point jacobian value
    static constexpr char const * jacobianString() { return "jacobian"; }

    /// string/key for element/particle length scale
    static constexpr char const * lengthScaleString() { return "lengthScale"; }

    /// string/key for tensile strength
    static constexpr char const * tensileStrengthString() { return "tensileStrength"; }

    /// string/key for compressive strength
    static constexpr char const * compressiveStrengthString() { return "compressiveStrength"; }

    /// string/key for maximum strength
    static constexpr char const * maximumStrengthString() { return "maximumStrength"; }

    /// string/key for crack speed
    static constexpr char const * crackSpeedString() { return "crackSpeed"; }

    /// string/key for third invariant dependence
    static constexpr char const * thirdInvariantDependenceString() { return "thirdInvariantDependence"; }

    /// string/key/ for damged material friction slope
    static constexpr char const * damagedMaterialFrictionSlopeString() { return "damagedMaterialFrictionSlope"; }

    //string/key for element/particle velocityGradient value
    static constexpr char const * velocityGradientString() { return "velocityGradient"; }

    /// string/key for quadrature point plasticStrain value 
    static constexpr char const * plasticStrainString() { return "plasticStrain"; }

    static constexpr char const * refStrainRateString() { return "refStrainRate"; } //added by SG

    static constexpr char const * strainRateString() { return "strainRate"; } // added by SG

    static constexpr char const * instTensileStrengthString() { return "instTensileStrength"; } // added by SG

    static constexpr char const * instCompressiveStrengthString() { return "instCompressiveStrength"; } // added by SG

    static constexpr char const * instStrengthString() { return "instStrength"; } // added by SG

    static constexpr char const * instPressureString() { return "instPressure"; } // added by SG

    static constexpr char const * rateSensitivityString() { return "rateSensitivity"; } //added by SG

    static constexpr char const * m2String() { return "m2"; } //added by SG

    


  };

  /**
   * @brief Create a instantiation of the CeramicDamageUpdate class that refers to the data in this.
   * @return An instantiation of CeramicDamageUpdate.
   */
  CeramicDamageUpdates createKernelUpdates() const//removed const by SG
  { 
    return CeramicDamageUpdates( m_damage,
                                 m_jacobian,
                                 m_lengthScale,
                                 m_strengthScale,
                                 m_porosity,
                                 m_referencePorosity,
                                 m_tensileStrength,
                                 m_compressiveStrength,
                                 m_maximumStrength,
                                 m_crackSpeed,
                                 m_damagedMaterialFrictionSlope,
                                 m_thirdInvariantDependence,
                                 m_velocityGradient,
                                 m_plasticStrain,
                                 m_bulkModulus, //added toView by SG     
                                 m_shearModulus,
                                 m_thermalExpansionCoefficient,
                                 m_newStress,
                                 m_oldStress,
                                 m_density,
                                 m_wavespeed,
                                 m_refStrainRate,
                                 m_strainRate,  // <-- Added this
                                 m_instTensileStrength,  // <-- Added this
                                 m_instCompressiveStrength,  // <-- Added this
                                 m_instStrength,  // <-- Added this
                                 m_instPressure,  // <-- Added this
                                 m_rateSensitivity, // <-- Added this
                                 m_m2, // <-- Added this
                                 m_disableInelasticity );
  }

  /**
   * @brief Construct an update kernel for a derived type.
   * @tparam UPDATE_KERNEL The type of update kernel from the derived type.
   * @tparam PARAMS The parameter pack to hold the constructor parameters for the derived update kernel.
   * @param constructorParams The constructor parameter for the derived type.
   * @return An @p UPDATE_KERNEL object.
   */
  template< typename UPDATE_KERNEL, typename ... PARAMS >
  UPDATE_KERNEL createDerivedKernelUpdates( PARAMS && ... constructorParams )
  {
    return UPDATE_KERNEL( std::forward< PARAMS >( constructorParams )...,
                          m_damage,
                          m_jacobian,
                          m_lengthScale,
                          m_strengthScale,
                          m_porosity,
                          m_referencePorosity,
                          m_tensileStrength,
                          m_compressiveStrength,
                          m_maximumStrength,
                          m_crackSpeed,
                          m_damagedMaterialFrictionSlope,
                          m_thirdInvariantDependence,
                          m_velocityGradient,
                          m_plasticStrain,                         
                          m_bulkModulus,
                          m_shearModulus,
                          m_thermalExpansionCoefficient,
                          m_newStress,
                          m_oldStress,
                          m_density,
                          m_wavespeed,           
                          m_refStrainRate,
                          m_strainRate,  // <-- Added this
                          m_instTensileStrength,  // <-- Added this
                          m_instCompressiveStrength,  // <-- Added this
                          m_instStrength,  // <-- Added this
                          m_instPressure,  // <-- Added this
                          m_rateSensitivity,  // <-- Added this
                          m_m2,  // <-- Added this
                          m_disableInelasticity );
  }


protected:
  virtual void postInputInitialization() override;

  /// State variable: The damage values for each quadrature point
  array2d< real64 > m_damage;

  /// State variable: The jacobian of the deformation
  array2d< real64 > m_jacobian;

  /// Discretization-sized variable: The length scale for each element/particle
  array1d< real64 > m_lengthScale;

  /// Material parameter: The strength scale values
  array1d< real64 > m_strengthScale;

  /// Material parameter: The porosity of particles
  array1d< real64 > m_porosity;

  /// Material parameter: The reference porosity
  array1d< real64 > m_referencePorosity;

  /// Material parameter: The value of unconfined tensile strength
  real64 m_tensileStrength;

  /// Material parameter: The value of unconfined compressive strength
  real64 m_compressiveStrength;

  /// Material parameter: The value of maximum theoretical strength
  real64 m_maximumStrength;

  /// Material parameter: The value of crack speed
  real64 m_crackSpeed;

  /// Material parameter: The damaged material friction slope
  real64 m_damagedMaterialFrictionSlope;

  /// Model parameter: Flag to enable third invariant dependence
  int m_thirdInvariantDependence;

  ///State variable: The velocity gradient for each element/particle
  array3d< real64 > m_velocityGradient;

  ///State variable: The plastic strain values for each quadrature point
  array3d< real64 > m_plasticStrain;
  
  /// Material parameter: The reference strain rate (scalar)
  real64 m_refStrainRate;

   /// Material parameter: The reference strain rate (scalar)
  array2d< real64 > m_strainRate; //added by SG

   /// Material parameter: The reference strain rate (scalar)
  array2d< real64 > m_instTensileStrength; //added by SG

   /// Material parameter: The reference strain rate (scalar)
  array2d< real64 > m_instCompressiveStrength; //added by SG

    /// Material parameter: The reference strain rate (scalar)
  array2d< real64 > m_instStrength; //added by SG

    /// Material parameter: The reference strain rate (scalar)
  array2d< real64 > m_instPressure; //added by SG

/// Material parameter: The rate sensitivity (scalar)
  real64 m_rateSensitivity;

  /// Material parameter: The second slope for pressure dependent strength scaling  
  real64 m_m2; //added by SG


};

} /* namespace constitutive */

} /* namespace geos */

#endif /* GEOS_CONSTITUTIVE_SOLID_CERAMICDAMAGE_HPP_ */
