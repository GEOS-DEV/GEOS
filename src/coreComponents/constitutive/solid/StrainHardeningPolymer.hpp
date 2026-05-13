/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file StrainHardeningPolymer.hpp
 * @brief Corotational strain-hardening / shear-softening polymer model for explicit MPM.
 *
 * @details
 * This material model implements a pressure-independent, J2-type flow-stress
 * law for a strain-hardening polymer with thermal property scaling, plastic
 * shear softening, stretch hardening, and stretch-based failure detection.
 *
 * The model is intended for explicit time integration, especially explicit MPM.
 * It does not provide a consistent implicit tangent. The stress update is
 * formulated in a corotational small-strain frame: finite rotations are handled
 * outside this class by the solver/wrapper, while the constitutive update itself
 * operates on unrotated stress, unrotated strain increment, and unrotated
 * tensorial history variables.
 *
 * Frame convention:
 * - On entry, the wrapper has rotated the beginning-of-step Cauchy stress and
 *   plastic strain into the material/corotational frame.
 * - The strain increment supplied to this update is also expressed in that
 *   corotational frame.
 * - This class updates stress and internal variables in the corotational frame.
 * - On exit, the wrapper rotates the updated stress and plastic strain back to
 *   the spatial frame for storage and force calculation.
 *
 * The elastic predictor uses isotropic linear hypoelasticity in the unrotated
 * frame,
 *
 *   sigma_trial = sigma_n
 *               + K(T) tr(Delta epsilon) I
 *               + 2 G(T) dev(Delta epsilon),
 *
 * where K(T) and G(T) are temperature-scaled bulk and shear moduli. The pressure
 * part of the trial stress is retained during plastic correction. Plastic flow
 * is governed by the von Mises invariant
 *
 *   q = sqrt(3/2 s:s),
 *
 * with radial return in deviatoric stress space when q_trial exceeds the current
 * flow strength.
 *
 * The flow strength is modeled as
 *
 *   sigma_y = sigma_0(T)
 *           + S(T) exp( - (gamma_p / r_1) ^ r_2 )
 *           + H(T) h(lambda),
 *
 * where:
 * - sigma_0(T) is the temperature-scaled base yield strength,
 * - S(T) is the temperature-scaled shear-softening magnitude,
 * - gamma_p is the accumulated magnitude of the plastic strain tensor,
 * - r_1 and r_2 control the decay rate and shape of the shear-softening term,
 * - H(T) is the temperature-scaled stretch-hardening slope,
 * - lambda is the tensile principal-stretch driver, obtained from the right
 *   stretch tensor associated with the deformation gradient, and
 * - h(lambda) is the stretch-hardening function used by the implementation.
 *
 * The shear-softening contribution is initially positive and decays with
 * increasing plastic strain magnitude. The stretch-hardening contribution
 * increases the flow strength as the tensile stretch grows.
 *
 * Temperature dependence is applied through logistic scaling factors of the
 * form
 *
 *   scale(T; T0, A, B) = 1 + A / (1 + exp( B (T - T0) )),
 *
 * unless A is zero, in which case the scale is one. Each temperature-dependent
 * material parameter has its own A, B, and T0 values.
 *
 * Failure is detected from the maximum principal stretch. If
 *
 *   lambda_max > lambda_failure(T),
 *
 * the damage variable is set to one. In MPM/DFG applications, this damage flag
 * may be consumed by the solver to partition material, release stress, or apply
 * other failure logic. The objective rotation of damage-related tensorial state
 * is not performed here.
 *
 * State variables updated by this model include:
 * - plastic strain,
 * - damage,
 * - Jacobian / volume change measure,
 * - current yield strength,
 * - temperature-dependent elastic properties when stored by the parent model.
 *
 * Voigt convention:
 * Strain-like quantities are stored using the GEOS symmetric Voigt convention,
 * with engineering shear components. Tensor rotations of strain-like quantities
 * must therefore convert shear components to tensor form before rotation and
 * convert them back afterward. In normal MPM use, this conversion and rotation
 * are handled by the wrapper.
 *
 * Main assumptions and limitations:
 * - explicit time integration,
 * - small elastic strain in the corotational frame,
 * - finite rotation handled externally,
 * - pressure-independent plastic flow,
 * - no rate dependence or viscoplasticity,
 * - no implicit consistent tangent,
 * - material parameters must remain finite and positive after temperature
 *   scaling.
 */

#ifndef GEOS_CONSTITUTIVE_SOLID_STRAINHARDENINGPOLYMER_HPP_
#define GEOS_CONSTITUTIVE_SOLID_STRAINHARDENINGPOLYMER_HPP_

#include "ElasticIsotropic.hpp"
#include "InvariantDecompositions.hpp"
#include "PropertyConversions.hpp"
#include "SolidModelDiscretizationOpsFullyAnisotropic.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geos
{

namespace constitutive
{

/**
 * @class StrainHardeningPolymerUpdates
 *
 * Class to provide material updates that may be
 * called from a kernel function.
 */
class StrainHardeningPolymerUpdates : public ElasticIsotropicUpdates
{
public:

  /**
   * @brief Constructor
   * @param[in] deformationGradient The ArrayView holding the deformation gradient for each element/particle.
   * @param[in] plasticStrain The ArrayView holding the plastic strain for each quadrature point
   * @param[in] damage The ArrayView holding the damage for each quadrature point.
   * @param[in] temperature The ArrayView holding the temperature for each element/particle.
   * @param[in] jacobian The ArrayView holding the jacobian for each quadrature point.
   * @param[in] yieldStrength The ArrayView holding the current yield strength
   * @param[in] defaultBulkModulus
   * @param[in] bulkModulusA
   * @param[in] bulkModulusB
   * @param[in] bulkModulusT0
   * @param[in] defaultShearModulus
   * @param[in] shearModulusA
   * @param[in] shearModulusB
   * @param[in] shearModulusT0
   * @param[in] defaultYieldStrength
   * @param[in] yieldStrengthA
   * @param[in] yieldStrengthB
   * @param[in] yieldStrengthT0
   * @param[in] strainHardeningSlope The strain hardening slope
   * @param[in] strainHardeningSlope,
   * @param[in] strainHardeningSlopeB
   * @param[in] strainHardeningSlopeT0
   * @param[in] shearSofteningMagnitude  The shear softening magnitude
   * @param[in] shearSofteningMagnitudeA
   * @param[in] shearSofteningMagnitudeB
   * @param[in] shearSofteningMagnitudeT0
   * @param[in] shearSofteningShapeParameter1 The shear softening shape parameter 1
   * @param[in] shearSofteningShapeParameter2 The shear softening shape parameter 2
   * @param[in] maximumStretch The maximum stretch
   * @param[in] maximumStretchA
   * @param[in] maximumStretchB
   * @param[in] maximumStretchT0
   * @param[in] bulkModulus The ArrayView holding the bulk modulus data for each element.
   * @param[in] shearModulus The ArrayView holding the shear modulus data for each element.
   * @param[in] thermalExpansionCoefficient The ArrayView holding the thermal expansion coefficient data for each element.
   * @param[in] newStress The ArrayView holding the new stress data for each quadrature point.
   * @param[in] oldStress The ArrayView holding the old stress data for each quadrature point.
   * @param[in] density
   * @param[in] wavespeed
   * @param[in] disableInelasticity
   */
  StrainHardeningPolymerUpdates( arrayView3d< real64 > const & deformationGradient,
                                 arrayView3d< real64 > const & plasticStrain,
                                 arrayView2d< real64 > const & damage,
                                 arrayView1d< real64 > const & temperature,
                                 arrayView2d< real64 > const & jacobian,
                                 arrayView1d< real64 > const & yieldStrength,
                                 real64 const & defaultBulkModulus,
                                 real64 const & bulkModulusA,
                                 real64 const & bulkModulusB,
                                 real64 const & bulkModulusT0,
                                 real64 const & defaultShearModulus,
                                 real64 const & shearModulusA,
                                 real64 const & shearModulusB,
                                 real64 const & shearModulusT0,
                                 real64 const & defaultYieldStrength,
                                 real64 const & yieldStrengthA,
                                 real64 const & yieldStrengthB,
                                 real64 const & yieldStrengthT0,
                                 real64 const & strainHardeningSlope,
                                 real64 const & strainHardeningSlopeA,
                                 real64 const & strainHardeningSlopeB,
                                 real64 const & strainHardeningSlopeT0,
                                 real64 const & shearSofteningMagnitude,
                                 real64 const & shearSofteningMagnitudeA,
                                 real64 const & shearSofteningMagnitudeB,
                                 real64 const & shearSofteningMagnitudeT0,
                                 real64 const & shearSofteningShapeParameter1,
                                 real64 const & shearSofteningShapeParameter2,
                                 real64 const & maximumStretch,
                                 real64 const & maximumStretchA,
                                 real64 const & maximumStretchB,
                                 real64 const & maximumStretchT0,
                                 arrayView1d< real64 > const & bulkModulus,
                                 arrayView1d< real64 > const & shearModulus,
                                 arrayView1d< real64 const > const & thermalExpansionCoefficient,
                                 arrayView3d< real64, solid::STRESS_USD > const & newStress,
                                 arrayView3d< real64, solid::STRESS_USD > const & oldStress,
                                 arrayView2d< real64 > const & density,
                                 arrayView2d< real64 > const & wavespeed,
                                 bool const & disableInelasticity ):
    ElasticIsotropicUpdates( bulkModulus,
                             shearModulus,
                             thermalExpansionCoefficient,
                             newStress,
                             oldStress,
                             density,
                             wavespeed,
                             disableInelasticity ),
    m_deformationGradient( deformationGradient ),
    m_plasticStrain( plasticStrain ),
    m_damage( damage ),
    m_temperature( temperature ),
    m_jacobian( jacobian ),
    m_yieldStrength( yieldStrength ),
    m_defaultBulkModulus( defaultBulkModulus ),
    m_bulkModulusA( bulkModulusA ),
    m_bulkModulusB( bulkModulusB ),
    m_bulkModulusT0( bulkModulusT0 ),
    m_defaultShearModulus( defaultShearModulus ),
    m_shearModulusA( shearModulusA ),
    m_shearModulusB( shearModulusB ),
    m_shearModulusT0( shearModulusT0 ),
    m_defaultYieldStrength( defaultYieldStrength ),
    m_yieldStrengthA( yieldStrengthA ),
    m_yieldStrengthB( yieldStrengthB ),
    m_yieldStrengthT0( yieldStrengthT0 ),
    m_strainHardeningSlope( strainHardeningSlope ),
    m_strainHardeningSlopeA( strainHardeningSlopeA ),
    m_strainHardeningSlopeB( strainHardeningSlopeB ),
    m_strainHardeningSlopeT0( strainHardeningSlopeT0 ),
    m_shearSofteningMagnitude( shearSofteningMagnitude ),
    m_shearSofteningMagnitudeA( shearSofteningMagnitudeA ),
    m_shearSofteningMagnitudeB( shearSofteningMagnitudeB ),
    m_shearSofteningMagnitudeT0( shearSofteningMagnitudeT0 ),
    m_shearSofteningShapeParameter1( shearSofteningShapeParameter1 ),
    m_shearSofteningShapeParameter2( shearSofteningShapeParameter2 ),
    m_maximumStretch( maximumStretch ),
    m_maximumStretchA( maximumStretchA ),
    m_maximumStretchB( maximumStretchB ),
    m_maximumStretchT0( maximumStretchT0 )
  {}

  /// Default copy constructor
  StrainHardeningPolymerUpdates( StrainHardeningPolymerUpdates const & ) = default;

  /// Default move constructor
  StrainHardeningPolymerUpdates( StrainHardeningPolymerUpdates && ) = default;

  /// Deleted default constructor
  StrainHardeningPolymerUpdates() = delete;

  /// Deleted copy assignment operator
  StrainHardeningPolymerUpdates & operator=( StrainHardeningPolymerUpdates const & ) = delete;

  /// Deleted move assignment operator
  StrainHardeningPolymerUpdates & operator=( StrainHardeningPolymerUpdates && ) =  delete;

  /// Use the uncompressed version of the stiffness bilinear form
  using DiscretizationOps = SolidModelDiscretizationOpsFullyAnisotropic; // TODO: typo in anistropic (fix in DiscOps PR)

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
                                             real64 const ( &beginningRotation )[3][3],
                                             real64 const ( &endRotation )[3][3],
                                             real64 const ( &strainIncrement )[6],
                                             real64 ( &stress )[6] ) const override;

  GEOS_HOST_DEVICE
  void smallStrainUpdateHelper( localIndex const k,
                                localIndex const q,
                                real64 const timeIncrement,
                                real64 const ( &beginningRotation )[3][3],
                                real64 const ( &endRotation )[3][3],
                                real64 const ( &strainIncrement )[6],
                                real64 ( &stress )[6] ) const;

  GEOS_HOST_DEVICE
  void computePlasticStrainIncrement ( localIndex const k,
                                       localIndex const q,
                                       const real64 timeIncrement,
                                       real64 const ( &strainIncrement )[6],
                                       real64 const ( &stressIncrement )[6],
                                       real64 ( &plasticStrainIncrement )[6] ) const;

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  virtual void saveConvergedState( localIndex const k,
                                   localIndex const q ) const override final
  {
    ElasticIsotropicUpdates::saveConvergedState( k, q );
  }

  GEOS_HOST_DEVICE
  real64 thermalScaling( const real64 & T,
                           const real64 & T0,
                           const real64 & A,
                           const real64 & B
                           ) const;

private:
  /// A reference to the ArrayView holding the deformation gradient for each element/particle.
  arrayView3d< real64 > const m_deformationGradient;

  /// A reference to the ArrayView holding the plastic strain for each quadrature point.
  arrayView3d< real64 > const m_plasticStrain;

  /// A reference to the ArrayView holding the damage for each quadrature point.
  arrayView2d< real64 > const m_damage;

  /// A reference to the ArrayView holding the temperature for each quadrature point.
  arrayView1d< real64 > const m_temperature;

  /// A reference to the ArrayView holding the jacobian for each quadrature point.
  arrayView2d< real64 > const m_jacobian;

  /// A reference to the ArrayView holding the yield strength for each element/particle
  arrayView1d< real64 > const m_yieldStrength;

  /// The temperature-independent bulk modulus value
  real64 const m_defaultBulkModulus;
  real64 const m_bulkModulusA;
  real64 const m_bulkModulusB;
  real64 const m_bulkModulusT0;

  /// The temperature-independent shear modulus value
  real64 const m_defaultShearModulus;
  real64 const m_shearModulusA;
  real64 const m_shearModulusB;
  real64 const m_shearModulusT0;

  // Yield strength before softening/hardening
  real64 const m_defaultYieldStrength;
  real64 const m_yieldStrengthA;
  real64 const m_yieldStrengthB;
  real64 const m_yieldStrengthT0;

  /// The strain hardening slope
  real64 const m_strainHardeningSlope;
  real64 const m_strainHardeningSlopeA;
  real64 const m_strainHardeningSlopeB;
  real64 const m_strainHardeningSlopeT0;

  /// The shear softening magnitude
  real64 m_shearSofteningMagnitude;
  real64 m_shearSofteningMagnitudeA;
  real64 m_shearSofteningMagnitudeB;
  real64 m_shearSofteningMagnitudeT0;

  /// The shear softening shape parameter 1
  real64 m_shearSofteningShapeParameter1;

  /// The shear softening shape parameter 2
  real64 m_shearSofteningShapeParameter2;

  /// The compressive strength
  real64 const m_maximumStretch;
  real64 const m_maximumStretchA;
  real64 const m_maximumStretchB;
  real64 const m_maximumStretchT0;
};


GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void StrainHardeningPolymerUpdates::smallStrainUpdate( localIndex const k,
                                                       localIndex const q,
                                                       real64 const & timeIncrement,
                                                       real64 const ( &strainIncrement )[6],
                                                       real64 ( & stress )[6],
                                                       real64 ( & stiffness )[6][6] ) const
{
  // CC: we don't call this for the MPM solver but other solvers may want to use this
  // Need to resolve if rotation matrix is always passed and update accordingly
  GEOS_UNUSED_VAR( k );
  GEOS_UNUSED_VAR( q );
  GEOS_UNUSED_VAR( timeIncrement );
  GEOS_UNUSED_VAR( strainIncrement );
  GEOS_UNUSED_VAR( stress );
  GEOS_UNUSED_VAR( stiffness );
  GEOS_ERROR( "smallStrainUpdate not implemented for StrainHardeningPolymer" );

  // // elastic predictor (assume strainIncrement is all elastic)
  // ElasticIsotropicUpdates::smallStrainUpdate( k,
  //                                             q,
  //                                             timeIncrement,
  //                                             strainIncrement,
  //                                             stress,
  //                                             stiffness );
  // m_jacobian[k][q] *= exp( strainIncrement[0] + strainIncrement[1] + strainIncrement[2] );

  // if( m_disableInelasticity )
  // {
  //   return;
  // }

  // // call the constitutive model
  // StrainHardeningPolymerUpdates::smallStrainUpdateHelper( k,
  //                                                         q,
  //                                                         timeIncrement,
  //                                                         strainIncrement,
  //                                                         stress );

  // // It doesn't make sense to modify stiffness with this model

  // // save new stress and return
  // saveStress( k, q, stress );
  // return;
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void StrainHardeningPolymerUpdates::smallStrainUpdate( localIndex const k,
                                                       localIndex const q,
                                                       real64 const & timeIncrement,
                                                       real64 const ( &strainIncrement )[6],
                                                       real64 ( & stress )[6],
                                                       DiscretizationOps & stiffness ) const
{
  // CC: we don't call this for the MPM solver but other solvers may want to use this
  // Need to resolve if rotation matrix is always passed and update accordingly
  GEOS_UNUSED_VAR( k );
  GEOS_UNUSED_VAR( q );
  GEOS_UNUSED_VAR( timeIncrement );
  GEOS_UNUSED_VAR( strainIncrement );
  GEOS_UNUSED_VAR( stress );
  GEOS_UNUSED_VAR( stiffness );
  GEOS_ERROR( "smallStrainUpdate not implemented for StrainHardeningPolymer" );

  // smallStrainUpdate( k,
  //                    q,
  //                    timeIncrement,
  //                    strainIncrement,
  //                    stress,
  //                    stiffness.m_c );
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void StrainHardeningPolymerUpdates::smallStrainUpdate_StressOnly( localIndex const k,
                                                                  localIndex const q,
                                                                  real64 const & timeIncrement,
                                                                  real64 const ( &strainIncrement )[6],
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
void StrainHardeningPolymerUpdates::smallStrainUpdate_StressOnly( localIndex const k,
                                                                  localIndex const q,
                                                                  real64 const & timeIncrement,
                                                                  real64 const ( &beginningRotation )[3][3],
                                                                  real64 const ( &endRotation )[3][3],
                                                                  real64 const ( &strainIncrement )[6],
                                                                  real64 ( & stress )[6] ) const
{
  // elastic predictor "trialStress" (assume strainIncrement is all elastic)
  // using current definitions of m_bulkModulus[k] and m_shearModulus[k]

  real64 bulkModulusScale = StrainHardeningPolymerUpdates::thermalScaling( m_temperature[k], m_bulkModulusT0, m_bulkModulusA, m_bulkModulusB );      // This
                                                                                                                                            // will
                                                                                                                                            // actually
                                                                                                                                            // be
                                                                                                                                            // some
                                                                                                                                            // function:
                                                                                                                                            // 
                                                                                                                                            // 
                                                                                                                                            // m_bulkModulus[k]
                                                                                                                                            // =
                                                                                                                                            // m_defaultBulkModulus
                                                                                                                                            // +
                                                                                                                                            // A*f(m_temperature[k]),
                                                                                                                                            // etc.
  m_bulkModulus[k] = m_defaultBulkModulus * bulkModulusScale;

  real64 shearModulusScale = StrainHardeningPolymerUpdates::thermalScaling( m_temperature[k], m_shearModulusT0, m_shearModulusA, m_shearModulusB );      // This
                                                                                                                                        // will
                                                                                                                                        // actually
                                                                                                                                        // be
                                                                                                                                        // some
                                                                                                                                        // function:
                                                                                                                                        //  
                                                                                                                                        // m_bulkModulus[k]
                                                                                                                                        // =
                                                                                                                                        // m_defaultBulkModulus
                                                                                                                                        // +
                                                                                                                                        // A*f(m_temperature[k]),
                                                                                                                                        // etc.
  m_shearModulus[k] = m_defaultShearModulus *  shearModulusScale; // This will actually be some function of m_temperature[k]

  ElasticIsotropicUpdates::smallStrainUpdate_StressOnly( k,
                                                         q,
                                                         timeIncrement,
                                                         strainIncrement,
                                                         stress );  // "stress" is overwritten to be trial stress
  m_jacobian[k][q] *= exp( strainIncrement[0] + strainIncrement[1] + strainIncrement[2] );

  if( m_disableInelasticity )
  {
    return;
  }

  // call the constitutive model
  StrainHardeningPolymerUpdates::smallStrainUpdateHelper( k,
                                                          q,
                                                          timeIncrement,
                                                          beginningRotation,
                                                          endRotation,
                                                          strainIncrement,
                                                          stress ); // This will update "stress" from trialStress to end-of-step stress

  // save new stress and return
  saveStress( k, q, stress );
  return;
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void StrainHardeningPolymerUpdates::smallStrainUpdateHelper( localIndex const k,
                                                             localIndex const q,
                                                             real64 const timeIncrement,
                                                             real64 const ( &beginningRotation )[3][3],
                                                             real64 const ( &endRotation )[3][3],
                                                             real64 const ( &strainIncrement )[6],
                                                             real64 ( & stress )[6] ) const // this is the trial stress, and will be
                                                                                            // overwritten.
{
  // Store trial stress for computing the plastic strain increment.
  real64 trialStress[6] = { };
  LvArray::tensorOps::copy< 6 >( trialStress, stress );

  // decompose into mean (P) and von Mises (Q) stress invariants
  real64 trialP;
  real64 trialQ;
  real64 deviator[6] = { };
  twoInvariant::stressDecomposition( trialStress,
                                     trialP,
                                     trialQ,
                                     deviator );

  // CC: model needs the unrotated deformation gradient
  // Right stretch tensor
  real64 rotationTranspose[3][3];
  LvArray::tensorOps::transpose< 3, 3 >( rotationTranspose, beginningRotation );

  real64 oldPlasticStrain[6] = { };
  LvArray::tensorOps::copy< 6 >( oldPlasticStrain, m_plasticStrain[k][q] );
  oldPlasticStrain[3] *= 0.5;
  oldPlasticStrain[4] *= 0.5;
  oldPlasticStrain[5] *= 0.5;

  real64 unrotatedOldPlasticStrain[6] = { };
  LvArray::tensorOps::Rij_eq_AikSymBklAjl< 3 >( unrotatedOldPlasticStrain, rotationTranspose, oldPlasticStrain );

  unrotatedOldPlasticStrain[3] *= 2.0;
  unrotatedOldPlasticStrain[4] *= 2.0;
  unrotatedOldPlasticStrain[5] *= 2.0;

  real64 unrotatedDeformationGradient[3][3] = { };
  LvArray::tensorOps::Rij_eq_AikBkj< 3, 3, 3 >( unrotatedDeformationGradient, rotationTranspose, m_deformationGradient[k] );

  real64 U[6] = { };
  LvArray::tensorOps::denseToSymmetric< 3 >( U, unrotatedDeformationGradient );

  real64 stretch[3] = { };
  real64 eigenVectors[3][3] = { };
  LvArray::tensorOps::symEigenvectors< 3 >( stretch, eigenVectors, U );

  // Find the largest eigenvalues and compare to max allowable failure stretch (which is temperature dependent)

  // Stretch to failure at current temperature:
  real64 failureStretch = m_maximumStretch * StrainHardeningPolymerUpdates::thermalScaling( m_temperature[k], m_maximumStretchT0, m_maximumStretchA, m_maximumStretchB );

  // Maximum principal stretch:
  real64 maximumStretch = 0.0;
  for( localIndex i = 0; i < 3; ++i )
  {
    maximumStretch = LvArray::math::max( stretch[i], maximumStretch );
  }

  if( maximumStretch > failureStretch )
  {
    m_damage[k][q] = 1.0;
  }

  // Return to yield surface requires iterative solution
  // Implemented fixed points, however a newton solver may be more efficient and applicable
  real64 relTol = 1e-10;   // CC: need to experiment with these for the best options
  int maxEvals = 100;   // Same as above

  // In initialization, yieldStrength is set to defaultYieldStrength, but we will generally want it to be modified by temp
  // Here we would update the m_bulkModulus[k] and m_shearModulus[k] with temperature dependent values:
  // These will be input paramters:
  real64 oldYieldStrength = m_yieldStrength[k];
  real64 yieldStrengthIter = oldYieldStrength;

  // Compute change in yield strength: yieldStrength = m_initialYield + plasticSoftening + stretchHardening;
  // Where each of the 3 terms on the right hace parameters modified by temperature with 3 parameters.
  real64 yield0 = m_defaultYieldStrength * StrainHardeningPolymerUpdates::thermalScaling( m_temperature[k], m_yieldStrengthT0, m_yieldStrengthA, m_yieldStrengthB );
  real64 strainHardeningSlope = m_strainHardeningSlope *
                                StrainHardeningPolymerUpdates::thermalScaling( m_temperature[k], m_strainHardeningSlopeT0, m_strainHardeningSlopeA, m_strainHardeningSlopeB );
  real64 shearSofteningMagnitude = m_shearSofteningMagnitude * StrainHardeningPolymerUpdates::thermalScaling( m_temperature[k], m_shearSofteningMagnitudeT0, m_shearSofteningMagnitudeA,
                                                                                                                m_shearSofteningMagnitudeB );

  real64 unrotatedTempPlasticStrain[6] = { };
  real64 plasticStrainIncrement[6] = { };

  // Fixed-point iteration to find plastic strain and consistent return to updated yield surface
  for( int iter=0; iter < maxEvals; ++iter )
  {
    LvArray::tensorOps::copy< 6 >( unrotatedTempPlasticStrain, unrotatedOldPlasticStrain );
    LvArray::tensorOps::add< 6 >( unrotatedTempPlasticStrain, plasticStrainIncrement );

    // Compute magnitude of plastic strain tensor
    real64 gamma_p = 0.0;
    for( int i = 0; i < 6; i++ )
    {
      gamma_p += 0.5*( 1 + (i < 3) ) * unrotatedTempPlasticStrain[i] * unrotatedTempPlasticStrain[i];
    }

    gamma_p = LvArray::math::sqrt( gamma_p );

    // This term starts at value r0 and decays with plastic shear strain to give plastic softening.
    // Put in a check to prevent roundoff error.
    real64 gamma_by_r1_to_r2 = LvArray::math::pow( gamma_p / m_shearSofteningShapeParameter1, m_shearSofteningShapeParameter2 );

    // Shear Softening:
    real64 plasticSoftening = shearSofteningMagnitude * LvArray::math::exp( LvArray::math::max( -1.0 * gamma_by_r1_to_r2, -16.0 ) );

    // Stretch hardening (only in tension)
    real64 const lambdaBar = LvArray::math::max( maximumStretch, 1.0 );
    real64 stretchHardening = strainHardeningSlope * ( lambdaBar * lambdaBar - 1.0 / lambdaBar );

    // Flow stress after temp, hardening, and softening modifications
    yieldStrengthIter = yield0 + plasticSoftening + stretchHardening;
    real64 newPlasticStrain[6] = { };

    // check yield function
    if( trialQ > yieldStrengthIter || iter > 0 )
    {
      // In case there is hardening, don't allow dev stress to exceed trial value or be negative (for the iter > 0 cases)
      real64 const returnedQ = LvArray::math::min( trialQ, LvArray::math::max( yieldStrengthIter, 0.0 ) );

      // re-construct stress = P*eye + sqrt(2/3)*Q*nhat
      real64 stressTemp[6] = {};
      twoInvariant::stressRecomposition( trialP,
                                         returnedQ,
                                         deviator,
                                         stressTemp );

      // Increment plastic strain
      real64 stressIncrement[6] = { };
      real64 oldStress[6] = { 0 };
      LvArray::tensorOps::copy< 6 >( oldStress, m_oldStress[k][q] );

      // We will compute:  plasticStrainIncrement = strainIncrement - elasticStrainIncrement = strainIncrement - C^inv*( newStress - oldStress )
      LvArray::tensorOps::copy< 6 >( stressIncrement, stressTemp );
      LvArray::tensorOps::subtract< 6 >( stressIncrement, oldStress );

      // increment plastic strain
      computePlasticStrainIncrement( k,
                                     q,
                                     timeIncrement,
                                     strainIncrement,
                                     stressIncrement,
                                     plasticStrainIncrement );

      real64 unrotatedNewPlasticStrain[6] = { };
      LvArray::tensorOps::copy< 6 >( unrotatedNewPlasticStrain, unrotatedOldPlasticStrain );
      LvArray::tensorOps::add< 6 >( unrotatedNewPlasticStrain, plasticStrainIncrement );

      real64 const yieldScale = LvArray::math::max( 1.0, LvArray::math::abs( yieldStrengthIter ) );
      if( LvArray::math::abs( yieldStrengthIter - oldYieldStrength ) <= relTol * yieldScale )
      {
        unrotatedNewPlasticStrain[3] *= 0.5;
        unrotatedNewPlasticStrain[4] *= 0.5;
        unrotatedNewPlasticStrain[5] *= 0.5;

        
        LvArray::tensorOps::Rij_eq_AikSymBklAjl< 3 >( newPlasticStrain, endRotation, unrotatedNewPlasticStrain );
        newPlasticStrain[3] *= 2.0;
        newPlasticStrain[4] *= 2.0;
        newPlasticStrain[5] *= 2.0;

        LvArray::tensorOps::copy< 6 >( stress, stressTemp );
        return;
      }
      else
      {
        oldYieldStrength = yieldStrengthIter;
      }
    }
    else
    {
      return;
    }

    // Store converged values:
    m_yieldStrength[k] = yieldStrengthIter;
    LvArray::tensorOps::copy< 6 >( m_plasticStrain[k][q], newPlasticStrain );

  }

  GEOS_ERROR( "Plastic strain of StrainHardeningPolymer model did not converge within max evals." );
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void StrainHardeningPolymerUpdates::computePlasticStrainIncrement ( localIndex const k,
                                                                    localIndex const q,
                                                                    const real64 timeIncrement,
                                                                    real64 const ( &strainIncrement )[6],
                                                                    real64 const ( &stressIncrement )[6],
                                                                    real64 ( & plasticStrainIncrement )[6] ) const
{
  GEOS_UNUSED_VAR( q );
  GEOS_UNUSED_VAR( timeIncrement );

  // For hypo-elastic models we compute the increment in plastic strain assuming
  // for some increment in total strain and stress and elastic properties.

  // Isotroptic-deviatoric decomposition;
  real64 trialP;
  real64 trialQ;
  real64 stressIncrementDeviator[6] = {};
  twoInvariant::stressDecomposition( stressIncrement,
                                     trialP,
                                     trialQ,
                                     stressIncrementDeviator );

  real64 stressIncrementIsostatic[6] = {};
  stressIncrementIsostatic[0] = trialP;
  stressIncrementIsostatic[1] = trialP;
  stressIncrementIsostatic[2] = trialP;

  // For damage or softening it there may be cases where bulk or shear are approx 0,
  // so we need to be careful that we compute this
  real64 elasticStrainIncrement[6] = {};
  for( int i = 0; i < 6; ++i )
  {
    if( m_bulkModulus[k] > 1.0e-12 )
    {
      // CC: off diagonal elements need x2 for strain
      elasticStrainIncrement[i] += ( 1 + (i >= 3) ) * stressIncrementIsostatic[i] * 1.0/3.0/m_bulkModulus[k];
    }
    if( m_shearModulus[k] > 1.0e-12 )
    {
      elasticStrainIncrement[i] += ( 1 + (i >= 3) ) * LvArray::math::sqrt( 2/3 ) * trialQ * stressIncrementDeviator[i] * 1.0/2.0/m_shearModulus[k];
    }
  }

  LvArray::tensorOps::copy< 6 >( plasticStrainIncrement, strainIncrement );
  LvArray::tensorOps::subtract< 6 >( plasticStrainIncrement, elasticStrainIncrement );
}

// Compute thermal scaling function.  
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 StrainHardeningPolymerUpdates::thermalScaling( const real64 & T,
                                                        const real64 & T0,
                                                        const real64 & A,
                                                        const real64 & B ) const
{
  // This is an empirical scaling function used to modify parameters based on available
  // data.  It is applied differently for stretch hardeng, shear softening, and
  // elasticity parameters.  Be advised that it is not a thermal softening model
  // with the expected response that scaling=1 when T=T0.  The model is disabled
  // by default and when A = 0.
  //
  // TODO: Reformulate this so we can have input values unchanged when T=T0
  // but still fit data for temp dependence of various parameters.  
  if( A > __DBL_EPSILON__ )
  {
    return 1.0 + A / (1.0 + LvArray::math::exp( B * (T-T0) ) );
  }
  else
  {
    return 1.0;
  }
}

/**
 * @class StrainHardeningPolymer
 *
 * Strain hardening polymer material model.
 */
class StrainHardeningPolymer : public ElasticIsotropic
{
public:

  /// @typedef Alias for StrainHardeningPolymerUpdates
  using KernelWrapper = StrainHardeningPolymerUpdates;

  /**
   * constructor
   * @param[in] name name of the instance in the catalog
   * @param[in] parent the group which contains this instance
   */
  StrainHardeningPolymer( string const & name, Group * const parent );

  /**
   * Default Destructor
   */
  virtual ~StrainHardeningPolymer() override;


  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  virtual void saveConvergedState() const override;

  /**
   * @name Static Factory Catalog members and functions
   */
  ///@{

  /// string name to use for this class in the catalog
  static constexpr auto m_catalogNameString = "StrainHardeningPolymer";

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
    // string/key for element/particle deformation gradient value
    static constexpr char const * deformationGradientString() { return "deformationGradient"; }

    // string/key for quadrature point plastic strain value
    static constexpr char const * plasticStrainString() { return "plasticStrain"; }

    /// string/key for quadrature point damage value
    static constexpr char const * damageString() { return "damage"; }

    /// string/key for quadrature point temperature value
    static constexpr char const * temperatureString() { return "temperature"; }

    /// string/key for quadrature point jacobian value
    static constexpr char const * jacobianString() { return "jacobian"; }

    /// string/key for yield strength
    static constexpr char const * yieldStrengthString() { return "yieldStrength"; }

    /// string/key for strain hardening slope
    static constexpr char const * strainHardeningSlopeString() { return "strainHardeningSlope"; }
    static constexpr char const * strainHardeningSlopeAString() { return "strainHardeningSlopeA"; }
    static constexpr char const * strainHardeningSlopeBString() { return "strainHardeningSlopeB"; }
    static constexpr char const * strainHardeningSlopeT0String() { return "strainHardeningSlopeT0"; }

    /// string/key for shear softening magnitude
    static constexpr char const * shearSofteningMagnitudeString() { return "shearSofteningMagnitude"; }
    static constexpr char const * shearSofteningMagnitudeAString() { return "shearSofteningMagnitudeA"; }
    static constexpr char const * shearSofteningMagnitudeBString() { return "shearSofteningMagnitudeB"; }
    static constexpr char const * shearSofteningMagnitudeT0String() { return "shearSofteningMagnitudeT0"; }

    /// string/key for shear softening shape parameter 1
    static constexpr char const * shearSofteningShapeParameter1String() { return "shearSofteningShapeParameter1"; }

    /// string/key for shear softening shape parameter 2
    static constexpr char const * shearSofteningShapeParameter2String() { return "shearSofteningShapeParameter2"; }

    /// string/key for default bulk modulus (temp independent value)
    static constexpr char const * bulkModulusAString() { return "bulkModulusA"; }
    static constexpr char const * bulkModulusBString() { return "bulkModulusB"; }
    static constexpr char const * bulkModulusT0String() { return "bulkModulusT0"; }

    /// string/key for default shear modulus (temp independent value)
    static constexpr char const * shearModulusAString() { return "shearModulusA"; }
    static constexpr char const * shearModulusBString() { return "shearModulusB"; }
    static constexpr char const * shearModulusT0String() { return "shearModulusT0"; }

    /// string/key for default yield strength
    static constexpr char const * defaultYieldStrengthString() { return "defaultYieldStrength"; }
    static constexpr char const * yieldStrengthAString() { return "yieldStrengthA"; }
    static constexpr char const * yieldStrengthBString() { return "yieldStrengthB"; }
    static constexpr char const * yieldStrengthT0String() { return "yieldStrengthT0"; }

    /// string/key for maximum stretch
    static constexpr char const * maximumStretchString() { return "maximumStretch"; }
    static constexpr char const * maximumStretchAString() { return "maximumStretchA"; }
    static constexpr char const * maximumStretchBString() { return "maximumStretchB"; }
    static constexpr char const * maximumStretchT0String() { return "maximumStretchT0"; }
  };

  /**
   * @brief Create a instantiation of the StrainHardeningPolymerUpdates class that refers to the data in this.
   * @return An instantiation of StrainHardeningPolymerUpdates.
   */
  StrainHardeningPolymerUpdates createKernelUpdates() const
  {
    return StrainHardeningPolymerUpdates( m_deformationGradient,
                                          m_plasticStrain,
                                          m_damage,
                                          m_temperature,
                                          m_jacobian,
                                          m_yieldStrength,
                                          m_defaultBulkModulus,
                                          m_bulkModulusA,
                                          m_bulkModulusB,
                                          m_bulkModulusT0,
                                          m_defaultShearModulus,
                                          m_shearModulusA,
                                          m_shearModulusB,
                                          m_shearModulusT0,
                                          m_defaultYieldStrength,
                                          m_yieldStrengthA,
                                          m_yieldStrengthB,
                                          m_yieldStrengthT0,
                                          m_strainHardeningSlope,
                                          m_strainHardeningSlopeA,
                                          m_strainHardeningSlopeB,
                                          m_strainHardeningSlopeT0,
                                          m_shearSofteningMagnitude,
                                          m_shearSofteningMagnitudeA,
                                          m_shearSofteningMagnitudeB,
                                          m_shearSofteningMagnitudeT0,
                                          m_shearSofteningShapeParameter1,
                                          m_shearSofteningShapeParameter2,
                                          m_maximumStretch,
                                          m_maximumStretchA,
                                          m_maximumStretchB,
                                          m_maximumStretchT0,
                                          m_bulkModulus,
                                          m_shearModulus,
                                          m_thermalExpansionCoefficient,
                                          m_newStress,
                                          m_oldStress,
                                          m_density,
                                          m_wavespeed,
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
                          m_deformationGradient,
                          m_plasticStrain,
                          m_damage,
                          m_temperature,
                          m_jacobian,
                          m_yieldStrength,
                          m_defaultBulkModulus,
                          m_bulkModulusA,
                          m_bulkModulusB,
                          m_bulkModulusT0,
                          m_defaultShearModulus,
                          m_shearModulusA,
                          m_shearModulusB,
                          m_shearModulusT0,
                          m_defaultYieldStrength,
                          m_yieldStrengthA,
                          m_yieldStrengthB,
                          m_yieldStrengthT0,
                          m_strainHardeningSlope,
                          m_strainHardeningSlopeA,
                          m_strainHardeningSlopeB,
                          m_strainHardeningSlopeT0,
                          m_shearSofteningMagnitude,
                          m_shearSofteningMagnitudeA,
                          m_shearSofteningMagnitudeB,
                          m_shearSofteningMagnitudeT0,
                          m_shearSofteningShapeParameter1,
                          m_shearSofteningShapeParameter2,
                          m_maximumStretch,
                          m_maximumStretchA,
                          m_maximumStretchB,
                          m_maximumStretchT0,
                          m_bulkModulus,
                          m_shearModulus,
                          m_thermalExpansionCoefficient,
                          m_newStress,
                          m_oldStress,
                          m_density,
                          m_wavespeed,
                          m_disableInelasticity );
  }


protected:
  virtual void postInputInitialization() override;
  // These are variables that aren't in the parent elasticIsotropic class, or that we wish
  // to over-write.  The m_bulkModulus array and m_defaultBulkModulus parameterm etc. should not
  // be here, since they are defined in the parent class.

  /// State variable: The bulkModulus values for each quadrature point
  // array1d< real64 > m_bulkModulus;

  /// State variable: The shear Modulus values for each quadrature point
  // array1d< real64 > m_shearModulus;

  /// State variable: The deformation gradient values for each element/particle.
  array3d< real64 > m_deformationGradient;

  /// State variable: The plastic strain values for each quadrature point
  array3d< real64 > m_plasticStrain;

  /// State variable: The damage values for each quadrature point
  array2d< real64 > m_damage;

  /// State variable: The temperature values for each element/particle
  array1d< real64 > m_temperature;

  /// State variable: The jacobian of the deformation gradient for each quadrature point
  array2d< real64 > m_jacobian;

  /// State variable: The yield strength
  array1d< real64 > m_yieldStrength;

  /// The default value of the first Lame constant for any new allocations.
  //real64 m_defaultBulkModulus;
  real64 m_bulkModulusA;
  real64 m_bulkModulusB;
  real64 m_bulkModulusT0;

  /// The default value of the second Lame constant for any new allocations.
  //real64 m_defaultShearModulus;
  real64 m_shearModulusA;
  real64 m_shearModulusB;
  real64 m_shearModulusT0;

  /// Material parameter: The value of default yield strength
  real64 m_defaultYieldStrength;
  real64 m_yieldStrengthA;
  real64 m_yieldStrengthB;
  real64 m_yieldStrengthT0;

  /// Material parameter: The value of strain hardening slope
  real64 m_strainHardeningSlope;
  real64 m_strainHardeningSlopeA;
  real64 m_strainHardeningSlopeB;
  real64 m_strainHardeningSlopeT0;

  /// Material parameter: The value of shear softening magnitude
  real64 m_shearSofteningMagnitude;
  real64 m_shearSofteningMagnitudeA;
  real64 m_shearSofteningMagnitudeB;
  real64 m_shearSofteningMagnitudeT0;

  /// Material parameter: The value of shear softening shape parameter 1
  real64 m_shearSofteningShapeParameter1;

  /// Material parameter: The value of shear softening shape parameter 2
  real64 m_shearSofteningShapeParameter2;

  /// Material parameter: The value of maximum theoretical strength
  real64 m_maximumStretch;
  real64 m_maximumStretchA;
  real64 m_maximumStretchB;
  real64 m_maximumStretchT0;
};

} /* namespace constitutive */

} /* namespace geos */

#endif /* GEOSX_CONSTITUTIVE_SOLID_STRAINHARDENINGPOLYMER_HPP_ */
