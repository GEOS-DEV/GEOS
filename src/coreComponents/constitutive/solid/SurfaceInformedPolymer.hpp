/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SurfaceInformedPolymer.hpp
 * @brief Reference-temperature polymer model with softening, stretch hardening, and pressure asymmetry.
 *
 * @details
 * SurfaceInformedPolymer is an explicit-MPM continuum polymer model that keeps the same integration
 * philosophy as StrainHardeningPolymer: the stress update is hypoelastic-plastic in a corotational
 * small-strain frame, finite rotation is handled by the MPM wrapper, and the model provides no
 * implicit tangent.  The new form is intended to make continuum and thin-film cohesive-zone polymer
 * models share one set of material functions.
 *
 * Material parameters are interpreted at a reference transition temperature, supplied through
 * glassTransitionTemperature.  A positive temperature scale S_T(T), normalized so that S_T(Tg)=1,
 * scales elastic moduli, yield strength, and the softening magnitude.  The stretch-hardening slope
 * has an independent exponent, H(T)=H_g S_T(T)^p_H, because inverted stress-strain surfaces show
 * that post-yield hardening often evolves more weakly with temperature than the initial flow stress.
 *
 * The scalar flow strength is
 *
 *   sigma_f^0 = Y(T,Xc) + S_soft(T) exp[-(kappa/r1)^r2]
 *             + H_g S_T(T)^p_H [ lambda_eff^2 - lambda_eff^{-1} ],
 *
 * where kappa is the accumulated equivalent plastic strain-like state and lambda_eff is a chain-stretch
 * measure based on the larger of the tensile maximum principal stretch and the isochoric stretch.  The
 * term S_soft initially raises the flow stress and then decays with plastic strain, representing
 * post-yield shear softening.  The stretch term supplies large-strain hardening when chains are
 * extended by tension, constrained compression, or shear.
 *
 * The yield envelope uses a weak pressure-sensitive correction,
 *
 *   Phi = q + eta(T) p_eff - sigma_f^0 <= 0,
 *
 * where q is the von Mises equivalent stress and p_eff is the mean stress with tensile sign
 * convention after clipping the compressive side to an optional cap.  The update deliberately uses a
 * non-associated radial return: the deviatoric stress is returned to the pressure-adjusted scalar
 * surface while the trial mean stress is retained.  This makes the term a bounded yield-strength
 * asymmetry for tension and compression without adding volumetric plastic flow.
 *
 * Crystallinity, if supplied, scales elastic stiffness and yield strength through bounded linear
 * factors that can be smoothly activated near the transition temperature.  Setting the crystallinity
 * coefficients to zero recovers an amorphous or crystallinity-insensitive model.
 */

#ifndef GEOS_CONSTITUTIVE_SOLID_SURFACEINFORMEDPOLYMER_HPP_
#define GEOS_CONSTITUTIVE_SOLID_SURFACEINFORMEDPOLYMER_HPP_

#include "ElasticIsotropic.hpp"
#include "InvariantDecompositions.hpp"
#include "PropertyConversions.hpp"
#include "SolidModelDiscretizationOpsFullyAnisotropic.hpp"
#include "SurfaceInformedPolymerHelpers.hpp"
#include "LvArray/src/tensorOps.hpp"

#include <cfloat>

namespace geos
{
namespace constitutive
{

class SurfaceInformedPolymerUpdates : public ElasticIsotropicUpdates
{
public:
  SurfaceInformedPolymerUpdates( arrayView3d< real64 > const & deformationGradient,
                                 arrayView3d< real64 > const & plasticStrain,
                                 arrayView2d< real64 > const & equivalentPlasticStrain,
                                 arrayView2d< real64 > const & damage,
                                 arrayView1d< real64 > const & temperature,
                                 arrayView2d< real64 > const & jacobian,
                                 arrayView1d< real64 > const & yieldStrength,
                                 real64 const & defaultBulkModulus,
                                 real64 const & defaultShearModulus,
                                 real64 const & defaultYieldStrength,
                                 real64 const & shearSofteningMagnitude,
                                 real64 const & shearSofteningShapeParameter1,
                                 real64 const & shearSofteningShapeParameter2,
                                 real64 const & strainHardeningSlope,
                                 real64 const & hardeningScaleExponent,
                                 real64 const & maximumStretch,
                                 real64 const & fractureStretchLambdaMin,
                                 real64 const & fractureStretchLambda0,
                                 real64 const & fractureStretchT0,
                                 real64 const & fractureStretchTemperatureScale,
                                 real64 const & glassTransitionTemperature,
                                 real64 const & temperatureColdSlope,
                                 real64 const & temperatureHotSlope,
                                 real64 const & temperatureTransitionMagnitude,
                                 real64 const & temperatureTransitionWidth,
                                 real64 const & crystallinity,
                                 real64 const & referenceCrystallinity,
                                 real64 const & crystallinityTransitionWidth,
                                 real64 const & elasticCrystallinityCoeff,
                                 real64 const & yieldStrengthCrystallinityCoeff,
                                 real64 const & pressureAsymmetryAmplitude,
                                 real64 const & pressureAsymmetryWidth,
                                 real64 const & compressivePressureStrengtheningCap,
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
    m_equivalentPlasticStrain( equivalentPlasticStrain ),
    m_damage( damage ),
    m_temperature( temperature ),
    m_jacobian( jacobian ),
    m_yieldStrength( yieldStrength ),
    m_defaultBulkModulus( defaultBulkModulus ),
    m_defaultShearModulus( defaultShearModulus ),
    m_defaultYieldStrength( defaultYieldStrength ),
    m_shearSofteningMagnitude( shearSofteningMagnitude ),
    m_shearSofteningShapeParameter1( shearSofteningShapeParameter1 ),
    m_shearSofteningShapeParameter2( shearSofteningShapeParameter2 ),
    m_strainHardeningSlope( strainHardeningSlope ),
    m_hardeningScaleExponent( hardeningScaleExponent ),
    m_maximumStretch( maximumStretch ),
    m_fractureStretchLambdaMin( fractureStretchLambdaMin ),
    m_fractureStretchLambda0( fractureStretchLambda0 ),
    m_fractureStretchT0( fractureStretchT0 ),
    m_fractureStretchTemperatureScale( fractureStretchTemperatureScale ),
    m_glassTransitionTemperature( glassTransitionTemperature ),
    m_temperatureColdSlope( temperatureColdSlope ),
    m_temperatureHotSlope( temperatureHotSlope ),
    m_temperatureTransitionMagnitude( temperatureTransitionMagnitude ),
    m_temperatureTransitionWidth( temperatureTransitionWidth ),
    m_crystallinity( crystallinity ),
    m_referenceCrystallinity( referenceCrystallinity ),
    m_crystallinityTransitionWidth( crystallinityTransitionWidth ),
    m_elasticCrystallinityCoeff( elasticCrystallinityCoeff ),
    m_yieldStrengthCrystallinityCoeff( yieldStrengthCrystallinityCoeff ),
    m_pressureAsymmetryAmplitude( pressureAsymmetryAmplitude ),
    m_pressureAsymmetryWidth( pressureAsymmetryWidth ),
    m_compressivePressureStrengtheningCap( compressivePressureStrengtheningCap )
  {}

  SurfaceInformedPolymerUpdates() = delete;
  SurfaceInformedPolymerUpdates( SurfaceInformedPolymerUpdates const & ) = default;
  SurfaceInformedPolymerUpdates( SurfaceInformedPolymerUpdates && ) = default;
  SurfaceInformedPolymerUpdates & operator=( SurfaceInformedPolymerUpdates const & ) = delete;
  SurfaceInformedPolymerUpdates & operator=( SurfaceInformedPolymerUpdates && ) = delete;

  using DiscretizationOps = SolidModelDiscretizationOpsFullyAnisotropic;
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
  GEOS_FORCE_INLINE
  void computePlasticStrainIncrement( localIndex const k,
                                      localIndex const q,
                                      real64 const timeIncrement,
                                      real64 const ( &strainIncrement )[6],
                                      real64 const ( &stressIncrement )[6],
                                      real64 ( &plasticStrainIncrement )[6] ) const;

private:
  arrayView3d< real64 > const m_deformationGradient;
  arrayView3d< real64 > const m_plasticStrain;
  arrayView2d< real64 > const m_equivalentPlasticStrain;
  arrayView2d< real64 > const m_damage;
  arrayView1d< real64 > const m_temperature;
  arrayView2d< real64 > const m_jacobian;
  arrayView1d< real64 > const m_yieldStrength;

  real64 const m_defaultBulkModulus;
  real64 const m_defaultShearModulus;
  real64 const m_defaultYieldStrength;
  real64 const m_shearSofteningMagnitude;
  real64 const m_shearSofteningShapeParameter1;
  real64 const m_shearSofteningShapeParameter2;
  real64 const m_strainHardeningSlope;
  real64 const m_hardeningScaleExponent;
  real64 const m_maximumStretch;
  real64 const m_fractureStretchLambdaMin;
  real64 const m_fractureStretchLambda0;
  real64 const m_fractureStretchT0;
  real64 const m_fractureStretchTemperatureScale;

  real64 const m_glassTransitionTemperature;
  real64 const m_temperatureColdSlope;
  real64 const m_temperatureHotSlope;
  real64 const m_temperatureTransitionMagnitude;
  real64 const m_temperatureTransitionWidth;

  real64 const m_crystallinity;
  real64 const m_referenceCrystallinity;
  real64 const m_crystallinityTransitionWidth;
  real64 const m_elasticCrystallinityCoeff;
  real64 const m_yieldStrengthCrystallinityCoeff;

  real64 const m_pressureAsymmetryAmplitude;
  real64 const m_pressureAsymmetryWidth;
  real64 const m_compressivePressureStrengtheningCap;
};

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void SurfaceInformedPolymerUpdates::smallStrainUpdate( localIndex const k,
                                                       localIndex const q,
                                                       real64 const & timeIncrement,
                                                       real64 const ( &strainIncrement )[6],
                                                       real64 ( &stress )[6],
                                                       real64 ( &stiffness )[6][6] ) const
{
  GEOS_UNUSED_VAR( k );
  GEOS_UNUSED_VAR( q );
  GEOS_UNUSED_VAR( timeIncrement );
  GEOS_UNUSED_VAR( strainIncrement );
  GEOS_UNUSED_VAR( stress );
  GEOS_UNUSED_VAR( stiffness );
  GEOS_ERROR( "smallStrainUpdate is not implemented for SurfaceInformedPolymer; use the stress-only MPM update." );
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void SurfaceInformedPolymerUpdates::smallStrainUpdate( localIndex const k,
                                                       localIndex const q,
                                                       real64 const & timeIncrement,
                                                       real64 const ( &strainIncrement )[6],
                                                       real64 ( &stress )[6],
                                                       DiscretizationOps & stiffness ) const
{
  GEOS_UNUSED_VAR( k );
  GEOS_UNUSED_VAR( q );
  GEOS_UNUSED_VAR( timeIncrement );
  GEOS_UNUSED_VAR( strainIncrement );
  GEOS_UNUSED_VAR( stress );
  GEOS_UNUSED_VAR( stiffness );
  GEOS_ERROR( "smallStrainUpdate is not implemented for SurfaceInformedPolymer; use the stress-only MPM update." );
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void SurfaceInformedPolymerUpdates::smallStrainUpdate_StressOnly( localIndex const k,
                                                                  localIndex const q,
                                                                  real64 const & timeIncrement,
                                                                  real64 const ( &strainIncrement )[6],
                                                                  real64 ( &stress )[6] ) const
{
  GEOS_UNUSED_VAR( k );
  GEOS_UNUSED_VAR( q );
  GEOS_UNUSED_VAR( timeIncrement );
  GEOS_UNUSED_VAR( strainIncrement );
  GEOS_UNUSED_VAR( stress );
  GEOS_ERROR( "The rotation-aware stress-only update must be used for SurfaceInformedPolymer." );
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void SurfaceInformedPolymerUpdates::smallStrainUpdate_StressOnly( localIndex const k,
                                                                  localIndex const q,
                                                                  real64 const & timeIncrement,
                                                                  real64 const ( &beginningRotation )[3][3],
                                                                  real64 const ( &endRotation )[3][3],
                                                                  real64 const ( &strainIncrement )[6],
                                                                  real64 ( &stress )[6] ) const
{
  real64 const T = m_temperature[k];
  real64 const temperatureScale = surfaceInformedPolymerHelpers::temperatureScale( T,
                                                                                  m_glassTransitionTemperature,
                                                                                  m_temperatureColdSlope,
                                                                                  m_temperatureHotSlope,
                                                                                  m_temperatureTransitionMagnitude,
                                                                                  m_temperatureTransitionWidth );
  real64 const elasticCrystallinityScale = surfaceInformedPolymerHelpers::crystallinityScale( T,
                                                                                             m_glassTransitionTemperature,
                                                                                             m_crystallinity,
                                                                                             m_referenceCrystallinity,
                                                                                             m_elasticCrystallinityCoeff,
                                                                                             m_crystallinityTransitionWidth );

  m_bulkModulus[k] = m_defaultBulkModulus * temperatureScale * elasticCrystallinityScale;
  m_shearModulus[k] = m_defaultShearModulus * temperatureScale * elasticCrystallinityScale;

  ElasticIsotropicUpdates::smallStrainUpdate_StressOnly( k,
                                                         q,
                                                         timeIncrement,
                                                         strainIncrement,
                                                         stress );

  m_jacobian[k][q] *= LvArray::math::exp( strainIncrement[0] + strainIncrement[1] + strainIncrement[2] );

  if( !m_disableInelasticity )
  {
    SurfaceInformedPolymerUpdates::smallStrainUpdateHelper( k,
                                                            q,
                                                            timeIncrement,
                                                            beginningRotation,
                                                            endRotation,
                                                            strainIncrement,
                                                            stress );
  }

  saveStress( k, q, stress );
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void SurfaceInformedPolymerUpdates::smallStrainUpdateHelper( localIndex const k,
                                                             localIndex const q,
                                                             real64 const timeIncrement,
                                                             real64 const ( &beginningRotation )[3][3],
                                                             real64 const ( &endRotation )[3][3],
                                                             real64 const ( &strainIncrement )[6],
                                                             real64 ( &stress )[6] ) const
{
  real64 trialStress[6] = { 0.0 };
  LvArray::tensorOps::copy< 6 >( trialStress, stress );

  real64 trialP = 0.0;
  real64 trialQ = 0.0;
  real64 deviator[6] = { 0.0 };
  twoInvariant::stressDecomposition( trialStress, trialP, trialQ, deviator );

  real64 rotationTranspose[3][3] = { { 0.0 } };
  LvArray::tensorOps::transpose< 3, 3 >( rotationTranspose, beginningRotation );

  real64 oldPlasticStrain[6] = { 0.0 };
  LvArray::tensorOps::copy< 6 >( oldPlasticStrain, m_plasticStrain[k][q] );
  oldPlasticStrain[3] *= 0.5;
  oldPlasticStrain[4] *= 0.5;
  oldPlasticStrain[5] *= 0.5;

  real64 unrotatedOldPlasticStrain[6] = { 0.0 };
  LvArray::tensorOps::Rij_eq_AikSymBklAjl< 3 >( unrotatedOldPlasticStrain, rotationTranspose, oldPlasticStrain );
  unrotatedOldPlasticStrain[3] *= 2.0;
  unrotatedOldPlasticStrain[4] *= 2.0;
  unrotatedOldPlasticStrain[5] *= 2.0;

  real64 unrotatedDeformationGradient[3][3] = { { 0.0 } };
  LvArray::tensorOps::Rij_eq_AikBkj< 3, 3, 3 >( unrotatedDeformationGradient, rotationTranspose, m_deformationGradient[k] );
  real64 const lambdaChain = surfaceInformedPolymerHelpers::chainStretch( unrotatedDeformationGradient );

  real64 const T = m_temperature[k];
  real64 const lambdaFailure = surfaceInformedPolymerHelpers::fractureStretch( T,
                                                                              m_maximumStretch,
                                                                              m_fractureStretchLambdaMin,
                                                                              m_fractureStretchLambda0,
                                                                              m_fractureStretchT0,
                                                                              m_fractureStretchTemperatureScale );
  if( lambdaChain > lambdaFailure )
  {
    m_damage[k][q] = 1.0;
  }
  real64 const ST = surfaceInformedPolymerHelpers::temperatureScale( T,
                                                                    m_glassTransitionTemperature,
                                                                    m_temperatureColdSlope,
                                                                    m_temperatureHotSlope,
                                                                    m_temperatureTransitionMagnitude,
                                                                    m_temperatureTransitionWidth );
  real64 const CY = surfaceInformedPolymerHelpers::crystallinityScale( T,
                                                                      m_glassTransitionTemperature,
                                                                      m_crystallinity,
                                                                      m_referenceCrystallinity,
                                                                      m_yieldStrengthCrystallinityCoeff,
                                                                      m_crystallinityTransitionWidth );
  real64 const yield0 = m_defaultYieldStrength * ST * CY;
  real64 const softeningMagnitude = m_shearSofteningMagnitude * ST;
  real64 const hardeningSlope = m_strainHardeningSlope * LvArray::math::pow( LvArray::math::max( ST, 1.0e-16 ), m_hardeningScaleExponent );
  real64 const pressureAsymmetry = surfaceInformedPolymerHelpers::pressureAsymmetry( T,
                                                                                    m_glassTransitionTemperature,
                                                                                    m_pressureAsymmetryAmplitude,
                                                                                    m_pressureAsymmetryWidth );
  real64 const stretchHardening = hardeningSlope * surfaceInformedPolymerHelpers::stretchHardeningMeasure( lambdaChain );
  real64 const pressureStrengtheningP = surfaceInformedPolymerHelpers::pressureStrengtheningMeanStress( trialP,
                                                                                                  m_compressivePressureStrengtheningCap );

  real64 plasticStrainIncrement[6] = { 0.0 };
  real64 yieldStrengthOld = m_yieldStrength[k] > 0.0 ? m_yieldStrength[k] : yield0 + softeningMagnitude + stretchHardening;
  real64 yieldStrengthIter = yieldStrengthOld;

  real64 const relTol = 1.0e-10;
  int const maxEvals = 100;

  for( int iter = 0; iter < maxEvals; ++iter )
  {
    real64 unrotatedTempPlasticStrain[6] = { 0.0 };
    LvArray::tensorOps::copy< 6 >( unrotatedTempPlasticStrain, unrotatedOldPlasticStrain );
    LvArray::tensorOps::add< 6 >( unrotatedTempPlasticStrain, plasticStrainIncrement );

    real64 const kappa = surfaceInformedPolymerHelpers::strainMagnitude( unrotatedTempPlasticStrain );
    real64 const plasticSoftening = surfaceInformedPolymerHelpers::softeningContribution( softeningMagnitude,
                                                                                          kappa,
                                                                                          m_shearSofteningShapeParameter1,
                                                                                          m_shearSofteningShapeParameter2 );
    yieldStrengthIter = yield0 + plasticSoftening + stretchHardening;

    real64 const yieldFunction = trialQ + pressureAsymmetry * pressureStrengtheningP - yieldStrengthIter;
    if( yieldFunction > 0.0 || iter > 0 )
    {
      real64 const returnedQ = LvArray::math::min( trialQ,
                                                  LvArray::math::max( yieldStrengthIter - pressureAsymmetry * pressureStrengtheningP, 0.0 ) );

      real64 stressTemp[6] = { 0.0 };
      twoInvariant::stressRecomposition( trialP, returnedQ, deviator, stressTemp );

      real64 stressIncrement[6] = { 0.0 };
      LvArray::tensorOps::copy< 6 >( stressIncrement, stressTemp );
      LvArray::tensorOps::subtract< 6 >( stressIncrement, m_oldStress[k][q] );

      computePlasticStrainIncrement( k,
                                     q,
                                     timeIncrement,
                                     strainIncrement,
                                     stressIncrement,
                                     plasticStrainIncrement );

      real64 unrotatedNewPlasticStrain[6] = { 0.0 };
      LvArray::tensorOps::copy< 6 >( unrotatedNewPlasticStrain, unrotatedOldPlasticStrain );
      LvArray::tensorOps::add< 6 >( unrotatedNewPlasticStrain, plasticStrainIncrement );

      real64 const yieldScale = LvArray::math::max( 1.0, LvArray::math::abs( yieldStrengthIter ) );
      if( LvArray::math::abs( yieldStrengthIter - yieldStrengthOld ) <= relTol * yieldScale )
      {
        real64 newPlasticStrainTensorVoigt[6] = { 0.0 };
        LvArray::tensorOps::copy< 6 >( newPlasticStrainTensorVoigt, unrotatedNewPlasticStrain );
        newPlasticStrainTensorVoigt[3] *= 0.5;
        newPlasticStrainTensorVoigt[4] *= 0.5;
        newPlasticStrainTensorVoigt[5] *= 0.5;

        real64 newPlasticStrain[6] = { 0.0 };
        LvArray::tensorOps::Rij_eq_AikSymBklAjl< 3 >( newPlasticStrain, endRotation, newPlasticStrainTensorVoigt );
        newPlasticStrain[3] *= 2.0;
        newPlasticStrain[4] *= 2.0;
        newPlasticStrain[5] *= 2.0;

        LvArray::tensorOps::copy< 6 >( m_plasticStrain[k][q], newPlasticStrain );
        m_equivalentPlasticStrain[k][q] = surfaceInformedPolymerHelpers::strainMagnitude( unrotatedNewPlasticStrain );
        m_yieldStrength[k] = yieldStrengthIter;
        LvArray::tensorOps::copy< 6 >( stress, stressTemp );
        return;
      }

      yieldStrengthOld = yieldStrengthIter;
    }
    else
    {
      m_equivalentPlasticStrain[k][q] = kappa;
      m_yieldStrength[k] = yieldStrengthIter;
      return;
    }
  }

  GEOS_ERROR( "Plastic correction for SurfaceInformedPolymer did not converge." );
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void SurfaceInformedPolymerUpdates::computePlasticStrainIncrement( localIndex const k,
                                                                   localIndex const q,
                                                                   real64 const timeIncrement,
                                                                   real64 const ( &strainIncrement )[6],
                                                                   real64 const ( &stressIncrement )[6],
                                                                   real64 ( &plasticStrainIncrement )[6] ) const
{
  GEOS_UNUSED_VAR( q );
  GEOS_UNUSED_VAR( timeIncrement );

  real64 deltaP = 0.0;
  real64 deltaQ = 0.0;
  real64 stressIncrementDeviator[6] = { 0.0 };
  twoInvariant::stressDecomposition( stressIncrement, deltaP, deltaQ, stressIncrementDeviator );

  real64 stressIncrementIsostatic[6] = { 0.0 };
  stressIncrementIsostatic[0] = deltaP;
  stressIncrementIsostatic[1] = deltaP;
  stressIncrementIsostatic[2] = deltaP;

  real64 elasticStrainIncrement[6] = { 0.0 };
  for( localIndex i = 0; i < 6; ++i )
  {
    if( m_bulkModulus[k] > 1.0e-12 )
    {
      elasticStrainIncrement[i] += ( 1 + ( i >= 3 ) ) * stressIncrementIsostatic[i] / ( 3.0 * m_bulkModulus[k] );
    }
    if( m_shearModulus[k] > 1.0e-12 )
    {
      elasticStrainIncrement[i] += ( 1 + ( i >= 3 ) ) * LvArray::math::sqrt( 2.0 / 3.0 ) *
                                   deltaQ * stressIncrementDeviator[i] / ( 2.0 * m_shearModulus[k] );
    }
  }

  LvArray::tensorOps::copy< 6 >( plasticStrainIncrement, strainIncrement );
  LvArray::tensorOps::subtract< 6 >( plasticStrainIncrement, elasticStrainIncrement );
}

class SurfaceInformedPolymer : public ElasticIsotropic
{
public:
  using KernelWrapper = SurfaceInformedPolymerUpdates;

  SurfaceInformedPolymer( string const & name, Group * const parent );
  virtual ~SurfaceInformedPolymer() override;

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  virtual void saveConvergedState() const override;

  static constexpr auto m_catalogNameString = "SurfaceInformedPolymer";
  static string catalogName() { return m_catalogNameString; }
  virtual string getCatalogName() const override { return catalogName(); }

  struct viewKeyStruct : public SolidBase::viewKeyStruct
  {
    static constexpr char const * deformationGradientString() { return "deformationGradient"; }
    static constexpr char const * plasticStrainString() { return "plasticStrain"; }
    static constexpr char const * equivalentPlasticStrainString() { return "equivalentPlasticStrain"; }
    static constexpr char const * damageString() { return "damage"; }
    static constexpr char const * temperatureString() { return "temperature"; }
    static constexpr char const * jacobianString() { return "jacobian"; }
    static constexpr char const * yieldStrengthString() { return "yieldStrength"; }

    static constexpr char const * defaultYieldStrengthString() { return "defaultYieldStrength"; }
    static constexpr char const * shearSofteningMagnitudeString() { return "shearSofteningMagnitude"; }
    static constexpr char const * shearSofteningShapeParameter1String() { return "shearSofteningShapeParameter1"; }
    static constexpr char const * shearSofteningShapeParameter2String() { return "shearSofteningShapeParameter2"; }
    static constexpr char const * strainHardeningSlopeString() { return "strainHardeningSlope"; }
    static constexpr char const * hardeningScaleExponentString() { return "hardeningScaleExponent"; }
    static constexpr char const * maximumStretchString() { return "maximumStretch"; }
    static constexpr char const * fractureStretchLambdaMinString() { return "fractureStretchLambdaMin"; }
    static constexpr char const * fractureStretchLambda0String() { return "fractureStretchLambda0"; }
    static constexpr char const * fractureStretchT0String() { return "fractureStretchT0"; }
    static constexpr char const * fractureStretchTemperatureScaleString() { return "fractureStretchTemperatureScale"; }

    static constexpr char const * glassTransitionTemperatureString() { return "glassTransitionTemperature"; }
    static constexpr char const * temperatureColdSlopeString() { return "temperatureColdSlope"; }
    static constexpr char const * temperatureHotSlopeString() { return "temperatureHotSlope"; }
    static constexpr char const * temperatureTransitionMagnitudeString() { return "temperatureTransitionMagnitude"; }
    static constexpr char const * temperatureTransitionWidthString() { return "temperatureTransitionWidth"; }

    static constexpr char const * crystallinityString() { return "crystallinity"; }
    static constexpr char const * referenceCrystallinityString() { return "referenceCrystallinity"; }
    static constexpr char const * crystallinityTransitionWidthString() { return "crystallinityTransitionWidth"; }
    static constexpr char const * elasticCrystallinityCoeffString() { return "elasticCrystallinityCoeff"; }
    static constexpr char const * yieldStrengthCrystallinityCoeffString() { return "yieldStrengthCrystallinityCoeff"; }

    static constexpr char const * pressureAsymmetryAmplitudeString() { return "pressureAsymmetryAmplitude"; }
    static constexpr char const * pressureAsymmetryWidthString() { return "pressureAsymmetryWidth"; }
    static constexpr char const * compressivePressureStrengtheningCapString() { return "compressivePressureStrengtheningCap"; }
  };

  SurfaceInformedPolymerUpdates createKernelUpdates() const
  {
    return SurfaceInformedPolymerUpdates( m_deformationGradient,
                                          m_plasticStrain,
                                          m_equivalentPlasticStrain,
                                          m_damage,
                                          m_temperature,
                                          m_jacobian,
                                          m_yieldStrength,
                                          m_defaultBulkModulus,
                                          m_defaultShearModulus,
                                          m_defaultYieldStrength,
                                          m_shearSofteningMagnitude,
                                          m_shearSofteningShapeParameter1,
                                          m_shearSofteningShapeParameter2,
                                          m_strainHardeningSlope,
                                          m_hardeningScaleExponent,
                                          m_maximumStretch,
                                          m_fractureStretchLambdaMin,
                                          m_fractureStretchLambda0,
                                          m_fractureStretchT0,
                                          m_fractureStretchTemperatureScale,
                                          m_glassTransitionTemperature,
                                          m_temperatureColdSlope,
                                          m_temperatureHotSlope,
                                          m_temperatureTransitionMagnitude,
                                          m_temperatureTransitionWidth,
                                          m_crystallinity,
                                          m_referenceCrystallinity,
                                          m_crystallinityTransitionWidth,
                                          m_elasticCrystallinityCoeff,
                                          m_yieldStrengthCrystallinityCoeff,
                                          m_pressureAsymmetryAmplitude,
                                          m_pressureAsymmetryWidth,
                                          m_compressivePressureStrengtheningCap,
                                          m_bulkModulus,
                                          m_shearModulus,
                                          m_thermalExpansionCoefficient,
                                          m_newStress,
                                          m_oldStress,
                                          m_density,
                                          m_wavespeed,
                                          m_disableInelasticity );
  }

  template< typename UPDATE_KERNEL, typename ... PARAMS >
  UPDATE_KERNEL createDerivedKernelUpdates( PARAMS && ... constructorParams )
  {
    return UPDATE_KERNEL( std::forward< PARAMS >( constructorParams )...,
                          m_deformationGradient,
                          m_plasticStrain,
                          m_equivalentPlasticStrain,
                          m_damage,
                          m_temperature,
                          m_jacobian,
                          m_yieldStrength,
                          m_defaultBulkModulus,
                          m_defaultShearModulus,
                          m_defaultYieldStrength,
                          m_shearSofteningMagnitude,
                          m_shearSofteningShapeParameter1,
                          m_shearSofteningShapeParameter2,
                          m_strainHardeningSlope,
                          m_hardeningScaleExponent,
                          m_maximumStretch,
                          m_fractureStretchLambdaMin,
                          m_fractureStretchLambda0,
                          m_fractureStretchT0,
                          m_fractureStretchTemperatureScale,
                          m_glassTransitionTemperature,
                          m_temperatureColdSlope,
                          m_temperatureHotSlope,
                          m_temperatureTransitionMagnitude,
                          m_temperatureTransitionWidth,
                          m_crystallinity,
                          m_referenceCrystallinity,
                          m_crystallinityTransitionWidth,
                          m_elasticCrystallinityCoeff,
                          m_yieldStrengthCrystallinityCoeff,
                          m_pressureAsymmetryAmplitude,
                          m_pressureAsymmetryWidth,
                          m_compressivePressureStrengtheningCap,
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

  array3d< real64 > m_deformationGradient;
  array3d< real64 > m_plasticStrain;
  array2d< real64 > m_equivalentPlasticStrain;
  array2d< real64 > m_damage;
  array1d< real64 > m_temperature;
  array2d< real64 > m_jacobian;
  array1d< real64 > m_yieldStrength;

  real64 m_defaultYieldStrength;
  real64 m_shearSofteningMagnitude;
  real64 m_shearSofteningShapeParameter1;
  real64 m_shearSofteningShapeParameter2;
  real64 m_strainHardeningSlope;
  real64 m_hardeningScaleExponent;
  real64 m_maximumStretch;
  real64 m_fractureStretchLambdaMin;
  real64 m_fractureStretchLambda0;
  real64 m_fractureStretchT0;
  real64 m_fractureStretchTemperatureScale;

  real64 m_glassTransitionTemperature;
  real64 m_temperatureColdSlope;
  real64 m_temperatureHotSlope;
  real64 m_temperatureTransitionMagnitude;
  real64 m_temperatureTransitionWidth;

  real64 m_crystallinity;
  real64 m_referenceCrystallinity;
  real64 m_crystallinityTransitionWidth;
  real64 m_elasticCrystallinityCoeff;
  real64 m_yieldStrengthCrystallinityCoeff;

  real64 m_pressureAsymmetryAmplitude;
  real64 m_pressureAsymmetryWidth;
  real64 m_compressivePressureStrengtheningCap;
};

} // namespace constitutive
} // namespace geos

#endif // GEOS_CONSTITUTIVE_SOLID_SURFACEINFORMEDPOLYMER_HPP_
