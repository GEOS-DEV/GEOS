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
 * @file SurfaceInformedPolymerCohesiveZone.hpp
 * @brief Thin-film cohesive-zone projection of SurfaceInformedPolymer.
 *
 * @details
 * This cohesive-zone law interprets the displacement jump as the deformation of a finite-thickness
 * polymer layer.  The normal and tangential jumps are converted to nominal film strain measures,
 *
 *   eps_n = delta_n / h0,  gamma = delta_t / h0,
 *
 * and elastic trial tractions are computed from the same temperature/crystallinity-scaled moduli used
 * by the continuum polymer.  The normal film response is decomposed into a retained volumetric mean
 * stress p_t=K eps_n and a return-mapped deviatoric normal stress s_n.  This is important for thin
 * layers in uniaxial-strain loading: the cohesive projection must preserve the large constrained
 * volumetric traction carried by the corresponding continuum layer instead of limiting the total
 * normal traction to the flow stress.
 *
 * The deviatoric film stress state is then returned to the same scalar surface used by the continuum
 * model,
 *
 *   q + eta(T) p_eff = Y(T,Xc) + S_soft(T) exp[-(kappa/r1)^r2]
 *                       + H_g S_T(T)^p_H [lambda_eff^2 - lambda_eff^{-1}],
 *
 * where q=sqrt((3 s_n/2)^2+3 tau^2) for the reduced normal-shear film state and p_eff is p_t after
 * optional compressive clipping.  The sign returned to the MPM cohesive infrastructure is the
 * existing GEOS cohesive sign convention: opening produces a negative normal stress value.  A
 * separate maximumStretch parameter is retained as a cohesive failure cutoff; once exceeded, the law
 * returns zero traction in that update.
 */

#ifndef GEOS_CONSTITUTIVE_SURFACEINFORMEDPOLYMERCOHESIVEZONE_HPP_
#define GEOS_CONSTITUTIVE_SURFACEINFORMEDPOLYMERCOHESIVEZONE_HPP_

#include "constitutive/cohesiveZone/CohesiveZoneBase.hpp"
#include "constitutive/solid/SurfaceInformedPolymerHelpers.hpp"
#include "LvArray/src/tensorOps.hpp"

#include <cfloat>

namespace geos
{
namespace constitutive
{

class SurfaceInformedPolymerCohesiveZoneUpdates : public CohesiveZoneBaseUpdates
{
public:
  SurfaceInformedPolymerCohesiveZoneUpdates( real64 const & thickness,
                                             real64 const & bulkModulus,
                                             real64 const & shearModulus,
                                             real64 const & defaultYieldStrength,
                                             real64 const & shearSofteningMagnitude,
                                             real64 const & shearSofteningShapeParameter1,
                                             real64 const & shearSofteningShapeParameter2,
                                             real64 const & strainHardeningSlope,
                                             real64 const & hardeningScaleExponent,
                                             real64 const & maximumStretch,
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
                                             arrayView1d< real64 > const & damage,
                                             arrayView1d< real64 > const & temperature,
                                             arrayView1d< real64 > const & previousLambda,
                                             arrayView1d< real64 > const & equivalentPlasticStrain,
                                             arrayView1d< real64 > const & plasticNormalStrain,
                                             arrayView1d< real64 > const & plasticTangentialStrain,
                                             arrayView2d< real64 > const & newNormalStress,
                                             arrayView2d< real64 > const & newShearStress,
                                             arrayView2d< real64 > const & oldNormalStress,
                                             arrayView2d< real64 > const & oldShearStress ):
    CohesiveZoneBaseUpdates( newNormalStress,
                             newShearStress,
                             oldNormalStress,
                             oldShearStress ),
    m_thickness( thickness ),
    m_bulkModulus( bulkModulus ),
    m_shearModulus( shearModulus ),
    m_defaultYieldStrength( defaultYieldStrength ),
    m_shearSofteningMagnitude( shearSofteningMagnitude ),
    m_shearSofteningShapeParameter1( shearSofteningShapeParameter1 ),
    m_shearSofteningShapeParameter2( shearSofteningShapeParameter2 ),
    m_strainHardeningSlope( strainHardeningSlope ),
    m_hardeningScaleExponent( hardeningScaleExponent ),
    m_maximumStretch( maximumStretch ),
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
    m_compressivePressureStrengtheningCap( compressivePressureStrengtheningCap ),
    m_damage( damage ),
    m_temperature( temperature ),
    m_previousLambda( previousLambda ),
    m_equivalentPlasticStrain( equivalentPlasticStrain ),
    m_plasticNormalStrain( plasticNormalStrain ),
    m_plasticTangentialStrain( plasticTangentialStrain )
  {}

  SurfaceInformedPolymerCohesiveZoneUpdates() = delete;
  SurfaceInformedPolymerCohesiveZoneUpdates( SurfaceInformedPolymerCohesiveZoneUpdates const & ) = default;
  SurfaceInformedPolymerCohesiveZoneUpdates( SurfaceInformedPolymerCohesiveZoneUpdates && ) = default;
  SurfaceInformedPolymerCohesiveZoneUpdates & operator=( SurfaceInformedPolymerCohesiveZoneUpdates const & ) = delete;
  SurfaceInformedPolymerCohesiveZoneUpdates & operator=( SurfaceInformedPolymerCohesiveZoneUpdates && ) = delete;

  GEOS_HOST_DEVICE
  void jumpDisplacementUpdate( localIndex const k,
                               real64 const & normalDisplacement,
                               real64 const & tangentialDisplacement,
                               real64 & normalStress,
                               real64 & shearStress ) const
  {
    if( m_damage[k] >= 1.0 )
    {
      normalStress = 0.0;
      shearStress = 0.0;
      return;
    }

    real64 const thickness = LvArray::math::max( m_thickness, 1.0e-16 );
    real64 const normalStrain = normalDisplacement / thickness;
    real64 const tangentialStrain = tangentialDisplacement / thickness;

    real64 filmDeformationGradient[3][3] = { { 0.0 } };
    filmDeformationGradient[0][0] = 1.0;
    filmDeformationGradient[0][1] = tangentialStrain;
    filmDeformationGradient[1][1] = 1.0 + normalStrain;
    filmDeformationGradient[2][2] = 1.0;

    real64 const lambdaChain = surfaceInformedPolymerHelpers::chainStretch( filmDeformationGradient );
    m_previousLambda[k] = lambdaChain;
    if( lambdaChain > m_maximumStretch )
    {
      m_damage[k] = 1.0;
      normalStress = 0.0;
      shearStress = 0.0;
      return;
    }

    real64 const T = m_temperature[k];
    real64 const ST = surfaceInformedPolymerHelpers::temperatureScale( T,
                                                                      m_glassTransitionTemperature,
                                                                      m_temperatureColdSlope,
                                                                      m_temperatureHotSlope,
                                                                      m_temperatureTransitionMagnitude,
                                                                      m_temperatureTransitionWidth );
    real64 const CE = surfaceInformedPolymerHelpers::crystallinityScale( T,
                                                                        m_glassTransitionTemperature,
                                                                        m_crystallinity,
                                                                        m_referenceCrystallinity,
                                                                        m_elasticCrystallinityCoeff,
                                                                        m_crystallinityTransitionWidth );
    real64 const CY = surfaceInformedPolymerHelpers::crystallinityScale( T,
                                                                        m_glassTransitionTemperature,
                                                                        m_crystallinity,
                                                                        m_referenceCrystallinity,
                                                                        m_yieldStrengthCrystallinityCoeff,
                                                                        m_crystallinityTransitionWidth );

    real64 const K = m_bulkModulus * ST * CE;
    real64 const G = m_shearModulus * ST * CE;

    real64 const yield0 = m_defaultYieldStrength * ST * CY;
    real64 const softeningMagnitude = m_shearSofteningMagnitude * ST;
    real64 const hardeningSlope = m_strainHardeningSlope * LvArray::math::pow( LvArray::math::max( ST, 1.0e-16 ), m_hardeningScaleExponent );
    real64 const stretchHardening = hardeningSlope * surfaceInformedPolymerHelpers::stretchHardeningMeasure( lambdaChain );
    real64 const pressureAsymmetry = surfaceInformedPolymerHelpers::pressureAsymmetry( T,
                                                                                      m_glassTransitionTemperature,
                                                                                      m_pressureAsymmetryAmplitude,
                                                                                      m_pressureAsymmetryWidth );

    real64 const plasticNormal0 = m_plasticNormalStrain[k];
    real64 const plasticTangential0 = m_plasticTangentialStrain[k];
    real64 const kappa0 = m_equivalentPlasticStrain[k];

    // The cohesive film is a reduced normal-shear representation of a finite-thickness continuum
    // layer.  For the normal part, use the uniaxial-strain split of the corresponding continuum
    // response: p_t=K eps_n is retained as a volumetric stress, while only the deviatoric normal
    // stress s_n=(4/3)G(eps_n-eps_n^p) participates in the radial return.  This keeps a nearly
    // incompressible film stiff in constrained tension/compression while still allowing plastic
    // flow of the deviatoric part.
    real64 const pTrial = K * normalStrain;
    real64 const pressureStrengtheningP = surfaceInformedPolymerHelpers::pressureStrengtheningMeanStress( pTrial,
                                                                                            m_compressivePressureStrengtheningCap );
    real64 const normalDeviatoricModulus = 4.0 * G / 3.0;
    real64 const sNTrial = normalDeviatoricModulus * ( normalStrain - plasticNormal0 );
    real64 const tauTrial = G * ( tangentialStrain - plasticTangential0 );
    real64 const qTrial = LvArray::math::sqrt( 2.25 * sNTrial * sNTrial + 3.0 * tauTrial * tauTrial );

    real64 kappa = kappa0;
    real64 sN = sNTrial;
    real64 tau = tauTrial;
    real64 sigmaN = pTrial + sNTrial;
    bool plastic = false;

    for( int iter = 0; iter < 16; ++iter )
    {
      real64 const flowStrength = yield0 +
                                  surfaceInformedPolymerHelpers::softeningContribution( softeningMagnitude,
                                                                                         kappa,
                                                                                         m_shearSofteningShapeParameter1,
                                                                                         m_shearSofteningShapeParameter2 ) +
                                  stretchHardening;
      real64 const yieldMeasure = qTrial + pressureAsymmetry * pressureStrengtheningP;
      if( yieldMeasure <= flowStrength || qTrial <= 1.0e-16 )
      {
        plastic = false;
        sN = sNTrial;
        tau = tauTrial;
        sigmaN = pTrial + sN;
        kappa = kappa0;
        break;
      }

      plastic = true;
      real64 const returnedQ = LvArray::math::min( qTrial,
                                                   LvArray::math::max( flowStrength - pressureAsymmetry * pressureStrengtheningP, 0.0 ) );
      real64 const scale = returnedQ / qTrial;
      real64 const sNNew = scale * sNTrial;
      real64 const tauNew = scale * tauTrial;
      real64 const deltaKappa = LvArray::math::max( 0.0, ( qTrial - returnedQ ) / ( 3.0 * LvArray::math::max( G, 1.0e-16 ) ) );
      real64 const kappaNew = kappa0 + deltaKappa;

      if( LvArray::math::abs( kappaNew - kappa ) <= 1.0e-10 * LvArray::math::max( 1.0, kappaNew ) )
      {
        sN = sNNew;
        tau = tauNew;
        sigmaN = pTrial + sN;
        kappa = kappaNew;
        break;
      }
      kappa = kappaNew;
      sN = sNNew;
      tau = tauNew;
      sigmaN = pTrial + sN;
    }

    if( plastic )
    {
      m_plasticNormalStrain[k] = normalStrain - sN / LvArray::math::max( normalDeviatoricModulus, 1.0e-16 );
      m_plasticTangentialStrain[k] = tangentialStrain - tau / LvArray::math::max( G, 1.0e-16 );
      m_equivalentPlasticStrain[k] = LvArray::math::max( m_equivalentPlasticStrain[k], kappa );
    }

    // GEOS cohesive normal stress uses the opposite sign from the tensile-positive film stress used
    // by the material surface above.  Thus opening, which has sigmaN>0, returns negative normalStress.
    normalStress = -sigmaN;
    shearStress = -tau;
  }

private:
  real64 m_thickness;
  real64 m_bulkModulus;
  real64 m_shearModulus;
  real64 m_defaultYieldStrength;
  real64 m_shearSofteningMagnitude;
  real64 m_shearSofteningShapeParameter1;
  real64 m_shearSofteningShapeParameter2;
  real64 m_strainHardeningSlope;
  real64 m_hardeningScaleExponent;
  real64 m_maximumStretch;
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

  arrayView1d< real64 > const m_damage;
  arrayView1d< real64 > const m_temperature;
  arrayView1d< real64 > const m_previousLambda;
  arrayView1d< real64 > const m_equivalentPlasticStrain;
  arrayView1d< real64 > const m_plasticNormalStrain;
  arrayView1d< real64 > const m_plasticTangentialStrain;
};

class SurfaceInformedPolymerCohesiveZone : public CohesiveZoneBase
{
public:
  using KernelWrapper = SurfaceInformedPolymerCohesiveZoneUpdates;

  SurfaceInformedPolymerCohesiveZone( string const & name, Group * const parent );
  virtual ~SurfaceInformedPolymerCohesiveZone() override;

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  static constexpr auto m_catalogNameString = "SurfaceInformedPolymerCohesiveZone";
  static std::string catalogName() { return m_catalogNameString; }
  virtual string getCatalogName() const override { return catalogName(); }

  struct viewKeyStruct : public CohesiveZoneBase::viewKeyStruct
  {
    static constexpr char const * thicknessString() { return "thickness"; }
    static constexpr char const * bulkModulusString() { return "bulkModulus"; }
    static constexpr char const * shearModulusString() { return "shearModulus"; }
    static constexpr char const * defaultYieldStrengthString() { return "defaultYieldStrength"; }
    static constexpr char const * shearSofteningMagnitudeString() { return "shearSofteningMagnitude"; }
    static constexpr char const * shearSofteningShapeParameter1String() { return "shearSofteningShapeParameter1"; }
    static constexpr char const * shearSofteningShapeParameter2String() { return "shearSofteningShapeParameter2"; }
    static constexpr char const * strainHardeningSlopeString() { return "strainHardeningSlope"; }
    static constexpr char const * hardeningScaleExponentString() { return "hardeningScaleExponent"; }
    static constexpr char const * maximumStretchString() { return "maximumStretch"; }

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

    static constexpr char const * damageString() { return "damage"; }
    static constexpr char const * temperatureString() { return "temperature"; }
    static constexpr char const * previousLambdaString() { return "previousLambda"; }
    static constexpr char const * equivalentPlasticStrainString() { return "equivalentPlasticStrain"; }
    static constexpr char const * plasticNormalStrainString() { return "plasticNormalStrain"; }
    static constexpr char const * plasticTangentialStrainString() { return "plasticTangentialStrain"; }
  };

  arrayView1d< real64 > const getDamage()
  {
    return m_damage;
  }

  arrayView1d< real64 const > const getDamage() const
  {
    return m_damage;
  }

  SurfaceInformedPolymerCohesiveZoneUpdates createKernelUpdates( bool const includeState = true ) const
  {
    GEOS_UNUSED_VAR( includeState );
    return SurfaceInformedPolymerCohesiveZoneUpdates( m_thickness,
                                                      m_bulkModulus,
                                                      m_shearModulus,
                                                      m_defaultYieldStrength,
                                                      m_shearSofteningMagnitude,
                                                      m_shearSofteningShapeParameter1,
                                                      m_shearSofteningShapeParameter2,
                                                      m_strainHardeningSlope,
                                                      m_hardeningScaleExponent,
                                                      m_maximumStretch,
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
                                                      m_damage,
                                                      m_temperature,
                                                      m_previousLambda,
                                                      m_equivalentPlasticStrain,
                                                      m_plasticNormalStrain,
                                                      m_plasticTangentialStrain,
                                                      m_newNormalStress,
                                                      m_newShearStress,
                                                      m_oldNormalStress,
                                                      m_oldShearStress );
  }

  template< typename UPDATE_KERNEL, typename ... PARAMS >
  UPDATE_KERNEL createDerivedKernelUpdates( PARAMS && ... constructorParams ) const
  {
    return UPDATE_KERNEL( std::forward< PARAMS >( constructorParams )...,
                          m_thickness,
                          m_bulkModulus,
                          m_shearModulus,
                          m_defaultYieldStrength,
                          m_shearSofteningMagnitude,
                          m_shearSofteningShapeParameter1,
                          m_shearSofteningShapeParameter2,
                          m_strainHardeningSlope,
                          m_hardeningScaleExponent,
                          m_maximumStretch,
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
                          m_damage,
                          m_temperature,
                          m_previousLambda,
                          m_equivalentPlasticStrain,
                          m_plasticNormalStrain,
                          m_plasticTangentialStrain,
                          m_newNormalStress,
                          m_newShearStress,
                          m_oldNormalStress,
                          m_oldShearStress );
  }

protected:
  virtual void postInputInitialization() override;

  real64 m_thickness;
  real64 m_bulkModulus;
  real64 m_shearModulus;
  real64 m_defaultYieldStrength;
  real64 m_shearSofteningMagnitude;
  real64 m_shearSofteningShapeParameter1;
  real64 m_shearSofteningShapeParameter2;
  real64 m_strainHardeningSlope;
  real64 m_hardeningScaleExponent;
  real64 m_maximumStretch;
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

  array1d< real64 > m_damage;
  array1d< real64 > m_temperature;
  array1d< real64 > m_previousLambda;
  array1d< real64 > m_equivalentPlasticStrain;
  array1d< real64 > m_plasticNormalStrain;
  array1d< real64 > m_plasticTangentialStrain;
};

} // namespace constitutive
} // namespace geos

#endif // GEOS_CONSTITUTIVE_SURFACEINFORMEDPOLYMERCOHESIVEZONE_HPP_
