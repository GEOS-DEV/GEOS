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
 * @file Graphite.hpp
 * @brief Transversely isotropic strength and damage model for graphitic materials.
 *
 * The model represents a graphite-like layered solid whose first material-direction
 * row is the normal to the basal/weak planes.  The elastic predictor is
 * transversely isotropic.  The inelastic corrector decomposes stress into
 * signed distortion, in-plane shear, and coupled basal-shear modes and
 * limits the combined modal stress with a quadratic TI strength index.  The
 * signed distortion invariant separates states where the basal-plane normal is
 * relatively more compressed from states where it is relatively less compressed
 * or closer to opening.  The scalar
 * @c damage field is retained as the host-visible DFG failure/separation
 * indicator; it may reach one for either a resolved fracture or a comminuted
 * powder cloud.  Internal basal-plane and comminution damage variables control
 * how the strength degrades mechanically so that a damaged material point can
 * still carry pressure-dependent residual powder strength under confinement.
 */

#ifndef GEOS_CONSTITUTIVE_SOLID_GRAPHITE_HPP_
#define GEOS_CONSTITUTIVE_SOLID_GRAPHITE_HPP_

#include "SolidBase.hpp"
#include "InvariantDecompositions.hpp"
#include "PropertyConversions.hpp"
#include "SolidModelDiscretizationOpsFullyAnisotropic.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geos
{

namespace constitutive
{

/**
 * @class GraphiteUpdates
 *
 * Class to provide material updates that may be
 * called from a kernel function.
 */
class GraphiteUpdates : public SolidBaseUpdates
{
public:

  /**
   * @brief Constructor
   * @param[in] velocityGradient The ArrayView holding the velocity gradient for each element/particle.
   * @param[in] plasticStrain The ArrayView holding the plastic strain for each quadrature point.
   * @param[in] relaxation The ArrayView holding the relaxation for each quadrature point.
   * @param[in] enableCrackTipStressConcentration Flag to activate crack-tip stress concentration
   * @param[in] crackTipStressConcentration The ArrayView holding the crack-tip stress concentration (plotting array not needed)
   * @param[in] distanceToCrackTip The ArrayView holding the crack tip distance computed in the solver and passed here.
   * @param[in] basalPlanePlasticWork The ArrayView holding the basal-plane plastic work for each quadrature point.
   * @param[in] plasticWork The ArrayView holding the total plastic work for each quadrature point.
   * @param[in] alphaL Thermal expansion logitudinal to crystal symmetry axis
   * @param[in] alphaT Thermal expansion transverse to crystal symmetry axis
   * @param[in] damage The ArrayView holding the host-visible DFG damage for each quadrature point.
   * @param[in] basalPlaneDamage The ArrayView holding basal-plane fracture damage for each quadrature point.
   * @param[in] comminutionDamage The ArrayView holding comminution/powder damage for each quadrature point.
   * @param[in] enableDistension Flag to activate stress-free closure of directional distension.
   * @param[in] basalNormalDistension The ArrayView holding basal-normal closeable distension for each quadrature point.
   * @param[in] transverseDistension The ArrayView holding transverse closeable distension for each quadrature point.
   * @param[in] porosity The ArrayView holding distension-derived porosity for each quadrature point.
   * @param[in] temperature The ArrayView holding the temperature for each element/particle.
   * @param[in] temperatureRate The ArrayView holding the temperature rate for each element/particle.
   * @param[in] jacobian The ArrayView holding the jacobian for each quadrature point.
   * @param[in] materialDirection The ArrayView holding the material direction for each element/particle.
   * @param[in] lengthScale The ArrayView holding the length scale for each element.
   * @param[in] strengthScale The ArrayView holding the strength scale for each element.
   * @param[in] failureStrength The failure strength.
   * @param[in] maxPrincipalStressDamage Flag to activate max tensile stress damage evolution
   * @param[in] crackSpeed Damage evolution rate
   * @param[in] scaleFractureEnergyReleaseRate legacy flag used to initialize fractureEnergyStrengthScaleExponent when the exponent is not supplied.
   * @param[in] fractureEnergyStrengthScaleExponent exponent eta_G in G_f^{eff}=G_f strengthScale^{eta_G}.
   * @param[in] basalPlaneFractureEnergyReleaseRate Fracture energy release rate for basal plane separation
   * @param[in] totalFractureEnergyReleaseRate Used to determine amount of plastic work needed to generate a damage surface
   * @param[in] damageEvolutionExponent exponent used by the work-based damage law.
   * @param[in] thermalExpansionCoefficient The ArrayView holding the thermal expansion coefficient data for each element.
   * @param[in] newStress The ArrayView holding the new stress data for each quadrature point.
   * @param[in] oldStress The ArrayView holding the old stress data for each quadrature point.
   */
  GraphiteUpdates( real64 const & defaultYoungModulusTransverse,
                   real64 const & defaultYoungModulusAxial,
                   real64 const & defaultPoissonRatioTransverse,
                   real64 const & defaultPoissonRatioAxialTransverse,
                   real64 const & defaultShearModulusAxialTransverse,
                   arrayView1d< real64 > const & effectiveBulkModulus,
                   arrayView1d< real64 > const & effectiveShearModulus,
                   arrayView3d< real64 > const & materialDirection,
                   real64 const defaultYoungModulusTransversePressureDerivative,
                   real64 const defaultYoungModulusAxialPressureDerivative,
                   real64 const defaultShearModulusAxialTransversePressureDerivative,
                   real64 const defaultYoungModulusTransversePressureScale,
                   real64 const defaultYoungModulusAxialPressureScale,
                   real64 const defaultShearModulusAxialTransversePressureScale,
                   arrayView3d< real64 > const & velocityGradient,
                   arrayView3d< real64 > const & plasticStrain,
                   arrayView2d< real64 > const & relaxation,
                   int const & enableCrackTipStressConcentration,
                   arrayView1d< real64 > const & crackTipStressConcentration,
                   arrayView1d< real64 > const & distanceToCrackTip,
                   arrayView2d< real64 > const & basalPlanePlasticWork,
                   arrayView2d< real64 > const & plasticWork,
                   real64 const & alphaL,
                   real64 const & alphaT,
                   arrayView2d< real64 > const & damage,
                   arrayView2d< real64 > const & basalPlaneDamage,
                   arrayView2d< real64 > const & comminutionDamage,
                   int const & enableDistension,
                   arrayView2d< real64 > const & basalNormalDistension,
                   arrayView2d< real64 > const & transverseDistension,
                   arrayView2d< real64 > const & porosity,
                   arrayView1d< real64 > const & temperature,
                   arrayView1d< real64 > const & temperatureRate,
                   arrayView2d< real64 > const & jacobian,
                   arrayView1d< real64 > const & lengthScale,
                   arrayView1d< real64 > const & strengthScale,
                   real64 const & failureStrength,
                   int const & maximumPrincipalStressDamage,
                   real64 const & crackSpeed,
                   int const & scaleFractureEnergyReleaseRate,
                   real64 const & fractureEnergyStrengthScaleExponent,
                   real64 const & basalPlaneFractureEnergyReleaseRate,
                   real64 const & totalFractureEnergyReleaseRate,
                   real64 const & damageEvolutionExponent,
                   real64 const & damagedMaterialFrictionalSlope,
                   real64 const & distortionShearResponseX2,
                   real64 const & distortionShearResponseY1,
                   real64 const & distortionShearResponseY2,
                   real64 const & distortionShearResponseM1,
                   real64 const & positiveDistortionShearResponseX2,
                   real64 const & positiveDistortionShearResponseY1,
                   real64 const & positiveDistortionShearResponseY2,
                   real64 const & positiveDistortionShearResponseM1,
                   real64 const & inPlaneShearResponseX2,
                   real64 const & inPlaneShearResponseY1,
                   real64 const & inPlaneShearResponseY2,
                   real64 const & inPlaneShearResponseM1,
                   real64 const & coupledShearResponseX2,
                   real64 const & coupledShearResponseY1,
                   real64 const & coupledShearResponseY2,
                   real64 const & coupledShearResponseM1,
                   real64 const & distortionStrainHardeningC0,
                   real64 const & inPlaneStrainHardeningC0,
                   real64 const & coupledStrainHardeningC0,
                   real64 const & maximumPlasticStrain,
                   arrayView1d< real64 > const & thermalExpansionCoefficient,
                   arrayView3d< real64, solid::STRESS_USD > const & newStress,
                   arrayView3d< real64, solid::STRESS_USD > const & oldStress,
                   arrayView2d< real64 > const & density,
                   arrayView2d< real64 > const & wavespeed,
                   bool const & disableInelasticity ):
    SolidBaseUpdates( newStress,
                      oldStress,
                      density,
                      wavespeed,
                      thermalExpansionCoefficient,
                      disableInelasticity ),
    m_defaultYoungModulusTransverse( defaultYoungModulusTransverse ),
    m_defaultYoungModulusAxial( defaultYoungModulusAxial ),
    m_defaultPoissonRatioTransverse( defaultPoissonRatioTransverse ),
    m_defaultPoissonRatioAxialTransverse( defaultPoissonRatioAxialTransverse ),
    m_defaultShearModulusAxialTransverse( defaultShearModulusAxialTransverse ),
    m_effectiveBulkModulus( effectiveBulkModulus ),
    m_effectiveShearModulus( effectiveShearModulus ),
    m_materialDirection( materialDirection ),
    m_defaultYoungModulusTransversePressureDerivative( defaultYoungModulusTransversePressureDerivative ),
    m_defaultYoungModulusAxialPressureDerivative( defaultYoungModulusAxialPressureDerivative ),
    m_defaultShearModulusAxialTransversePressureDerivative( defaultShearModulusAxialTransversePressureDerivative ),
    m_defaultYoungModulusTransversePressureScale( defaultYoungModulusTransversePressureScale ),
    m_defaultYoungModulusAxialPressureScale( defaultYoungModulusAxialPressureScale ),
    m_defaultShearModulusAxialTransversePressureScale( defaultShearModulusAxialTransversePressureScale ),
    m_velocityGradient( velocityGradient ),
    m_plasticStrain( plasticStrain ),
    m_relaxation( relaxation ),
    m_enableCrackTipStressConcentration( enableCrackTipStressConcentration ),
    m_crackTipStressConcentration( crackTipStressConcentration ),
    m_distanceToCrackTip( distanceToCrackTip ),
    m_basalPlanePlasticWork( basalPlanePlasticWork ),
    m_plasticWork( plasticWork ),
    m_alphaL( alphaL ),
    m_alphaT( alphaT ),
    m_damage( damage ),
    m_basalPlaneDamage( basalPlaneDamage ),
    m_comminutionDamage( comminutionDamage ),
    m_enableDistension( enableDistension ),
    m_basalNormalDistension( basalNormalDistension ),
    m_transverseDistension( transverseDistension ),
    m_porosity( porosity ),
    m_temperature( temperature ),
    m_temperatureRate( temperatureRate ),
    m_jacobian( jacobian ),
    m_lengthScale( lengthScale ),
    m_strengthScale( strengthScale ),
    m_failureStrength( failureStrength ),
    m_maximumPrincipalStressDamage( maximumPrincipalStressDamage ),
    m_crackSpeed( crackSpeed ),
    m_scaleFractureEnergyReleaseRate( scaleFractureEnergyReleaseRate ),
    m_fractureEnergyStrengthScaleExponent( fractureEnergyStrengthScaleExponent ),
    m_basalPlaneFractureEnergyReleaseRate( basalPlaneFractureEnergyReleaseRate ),
    m_totalFractureEnergyReleaseRate( totalFractureEnergyReleaseRate ),
    m_damageEvolutionExponent( damageEvolutionExponent ),
    m_damagedMaterialFrictionalSlope( damagedMaterialFrictionalSlope ),
    m_distortionShearResponseX2( distortionShearResponseX2 ),
    m_distortionShearResponseY1( distortionShearResponseY1 ),
    m_distortionShearResponseY2( distortionShearResponseY2 ),
    m_distortionShearResponseM1( distortionShearResponseM1 ),
    m_positiveDistortionShearResponseX2( positiveDistortionShearResponseX2 ),
    m_positiveDistortionShearResponseY1( positiveDistortionShearResponseY1 ),
    m_positiveDistortionShearResponseY2( positiveDistortionShearResponseY2 ),
    m_positiveDistortionShearResponseM1( positiveDistortionShearResponseM1 ),
    m_inPlaneShearResponseX2( inPlaneShearResponseX2 ),
    m_inPlaneShearResponseY1( inPlaneShearResponseY1 ),
    m_inPlaneShearResponseY2( inPlaneShearResponseY2 ),
    m_inPlaneShearResponseM1( inPlaneShearResponseM1 ),
    m_coupledShearResponseX2( coupledShearResponseX2 ),
    m_coupledShearResponseY1( coupledShearResponseY1 ),
    m_coupledShearResponseY2( coupledShearResponseY2 ),
    m_coupledShearResponseM1( coupledShearResponseM1 ),
    m_distortionStrainHardeningC0( distortionStrainHardeningC0 ),
    m_inPlaneStrainHardeningC0( inPlaneStrainHardeningC0 ),
    m_coupledStrainHardeningC0( coupledStrainHardeningC0 ),
    m_maximumPlasticStrain( maximumPlasticStrain )
  {}

  /// Default copy constructor
  GraphiteUpdates( GraphiteUpdates const & ) = default;

  /// Default move constructor
  GraphiteUpdates( GraphiteUpdates && ) = default;

  /// Deleted default constructor
  GraphiteUpdates() = delete;

  /// Deleted copy assignment operator
  GraphiteUpdates & operator=( GraphiteUpdates const & ) = delete;

  /// Deleted move assignment operator
  GraphiteUpdates & operator=( GraphiteUpdates && ) =  delete;

  /// Use the uncompressed version of the stiffness bilinear form
  using DiscretizationOps = SolidModelDiscretizationOpsFullyAnisotropic; // TODO: typo in anistropic (fix in DiscOps PR)

  // Bring in base implementations to prevent hiding warnings
  using SolidBaseUpdates::smallStrainUpdate;

  GEOS_HOST_DEVICE
  real64 smoothStep( const real64 x,
                     const real64 xmin,
                     const real64 xmax ) const;

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
                                real64 ( &stress )[6] ) const;

  GEOS_HOST_DEVICE
  void computeTransverselyIsotropicTrialStress( const real64 timeIncrement,
                                                const real64 Ez,
                                                const real64 Ep,
                                                const real64 nuzp,
                                                const real64 nup,
                                                const real64 Gzp,
                                                real64 const (&materialDirection)[3],
                                                real64 const (&oldStress)[6],
                                                real64 const (&D)[6],
                                                real64 ( &newStress ) [6],
                                                localIndex const k ) const;

  GEOS_HOST_DEVICE
  void computeTransverselyIsotropicPlasticStrainIncrement( real64 const ( &velocityGradient )[3][3],
                                                           real64 const ( &oldStress )[6],
                                                           real64 const ( &newStress )[6],
                                                           const real64 Ez,
                                                           const real64 Ep,
                                                           const real64 nuzp,
                                                           const real64 nup,
                                                           const real64 Gzp,
                                                           real64 const ( &materialDirection )[3],
                                                           const real64 timeIncrement,
                                                           real64 ( &plasticStrainIncrement )[6] ) const;

  GEOS_HOST_DEVICE
  void computeTransverselyIsotropicElasticStrainIncrementFromStressIncrement( real64 const ( &stressIncrement )[6],
                                                                              const real64 Ez,
                                                                              const real64 Ep,
                                                                              const real64 nuzp,
                                                                              const real64 nup,
                                                                              const real64 Gzp,
                                                                              real64 const ( &materialDirection )[3],
                                                                              real64 ( &elasticStrainIncrement )[6] ) const;

  GEOS_HOST_DEVICE
  real64 computeTransverselyIsotropicNormalCompliance( const real64 Ez,
                                                       const real64 Ep,
                                                       const real64 nuzp,
                                                       const real64 nup,
                                                       const real64 Gzp,
                                                       real64 const ( &materialDirection )[3],
                                                       real64 const ( &normal )[3] ) const;

  GEOS_HOST_DEVICE
  real64 computeUniaxialPlaneNormalTensileStrength( const real64 strengthScale,
                                                    const real64 distortionHardeningMultiplier ) const;

  GEOS_HOST_DEVICE
  real64 symmetricStressAlongNormal( real64 const ( &stress )[6],
                                     real64 const ( &normal )[3] ) const;

  GEOS_HOST_DEVICE
  real64 computeMaximumPrincipalStressAndDirection( real64 const ( &stress )[6],
                                                    real64 ( &normal )[3] ) const;

  GEOS_HOST_DEVICE
  void subtractNormalStressComponent( const real64 stressDrop,
                                      real64 const ( &normal )[3],
                                      real64 ( &stress )[6] ) const;

  GEOS_HOST_DEVICE
  real64 limitDamageIncrementByCrackSpeed( localIndex const k,
                                           const real64 timeIncrement,
                                           const real64 oldDamage,
                                           const real64 computedDamage ) const;

  GEOS_HOST_DEVICE
  real64 computeTransverseSpectralStrainPart( real64 const ( &strainDense )[3][3],
                                              real64 const ( &materialDirection )[3],
                                              const bool compressivePart,
                                              real64 ( &strainPartDense )[3][3] ) const;

  GEOS_HOST_DEVICE
  bool applyTransverselyIsotropicDistensionClosure( localIndex const k,
                                                    localIndex const q,
                                                    const real64 timeIncrement,
                                                    real64 const ( &materialDirection )[3],
                                                    real64 ( &DForStress )[6],
                                                    real64 ( &velocityGradientForStress )[3][3],
                                                    real64 ( &closurePlasticStrainIncrement )[6] ) const;

  GEOS_HOST_DEVICE
  void addUnrotatedPlasticStrainStateIncrement( localIndex const k,
                                                localIndex const q,
                                                real64 const ( &rotationTranspose )[3][3],
                                                real64 const ( &endRotation )[3][3],
                                                real64 const ( &plasticStrainStateIncrement )[6] ) const;

  GEOS_HOST_DEVICE
  void updateDistensionPorosity( localIndex const k,
                                 localIndex const q ) const;

  GEOS_HOST_DEVICE
  void addDistensionFromTensileOpening( localIndex const k,
                                        localIndex const q,
                                        const real64 openingStrain,
                                        real64 const ( &openingNormal )[3],
                                        real64 const ( &materialDirection )[3],
                                        const bool basalOpening ) const;

  GEOS_HOST_DEVICE
  void addDistensionFromInelasticStrain( localIndex const k,
                                         localIndex const q,
                                         real64 const ( &materialDirection )[3],
                                         real64 const ( &inelasticStrainIncrement )[6] ) const;

  GEOS_HOST_DEVICE
  bool applyEnergyRegularizedBrittleTensileReturn( localIndex const k,
                                                   localIndex const q,
                                                   const real64 timeIncrement,
                                                   const real64 stressMeasure,
                                                   const real64 oldDamage,
                                                   const real64 effectiveStrength,
                                                   const real64 fractureEnergyReleaseRate,
                                                   const real64 normalCompliance,
                                                   real64 const ( &normal )[3],
                                                   real64 const ( &materialDirection )[3],
                                                   const bool updateBasalPlaneDamage,
                                                   real64 ( &stress )[6] ) const;

  GEOS_HOST_DEVICE
  bool applyEnergyRegularizedBrittleTensileReturns( localIndex const k,
                                                    localIndex const q,
                                                    const real64 timeIncrement,
                                                    const real64 principalTensileStrength,
                                                    const real64 basalNormalTensileStrength,
                                                    const real64 basalPlaneFractureEnergyReleaseRate,
                                                    const real64 totalFractureEnergyReleaseRate,
                                                    const real64 Ez,
                                                    const real64 Ep,
                                                    const real64 nuzp,
                                                    const real64 nup,
                                                    const real64 Gzp,
                                                    real64 const ( &materialDirection )[3],
                                                    real64 ( &stress )[6] ) const;

  GEOS_HOST_DEVICE
  real64 transverselyIsotropicB1( real64 const (&materialDirection)[3],
                                  const int i,
                                  const int j,
                                  const int p,
                                  const int w ) const;


  GEOS_HOST_DEVICE
  real64 transverselyIsotropicB2( real64 const (&materialDirection)[3],
                                  const int i,
                                  const int j,
                                  const int p,
                                  const int w ) const;


  GEOS_HOST_DEVICE
  real64 transverselyIsotropicB3( real64 const (&materialDirection)[3],
                                  const int i,
                                  const int j,
                                  const int p,
                                  const int w ) const;


  GEOS_HOST_DEVICE
  real64 transverselyIsotropicB4( real64 const (&materialDirection)[3],
                                  const int i,
                                  const int j,
                                  const int p,
                                  const int w ) const;


  GEOS_HOST_DEVICE
  real64 transverselyIsotropicB5( real64 const (&materialDirection)[3],
                                  const int i,
                                  const int j,
                                  const int p,
                                  const int w ) const;

  GEOS_HOST_DEVICE
  real64 delta( const int i,
                const int j )  const;

  GEOS_HOST_DEVICE
  real64 slopePoint0( const real64 x,
                      const real64 x1,
                      const real64 x2,
                      const real64 y1,
                      const real64 y2,
                      const real64 m1 ) const;


  GEOS_HOST_DEVICE
  real64 evaluateSaturatingPressureArgument( const real64 pressureArgument,
                                             const real64 pressureScale ) const;

  GEOS_HOST_DEVICE
  real64 evaluatePressureDependentStrength( const real64 pressure,
                                            const real64 x2,
                                            const real64 y1,
                                            const real64 y2,
                                            const real64 m1 ) const;

  GEOS_HOST_DEVICE
  real64 evaluateGraphiteModeStrength( const real64 pressure,
                                        const real64 damage,
                                        const real64 strengthScale,
                                        const real64 x2,
                                        const real64 y1,
                                        const real64 y2,
                                        const real64 m1,
                                        const real64 strainHardeningMultiplier ) const;

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  virtual void saveConvergedState( localIndex const k,
                                   localIndex const q ) const override final
  {
    SolidBaseUpdates::saveConvergedState( k, q );
  }

private:
  real64 const m_defaultYoungModulusTransverse;

  real64 const m_defaultYoungModulusAxial;

  real64 const m_defaultPoissonRatioTransverse;

  real64 const m_defaultPoissonRatioAxialTransverse;

  real64 const m_defaultShearModulusAxialTransverse;

  arrayView1d< real64 > const m_effectiveBulkModulus;

  arrayView1d< real64 > const m_effectiveShearModulus;

  arrayView3d< real64 > const m_materialDirection;

  real64 const m_defaultYoungModulusTransversePressureDerivative;

  real64 const m_defaultYoungModulusAxialPressureDerivative;

  real64 const m_defaultShearModulusAxialTransversePressureDerivative;

  real64 const m_defaultYoungModulusTransversePressureScale;

  real64 const m_defaultYoungModulusAxialPressureScale;

  real64 const m_defaultShearModulusAxialTransversePressureScale;

  /// A reference to the ArrayView holding the velocity gradient for each element/particle.
  arrayView3d< real64 > const m_velocityGradient;

  /// A reference to the ArrayView holding the damage for each quadrature point.
  arrayView3d< real64 > const m_plasticStrain;

  /// A reference to the ArrayView holding the damage for each quadrature point.
  arrayView2d< real64 > const m_relaxation;

  /// Model parameter: Flag to enable stress concenration for crack-tip
  int m_enableCrackTipStressConcentration;
  ///State variable: Lets us plot the stress concentration, not actually needed by the algorithm.
  arrayView1d< real64 > const m_crackTipStressConcentration;
  /// A reference to the ArrayView holding the distance to crack tip.
  arrayView1d< real64 > const m_distanceToCrackTip;


  // A reference to the ArrayView holding the accumulated plastic work for each quadrature point.
  arrayView2d< real64 > const m_basalPlanePlasticWork;
  arrayView2d< real64 > const m_plasticWork;

  // thermal expansion in direction lateral to crystal symmetry
  real64 const m_alphaL;

  // thermal expansion in direction transverse to crystal symmetry
  real64 const m_alphaT;

  /// Host-visible DFG damage for each quadrature point. This scalar is intentionally
  /// allowed to reach one for either fracture or comminution so the DFG solver can
  /// split velocity fields, open/slip surfaces, and separate powder clouds.
  arrayView2d< real64 > const m_damage;

  /// Mode-resolved basal-plane fracture damage. This internal variable degrades
  /// basal opening/sliding-related strength without forcing the whole material
  /// point to behave as a fully comminuted powder.
  arrayView2d< real64 > const m_basalPlaneDamage;

  /// Internal comminution/powder damage. This variable blends the crystal strength
  /// to the residual pressure-dependent powder envelope for all shear modes.
  arrayView2d< real64 > const m_comminutionDamage;

  /// Flag to activate stress-free closure of transversely isotropic distension reservoirs.
  int const m_enableDistension;

  /// Closeable distension strain normal to the basal planes.
  arrayView2d< real64 > const m_basalNormalDistension;

  /// Closeable distension strain in the basal/transverse plane.
  arrayView2d< real64 > const m_transverseDistension;

  /// Distension-derived porosity for plotting/coupling, phi = 1 - exp[-(zeta_N+zeta_T)].
  arrayView2d< real64 > const m_porosity;

  /// A reference to the ArrayView holding the temperature for each quadrature point.
  arrayView1d< real64 > const m_temperature;

  /// A reference to the ArrayView holding the temperature for each quadrature point.
  arrayView1d< real64 > const m_temperatureRate;

  /// A reference to the ArrayView holding the jacobian for each quadrature point.
  arrayView2d< real64 > const m_jacobian;

  /// A reference to the ArrayView holding the length scale for each element/particle.
  arrayView1d< real64 > const m_lengthScale;

  /// Discretization-sized variable: The strength scale for each element/particle
  arrayView1d< real64 > const m_strengthScale;

  /// The maximum theoretical strength
  real64 const m_failureStrength;

  // Energy regularization of fracture evolution
  int const m_maximumPrincipalStressDamage;
  real64 const m_crackSpeed;

  int const m_scaleFractureEnergyReleaseRate;   // Legacy flag used to initialize the fracture-energy scale exponent.
  real64 const m_fractureEnergyStrengthScaleExponent;
  real64 const m_basalPlaneFractureEnergyReleaseRate;
  real64 const m_totalFractureEnergyReleaseRate;
  real64 const m_damageEvolutionExponent;

  // The damaged material frictional slope
  real64 const m_damagedMaterialFrictionalSlope;

  /// Material parameters: The values controlling the failure envelope.
  /// The distortionShearResponse* inputs define the negative signed
  /// distortion branch I_d^- = max(-I_d,0), where
  /// I_d = sigma_aa - (sigma_bb + sigma_cc)/2.
  real64 const m_distortionShearResponseX2;
  real64 const m_distortionShearResponseY1;
  real64 const m_distortionShearResponseY2;
  real64 const m_distortionShearResponseM1;

  /// Optional positive signed-distortion branch I_d^+ = max(I_d,0).
  /// When omitted, these parameters default to the negative branch values,
  /// which gives the unsigned distortion response |I_d|/Y_d.
  real64 const m_positiveDistortionShearResponseX2;
  real64 const m_positiveDistortionShearResponseY1;
  real64 const m_positiveDistortionShearResponseY2;
  real64 const m_positiveDistortionShearResponseM1;

  real64 const m_inPlaneShearResponseX2;
  real64 const m_inPlaneShearResponseY1;
  real64 const m_inPlaneShearResponseY2;
  real64 const m_inPlaneShearResponseM1;

  real64 const m_coupledShearResponseX2;
  real64 const m_coupledShearResponseY1;
  real64 const m_coupledShearResponseY2;
  real64 const m_coupledShearResponseM1;

  real64 const m_distortionStrainHardeningC0;
  real64 const m_inPlaneStrainHardeningC0;
  real64 const m_coupledStrainHardeningC0;

  real64 const m_maximumPlasticStrain;

};

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 GraphiteUpdates::smoothStep( const real64 x,
                                    const real64 xmin,
                                    const real64 xmax ) const
{
  // Smooth blending function from 0 to 1 as
  // x goes from xmin to xmax.
  //
  // will fail if xmax=xmin, so don't do that.

  if( x <= xmin )
  {
    return 0.0;
  }
  else if( x >= xmax )
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
void GraphiteUpdates::smallStrainUpdate( localIndex const k,
                                         localIndex const q,
                                         real64 const & timeIncrement,
                                         real64 const ( &strainIncrement )[6],
                                         real64 ( & stress )[6],
                                         real64 ( & stiffness )[6][6] ) const
{
  GEOS_UNUSED_VAR( k );
  GEOS_UNUSED_VAR( q );
  GEOS_UNUSED_VAR( timeIncrement );
  GEOS_UNUSED_VAR( strainIncrement );
  GEOS_UNUSED_VAR( stress );
  GEOS_UNUSED_VAR( stiffness );
  GEOS_ERROR( "smallStrainUpdate not implemented for Graphite model" );
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void GraphiteUpdates::smallStrainUpdate( localIndex const k,
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
void GraphiteUpdates::smallStrainUpdate_StressOnly( localIndex const k,
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
  GEOS_ERROR( "smallStrainUpdateStressOnly overload not implemented for Graphite model" );
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void GraphiteUpdates::smallStrainUpdate_StressOnly( localIndex const k,
                                                    localIndex const q,
                                                    real64 const & timeIncrement,
                                                    real64 const ( &beginningRotation )[3][3],
                                                    real64 const ( &endRotation )[3][3],
                                                    real64 const ( &strainIncrement )[6],
                                                    real64 ( & stress )[6] ) const
{
  // elastic predictor (assume strainIncrement is all elastic)
  // ElasticTransverseIsotropicPressureDependentUpdates::smallStrainUpdate_StressOnly( k,
  //                                                                                   q,
  //                                                                                   timeIncrement,
  //                                                                                   beginningRotation,
  //                                                                                   endRotation,
  //                                                                                   strainIncrement,
  //                                                                                   stress );

  m_jacobian[k][q] *= LvArray::math::exp( strainIncrement[0] + strainIncrement[1] + strainIncrement[2] );

  // if( m_disableInelasticity )
  // {
  //   return;
  // }

  // call the constitutive model
  GraphiteUpdates::smallStrainUpdateHelper( k,
                                            q,
                                            timeIncrement,
                                            beginningRotation,
                                            endRotation,
                                            stress );

  // save new stress and return
  saveStress( k, q, stress );

  // CC: debug
  // GEOS_LOG_RANK( "Particle " << k << ", Saved stress: {" << stress[0] << ", " << stress[1] << ", " << stress[2] << ", " << stress[3] <<
  // ", " << stress[4] << ", " << stress[5] << "}" );
}


GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void GraphiteUpdates::smallStrainUpdateHelper( localIndex const k,
                                               localIndex const q,
                                               real64 const timeIncrement,
                                               real64 const ( &beginningRotation )[3][3],
                                               real64 const ( &endRotation )[3][3],
                                               real64 ( & stress )[6] ) const
{
  // CC: debug
  // GEOS_LOG_RANK( "Particle " << k << ", Stress in: {" << stress[0] << ", " << stress[1] << ", " << stress[2] << ", " << stress[3] << ", "
  // << stress[4] << ", " << stress[5] << "}" );
  real64 oldStress[6] = {};
  LvArray::tensorOps::copy< 6 >( oldStress, m_oldStress[k][q] );  //stress );
  // CC: debug
  // GEOS_LOG_RANK( "Particle " << k << ", Old stress copy: {" << oldStress[0] << ", " << oldStress[1] << ", " << oldStress[2] << ", " <<
  // oldStress[3] << ", " << oldStress[4] << ", " << oldStress[5] << "}" );

  real64 rotationTranspose[3][3] = {};
  LvArray::tensorOps::transpose< 3, 3 >( rotationTranspose, beginningRotation );

  real64 tempMat[ 3 ][ 3 ]= {};

  real64 unrotatedVelocityGradient[3][3] = {};
  LvArray::tensorOps::Rij_eq_AkiBkj< 3, 3, 3 >( tempMat, beginningRotation, m_velocityGradient[k] );
  LvArray::tensorOps::Rij_eq_AikBkj< 3, 3, 3 >( unrotatedVelocityGradient, tempMat, beginningRotation );

  real64 unrotatedVelocityGradientTranspose[3][3] = {};
  LvArray::tensorOps::transpose< 3, 3 >( unrotatedVelocityGradientTranspose, unrotatedVelocityGradient );

  // CC: Is there an LvArray operation to get the symmetric part of a matrix?
  real64 denseD[3][3] = {};
  LvArray::tensorOps::copy< 3, 3 >( denseD, unrotatedVelocityGradient );
  LvArray::tensorOps::add< 3, 3 >( denseD, unrotatedVelocityGradientTranspose );
  LvArray::tensorOps::scale< 3, 3 >( denseD, 0.5 );

  real64 D[6] = {};
  LvArray::tensorOps::denseToSymmetric< 3 >( D, denseD );

  // make sure material direction is normalized.
  // Check updates to material direction are correct, this model only utilizes the x basis to denote the anisotropic c-axis of the crystal
  real64 materialDirection[3] = {m_materialDirection[k][0][0],
                                 m_materialDirection[k][0][1],
                                 m_materialDirection[k][0][2]};

  // Should be normalized already:
  LvArray::tensorOps::normalize< 3 >( materialDirection );

  // Unrotate material direction
  real64 unrotatedMaterialDirection[3] = {};
  LvArray::tensorOps::Ri_eq_AijBj< 3, 3 >( unrotatedMaterialDirection, rotationTranspose, materialDirection );

  real64 unrotatedVelocityGradientForStress[3][3] = {};
  LvArray::tensorOps::copy< 3, 3 >( unrotatedVelocityGradientForStress, unrotatedVelocityGradient );

  real64 DForStress[6] = {};
  LvArray::tensorOps::copy< 6 >( DForStress, D );

  real64 distensionClosurePlasticStrainIncrement[6] = {};
  bool distensionClosureOccurred = false;
  bool const distensionClosurePossible =
    ( m_enableDistension == 1 ) &&
    ( !m_disableInelasticity ) &&
    ( ( m_basalNormalDistension[k][q] > 0.0 ) ||
      ( m_transverseDistension[k][q] > 0.0 ) );
  if( distensionClosurePossible )
  {
    distensionClosureOccurred = applyTransverselyIsotropicDistensionClosure( k,
                                                                             q,
                                                                             timeIncrement,
                                                                             unrotatedMaterialDirection,
                                                                             DForStress,
                                                                             unrotatedVelocityGradientForStress,
                                                                             distensionClosurePlasticStrainIncrement );
  }

  // Use beginning of step normal stress to compute stress dependence of Ez
  real64 temp[3] = {};
  int voigtMap[3][3] = { {0, 5, 4}, {5, 1, 3}, {4, 3, 2} };
  LvArray::tensorOps::Ri_eq_symAijBj< 3 >( temp, oldStress, unrotatedMaterialDirection );
  real64 oldPlaneNormalStress = LvArray::tensorOps::AiBi< 3 >( unrotatedMaterialDirection, temp );    // Used (sometimes) to integrate plastic work for damage criteria

  // Beginning of step pressure to compute pressure-dependence of elastic moduli
  real64 oldPressure = (-1.0/3.0)*( oldStress[0] + oldStress[1] + oldStress[2] );   // CC: Unused?

  // This is a transversely isotropic material for graphite-like crystals having
  // some weak plane with plane-normal-stress- and pressure-dependent elastic properties:

  // The elastic moduli use non-negative pressure-like arguments from the
  // beginning-of-step stress.  A rational saturation function can be applied
  // independently to each modulus,
  //
  //     f(xi; P) = xi / (1 + xi/P),
  //
  // so that finite P gives monotone pressure stiffening with a bounded
  // derivative, while P=DBL_MAX recovers the original linear dependence.
  // The axial modulus uses a mixed mean-pressure/plane-normal argument so
  // that compression normal to the graphitic planes can stiffen the axial
  // response more directly than purely in-plane confinement.
  real64 const axialPressureArgument = LvArray::math::max( 0.0, -0.5*oldPlaneNormalStress + 0.5*oldPressure );
  real64 const meanPressureArgument = LvArray::math::max( 0.0, oldPressure );

  real64 Ez = LvArray::math::max( 0.001*m_defaultYoungModulusAxial,
                                  m_defaultYoungModulusAxial +
                                  m_defaultYoungModulusAxialPressureDerivative *
                                  evaluateSaturatingPressureArgument( axialPressureArgument,
                                                                      m_defaultYoungModulusAxialPressureScale ) );
  real64 Ep = LvArray::math::max( 0.001*m_defaultYoungModulusTransverse,
                                  m_defaultYoungModulusTransverse +
                                  m_defaultYoungModulusTransversePressureDerivative *
                                  evaluateSaturatingPressureArgument( meanPressureArgument,
                                                                      m_defaultYoungModulusTransversePressureScale ) );
  real64 Gzp  = LvArray::math::max( 0.001*m_defaultShearModulusAxialTransverse,
                                    m_defaultShearModulusAxialTransverse +
                                    m_defaultShearModulusAxialTransversePressureDerivative *
                                    evaluateSaturatingPressureArgument( meanPressureArgument,
                                                                        m_defaultShearModulusAxialTransversePressureScale ) );
  real64 nuzp = LvArray::math::min( 0.4999, m_defaultPoissonRatioAxialTransverse );
  real64 nup = LvArray::math::min( 0.4999, m_defaultPoissonRatioTransverse );

  // Maintain positive definiteness of the actual pressure-dependent TI
  // stiffness used in the update.  The input checks guarantee that the
  // reference stiffness is admissible, but independent pressure dependence of
  // Ez and Ep can otherwise violate
  //
  //     1 - nu_p - 2*(Ep/Ez)*nu_zp^2 > 0
  //
  // at high pressure.  If needed, cap Ep just below the admissible limit while
  // preserving the monotone pressure-stiffening trends of the other moduli.
  real64 const stabilityTolerance = 1.0e-12;
  real64 const nuzpSquared = nuzp*nuzp;
  if( nuzpSquared > 0.0 )
  {
    real64 const maximumStableEp = LvArray::math::max( 1.0e-12*m_defaultYoungModulusTransverse,
                                                       ( 1.0 - nup - stabilityTolerance ) * Ez /
                                                       ( 2.0*nuzpSquared ) );
    Ep = LvArray::math::min( Ep, maximumStableEp );
  }

  // Update effective elastic properties.  The denominator should be negative
  // for a positive-definite TI stiffness; the guard avoids a non-positive or
  // singular effective bulk modulus in wave-speed diagnostics if a parameter
  // set is very near the stability boundary.
  real64 const bulkDenominator = 2*Ez*( nup+nuzp-1 ) + Ep*( 2*nuzp-1 );
  real64 const safeBulkDenominator = LvArray::math::min( bulkDenominator, -1.0e-20 );
  m_effectiveBulkModulus[k] = -Ep*Ez/safeBulkDenominator;
  m_effectiveShearModulus[k] = 0.6*m_effectiveBulkModulus[k];

  // TODO: make elastic properties safe to avoid floating point errors.
  m_wavespeed[k][0] = LvArray::math::sqrt( ( m_effectiveBulkModulus[k] + (4.0/3.0) * m_effectiveShearModulus[k] ) / m_density[k][0] );


  // CC: debug
  // GEOS_LOG_RANK( "Particle " << k << ":\n"
  //                 "Old pressure: " << oldPressure << "\n"
  //                 "\tEz=" << Ez << "\n" <<
  //                 "\tEp=" << Ep << "\n" <<
  //                 "\tGzp=" << Gzp << "\n" <<
  //                 "\tnuzp=" << nuzp << "\n" <<
  //                 "\tnup=" << nup << "\n" <<
  //                 "\teffK=" << m_effectiveBulkModulus[k] << "\n" <<
  //                 "\teffG=" << m_effectiveShearModulus[k] << "\n");

  // real64 Ez = c33 - c13 * c13 / ( c11 - c66 );
  // real64 Ep = 4 * c66 * ( c11 * c33 - c66 * c33 - c13 * c13 ) / ( c11 * c33 - c13 * c13 );
  // real64 Gzp = c44 / 2.0;
  // real64 nuzp = c13 / ( 2 * ( c11 - c66 ) );
  // real64 nup = Ep / c66 - 1; //4 * ( c11 * c33 - c66 * c33 - c13 * c13 ) / ( c11 * c33 - c13 * c13 ) - 1;

  // Hypoelastic trial stress update.
  computeTransverselyIsotropicTrialStress( timeIncrement,        // time step
                                           Ez,                   // Elastic modulus preferred direction
                                           Ep,                   // Elastic modulus transverse plane
                                           nuzp,                 // Poisson ratio coupled
                                           nup,                  // Poisson ratio transverse plane
                                           Gzp,                  // Shear modulus coupled plane
                                           unrotatedMaterialDirection,     // preferred direction
                                           oldStress,            // stress at start of step
                                           DForStress,           // D=sym(L) with any stress-free distension closure removed
                                           stress,               // stress at end of step
                                           k );


  // CC: debug
  // GEOS_LOG_RANK( "Particle " << k << ", Trial stress: {" << stress[0] << ", " << stress[1] << ", " << stress[2] << ", " << stress[3] <<
  // ", " << stress[4] << ", " << stress[5] << "}" );

  // m_jacobian[k][q] *= LvArray::math::exp( strainIncrement[0] + strainIncrement[1] + strainIncrement[2] );

  if( m_disableInelasticity )
  {
    return;
  }

  // strengthScale is a particle-scale strength multiplier supplied by the
  // host code.  It modifies low-pressure flaw-controlled strengths and slopes,
  // but it does not modify the high-pressure Y2 plateaus.  The Y2 values
  // represent confined strength after defects and microcracks close and are
  // shared by the intact and comminuted-powder branches.  Fracture energy
  // scaling is controlled independently through
  // G_f^eff = G_f * strengthScale^eta_G.
  real64 const strengthScale = LvArray::math::max( m_strengthScale[k], 1.0e-12 );
  real64 const fractureEnergyScale = LvArray::math::pow( strengthScale, m_fractureEnergyStrengthScaleExponent );
  real64 const basalPlaneFractureEnergyReleaseRate =
    ( m_basalPlaneFractureEnergyReleaseRate < DBL_MAX ) ?
    m_basalPlaneFractureEnergyReleaseRate * fractureEnergyScale : DBL_MAX;
  real64 const totalFractureEnergyReleaseRate =
    ( m_totalFractureEnergyReleaseRate < DBL_MAX ) ?
    m_totalFractureEnergyReleaseRate * fractureEnergyScale : DBL_MAX;
  real64 const failureStrength = m_failureStrength * strengthScale;

  real64 const distortionHardeningMultiplier =
    1.0 + ( m_distortionStrainHardeningC0 - 1.0 ) * smoothStep( m_relaxation[k][q], 0.0, 1.0 );

  bool const basalEnergyRegularizedBrittle = basalPlaneFractureEnergyReleaseRate < DBL_MAX;
  bool const principalEnergyRegularizedBrittle =
    ( m_maximumPrincipalStressDamage == 1 ) && ( totalFractureEnergyReleaseRate < DBL_MAX );

  real64 const basalNormalTensileStrength = basalEnergyRegularizedBrittle ?
    computeUniaxialPlaneNormalTensileStrength( strengthScale, distortionHardeningMultiplier ) : 0.0;

  real64 brittleStrainIncrement[6] = {};

  // If the particle is a crack-tip particle, the distanceToCrackTip will be greater than 0, and we compute the
  // stress concentration.  We don't actually need to store this as a state variable, it is sufficient to store
  // the distanceToCrackTip, but we've added this field to allow plotting of the stress concentration.
  m_crackTipStressConcentration[k] = 1.0;
  real64 referenceFractureEnergyReleaseRate = DBL_MAX;
  if( ( basalPlaneFractureEnergyReleaseRate < DBL_MAX ) && ( totalFractureEnergyReleaseRate < DBL_MAX ) )
  {
    referenceFractureEnergyReleaseRate = 0.5*( basalPlaneFractureEnergyReleaseRate + totalFractureEnergyReleaseRate );
  }
  else if( totalFractureEnergyReleaseRate < DBL_MAX )
  {
    referenceFractureEnergyReleaseRate = totalFractureEnergyReleaseRate;
  }
  else if( basalPlaneFractureEnergyReleaseRate < DBL_MAX )
  {
    referenceFractureEnergyReleaseRate = basalPlaneFractureEnergyReleaseRate;
  }

  if( ( m_enableCrackTipStressConcentration == 1 ) &&
      ( m_distanceToCrackTip[k] > 0.0 ) &&
      ( referenceFractureEnergyReleaseRate < DBL_MAX ) )
  {
    real64 const constrainedModulus =  m_effectiveBulkModulus[k] + (4.0/3.0) * m_effectiveShearModulus[k];
    real64 const nominalIntactStrength = failureStrength;

    real64 const fractureProcessZoneRadius =
      LvArray::math::max( 1.e-12, constrainedModulus * referenceFractureEnergyReleaseRate /( 6.283185307179586 * LvArray::math::max( 1.e-12, nominalIntactStrength * nominalIntactStrength ) ) );
    // distanceToCrackTip is the inverse-kernel DFG estimate of the
    // unresolved crack-tip distance. If that distance is smaller than the
    // fracture-process-zone radius, the resolved strength is unchanged. If it
    // is larger than the fracture-process-zone radius, reduce the effective
    // strength by the unresolved LEFM stress-concentration factor.
    m_crackTipStressConcentration[k] =
      LvArray::math::max( 1.0,
                          LvArray::math::sqrt( m_distanceToCrackTip[k] /
                                                fractureProcessZoneRadius ) );
  }

  if( basalEnergyRegularizedBrittle || principalEnergyRegularizedBrittle )
  {
    real64 trialStress[6] = {};
    LvArray::tensorOps::copy< 6 >( trialStress, stress );

    bool const brittleReturnOccurred =
      applyEnergyRegularizedBrittleTensileReturns( k,
                                                   q,
                                                   timeIncrement,
                                                   failureStrength,
                                                   basalNormalTensileStrength,
                                                   basalPlaneFractureEnergyReleaseRate,
                                                   totalFractureEnergyReleaseRate,
                                                   Ez,
                                                   Ep,
                                                   nuzp,
                                                   nup,
                                                   Gzp,
                                                   unrotatedMaterialDirection,
                                                   stress );
    if( brittleReturnOccurred )
    {
      real64 brittleStressDrop[6] = {};
      LvArray::tensorOps::copy< 6 >( brittleStressDrop, trialStress );
      LvArray::tensorOps::subtract< 6 >( brittleStressDrop, stress );
      computeTransverselyIsotropicElasticStrainIncrementFromStressIncrement( brittleStressDrop,
                                                                             Ez,
                                                                             Ep,
                                                                             nuzp,
                                                                             nup,
                                                                             Gzp,
                                                                             unrotatedMaterialDirection,
                                                                             brittleStrainIncrement );
    }
  }

  // Check for legacy tensile failure in the preferred direction only when the
  // basal normal-opening fracture energy is disabled.  With finite basal
  // fracture energy, normal opening is handled by the pre-plastic brittle
  // tensile return whose peak strength is derived from the positive-distortion
  // yield surface.
  LvArray::tensorOps::Ri_eq_symAijBj< 3 >( temp, stress, unrotatedMaterialDirection );
  real64 planeNormalStress = LvArray::tensorOps::AiBi< 3 >( unrotatedMaterialDirection, temp );
  int basalModeIFracture = 0;
  if( ( !basalEnergyRegularizedBrittle ) && ( planeNormalStress > failureStrength ) )
  {
    basalModeIFracture = 1;
  }

  bool const legacyPlasticRateDamage =
    ( m_maximumPrincipalStressDamage == 0 ) &&
    !( basalPlaneFractureEnergyReleaseRate < DBL_MAX ) &&
    !( totalFractureEnergyReleaseRate < DBL_MAX ) &&
    ( m_crackSpeed < DBL_MAX );

  // Test if damage is exactly 0 but there is significant accumulated plastic work.
  // This should only occur if the damage has been healed by an event in the solver,
  // in which case we need to also reset the accumulated plastic work so it
  // doesn't immediately re-damage.
  //
  if ( isZero( m_damage[k][q] ) && ( (m_plasticWork[k][q] > ( 1.e-3*totalFractureEnergyReleaseRate / m_lengthScale[k] ) ) ||
                                    (m_basalPlanePlasticWork[k][q] > ( 1.e-3*basalPlaneFractureEnergyReleaseRate / m_lengthScale[k] ) ) ) )
                                    {
                                      m_plasticWork[k][q] = 0.0;
                                      m_basalPlanePlasticWork[k][q] = 0.0;
                                      m_basalPlaneDamage[k][q] = 0.0;
                                      m_comminutionDamage[k][q] = 0.0;
                                    }

  // The public damage field is the host-visible failure variable consumed by
  // DFG.  It is allowed to reach one for either a resolved fracture surface or a
  // comminuted powder cloud.  The constitutive strength degradation below uses
  // the mode-resolved internal damage variables instead: basalPlaneDamage
  // weakens basal opening/sliding modes, while comminutionDamage blends all
  // shear modes toward a pressure-dependent powder residual.  This separation
  // lets DFG split velocity fields and prevent damaged/intact smearing without
  // forcing every fractured particle to become a zero-strength material.
  real64 const basalPlaneDamage = LvArray::math::max( 0.0, LvArray::math::min( 1.0, m_basalPlaneDamage[k][q] ) );
  real64 const comminutionDamage = LvArray::math::max( 0.0, LvArray::math::min( 1.0, m_comminutionDamage[k][q] ) );
  real64 const dfgDamage = LvArray::math::max( LvArray::math::max( basalPlaneDamage, comminutionDamage ),
                                               LvArray::math::max( 0.0, LvArray::math::min( 1.0, m_damage[k][q] ) ) );
  m_damage[k][q] = dfgDamage;

  real64 const tensileStressMultiplier = 1.0 - smoothStep( dfgDamage, 0.0, 1.0 );

  // Decompose the stress tensor after any pre-plastic brittle tensile return.
  real64 sigma1Dense[3][3] = {};
  real64 sigma2Dense[3][3] = {};
  real64 sigma4Dense[3][3] = {};
  real64 sigma5Dense[3][3] = {};
  for( int i=0; i<3; i++ )
  {
    for( int j=0; j<3; j++ )
    {
      for( int p=0; p<3; p++ )
      {
        for( int w=0; w<3; w++ )
        {
          sigma1Dense[i][j] += transverselyIsotropicB1( unrotatedMaterialDirection, i, j, p, w )*stress[voigtMap[p][w]];
          sigma2Dense[i][j] += transverselyIsotropicB2( unrotatedMaterialDirection, i, j, p, w )*stress[voigtMap[p][w]];
          sigma4Dense[i][j] += transverselyIsotropicB4( unrotatedMaterialDirection, i, j, p, w )*stress[voigtMap[p][w]];
          sigma5Dense[i][j] += transverselyIsotropicB5( unrotatedMaterialDirection, i, j, p, w )*stress[voigtMap[p][w]];
        }
      }
    }
  }

  real64 sigma1[6] = {0};   // axial
  real64 sigma2[6] = {0};   // in-plane normal
  real64 sigma4[6] = {0};   // in-plane total stress
  real64 sigma5[6] = {0};   // weak plane - shear
  LvArray::tensorOps::denseToSymmetric< 3 >( sigma1, sigma1Dense );
  LvArray::tensorOps::denseToSymmetric< 3 >( sigma2, sigma2Dense );
  LvArray::tensorOps::denseToSymmetric< 3 >( sigma4, sigma4Dense );
  LvArray::tensorOps::denseToSymmetric< 3 >( sigma5, sigma5Dense );

  // Post-brittle trial pressure to compute pressure-dependence of strength.
  real64 pressure = (-1.0/3.0)*( stress[0] + stress[1] + stress[2] );

  // flags indicating whether plastic strain needs to be updated.
  bool plastic = false;

  // find in-plane isotropic and deviatoric stress
  real64 inPlaneIso[6];
  LvArray::tensorOps::copy< 6 >( inPlaneIso, sigma2 );
  LvArray::tensorOps::scale< 6 >( inPlaneIso, 0.5 );

  real64 inPlaneDev[6];
  LvArray::tensorOps::copy< 6 >( inPlaneDev, sigma4 );
  LvArray::tensorOps::subtract< 6 >( inPlaneDev, inPlaneIso );

  // Plane-normal tension is limited by the intact Weibull-scaled tensile
  // strength and is then continuously degraded by damage.  Compression is left
  // intact so that damaged material can compact and develop frictional strength.
  if( ( !basalEnergyRegularizedBrittle ) && ( planeNormalStress > 0.0 ) )
  {
    real64 const tensileCapMultiplier =
      LvArray::math::min( 1.0, failureStrength / LvArray::math::max( planeNormalStress, 1.0e-20 ) );
    real64 const planeNormalMultiplier = tensileStressMultiplier * tensileCapMultiplier;
    if( planeNormalMultiplier < 0.999999999999 )
    {
      LvArray::tensorOps::scale< 6 >( sigma1, planeNormalMultiplier );
      plastic = true;
    }
  }

  // find distortion stress from in-plane isotropic and plane-normal stress
  real64 distortion[6] = {};
  LvArray::tensorOps::copy< 6 >( distortion, inPlaneIso );
  LvArray::tensorOps::add< 6 >( distortion, sigma1 );

  real64 distortionMeanStress = ( distortion[0] + distortion[1] + distortion[2] ) / 3.0;

  real64 distortion_iso[6] = {0};
  distortion_iso[0] = distortionMeanStress;
  distortion_iso[1] = distortionMeanStress;
  distortion_iso[2] = distortionMeanStress;

  real64 distortion_dev[6] = {};
  LvArray::tensorOps::copy< 6 >( distortion_dev, distortion );
  LvArray::tensorOps::subtract< 6 >( distortion_dev, distortion_iso );

  // Remove spherical tensile strength smoothly as damage grows, while leaving
  // compressive pressure untouched.  This gives a fully damaged powder no
  // cohesive opening resistance but still allows pressure-dependent granular
  // shear strength under confinement.
  if( distortionMeanStress > 0.0 && tensileStressMultiplier < 1.0 )
  {
    LvArray::tensorOps::scale< 6 >( distortion_iso, tensileStressMultiplier );
    plastic = true;
  }

  // Define pressure-dependent mode strengths.  A basal fracture degrades the
  // plane-normal distortion and coupled weak-plane shear modes.  A comminution
  // damage state degrades all modes toward residual powder strength.  In the
  // fully comminuted limit each mode has zero cohesion, uses
  // damagedMaterialFrictionalSlope, and retains its unscaled Y2 plateau.
  //
  // The distortion mode is signed.  I_d = sigma_aa - (sigma_bb+sigma_cc)/2
  // distinguishes basal-normal relative extension/weak opening-like distortion
  // from basal-normal relative compression.  The
  // distortionShearResponse* inputs define I_d^- = max(-I_d,0).  The optional
  // positiveDistortionShearResponse* inputs define I_d^+ = max(I_d,0).  When
  // both branches use the same parameters, the surface reduces to |I_d|/Y_d.
  real64 const distortionDamage = LvArray::math::max( basalPlaneDamage, comminutionDamage );
  real64 const coupledDamage = LvArray::math::max( basalPlaneDamage, comminutionDamage );
  real64 const inPlaneDamage = comminutionDamage;

  real64 const negativeDistortionStrength =
    evaluateGraphiteModeStrength( pressure, distortionDamage, strengthScale,
                                  m_distortionShearResponseX2,
                                  m_distortionShearResponseY1,
                                  m_distortionShearResponseY2,
                                  m_distortionShearResponseM1,
                                  distortionHardeningMultiplier );
  real64 const positiveDistortionStrength =
    evaluateGraphiteModeStrength( pressure, distortionDamage, strengthScale,
                                  m_positiveDistortionShearResponseX2,
                                  m_positiveDistortionShearResponseY1,
                                  m_positiveDistortionShearResponseY2,
                                  m_positiveDistortionShearResponseM1,
                                  distortionHardeningMultiplier );

  real64 const coupledHardeningMultiplier =
    1.0 + ( m_coupledStrainHardeningC0 - 1.0 ) * smoothStep( m_relaxation[k][q], 0.0, 1.0 );
  real64 const coupledYieldStrength =
    evaluateGraphiteModeStrength( pressure, coupledDamage, strengthScale,
                                  m_coupledShearResponseX2,
                                  m_coupledShearResponseY1,
                                  m_coupledShearResponseY2,
                                  m_coupledShearResponseM1,
                                  coupledHardeningMultiplier );

  real64 const inPlaneHardeningMultiplier =
    1.0 + ( m_inPlaneStrainHardeningC0 - 1.0 ) * smoothStep( m_relaxation[k][q], 0.0, 1.0 );
  real64 const inPlaneShearStrength =
    evaluateGraphiteModeStrength( pressure, inPlaneDamage, strengthScale,
                                  m_inPlaneShearResponseX2,
                                  m_inPlaneShearResponseY1,
                                  m_inPlaneShearResponseY2,
                                  m_inPlaneShearResponseM1,
                                  inPlaneHardeningMultiplier );

  // Enforce strength ---------------------------------------------

  // Compute the TI modal stress measures.  The in-plane and coupled shear
  // scaling constants are chosen so that those measures coincide with
  // von-Mises-equivalent stress for their deviatoric subspaces when all mode
  // strengths are equal.  The distortion mode is one-dimensional; its
  // von-Mises-equivalent measure is the signed scalar
  // I_d = sigma_aa - (sigma_bb+sigma_cc)/2.  Splitting I_d into positive and
  // negative parts gives tension/compression-like asymmetry for normal
  // distortion and reduces to an unsigned distortion measure when both branch
  // strengths are equal.
  real64 const currentPlaneNormalStress = sigma1[0] + sigma1[1] + sigma1[2];
  real64 const currentInPlaneMeanStress = 0.5 * ( inPlaneIso[0] + inPlaneIso[1] + inPlaneIso[2] );
  real64 const signedDistortionStress = currentPlaneNormalStress - currentInPlaneMeanStress;
  real64 const positiveDistortionStress = LvArray::math::max( signedDistortionStress, 0.0 );
  real64 const negativeDistortionStress = LvArray::math::max( -signedDistortionStress, 0.0 );
  inPlaneDev[3] *= 1.41421356237;
  inPlaneDev[4] *= 1.41421356237;
  inPlaneDev[5] *= 1.41421356237;
  real64 inPlaneShearStress = 1.224744871391589 * LvArray::tensorOps::l2Norm< 6 >( inPlaneDev );
  inPlaneDev[3] /= 1.41421356237;
  inPlaneDev[4] /= 1.41421356237;
  inPlaneDev[5] /= 1.41421356237;

  sigma5[3] *= 1.41421356237;
  sigma5[4] *= 1.41421356237;
  sigma5[5] *= 1.41421356237;
  real64 coupledShearStress = 1.224744871391589 * LvArray::tensorOps::l2Norm< 6 >( sigma5 );
  sigma5[3] /= 1.41421356237;
  sigma5[4] /= 1.41421356237;
  sigma5[5] /= 1.41421356237;

  real64 const equivalentStrengthRatioSquared =
    ( positiveDistortionStress * positiveDistortionStress ) / ( positiveDistortionStrength * positiveDistortionStrength ) +
    ( negativeDistortionStress * negativeDistortionStress ) / ( negativeDistortionStrength * negativeDistortionStrength ) +
    ( inPlaneShearStress * inPlaneShearStress ) / ( inPlaneShearStrength * inPlaneShearStrength ) +
    ( coupledShearStress * coupledShearStress ) / ( coupledYieldStrength * coupledYieldStrength );

  if( equivalentStrengthRatioSquared > 1.0 )
  {
    real64 const radialReturnScale = 1.0 / LvArray::math::sqrt( equivalentStrengthRatioSquared );
    LvArray::tensorOps::scale< 6 >( distortion_dev, radialReturnScale );
    LvArray::tensorOps::scale< 6 >( inPlaneDev, radialReturnScale );
    LvArray::tensorOps::scale< 6 >( sigma5, radialReturnScale );
    plastic = true;
  }

  // reassemble end-of-step stress
  real64 newStress[6] = {0};
  LvArray::tensorOps::add< 6 >( newStress, distortion_iso );
  LvArray::tensorOps::add< 6 >( newStress, distortion_dev );
  LvArray::tensorOps::add< 6 >( newStress, inPlaneDev );
  LvArray::tensorOps::add< 6 >( newStress, sigma5 );

  // Copy the new stress
  LvArray::tensorOps::copy< 6 >( stress, newStress );

  //////////////////////////////////////////////////////
  // Evolve state variables.
  //
  // For inelastic response evolve damage if pressure is below the
  // brittle ductile transition pressure and evolve plastic strain
  // if pressure is above brittle ductile transition pressure.
  // plastic strain doesn't currently affect material response, but
  // it may be useful for plotting regions of high-pressure yield.
  real64 plasticStrainStateIncrement[6] = {};
  bool updatePlasticStrainState = false;
  if( distensionClosureOccurred )
  {
    LvArray::tensorOps::copy< 6 >( plasticStrainStateIncrement, distensionClosurePlasticStrainIncrement );
    LvArray::tensorOps::scale< 6 >( plasticStrainStateIncrement, -1.0 );
    updatePlasticStrainState = true;
  }

  // When directional distension is enabled, store the brittle tensile opening
  // strain in the accumulated inelastic/plastic-strain state as a closeable
  // contribution.  Subsequent distension closure subtracts the corresponding
  // strain from this state, while hardening, plastic work, and damage remain
  // driven only by the stress-producing plastic increment below.
  bool brittleStrainOccurred = false;
  for( int i=0; i<6; ++i )
  {
    brittleStrainOccurred = brittleStrainOccurred ||
                            ( LvArray::math::abs( brittleStrainIncrement[i] ) > 0.0 );
  }
  if( ( m_enableDistension == 1 ) && brittleStrainOccurred )
  {
    LvArray::tensorOps::add< 6 >( plasticStrainStateIncrement, brittleStrainIncrement );
    updatePlasticStrainState = true;
  }

  if( plastic )
  {
    // CC: debug
    // GEOS_LOG_RANK( "Particle " << k << ": PlaneNormalStress: " << planeNormalStress << ", " <<
    //                                      "SignedDistortionStress: " << signedDistortionStress << ", " <<
    //                                      "PositiveDistortionStrength: " << positiveDistortionStrength << ", " <<
    //                                      "NegativeDistortionStrength: " << negativeDistortionStrength << ", " <<
    //                                      "inPlaneShearStress: " << inPlaneShearStress << " (" << inPlaneShearStrength << ")" << ", " <<
    //                                      "coupledShearStress: " << coupledShearStress << " (" << coupledYieldStrength << ")");

    // Increment stress-producing plastic strain. This excludes stress-free
    // distension closure so hardening, plastic work, and damage remain driven
    // by irreversible inelastic deformation. The plotted plastic strain state
    // is updated separately below and includes the negative closure increment
    // so dilatant plastic volume can be recovered as gaps/pores close.
    real64 plasticStrainIncrement[6] = {0};
    computeTransverselyIsotropicPlasticStrainIncrement( unrotatedVelocityGradientForStress, // Velocity gradient tensor with stress-free distension closure removed
                                                        oldStress,                     // Stress at start of step
                                                        stress,                        // Stress at end of step
                                                        Ez,                            // Elastic modulus preferred direction
                                                        Ep,                            // Elastic modulus transverse plane
                                                        nuzp,                          // Poisson ratio coupled
                                                        nup,                           // Poisson ratio transverse plane
                                                        Gzp,                           // Shear modulus coupled plane
                                                        unrotatedMaterialDirection,            // material direction (unit vector, plane
                                                                                               // normal)
                                                        timeIncrement,                 // timeStep
                                                        plasticStrainIncrement );
    LvArray::tensorOps::subtract< 6 >( plasticStrainIncrement, brittleStrainIncrement );
    LvArray::tensorOps::add< 6 >( plasticStrainStateIncrement, plasticStrainIncrement );
    updatePlasticStrainState = true;

    // GEOS_LOG_RANK( "Particle: " << k << "\n " <<
    //               "\tOld Stress: {" << oldStress[0] << ", " << oldStress[1] << ", "<< oldStress[2] << ", "<< oldStress[3] << ", "<<
    // oldStress[4] << ", "<< oldStress[5] << "}\n " <<
    //               "\tNew Stress: {" << stress[0] << ", " << stress[1] << ", "<< stress[2] << ", "<< stress[3] << ", "<< stress[4] << ",
    // "<< stress[5] << "}\n " <<
    //               "\tStress Incr: {" << stress[0] - oldStress[0] << ", " << stress[1] - oldStress[1] << ", "<< stress[2] - oldStress[2]
    // << ", "<< stress[3] - oldStress[3] << ", "<< stress[4] - oldStress[4] << ", "<< stress[5] - oldStress[5] << "}\n " <<
    //               "\tOld plastic strain: {" << unrotatedOldPlasticStrain[0] << ", " << unrotatedOldPlasticStrain[1] << ", "<<
    // unrotatedOldPlasticStrain[2] << ", "<< unrotatedOldPlasticStrain[3] << ", "<< unrotatedOldPlasticStrain[4]<< ", "<<
    // unrotatedOldPlasticStrain[5] << "}"
    //               "\tPlastic strain increment: {" << plasticStrainIncrement[0] << ", " << plasticStrainIncrement[1] << ", " <<
    // plasticStrainIncrement[2] << ", " << plasticStrainIncrement[3] << ", "<< plasticStrainIncrement[4]<< ", " <<
    // plasticStrainIncrement[5] << "}" );

    // increment relaxation
    // For symmetric matrix need to double off diagonal elements for l2Norm?
    real64 plasticStrainIncrementForNorm[6] = { plasticStrainIncrement[0],
                                                plasticStrainIncrement[1],
                                                plasticStrainIncrement[2],
                                                plasticStrainIncrement[3] / 1.41421356237,
                                                plasticStrainIncrement[4] / 1.41421356237,
                                                plasticStrainIncrement[5] / 1.41421356237 };
    m_relaxation[k][q] += LvArray::tensorOps::l2Norm< 6 >( plasticStrainIncrementForNorm ) / m_maximumPlasticStrain;
    m_relaxation[k][q] = LvArray::math::min( 1.0, m_relaxation[k][q] );

    if( legacyPlasticRateDamage )
    {
      real64 const timeToFailure = m_lengthScale[k] / m_crackSpeed;
      if( basalModeIFracture == 1 )
      {
        real64 const rateDamage = LvArray::math::min( m_basalPlaneDamage[k][q] + timeIncrement / timeToFailure, 1.0 );
        m_basalPlaneDamage[k][q] = LvArray::math::max( m_basalPlaneDamage[k][q], rateDamage );
        m_damage[k][q] = LvArray::math::max( m_damage[k][q], rateDamage );
      }
      else
      {
        real64 const rateDamage = LvArray::math::min( m_comminutionDamage[k][q] + timeIncrement / timeToFailure, 1.0 );
        m_comminutionDamage[k][q] = LvArray::math::max( m_comminutionDamage[k][q], rateDamage );
        m_damage[k][q] = LvArray::math::max( m_damage[k][q], rateDamage );
      }
    }

    // If we have exceeded the tensile strength for normal stress, incement the plane-normal plastic
    // work that is used to compute damage.
    if ( m_basalPlaneFractureEnergyReleaseRate < DBL_MAX )
    {
      if( ( basalModeIFracture == 1 ) and ( pressure < m_distortionShearResponseX2 ) )
      {
      // Mode I basal plane fracture:
      real64 dEpsilonIJnJ[3] = {};
      real64 plasticStrainIncrementForTensorOp[6] = { plasticStrainIncrement[0],
                                                      plasticStrainIncrement[1],
                                                      plasticStrainIncrement[2],
                                                      plasticStrainIncrement[3] * 0.5,
                                                      plasticStrainIncrement[4] * 0.5,
                                                      plasticStrainIncrement[5] * 0.5 };
      LvArray::tensorOps::Ri_eq_symAijBj< 3 >(dEpsilonIJnJ, plasticStrainIncrementForTensorOp, unrotatedMaterialDirection);

      real64 plasticNormalStrainIncrement = LvArray::tensorOps::AiBi< 3 >(unrotatedMaterialDirection, dEpsilonIJnJ );
      real64 avePlaneNormalStress = 0.5*( planeNormalStress + oldPlaneNormalStress );
      m_basalPlanePlasticWork[k][q] += LvArray::math::max( 0.0 , plasticNormalStrainIncrement*avePlaneNormalStress );  // Should be non-negative, max may not be needed.
      }

      if ( pressure < m_distortionShearResponseX2 )
      { // Increment in basal plane shear work:
        for ( int i = 0; i < 6; ++i )
        {
          m_basalPlanePlasticWork[k][q] += plasticStrainIncrement[i]*sigma5[i];
        }
      }

      // Compare accumulated work to the effective basal fracture energy
      // normalized by the host-provided length scale.  The exponent controls how
      // sharply the work ratio approaches the DFG-visible failed state.
      real64 const normalizedWork = LvArray::math::max( 0.0, LvArray::math::min( 1.0, m_basalPlanePlasticWork[k][q] / ( basalPlaneFractureEnergyReleaseRate / m_lengthScale[k] ) ) );
      real64 const computedBasalPlaneDamage = LvArray::math::pow( normalizedWork, m_damageEvolutionExponent );
      real64 const newBasalPlaneDamage = limitDamageIncrementByCrackSpeed( k,
                                                                           timeIncrement,
                                                                           m_basalPlaneDamage[k][q],
                                                                           computedBasalPlaneDamage );
      m_basalPlaneDamage[k][q] = LvArray::math::max( m_basalPlaneDamage[k][q], newBasalPlaneDamage );
      m_damage[k][q] = LvArray::math::max( m_damage[k][q], m_basalPlaneDamage[k][q] );
    }

    // Compute total plastic work and compare to the effective value for the fracture energy release
    // rate (which should be defined to account for the fraction of inelastic dissipation that creates
    // new fracture surfaces).  Only evolve damage if pressure is below a brittle-ductile limit
    // defined as the X2 parameters for the distortion shear response.
    if( ( m_totalFractureEnergyReleaseRate < DBL_MAX ) and ( pressure < m_distortionShearResponseX2 ) )
    {
      // Increment in total plastic work.
      for ( int i = 0; i < 6; ++i )
      {
        // Here the stress is what's left after subtracting off the basal stress plane normal
        // and shear stress
        real64 dWp = plasticStrainIncrement[i]*(newStress[i] - sigma5[i] - sigma1[i]);
        m_plasticWork[k][q] += LvArray::math::max( 0.0, dWp );
      }

      // Comminution damage is driven by the non-basal plastic work budget and
      // blends all modal strengths toward the pressure-dependent residual powder
      // branch.  The public damage is still updated so DFG can introduce strong
      // discontinuities and separate powder clouds.
      real64 const normalizedWork = LvArray::math::max( 0.0, LvArray::math::min( 1.0, m_plasticWork[k][q] / ( totalFractureEnergyReleaseRate / m_lengthScale[k] ) ) );
      real64 const computedComminutionDamage = LvArray::math::pow( normalizedWork, m_damageEvolutionExponent );
      real64 const newComminutionDamage = limitDamageIncrementByCrackSpeed( k,
                                                                            timeIncrement,
                                                                            m_comminutionDamage[k][q],
                                                                            computedComminutionDamage );
      m_comminutionDamage[k][q] = LvArray::math::max( m_comminutionDamage[k][q], newComminutionDamage );
      m_damage[k][q] = LvArray::math::max( m_damage[k][q], m_comminutionDamage[k][q] );
    }

    if( m_enableDistension == 1 )
    {
      addDistensionFromInelasticStrain( k,
                                        q,
                                        unrotatedMaterialDirection,
                                        plasticStrainIncrement );
    }

    m_damage[k][q] = LvArray::math::max( m_damage[k][q],
                                         LvArray::math::max( m_basalPlaneDamage[k][q],
                                                             m_comminutionDamage[k][q] ) );
  }

  if( updatePlasticStrainState )
  {
    addUnrotatedPlasticStrainStateIncrement( k,
                                             q,
                                             rotationTranspose,
                                             endRotation,
                                             plasticStrainStateIncrement );
  }
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void GraphiteUpdates::computeTransverselyIsotropicTrialStress( const real64 timeIncrement,      // time step
                                                               const real64 Ez,                // Elastic modulus preferred direction
                                                               const real64 Ep,                // Elastic modulus transverse plane
                                                               const real64 nuzp,              // Poisson ratio coupled
                                                               const real64 nup,               // Poisson ratio transverse plane
                                                               const real64 Gzp,               // Shear modulus coupled plane
                                                               real64 const (&materialDirection)[3], // preferred direction
                                                               real64 const (&oldStress)[6],   // stress at start of step.
                                                               real64 const (&D)[6],           // D=sym(L)
                                                               real64 (& newStress) [6],        // stress at end of step
                                                               const localIndex k ) const
{
  // These are the TI elastic stiffness coefficients using Brannon's TI basis tensors:
  real64 h1 = ( Ez*Ez*(-1 + nup) ) / ( Ez*(-1 + nup) + 2*Ep*nuzp*nuzp );
  real64 h2 = -( ( Ep*(Ez*nup + Ep*nuzp*nuzp) ) / ( (1 + nup)*(Ez*(-1 + nup) + 2*Ep*nuzp*nuzp ) ) );
  real64 h3 = ( Ep*Ez*nuzp ) / ( Ez - Ez*nup - 2*Ep*nuzp*nuzp );
  real64 h4 = Ep/( 1 + nup );
  real64 h5 = 2*Gzp;

  // Construct the dense 3x3 transversely isotropic thermal expansion tensor.
  // We do this in a loop to easily implement the indicial expression, but
  // this could be more efficient since we only need 6 of the components.
  real64 alphaDense[3][3] = { };
  for( int i = 0; i < 3; ++i )
  {
    for( int j = 0; j < 3; ++j )
    {
      real64 delta = (i == j) ? 1.0 : 0.0;
      alphaDense[i][j] = (m_alphaL - m_alphaT) * materialDirection[i] * materialDirection[j] + delta * m_alphaT;
    }
  }

  // The trial stress increment is computed using the transversely isotropic basis construction of the stiffness
  // tensor, and we subtract a thermal strain rate from the symmetric portion of the velocity gradient too
  // allow for thermal expansion effects:
  real64 stressIncrementDense[3][3] = { };
  int voigtMap[3][3] = { {0, 5, 4}, {5, 1, 3}, {4, 3, 2} };
  for( int i=0; i<3; i++ )
  {
    for( int j=0; j<3; j++ )
    {
      for( int p=0; p<3; p++ )
      {
        for( int w=0; w<3; w++ )
        {
          // newStress[voigtMap[i][j]]
          stressIncrementDense[i][j] += ( h1*transverselyIsotropicB1( materialDirection, i, j, p, w ) +
                                          h2*transverselyIsotropicB2( materialDirection, i, j, p, w ) +
                                          h3*transverselyIsotropicB3( materialDirection, i, j, p, w ) +
                                          h4*transverselyIsotropicB4( materialDirection, i, j, p, w ) +
                                          h5*transverselyIsotropicB5( materialDirection, i, j, p, w ))*( D[voigtMap[p][w]] - alphaDense[p][w]*m_temperatureRate[k] )*timeIncrement;
        }
      }
    }
  }

  real64 stressIncrement[6] = {};
  LvArray::tensorOps::denseToSymmetric< 3 >( stressIncrement, stressIncrementDense );

  LvArray::tensorOps::copy< 6 >( newStress, oldStress );
  LvArray::tensorOps::add< 6 >( newStress, stressIncrement );

// CC: debug
// GEOS_LOG_RANK( "h constants: {" << h1 << ", " << h2 << ", " << h3 << ", " << h4 << ", " << h5 << "}\n" <<
//                "Mat dir: {" << materialDirection[0] << ", " << materialDirection[1] << ", " << materialDirection[2] << "}\n" <<
//                "D: {" << D[0] << ", " << D[1] << ", " << D[2] << ", " << D[3] << ", " << D[4] << ", " << D[5] << "}\n" <<
//                "Ddt: {" << D[0]*timeIncrement << ", " << D[1]*timeIncrement << ", " << D[2]*timeIncrement << ", " << D[3]*timeIncrement
// << ", " << D[4]*timeIncrement << ", " << D[5]*timeIncrement << "}\n" <<
//                "Old stress: {" << oldStress[0] << ", " << oldStress[1] << ", " << oldStress[2] << ", " << oldStress[3] << ", " <<
// oldStress[4] << ", " << oldStress[5] << "}\n" <<
//                "Stress Incr: {" << stressIncrement[0] << ", " << stressIncrement[1] << ", " << stressIncrement[2] << ", " <<
// stressIncrement[3] << ", " << stressIncrement[4] << ", " << stressIncrement[5] << "}\n" <<
//                "New stress: {" << newStress[0] << ", " << newStress[1] << ", " << newStress[2] << ", " << newStress[3] << ", " <<
// newStress[4] << ", " << newStress[5] << "}\n");
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void GraphiteUpdates::computeTransverselyIsotropicPlasticStrainIncrement( real64 const ( &velocityGradient )[3][3],           // Velocity
                                                                                                                              // gradient
                                                                                                                              // tensor
                                                                          real64 const ( &oldStress )[6],  // Stress at start of step
                                                                          real64 const ( &newStress )[6],  // Stress at end of step
                                                                          const real64 Ez,             // Elastic modulus preferred
                                                                                                       // direction
                                                                          const real64 Ep,             // Elastic modulus transverse plane
                                                                          const real64 nuzp,           // Poisson ratio coupled
                                                                          const real64 nup,            // Poisson ratio transverse plane
                                                                          const real64 Gzp,            // Shear modulus coupled plane
                                                                          real64 const ( &materialDirection )[3],       // material
                                                                                                                        // direction (unit
                                                                                                                        // vector)
                                                                          const real64 timeIncrement,
                                                                          real64 ( & plasticStrainIncrement )[6] ) const
{  // For hypo-elastic transversely isotropic models we compute the increment in plastic strain assuming
   // for some increment in total strain and stress and elastic properties.

  // New stress minus old stress
  real64 stressIncrement[6] = {};
  LvArray::tensorOps::copy< 6 >( stressIncrement, newStress );
  LvArray::tensorOps::subtract< 6 >( stressIncrement, oldStress );

  // These are the TI elastic compliance coefficients using Brannon's TI basis tensors:
  real64 s1 = 1.0/Ez;
  real64 s2 = -nup/Ep;
  real64 s3 = -nuzp/Ez;
  real64 s4 = (1+nup)/Ep;
  real64 s5 = 1.0/(2.0*Gzp);

//Could be issue with using symmetric voigt vector for strain here, try using dense tensor
  // real64 elasticStrainIncrement[6] = {};
  real64 elasticStrainIncrementDense[3][3] = {};
  int voigtMap[3][3] = { {0, 5, 4}, {5, 1, 3}, {4, 3, 2} };
  for( int i=0; i<3; i++ )
  {
    for( int j=0; j<3; j++ )
    {
      for( int p=0; p<3; p++ )
      {
        for( int w=0; w<3; w++ )
        {
          // elasticStrainIncrement[voigtMap[i][j]]
          elasticStrainIncrementDense[i][j] += ( s1 * transverselyIsotropicB1( materialDirection, i, j, p, w ) +
                                                 s2 * transverselyIsotropicB2( materialDirection, i, j, p, w ) +
                                                 s3 * transverselyIsotropicB3( materialDirection, i, j, p, w ) +
                                                 s4 * transverselyIsotropicB4( materialDirection, i, j, p, w ) +
                                                 s5 * transverselyIsotropicB5( materialDirection, i, j, p, w ) ) * stressIncrement[voigtMap[p][w]];
        }
      }
    }
  }

  // GEOS_LOG_RANK( "Elastic Strain Incr: {" << elasticStrainIncrement[0] << ", " << elasticStrainIncrement[1] << ", " <<
  // elasticStrainIncrement[2] << ", " << elasticStrainIncrement[3] << ", " << elasticStrainIncrement[4] << ", " <<
  // elasticStrainIncrement[5] << "}" );

  // plastic strain increment = Total strain increment - elastic strain increment
  real64 plasticStrainIncrementDense[3][3] = {};
  for( int i = 0; i < 3; ++i )
  {
    for( int j = 0; j < 3; ++j )
    {
      //Old GEOS code did not index plasticStrainIncrement, how did this not throw an error
      //does the voigt notation include the factore of 2 for off axis entries of the tensor?
      // plasticStrainIncrement[voigtMap[i][j]] += 0.5 * ( velocityGradient[i][j] + velocityGradient[j][i] ) * timeIncrement -
      // elasticStrainIncrement[voigtMap[i][j]];
      plasticStrainIncrementDense[i][j] += 0.5 * ( velocityGradient[i][j] + velocityGradient[j][i] ) * timeIncrement - elasticStrainIncrementDense[i][j];
    }
  }

  LvArray::tensorOps::denseToSymmetric< 3 >( plasticStrainIncrement, plasticStrainIncrementDense );
  // Strain off diagonal elements stored in voigt notation (2x symmetric elements)
  plasticStrainIncrement[3] *= 2.0;
  plasticStrainIncrement[4] *= 2.0;
  plasticStrainIncrement[5] *= 2.0;
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void GraphiteUpdates::computeTransverselyIsotropicElasticStrainIncrementFromStressIncrement( real64 const ( &stressIncrement )[6],
                                                                                              const real64 Ez,
                                                                                              const real64 Ep,
                                                                                              const real64 nuzp,
                                                                                              const real64 nup,
                                                                                              const real64 Gzp,
                                                                                              real64 const ( &materialDirection )[3],
                                                                                              real64 ( &elasticStrainIncrement )[6] ) const
{
  real64 const s1 = 1.0/Ez;
  real64 const s2 = -nup/Ep;
  real64 const s3 = -nuzp/Ez;
  real64 const s4 = (1.0+nup)/Ep;
  real64 const s5 = 1.0/(2.0*Gzp);

  real64 elasticStrainIncrementDense[3][3] = {};
  int const voigtMap[3][3] = { {0, 5, 4}, {5, 1, 3}, {4, 3, 2} };
  for( int i=0; i<3; ++i )
  {
    for( int j=0; j<3; ++j )
    {
      for( int p=0; p<3; ++p )
      {
        for( int w=0; w<3; ++w )
        {
          elasticStrainIncrementDense[i][j] +=
            ( s1 * transverselyIsotropicB1( materialDirection, i, j, p, w ) +
              s2 * transverselyIsotropicB2( materialDirection, i, j, p, w ) +
              s3 * transverselyIsotropicB3( materialDirection, i, j, p, w ) +
              s4 * transverselyIsotropicB4( materialDirection, i, j, p, w ) +
              s5 * transverselyIsotropicB5( materialDirection, i, j, p, w ) ) * stressIncrement[voigtMap[p][w]];
        }
      }
    }
  }

  LvArray::tensorOps::denseToSymmetric< 3 >( elasticStrainIncrement, elasticStrainIncrementDense );
  elasticStrainIncrement[3] *= 2.0;
  elasticStrainIncrement[4] *= 2.0;
  elasticStrainIncrement[5] *= 2.0;
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 GraphiteUpdates::computeTransverselyIsotropicNormalCompliance( const real64 Ez,
                                                                       const real64 Ep,
                                                                       const real64 nuzp,
                                                                       const real64 nup,
                                                                       const real64 Gzp,
                                                                       real64 const ( &materialDirection )[3],
                                                                       real64 const ( &normal )[3] ) const
{
  real64 const s1 = 1.0/Ez;
  real64 const s2 = -nup/Ep;
  real64 const s3 = -nuzp/Ez;
  real64 const s4 = (1.0+nup)/Ep;
  real64 const s5 = 1.0/(2.0*Gzp);

  real64 compliance = 0.0;
  for( int i=0; i<3; ++i )
  {
    for( int j=0; j<3; ++j )
    {
      real64 const Nij = normal[i]*normal[j];
      for( int p=0; p<3; ++p )
      {
        for( int w=0; w<3; ++w )
        {
          real64 const Npw = normal[p]*normal[w];
          compliance += Nij *
            ( s1 * transverselyIsotropicB1( materialDirection, i, j, p, w ) +
              s2 * transverselyIsotropicB2( materialDirection, i, j, p, w ) +
              s3 * transverselyIsotropicB3( materialDirection, i, j, p, w ) +
              s4 * transverselyIsotropicB4( materialDirection, i, j, p, w ) +
              s5 * transverselyIsotropicB5( materialDirection, i, j, p, w ) ) * Npw;
        }
      }
    }
  }

  return LvArray::math::max( compliance, 1.0e-30 );
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 GraphiteUpdates::computeUniaxialPlaneNormalTensileStrength( const real64 strengthScale,
                                                                    const real64 distortionHardeningMultiplier ) const
{
  real64 const safeStrengthScale = LvArray::math::max( strengthScale, 1.0e-12 );
  real64 const safeHardeningMultiplier = LvArray::math::max( distortionHardeningMultiplier, 0.0 );
  real64 const safeY2 = LvArray::math::max( m_positiveDistortionShearResponseY2, 1.0e-20 );
  real64 const safeX2 = LvArray::math::max( m_positiveDistortionShearResponseX2, 1.0e-20 );

  real64 const maxY1 = safeY2 * ( 1.0 - 1.0e-12 );
  real64 const intactY1 = LvArray::math::max( 0.0,
                                              LvArray::math::min( safeStrengthScale * safeHardeningMultiplier *
                                                                  m_positiveDistortionShearResponseY1,
                                                                  maxY1 ) );
  real64 const intactSecantSlope = ( safeY2 - intactY1 ) / safeX2;
  real64 const intactM1 = LvArray::math::max( safeStrengthScale * safeHardeningMultiplier *
                                              m_positiveDistortionShearResponseM1,
                                              intactSecantSlope );

  // For uniaxial basal-normal tension, sigma = T a x a, the positive signed-distortion
  // strength branch is evaluated at p = -T/3 and the tensile intersection satisfies
  // T = Y_1^eff - M_1^eff T/3 on the tensile side of the pressure-dependent envelope.
  return LvArray::math::max( 1.0e-20, intactY1 / ( 1.0 + intactM1 / 3.0 ) );
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 GraphiteUpdates::symmetricStressAlongNormal( real64 const ( &stress )[6],
                                                    real64 const ( &normal )[3] ) const
{
  real64 stressNormal[3] = {};
  LvArray::tensorOps::Ri_eq_symAijBj< 3 >( stressNormal, stress, normal );
  return LvArray::tensorOps::AiBi< 3 >( normal, stressNormal );
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 GraphiteUpdates::computeMaximumPrincipalStressAndDirection( real64 const ( &stress )[6],
                                                                    real64 ( &normal )[3] ) const
{
  real64 principalStresses[3] = { 0.0, 0.0, 0.0 };
  real64 eigenVectors[3][3] = { { 0.0 } };
  LvArray::tensorOps::symEigenvectors< 3 >( principalStresses, eigenVectors, stress );

  int maximumIndex = 0;
  for( int i=1; i<3; ++i )
  {
    if( principalStresses[i] > principalStresses[maximumIndex] )
    {
      maximumIndex = i;
    }
  }

  normal[0] = eigenVectors[0][maximumIndex];
  normal[1] = eigenVectors[1][maximumIndex];
  normal[2] = eigenVectors[2][maximumIndex];
  LvArray::tensorOps::normalize< 3 >( normal );

  return principalStresses[maximumIndex];
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void GraphiteUpdates::subtractNormalStressComponent( const real64 stressDrop,
                                                     real64 const ( &normal )[3],
                                                     real64 ( &stress )[6] ) const
{
  if( stressDrop <= 0.0 )
  {
    return;
  }

  real64 correctionDense[3][3] = {};
  for( int i=0; i<3; ++i )
  {
    for( int j=0; j<3; ++j )
    {
      correctionDense[i][j] = stressDrop * normal[i] * normal[j];
    }
  }

  real64 correction[6] = {};
  LvArray::tensorOps::denseToSymmetric< 3 >( correction, correctionDense );
  LvArray::tensorOps::subtract< 6 >( stress, correction );
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 GraphiteUpdates::limitDamageIncrementByCrackSpeed( localIndex const k,
                                                          const real64 timeIncrement,
                                                          const real64 oldDamage,
                                                          const real64 computedDamage ) const
{
  real64 const boundedOldDamage = LvArray::math::max( 0.0, LvArray::math::min( 1.0, oldDamage ) );
  real64 const boundedComputedDamage = LvArray::math::max( boundedOldDamage,
                                                           LvArray::math::min( 1.0, computedDamage ) );

  if( m_crackSpeed >= DBL_MAX )
  {
    return boundedComputedDamage;
  }

  real64 const rateLimitedDamage =
    LvArray::math::min( 1.0,
                        boundedOldDamage +
                        m_crackSpeed * timeIncrement / LvArray::math::max( m_lengthScale[k], 1.0e-20 ) );

  return LvArray::math::max( boundedOldDamage,
                             LvArray::math::min( boundedComputedDamage, rateLimitedDamage ) );
}


GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 GraphiteUpdates::computeTransverseSpectralStrainPart( real64 const ( &strainDense )[3][3],
                                                             real64 const ( &materialDirection )[3],
                                                             const bool compressivePart,
                                                             real64 ( &strainPartDense )[3][3] ) const
{
  for( int i=0; i<3; ++i )
  {
    for( int j=0; j<3; ++j )
    {
      strainPartDense[i][j] = 0.0;
    }
  }

  real64 transverseProjector[3][3] = {};
  for( int i=0; i<3; ++i )
  {
    for( int j=0; j<3; ++j )
    {
      transverseProjector[i][j] = delta( i, j ) - materialDirection[i]*materialDirection[j];
    }
  }

  real64 transverseStrainDense[3][3] = {};
  for( int i=0; i<3; ++i )
  {
    for( int j=0; j<3; ++j )
    {
      for( int p=0; p<3; ++p )
      {
        for( int w=0; w<3; ++w )
        {
          transverseStrainDense[i][j] += transverseProjector[i][p] * strainDense[p][w] * transverseProjector[w][j];
        }
      }
    }
  }

  real64 transverseTrace = 0.0;
  real64 transverseFrobeniusSquared = 0.0;
  for( int i=0; i<3; ++i )
  {
    transverseTrace += transverseStrainDense[i][i];
    for( int j=0; j<3; ++j )
    {
      transverseFrobeniusSquared += transverseStrainDense[i][j] * transverseStrainDense[i][j];
    }
  }

  real64 const transverseDiscriminant = LvArray::math::max( 0.0,
                                                            2.0*transverseFrobeniusSquared -
                                                            transverseTrace*transverseTrace );
  real64 const transverseSpread = LvArray::math::sqrt( transverseDiscriminant );
  real64 const transverseMinimum = 0.5 * ( transverseTrace - transverseSpread );
  real64 const transverseMaximum = 0.5 * ( transverseTrace + transverseSpread );
  real64 const spectralTolerance = 1.0e-14 * ( LvArray::math::abs( transverseTrace ) + transverseSpread + 1.0 );

  if( compressivePart )
  {
    if( transverseMinimum >= -spectralTolerance )
    {
      return 0.0;
    }
  }
  else if( transverseMaximum <= spectralTolerance )
  {
    return 0.0;
  }

  real64 transverseStrain[6] = {};
  LvArray::tensorOps::denseToSymmetric< 3 >( transverseStrain, transverseStrainDense );

  real64 principalStrains[3] = { 0.0, 0.0, 0.0 };
  real64 eigenVectors[3][3] = { { 0.0 } };
  LvArray::tensorOps::symEigenvectors< 3 >( principalStrains, eigenVectors, transverseStrain );

  real64 tracePart = 0.0;
  for( int a=0; a<3; ++a )
  {
    real64 const spectralMagnitude = compressivePart ?
      LvArray::math::max( 0.0, -principalStrains[a] ) :
      LvArray::math::max( 0.0, principalStrains[a] );

    if( spectralMagnitude > 0.0 )
    {
      real64 projectedDirection[3] = { eigenVectors[0][a], eigenVectors[1][a], eigenVectors[2][a] };
      real64 const normalComponent = LvArray::tensorOps::AiBi< 3 >( projectedDirection, materialDirection );
      for( int i=0; i<3; ++i )
      {
        projectedDirection[i] -= normalComponent * materialDirection[i];
      }

      real64 const projectedNorm = LvArray::math::sqrt( LvArray::math::max( 0.0,
        projectedDirection[0]*projectedDirection[0] +
        projectedDirection[1]*projectedDirection[1] +
        projectedDirection[2]*projectedDirection[2] ) );

      if( projectedNorm > 1.0e-10 )
      {
        real64 const inverseProjectedNorm = 1.0 / projectedNorm;
        for( int i=0; i<3; ++i )
        {
          projectedDirection[i] *= inverseProjectedNorm;
        }

        for( int i=0; i<3; ++i )
        {
          for( int j=0; j<3; ++j )
          {
            strainPartDense[i][j] += spectralMagnitude * projectedDirection[i] * projectedDirection[j];
          }
        }
        tracePart += spectralMagnitude;
      }
    }
  }

  return tracePart;
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
bool GraphiteUpdates::applyTransverselyIsotropicDistensionClosure( localIndex const k,
                                                                   localIndex const q,
                                                                   const real64 timeIncrement,
                                                                   real64 const ( &materialDirection )[3],
                                                                   real64 ( &DForStress )[6],
                                                                   real64 ( &velocityGradientForStress )[3][3],
                                                                   real64 ( &closurePlasticStrainIncrement )[6] ) const
{
  if( timeIncrement <= 0.0 )
  {
    return false;
  }

  real64 const normalReservoir = LvArray::math::max( 0.0, m_basalNormalDistension[k][q] );
  real64 const transverseReservoir = LvArray::math::max( 0.0, m_transverseDistension[k][q] );
  if( ( normalReservoir <= 0.0 ) && ( transverseReservoir <= 0.0 ) )
  {
    return false;
  }

  int const voigtMap[3][3] = { {0, 5, 4}, {5, 1, 3}, {4, 3, 2} };
  real64 strainIncrementDense[3][3] = {};
  for( int i=0; i<3; ++i )
  {
    for( int j=0; j<3; ++j )
    {
      strainIncrementDense[i][j] = DForStress[voigtMap[i][j]] * timeIncrement;
    }
  }

  real64 closureStrainDense[3][3] = {};
  bool closureOccurred = false;

  real64 normalStrainIncrement = 0.0;
  for( int i=0; i<3; ++i )
  {
    for( int j=0; j<3; ++j )
    {
      normalStrainIncrement += materialDirection[i] * strainIncrementDense[i][j] * materialDirection[j];
    }
  }

  if( ( normalStrainIncrement < 0.0 ) && ( normalReservoir > 0.0 ) )
  {
    real64 const normalClosure = LvArray::math::min( normalReservoir, -normalStrainIncrement );
    closureOccurred = closureOccurred || ( normalClosure > 0.0 );
    m_basalNormalDistension[k][q] = LvArray::math::max( 0.0, m_basalNormalDistension[k][q] - normalClosure );
    for( int i=0; i<3; ++i )
    {
      for( int j=0; j<3; ++j )
      {
        closureStrainDense[i][j] += normalClosure * materialDirection[i] * materialDirection[j];
      }
    }
  }

  if( transverseReservoir > 0.0 )
  {
    real64 transverseCompressionDense[3][3] = {};
    real64 const transverseCompression = computeTransverseSpectralStrainPart( strainIncrementDense,
                                                                              materialDirection,
                                                                              true,
                                                                              transverseCompressionDense );

    if( transverseCompression > 0.0 )
    {
      real64 const transverseClosure = LvArray::math::min( transverseReservoir, transverseCompression );
      closureOccurred = closureOccurred || ( transverseClosure > 0.0 );
      real64 const transverseScale = transverseClosure / LvArray::math::max( transverseCompression, 1.0e-20 );
      m_transverseDistension[k][q] = LvArray::math::max( 0.0, m_transverseDistension[k][q] - transverseClosure );
      for( int i=0; i<3; ++i )
      {
        for( int j=0; j<3; ++j )
        {
          closureStrainDense[i][j] += transverseScale * transverseCompressionDense[i][j];
        }
      }
    }
  }

  if( !closureOccurred )
  {
    return false;
  }

  real64 closureStrain[6] = {};
  LvArray::tensorOps::denseToSymmetric< 3 >( closureStrain, closureStrainDense );
  real64 const inverseTimeIncrement = 1.0 / timeIncrement;
  for( int i=0; i<6; ++i )
  {
    DForStress[i] += closureStrain[i] * inverseTimeIncrement;
  }

  closurePlasticStrainIncrement[0] += closureStrain[0];
  closurePlasticStrainIncrement[1] += closureStrain[1];
  closurePlasticStrainIncrement[2] += closureStrain[2];
  closurePlasticStrainIncrement[3] += 2.0 * closureStrain[3];
  closurePlasticStrainIncrement[4] += 2.0 * closureStrain[4];
  closurePlasticStrainIncrement[5] += 2.0 * closureStrain[5];

  for( int i=0; i<3; ++i )
  {
    for( int j=0; j<3; ++j )
    {
      velocityGradientForStress[i][j] += closureStrainDense[i][j] * inverseTimeIncrement;
    }
  }

  updateDistensionPorosity( k, q );
  return true;
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void GraphiteUpdates::addUnrotatedPlasticStrainStateIncrement( localIndex const k,
                                                               localIndex const q,
                                                               real64 const ( &rotationTranspose )[3][3],
                                                               real64 const ( &endRotation )[3][3],
                                                               real64 const ( &plasticStrainStateIncrement )[6] ) const
{
  real64 oldPlasticStrain[6] = {};
  LvArray::tensorOps::copy< 6 >( oldPlasticStrain, m_plasticStrain[k][q] );
  oldPlasticStrain[3] *= 0.5;
  oldPlasticStrain[4] *= 0.5;
  oldPlasticStrain[5] *= 0.5;

  real64 unrotatedOldPlasticStrain[6] = {};
  LvArray::tensorOps::Rij_eq_AikSymBklAjl< 3 >( unrotatedOldPlasticStrain, rotationTranspose, oldPlasticStrain );
  unrotatedOldPlasticStrain[3] *= 2.0;
  unrotatedOldPlasticStrain[4] *= 2.0;
  unrotatedOldPlasticStrain[5] *= 2.0;

  real64 unrotatedNewPlasticStrain[6] = {};
  LvArray::tensorOps::copy< 6 >( unrotatedNewPlasticStrain, unrotatedOldPlasticStrain );
  LvArray::tensorOps::add< 6 >( unrotatedNewPlasticStrain, plasticStrainStateIncrement );

  unrotatedNewPlasticStrain[3] *= 0.5;
  unrotatedNewPlasticStrain[4] *= 0.5;
  unrotatedNewPlasticStrain[5] *= 0.5;

  real64 newPlasticStrain[6] = {};
  LvArray::tensorOps::Rij_eq_AikSymBklAjl< 3 >( newPlasticStrain, endRotation, unrotatedNewPlasticStrain );
  newPlasticStrain[3] *= 2.0;
  newPlasticStrain[4] *= 2.0;
  newPlasticStrain[5] *= 2.0;

  LvArray::tensorOps::copy< 6 >( m_plasticStrain[k][q], newPlasticStrain );
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void GraphiteUpdates::updateDistensionPorosity( localIndex const k,
                                                localIndex const q ) const
{
  real64 const combinedDistension =
    LvArray::math::max( 0.0, m_basalNormalDistension[k][q] ) +
    LvArray::math::max( 0.0, m_transverseDistension[k][q] );
  m_porosity[k][q] = LvArray::math::max( 0.0,
                                         LvArray::math::min( 1.0,
                                                             1.0 - LvArray::math::exp( -combinedDistension ) ) );
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void GraphiteUpdates::addDistensionFromTensileOpening( localIndex const k,
                                                       localIndex const q,
                                                       const real64 openingStrain,
                                                       real64 const ( &openingNormal )[3],
                                                       real64 const ( &materialDirection )[3],
                                                       const bool basalOpening ) const
{
  if( ( m_enableDistension != 1 ) || ( openingStrain <= 0.0 ) )
  {
    return;
  }

  real64 const normalProjection = LvArray::tensorOps::AiBi< 3 >( openingNormal, materialDirection );
  real64 const basalWeight = basalOpening ? 1.0 :
    LvArray::math::max( 0.0, LvArray::math::min( 1.0, normalProjection * normalProjection ) );
  real64 const transverseWeight = basalOpening ? 0.0 : LvArray::math::max( 0.0, 1.0 - basalWeight );

  real64 const normalDistensionIncrement = basalWeight * openingStrain;
  real64 const transverseDistensionIncrement = transverseWeight * openingStrain;
  if( ( normalDistensionIncrement <= 0.0 ) && ( transverseDistensionIncrement <= 0.0 ) )
  {
    return;
  }

  m_basalNormalDistension[k][q] += normalDistensionIncrement;
  m_transverseDistension[k][q] += transverseDistensionIncrement;
  updateDistensionPorosity( k, q );
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void GraphiteUpdates::addDistensionFromInelasticStrain( localIndex const k,
                                                        localIndex const q,
                                                        real64 const ( &materialDirection )[3],
                                                        real64 const ( &inelasticStrainIncrement )[6] ) const
{
  real64 const boundedBasalDamage = LvArray::math::max( 0.0, LvArray::math::min( 1.0, m_basalPlaneDamage[k][q] ) );
  real64 const boundedComminutionDamage = LvArray::math::max( 0.0, LvArray::math::min( 1.0, m_comminutionDamage[k][q] ) );
  real64 const normalWeight = smoothStep( LvArray::math::max( boundedBasalDamage, boundedComminutionDamage ), 0.0, 1.0 );
  real64 const transverseWeight = smoothStep( boundedComminutionDamage, 0.0, 1.0 );

  if( ( normalWeight <= 0.0 ) && ( transverseWeight <= 0.0 ) )
  {
    return;
  }

  real64 strainDense[3][3] = { { inelasticStrainIncrement[0], 0.5*inelasticStrainIncrement[5], 0.5*inelasticStrainIncrement[4] },
                               { 0.5*inelasticStrainIncrement[5], inelasticStrainIncrement[1], 0.5*inelasticStrainIncrement[3] },
                               { 0.5*inelasticStrainIncrement[4], 0.5*inelasticStrainIncrement[3], inelasticStrainIncrement[2] } };

  bool distensionChanged = false;
  if( normalWeight > 0.0 )
  {
    real64 normalOpening = 0.0;
    for( int i=0; i<3; ++i )
    {
      for( int j=0; j<3; ++j )
      {
        normalOpening += materialDirection[i] * strainDense[i][j] * materialDirection[j];
      }
    }

    real64 const normalIncrement = normalWeight * LvArray::math::max( 0.0, normalOpening );
    if( normalIncrement > 0.0 )
    {
      m_basalNormalDistension[k][q] += normalIncrement;
      distensionChanged = true;
    }
  }

  if( transverseWeight > 0.0 )
  {
    real64 transverseOpeningDense[3][3] = {};
    real64 const transverseOpening = computeTransverseSpectralStrainPart( strainDense,
                                                                          materialDirection,
                                                                          false,
                                                                          transverseOpeningDense );
    real64 const transverseIncrement = transverseWeight * LvArray::math::max( 0.0, transverseOpening );
    if( transverseIncrement > 0.0 )
    {
      m_transverseDistension[k][q] += transverseIncrement;
      distensionChanged = true;
    }
  }

  if( distensionChanged )
  {
    updateDistensionPorosity( k, q );
  }
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
bool GraphiteUpdates::applyEnergyRegularizedBrittleTensileReturn( localIndex const k,
                                                                   localIndex const q,
                                                                   const real64 timeIncrement,
                                                                   const real64 stressMeasure,
                                                                   const real64 oldDamage,
                                                                   const real64 effectiveStrength,
                                                                   const real64 fractureEnergyReleaseRate,
                                                                   const real64 normalCompliance,
                                                                   real64 const ( &normal )[3],
                                                                   real64 const ( &materialDirection )[3],
                                                                   const bool updateBasalPlaneDamage,
                                                                   real64 ( &stress )[6] ) const
{
  if( ( fractureEnergyReleaseRate >= DBL_MAX ) || ( stressMeasure <= 0.0 ) )
  {
    return false;
  }

  real64 const boundedOldDamage = LvArray::math::max( 0.0, LvArray::math::min( 1.0, oldDamage ) );
  real64 const safeStrength = LvArray::math::max( effectiveStrength, 1.0e-20 );
  real64 const safeCompliance = LvArray::math::max( normalCompliance, 1.0e-30 );
  real64 const oldCap = safeStrength * ( 1.0 - smoothStep( boundedOldDamage, 0.0, 1.0 ) );

  if( stressMeasure <= oldCap * ( 1.0 + 1.0e-12 ) + 1.0e-20 )
  {
    return false;
  }

  real64 const regularizedFractureEnergy = fractureEnergyReleaseRate / LvArray::math::max( m_lengthScale[k], 1.0e-20 );
  real64 const trialEnergy = 0.5 * stressMeasure * stressMeasure * safeCompliance;
  real64 const peakEnergy = 0.5 * safeStrength * safeStrength * safeCompliance;

  real64 const maximumStepDamage =
    limitDamageIncrementByCrackSpeed( k, timeIncrement, boundedOldDamage, 1.0 );

  real64 newDamage = boundedOldDamage;
  if( maximumStepDamage <= boundedOldDamage )
  {
    newDamage = boundedOldDamage;
  }
  else if( ( trialEnergy >= regularizedFractureEnergy ) || ( peakEnergy >= regularizedFractureEnergy ) )
  {
    newDamage = maximumStepDamage;
  }
  else
  {
    real64 const softeningStrain =
      2.0 * ( regularizedFractureEnergy - peakEnergy ) / safeStrength;

    real64 lowerDamage = boundedOldDamage;
    real64 upperDamage = maximumStepDamage;
    real64 const upperCap = safeStrength * ( 1.0 - smoothStep( upperDamage, 0.0, 1.0 ) );
    real64 const upperResidual =
      upperDamage - boundedOldDamage -
      safeCompliance * LvArray::math::max( 0.0, stressMeasure - upperCap ) /
      LvArray::math::max( softeningStrain, 1.0e-20 );

    if( upperResidual < 0.0 )
    {
      newDamage = upperDamage;
    }
    else
    {
      for( int iter=0; iter<32; ++iter )
      {
        real64 const middleDamage = 0.5 * ( lowerDamage + upperDamage );
        real64 const middleCap = safeStrength * ( 1.0 - smoothStep( middleDamage, 0.0, 1.0 ) );
        real64 const middleResidual =
          middleDamage - boundedOldDamage -
          safeCompliance * LvArray::math::max( 0.0, stressMeasure - middleCap ) /
          LvArray::math::max( softeningStrain, 1.0e-20 );

        if( middleResidual >= 0.0 )
        {
          upperDamage = middleDamage;
        }
        else
        {
          lowerDamage = middleDamage;
        }
      }
      newDamage = upperDamage;
    }
  }

  real64 const newCap = safeStrength * ( 1.0 - smoothStep( newDamage, 0.0, 1.0 ) );
  real64 const stressDrop = LvArray::math::max( 0.0, stressMeasure - newCap );
  subtractNormalStressComponent( stressDrop, normal, stress );

  addDistensionFromTensileOpening( k,
                                   q,
                                   safeCompliance * stressDrop,
                                   normal,
                                   materialDirection,
                                   updateBasalPlaneDamage );

  if( updateBasalPlaneDamage )
  {
    m_basalPlaneDamage[k][q] = LvArray::math::max( m_basalPlaneDamage[k][q], newDamage );
  }
  m_damage[k][q] = LvArray::math::max( m_damage[k][q], newDamage );

  return stressDrop > 0.0;
}


GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
bool GraphiteUpdates::applyEnergyRegularizedBrittleTensileReturns( localIndex const k,
                                                                    localIndex const q,
                                                                    const real64 timeIncrement,
                                                                    const real64 principalTensileStrength,
                                                                    const real64 basalNormalTensileStrength,
                                                                    const real64 basalPlaneFractureEnergyReleaseRate,
                                                                    const real64 totalFractureEnergyReleaseRate,
                                                                    const real64 Ez,
                                                                    const real64 Ep,
                                                                    const real64 nuzp,
                                                                    const real64 nup,
                                                                    const real64 Gzp,
                                                                    real64 const ( &materialDirection )[3],
                                                                    real64 ( &stress )[6] ) const
{
  bool anyReturn = false;
  real64 const crackTipFactor = LvArray::math::max( m_crackTipStressConcentration[k], 1.0 );
  bool const basalEnabled = basalPlaneFractureEnergyReleaseRate < DBL_MAX;
  bool const principalEnabled = ( m_maximumPrincipalStressDamage == 1 ) &&
                                ( totalFractureEnergyReleaseRate < DBL_MAX );

  for( int activeSetIteration=0; activeSetIteration<2; ++activeSetIteration )
  {
    real64 basalNormal[3] = { materialDirection[0], materialDirection[1], materialDirection[2] };
    real64 basalStress = 0.0;
    real64 basalOldDamage = 0.0;
    real64 basalEffectiveStrength = 1.0e-20;
    real64 basalViolationRatio = -1.0;

    if( basalEnabled )
    {
      basalStress = symmetricStressAlongNormal( stress, basalNormal );
      basalOldDamage = LvArray::math::max( LvArray::math::max( 0.0, LvArray::math::min( 1.0, m_basalPlaneDamage[k][q] ) ),
                                           LvArray::math::max( 0.0, LvArray::math::min( 1.0, m_comminutionDamage[k][q] ) ) );
      basalEffectiveStrength = LvArray::math::max( basalNormalTensileStrength / crackTipFactor, 1.0e-20 );
      real64 const basalCap = basalEffectiveStrength * ( 1.0 - smoothStep( basalOldDamage, 0.0, 1.0 ) );
      if( basalStress > basalCap * ( 1.0 + 1.0e-12 ) + 1.0e-20 )
      {
        basalViolationRatio = basalStress / LvArray::math::max( basalCap, 1.0e-20 );
      }
    }

    real64 principalNormal[3] = { 1.0, 0.0, 0.0 };
    real64 principalStress = 0.0;
    real64 principalOldDamage = 0.0;
    real64 principalEffectiveStrength = 1.0e-20;
    real64 principalViolationRatio = -1.0;

    if( principalEnabled )
    {
      principalStress = computeMaximumPrincipalStressAndDirection( stress, principalNormal );
      principalOldDamage = LvArray::math::max( LvArray::math::max( 0.0, LvArray::math::min( 1.0, m_damage[k][q] ) ),
                                               LvArray::math::max( LvArray::math::max( 0.0, LvArray::math::min( 1.0, m_basalPlaneDamage[k][q] ) ),
                                                                    LvArray::math::max( 0.0, LvArray::math::min( 1.0, m_comminutionDamage[k][q] ) ) ) );
      principalEffectiveStrength = LvArray::math::max( principalTensileStrength / crackTipFactor, 1.0e-20 );
      real64 const principalCap = principalEffectiveStrength * ( 1.0 - smoothStep( principalOldDamage, 0.0, 1.0 ) );
      if( principalStress > principalCap * ( 1.0 + 1.0e-12 ) + 1.0e-20 )
      {
        principalViolationRatio = principalStress / LvArray::math::max( principalCap, 1.0e-20 );
      }
    }

    bool returnedThisIteration = false;
    if( ( principalViolationRatio < 0.0 ) && ( basalViolationRatio < 0.0 ) )
    {
      break;
    }
    else if( principalViolationRatio >= basalViolationRatio )
    {
      real64 const principalCompliance =
        computeTransverselyIsotropicNormalCompliance( Ez, Ep, nuzp, nup, Gzp, materialDirection, principalNormal );
      returnedThisIteration = applyEnergyRegularizedBrittleTensileReturn( k,
                                                                          q,
                                                                          timeIncrement,
                                                                          principalStress,
                                                                          principalOldDamage,
                                                                          principalEffectiveStrength,
                                                                          totalFractureEnergyReleaseRate,
                                                                          principalCompliance,
                                                                          principalNormal,
                                                                          materialDirection,
                                                                          false,
                                                                          stress );
    }
    else
    {
      real64 const basalCompliance = LvArray::math::max( 1.0 / Ez, 1.0e-30 );
      returnedThisIteration = applyEnergyRegularizedBrittleTensileReturn( k,
                                                                          q,
                                                                          timeIncrement,
                                                                          basalStress,
                                                                          basalOldDamage,
                                                                          basalEffectiveStrength,
                                                                          basalPlaneFractureEnergyReleaseRate,
                                                                          basalCompliance,
                                                                          basalNormal,
                                                                          materialDirection,
                                                                          true,
                                                                          stress );
    }

    anyReturn = anyReturn || returnedThisIteration;
    if( !returnedThisIteration )
    {
      break;
    }
  }

  return anyReturn;
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 GraphiteUpdates::transverselyIsotropicB1( real64 const (&materialDirection)[3],
                                                 const int i,
                                                 const int j,
                                                 const int p,
                                                 const int w ) const
{ // Return B1_ijkl component, of the transversely isotropic basis tensor B1
  // described in Brannon's rotation/tensor book.

  real64 B1 = materialDirection[i]*materialDirection[j]*materialDirection[p]*materialDirection[w];
  return B1;
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 GraphiteUpdates::transverselyIsotropicB2( real64 const (&materialDirection)[3],
                                                 const int i,
                                                 const int j,
                                                 const int p,
                                                 const int w ) const
{ // Return B2_ijkl component, of the transversely isotropic basis tensor B2
  // described in Brannon's rotation/tensor book.

  real64 B2 = delta( i, j )*delta( p, w ) -
              ( materialDirection[i]*materialDirection[j]*delta( p, w ) +
                delta( i, j )*materialDirection[p]*materialDirection[w] ) +
              materialDirection[i]*materialDirection[j]*materialDirection[p]*materialDirection[w];
  return B2;
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 GraphiteUpdates::transverselyIsotropicB3( real64 const (&materialDirection)[3],
                                                 const int i,
                                                 const int j,
                                                 const int p,
                                                 const int w ) const
{ // Return B3_ijkl component, of the transversely isotropic basis tensor B3
  // described in Brannon's rotation/tensor book.

  real64 B3 = materialDirection[i]*materialDirection[j]*delta( p, w ) +
              delta( i, j )*materialDirection[p]*materialDirection[w] -
              2.0*materialDirection[i]*materialDirection[j]*materialDirection[p]*materialDirection[w];
  return B3;
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 GraphiteUpdates::transverselyIsotropicB4( real64 const (&materialDirection)[3],
                                                 const int i,
                                                 const int j,
                                                 const int p,
                                                 const int w ) const
{ // Return B4_ijkl component, of the transversely isotropic basis tensor B4
  // described in Brannon's rotation/tensor book.

  real64 B4 = 0.5*( delta( i, p )*delta( j, w ) + delta( i, w )*delta( j, p ) ) -
              0.5*( delta( i, w )*materialDirection[j]*materialDirection[p] +
                    materialDirection[i]*delta( j, p )*materialDirection[w] +
                    materialDirection[i]*delta( j, w )*materialDirection[p] +
                    delta( i, p )*materialDirection[j]*materialDirection[w] ) +
              materialDirection[i]*materialDirection[j]*materialDirection[p]*materialDirection[w];
  return B4;
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 GraphiteUpdates::transverselyIsotropicB5( real64 const (&materialDirection)[3],
                                                 const int i,
                                                 const int j,
                                                 const int p,
                                                 const int w ) const
{ // Return B5_ijkl component, of the transversely isotropic basis tensor B5
  // described in Brannon's rotation/tensor book.

  real64 B5 = 0.5*( delta( i, p )*materialDirection[j]*materialDirection[w] +
                    delta( i, w )*materialDirection[j]*materialDirection[p] +
                    materialDirection[i]*delta( j, w )*materialDirection[p] +
                    materialDirection[i]*delta( j, p )*materialDirection[w] ) -
              2.0*materialDirection[i]*materialDirection[p]*materialDirection[j]*materialDirection[w];
  return B5;
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 GraphiteUpdates::delta( const int i,
                               const int j ) const
{ // Dirac delta function that returns a real
  return i == j ? 1.0 : 0.0;
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 GraphiteUpdates::slopePoint0( const real64 x,
                                     const real64 x1,
                                     const real64 x2,
                                     const real64 y1,
                                     const real64 y2,
                                     const real64 m1 ) const
{
  // returns a function h(x), for x1 <= x <= x2, with h(x1)=y1, h(x2)=y2, h'(x1)=m1, x'(x2)=0
  //realT beta = (m1*x1 - m1*x2 - y1 + y2)/(y1 - y2);
  //return -(((m1*(x1 - x2) + m1*powf((x - x2)/(x1 - x2),beta)*(-x + x2) + (m1*(-x1 + x2)*y1)/(y1 - y2))*(y1 - y2))/(m1*(x1 - x2)));
  // This should be the same, but simplified:
  return y2 + (y1-y2)*LvArray::math::pow( (x-x2)/(x1-x2), m1*(x1-x2)/(y1-y2) );
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 GraphiteUpdates::evaluateSaturatingPressureArgument( const real64 pressureArgument,
                                                            const real64 pressureScale ) const
{
  real64 const positivePressureArgument = LvArray::math::max( 0.0, pressureArgument );

  // A very large scale is the default and recovers the original linear law.
  // The branch also avoids overflow in the product form below for DBL_MAX.
  if( pressureScale > 1.0e100 )
  {
    return positivePressureArgument;
  }

  real64 const safePressureScale = LvArray::math::max( pressureScale, 1.0e-20 );
  return positivePressureArgument * safePressureScale /
         ( positivePressureArgument + safePressureScale );
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 GraphiteUpdates::evaluatePressureDependentStrength( const real64 pressure,
                                                           const real64 x2,
                                                           const real64 y1,
                                                           const real64 y2,
                                                           const real64 m1 ) const
{
  // The graphite strength envelopes use a point-slope curve with x1=0.
  // This helper is used for both the intact crystal branch and the
  // comminuted-powder branch.  The high-pressure strength y2 is treated as
  // the fixed pressure-dominated limit.  The low-pressure intercept y1 and
  // initial slope m1 may be scaled before entering this helper; they are then
  // limited here so the envelope remains monotone and concave up to the
  // plateau.  The small slope margin makes the curve tangent to the plateau
  // rather than exactly linear to it.
  real64 const x1 = 0.0;
  real64 const safeDx = LvArray::math::max( x2 - x1, 1.0e-20 );
  real64 const safeY2 = LvArray::math::max( y2, 1.0e-20 );
  real64 const maxY1 = ( 1.0 - 1.0e-12 ) * safeY2;
  real64 const safeY1 = LvArray::math::max( 0.0, LvArray::math::min( y1, maxY1 ) );
  real64 const minimumConcaveSlope = ( 1.0 + 1.0e-12 ) * ( safeY2 - safeY1 ) / safeDx;
  real64 const safeSlope = LvArray::math::max( m1, minimumConcaveSlope );

  if( pressure < x1 )
  {
    return LvArray::math::max( 0.0, safeY1 - ( x1 - pressure ) * safeSlope );
  }
  else if( pressure < x2 )
  {
    return slopePoint0( pressure, x1, x2, safeY1, safeY2, safeSlope );
  }
  else
  {
    return safeY2;
  }
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 GraphiteUpdates::evaluateGraphiteModeStrength( const real64 pressure,
                                                       const real64 damage,
                                                       const real64 strengthScale,
                                                       const real64 x2,
                                                       const real64 y1,
                                                       const real64 y2,
                                                       const real64 m1,
                                                       const real64 strainHardeningMultiplier ) const
{
  // strengthScale is a particle-scale Weibull/aleatory multiplier for the
  // low-pressure, flaw-controlled part of the intact crystal strength.  The
  // high-pressure plateau Y2 is not scaled: at sufficiently high confinement
  // the effect of pre-existing microcracks and weak flaws is assumed to close
  // out, so the intact branch and the comminuted powder branch share the same
  // mode-specific Y2.  This prevents a weak particle with strengthScale < 1
  // from hardening merely because damage blends it into the powder residual.
  //
  // The low-pressure ordinate and tangent slope are scaled by strengthScale
  // and the relaxation hardening multiplier, but are then limited so the
  // point-slope envelope remains a concave, monotone cap in p-q space.  The
  // scaled Y1 cannot exceed the fixed Y2 plateau, and the scaled M1 is raised
  // to at least the secant slope (Y2-Y1)/X2 when required.
  //
  // The fully damaged branch represents a comminuted granular residual.  It
  // has zero cohesion at p=0, uses damagedMaterialFrictionalSlope as the
  // low-pressure frictional slope, and retains the same unscaled mode-specific
  // Y2 plateau.
  real64 const hDamage = smoothStep( LvArray::math::max( 0.0, LvArray::math::min( 1.0, damage ) ), 0.0, 1.0 );
  real64 const safeStrengthScale = LvArray::math::max( strengthScale, 1.0e-12 );
  real64 const safeHardeningMultiplier = LvArray::math::max( strainHardeningMultiplier, 0.0 );
  real64 const safeY2 = LvArray::math::max( y2, 1.0e-20 );
  real64 const safeX2 = LvArray::math::max( x2, 1.0e-20 );

  real64 const maxY1 = safeY2 * ( 1.0 - 1.0e-12 );
  real64 const intactY1 = LvArray::math::max( 0.0,
                                              LvArray::math::min( safeStrengthScale * safeHardeningMultiplier * y1, maxY1 ) );
  real64 const intactSecantSlope = ( safeY2 - intactY1 ) / safeX2;
  real64 const intactM1 = LvArray::math::max( safeStrengthScale * safeHardeningMultiplier * m1, intactSecantSlope );

  real64 const intactStrength = evaluatePressureDependentStrength( pressure,
                                                                   x2,
                                                                   intactY1,
                                                                   safeY2,
                                                                   intactM1 );

  // The powder branch is also limited from above by the intact branch at
  // the same pressure.  Without this upper bound a weak particle with
  // strengthScale < 1, or a large damagedMaterialFrictionalSlope, could gain
  // shear strength as damage blends the response toward the powder residual.
  // For the point-slope curve used here, Y_powder(p) <= Y_intact(p) on
  // 0 <= p <= X2 is guaranteed when the powder exponent is no larger than
  // the intact exponent.  The resulting upper bound on the powder tangent is
  // M_powder <= M_intact * Y2 / (Y2 - Y1_intact), while the lower bound
  // M_powder >= Y2 / X2 preserves the monotone/concave cap to the plateau.
  real64 const powderMinimumSlope = safeY2 / safeX2;
  real64 const powderMaximumSlope = intactM1 * safeY2 / LvArray::math::max( safeY2 - intactY1, 1.0e-20 );
  real64 const powderM1 = LvArray::math::min( LvArray::math::max( m_damagedMaterialFrictionalSlope, powderMinimumSlope ),
                                              LvArray::math::max( powderMinimumSlope, powderMaximumSlope ) );

  real64 const powderStrength = evaluatePressureDependentStrength( pressure,
                                                                   x2,
                                                                   0.0,
                                                                   safeY2,
                                                                   powderM1 );

  return LvArray::math::max( 1.0e-20, ( 1.0 - hDamage ) * intactStrength + hDamage * powderStrength );
}


/**
 * @class Graphite
 *
 * Graphite material model.
 */
class Graphite : public SolidBase
{
public:

  /// @typedef Alias for GraphiteUpdates
  using KernelWrapper = GraphiteUpdates;

  /**
   * constructor
   * @param[in] name name of the instance in the catalog
   * @param[in] parent the group which contains this instance
   */
  Graphite( string const & name, Group * const parent );

  /**
   * Default Destructor
   */
  virtual ~Graphite() override;


  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  virtual void saveConvergedState() const override;

  /**
   * @name Static Factory Catalog members and functions
   */
  ///@{

  /// string name to use for this class in the catalog
  static constexpr auto m_catalogNameString = "Graphite";

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
    /// string/key for transverse Young's modulus
    static constexpr char const * defaultYoungModulusTransverseString() { return "defaultYoungModulusTransverse"; }

    /// string/key for axial Young's modulus
    static constexpr char const * defaultYoungModulusAxialString() { return "defaultYoungModulusAxial"; }

    /// string/key for transverse Poisson's Ratio
    static constexpr char const * defaultPoissonRatioTransverseString() { return "defaultPoissonRatioTransverse"; }

    /// string/key for axial Poisson's Ratio
    static constexpr char const * defaultPoissonRatioAxialTransverseString() { return "defaultPoissonRatioAxialTransverse"; }

    /// string/key for transverse shear modulus
    static constexpr char const * defaultShearModulusAxialTransverseString() { return "defaultShearModulusAxialTransverse"; }

    /// string/key for default transverse Young's modulus presssure derivative
    static constexpr char const * defaultYoungModulusTransversePressureDerivativeString() { return "defaultYoungModulusTransversePressureDerivative"; }

    /// string/key for default axial Young's modulus presssure derivative
    static constexpr char const * defaultYoungModulusAxialPressureDerivativeString() { return "defaultYoungModulusAxialPressureDerivative"; }

    /// string/key for default axial-transverse shear modulus presssure derivative
    static constexpr char const * defaultShearModulusAxialTransversePressureDerivativeString() { return "defaultShearModulusAxialTransversePressureDerivative"; }

    /// string/key for transverse Young's modulus pressure scale
    static constexpr char const * defaultYoungModulusTransversePressureScaleString() { return "defaultYoungModulusTransversePressureScale"; }

    /// string/key for axial Young's modulus pressure scale
    static constexpr char const * defaultYoungModulusAxialPressureScaleString() { return "defaultYoungModulusAxialPressureScale"; }

    /// string/key for axial-transverse shear modulus pressure scale
    static constexpr char const * defaultShearModulusAxialTransversePressureScaleString() { return "defaultShearModulusAxialTransversePressureScale"; }

    //string/key for element/particle velocityGradient value
    static constexpr char const * velocityGradientString() { return "velocityGradient"; }

    /// string/key for quadrature point plasticStrain value
    static constexpr char const * plasticStrainString() { return "plasticStrain"; }

    /// string/key for quadrature point relaxation value
    static constexpr char const * relaxationString() { return "relaxation"; }

    /// string/key for Flag to enable stress concenration for crack-tip
    static constexpr char const * enableCrackTipStressConcentrationString() { return "enableCrackTipStressConcentration"; }
    /// string/key for crackTipStressConcentration plotting variable
    static constexpr char const * crackTipStressConcentrationString() { return "crackTipStressConcentration"; }
    /// string/key for strength scale value
    static constexpr char const * distanceToCrackTipString() { return "distanceToCrackTip"; }


    /// string/key for quadrature point accumulated plastic work value
    static constexpr char const * basalPlanePlasticWorkString() { return "basalPlanePlasticWork"; }
    static constexpr char const * plasticWorkString() { return "plasticWork"; }

    /// constant for thermal expansion lateral to symmetry axis
    static constexpr char const * alphaLString() { return "alphaL"; }

    /// constant for thermal expansion transverse to symmetry axis
    static constexpr char const * alphaTString() { return "alphaT"; }

    /// string/key for host-visible quadrature point damage value used by DFG.
    static constexpr char const * damageString() { return "damage"; }

    /// string/key for internal basal-plane damage state.
    static constexpr char const * basalPlaneDamageString() { return "basalPlaneDamage"; }

    /// string/key for internal comminution damage state.
    static constexpr char const * comminutionDamageString() { return "comminutionDamage"; }

    /// string/key for flag to enable stress-free closure of directional distension.
    static constexpr char const * enableDistensionString() { return "enableDistension"; }

    /// string/key for basal-normal closeable distension.
    static constexpr char const * basalNormalDistensionString() { return "basalNormalDistension"; }

    /// string/key for transverse closeable distension.
    static constexpr char const * transverseDistensionString() { return "transverseDistension"; }

    /// string/key for distension-derived porosity.
    static constexpr char const * porosityString() { return "porosity"; }

    /// string/key for quadrature point temperature value
    static constexpr char const * temperatureString() { return "temperature"; }

    /// string/key for quadrature point temperature value
    static constexpr char const * temperatureRateString() { return "temperatureRate"; }

    /// string/key for quadrature point jacobian value
    static constexpr char const * jacobianString() { return "jacobian"; }

    /// string/key for element/particle length scale
    static constexpr char const * lengthScaleString() { return "lengthScale"; }

    /// string/key for element/particle length scale
    static constexpr char const * strengthScaleString() { return "strengthScale"; }

    /// string/key for maximum strength
    static constexpr char const * failureStrengthString() { return "failureStrength"; }

    /// string/key for flag to apply max principal stress failure criterion
    static constexpr char const * maximumPrincipalStressDamageString() { return "maximumPrincipalStressDamage"; }

    /// string/key for damage evolution rate parameters
    static constexpr char const * crackSpeedString() { return "crackSpeed"; }

    /// string/key for legacy flag to apply strength scale to fracture energy release rate.
    static constexpr char const * scaleFractureEnergyReleaseRateString() { return "scaleFractureEnergyReleaseRate"; }

    /// string/key for exponent used to scale fracture energy with strengthScale.
    static constexpr char const * fractureEnergyStrengthScaleExponentString() { return "fractureEnergyStrengthScaleExponent"; }

    /// string/key for fracture energy release rate basal plane shear and normal modes
    static constexpr char const * basalPlaneFractureEnergyReleaseRateString() { return "basalPlaneFractureEnergyReleaseRate"; }

    /// string/key for totalFractureEnergyReleaseRate
    static constexpr char const * totalFractureEnergyReleaseRateString() { return "totalFractureEnergyReleaseRate"; }

    /// string/key for the damage-evolution exponent applied to regularized plastic work.
    static constexpr char const * damageEvolutionExponentString() { return "damageEvolutionExponent"; }

    /// string/key for damaged material frictional slope
    static constexpr char const * damagedMaterialFrictionalSlopeString() { return "damagedMaterialFrictionalSlope"; }

    // string/key for distortion shear response parameter x2
    static constexpr char const * distortionShearResponseX2String() { return "distortionShearResponseX2"; }

    // string/key for distortion shear response parameter y1
    static constexpr char const * distortionShearResponseY1String() { return "distortionShearResponseY1"; }

    // string/key for distortion shear response parameter y2
    static constexpr char const * distortionShearResponseY2String() { return "distortionShearResponseY2"; }

    // string/key for distortion shear response parameter m1
    static constexpr char const * distortionShearResponseM1String() { return "distortionShearResponseM1"; }

    // string/key for positive signed-distortion response parameter x2
    static constexpr char const * positiveDistortionShearResponseX2String() { return "positiveDistortionShearResponseX2"; }

    // string/key for positive signed-distortion response parameter y1
    static constexpr char const * positiveDistortionShearResponseY1String() { return "positiveDistortionShearResponseY1"; }

    // string/key for positive signed-distortion response parameter y2
    static constexpr char const * positiveDistortionShearResponseY2String() { return "positiveDistortionShearResponseY2"; }

    // string/key for positive signed-distortion response parameter m1
    static constexpr char const * positiveDistortionShearResponseM1String() { return "positiveDistortionShearResponseM1"; }

    // string/key for in plane shear response x2
    static constexpr char const * inPlaneShearResponseX2String() { return "inPlaneShearResponseX2"; }

    // string/key for in plane shear response y1
    static constexpr char const * inPlaneShearResponseY1String() { return "inPlaneShearResponseY1"; }

    // string/key for in plane shear response y2
    static constexpr char const * inPlaneShearResponseY2String() { return "inPlaneShearResponseY2"; }

    // string/key for in plane shear response m1
    static constexpr char const * inPlaneShearResponseM1String() { return "inPlaneShearResponseM1"; }

    // string/key for coupled shear response x2
    static constexpr char const * coupledShearResponseX2String() { return "coupledShearResponseX2"; }

    // string/key for coupled shear response y1
    static constexpr char const * coupledShearResponseY1String() { return "coupledShearResponseY1"; }

    // string/key for coupled shear response y2
    static constexpr char const * coupledShearResponseY2String() { return "coupledShearResponseY2"; }

    // string/key for coupled shear response m1
    static constexpr char const * coupledShearResponseM1String() { return "coupledShearResponseM1"; }

    // string/key for strain hardening constant C0
    static constexpr char const * distortionStrainHardeningC0() { return "distortionStrainHardeningC0"; }
    static constexpr char const * inPlaneStrainHardeningC0() { return "inPlaneStrainHardeningC0"; }
    static constexpr char const * coupledStrainHardeningC0() { return "coupledStrainHardeningC0"; }

    // string/key for maximum plastic strain
    static constexpr char const * maximumPlasticStrainString() { return "maximumPlasticStrain"; }

    /// string/key for effective bulk modulus
    static constexpr char const * effectiveBulkModulusString() { return "effectiveBulkModulus"; }

    /// string/key for effective shear modulus
    static constexpr char const * effectiveShearModulusString() { return "effectiveShearModulus"; }

    /// string/key for material direction value
    static constexpr char const * materialDirectionString() { return "materialDirection"; }
  };

  /**
   * @brief Create a instantiation of the GraphiteUpdate class that refers to the data in this.
   * @return An instantiation of GraphiteUpdate.
   */
  GraphiteUpdates createKernelUpdates() const
  {
    return GraphiteUpdates( m_defaultYoungModulusTransverse,
                            m_defaultYoungModulusAxial,
                            m_defaultPoissonRatioTransverse,
                            m_defaultPoissonRatioAxialTransverse,
                            m_defaultShearModulusAxialTransverse,
                            m_effectiveBulkModulus,
                            m_effectiveShearModulus,
                            m_materialDirection,
                            m_defaultYoungModulusTransversePressureDerivative,
                            m_defaultYoungModulusAxialPressureDerivative,
                            m_defaultShearModulusAxialTransversePressureDerivative,
                            m_defaultYoungModulusTransversePressureScale,
                            m_defaultYoungModulusAxialPressureScale,
                            m_defaultShearModulusAxialTransversePressureScale,
                            m_velocityGradient,
                            m_plasticStrain,
                            m_relaxation,
                            m_enableCrackTipStressConcentration,
                            m_crackTipStressConcentration,
                            m_distanceToCrackTip,
                            m_basalPlanePlasticWork,
                            m_plasticWork,
                            m_alphaL,
                            m_alphaT,
                            m_damage,
                            m_basalPlaneDamage,
                            m_comminutionDamage,
                            m_enableDistension,
                            m_basalNormalDistension,
                            m_transverseDistension,
                            m_porosity,
                            m_temperature,
                            m_temperatureRate,
                            m_jacobian,
                            m_lengthScale,
                            m_strengthScale,
                            m_failureStrength,
                            m_maximumPrincipalStressDamage,
                            m_crackSpeed,
                            m_scaleFractureEnergyReleaseRate,
                            m_fractureEnergyStrengthScaleExponent,
                            m_basalPlaneFractureEnergyReleaseRate,
                            m_totalFractureEnergyReleaseRate,
                            m_damageEvolutionExponent,
                            m_damagedMaterialFrictionalSlope,
                            m_distortionShearResponseX2,
                            m_distortionShearResponseY1,
                            m_distortionShearResponseY2,
                            m_distortionShearResponseM1,
                            m_positiveDistortionShearResponseX2,
                            m_positiveDistortionShearResponseY1,
                            m_positiveDistortionShearResponseY2,
                            m_positiveDistortionShearResponseM1,
                            m_inPlaneShearResponseX2,
                            m_inPlaneShearResponseY1,
                            m_inPlaneShearResponseY2,
                            m_inPlaneShearResponseM1,
                            m_coupledShearResponseX2,
                            m_coupledShearResponseY1,
                            m_coupledShearResponseY2,
                            m_coupledShearResponseM1,
                            m_distortionStrainHardeningC0,
                            m_inPlaneStrainHardeningC0,
                            m_coupledStrainHardeningC0,
                            m_maximumPlasticStrain,
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
                          m_defaultYoungModulusTransverse,
                          m_defaultYoungModulusAxial,
                          m_defaultPoissonRatioTransverse,
                          m_defaultPoissonRatioAxialTransverse,
                          m_defaultShearModulusAxialTransverse,
                          m_effectiveBulkModulus,
                          m_effectiveShearModulus,
                          m_materialDirection,
                          m_defaultYoungModulusTransversePressureDerivative,
                          m_defaultYoungModulusAxialPressureDerivative,
                          m_defaultShearModulusAxialTransversePressureDerivative,
                          m_defaultYoungModulusTransversePressureScale,
                          m_defaultYoungModulusAxialPressureScale,
                          m_defaultShearModulusAxialTransversePressureScale,
                          m_velocityGradient,
                          m_plasticStrain,
                          m_relaxation,
                          m_enableCrackTipStressConcentration,
                          m_crackTipStressConcentration,
                          m_distanceToCrackTip,
                          m_basalPlanePlasticWork,
                          m_plasticWork,
                          m_alphaL,
                          m_alphaT,
                          m_damage,
                          m_basalPlaneDamage,
                          m_comminutionDamage,
                          m_enableDistension,
                          m_basalNormalDistension,
                          m_transverseDistension,
                          m_porosity,
                          m_temperature,
                          m_temperatureRate,
                          m_jacobian,
                          m_lengthScale,
                          m_strengthScale,
                          m_failureStrength,
                          m_maximumPrincipalStressDamage,
                          m_crackSpeed,
                          m_scaleFractureEnergyReleaseRate,
                          m_fractureEnergyStrengthScaleExponent,
                          m_basalPlaneFractureEnergyReleaseRate,
                          m_totalFractureEnergyReleaseRate,
                          m_damageEvolutionExponent,
                          m_damagedMaterialFrictionalSlope,
                          m_distortionShearResponseX2,
                          m_distortionShearResponseY1,
                          m_distortionShearResponseY2,
                          m_distortionShearResponseM1,
                          m_positiveDistortionShearResponseX2,
                          m_positiveDistortionShearResponseY1,
                          m_positiveDistortionShearResponseY2,
                          m_positiveDistortionShearResponseM1,
                          m_inPlaneShearResponseX2,
                          m_inPlaneShearResponseY1,
                          m_inPlaneShearResponseY2,
                          m_inPlaneShearResponseM1,
                          m_coupledShearResponseX2,
                          m_coupledShearResponseY1,
                          m_coupledShearResponseY2,
                          m_coupledShearResponseM1,
                          m_distortionStrainHardeningC0,
                          m_inPlaneStrainHardeningC0,
                          m_coupledStrainHardeningC0,
                          m_maximumPlasticStrain,
                          m_thermalExpansionCoefficient,
                          m_newStress,
                          m_oldStress,
                          m_density,
                          m_wavespeed,
                          m_disableInelasticity );
  }

  /**
   * @brief Getter for default transverse Young's modulus
   * @return The value of the default transverse Young's modulus.
   */
  real64 getDefaultYoungModulusTransverse() const
  {
    return m_defaultYoungModulusTransverse;
  }

  /**
   * @brief Setter for the default transverse Young's modulus.
   * @param[in] input New value for the default transverse Young's modulus
   */
  void setDefaultYoungModulusTransverse( real64 const input )
  {
    m_defaultYoungModulusTransverse = input;
  }

  /**
   * @brief Getter for default axial Young's modulus
   * @return The value of the default axial Young's modulus.
   */
  real64 getDefaultYoungModulusAxial() const
  {
    return m_defaultYoungModulusAxial;
  }

  /**
   * @brief Setter for the default axial Young's modulus.
   * @param[in] input New value for the default axial Young's modulus
   */
  void setDefaultYoungModulusAxial( real64 const input )
  {
    m_defaultYoungModulusAxial = input;
  }

  /**
   * @brief Getter for default transverse Poisson's ratio
   * @return The value of the default transverse Poisson's ratio.
   */
  real64 getDefaultPoissonRatioTransverse() const
  {
    return m_defaultPoissonRatioTransverse;
  }

  /**
   * @brief Setter for the default transverse Poisson's ratio.
   * @param[in] input New value for the default transverse Poisson's ratio
   */
  void setDefaultPoissonRatioTransverse( real64 const input )
  {
    m_defaultPoissonRatioTransverse = input;
  }

  /**
   * @brief Getter for default axial Poisson's ratio
   * @return The value of the default axial/transverse Poisson's modulus.
   */
  real64 getDefaultPoissonRatioAxialTransverse() const
  {
    return m_defaultPoissonRatioAxialTransverse;
  }

  /**
   * @brief Setter for the default axial Poisson's modulus.
   * @param[in] input New value for the default axial/transverse Poisson's
   *             modulus
   */
  void setDefaultPoissonRatioAxialTransverse( real64 const input )
  {
    m_defaultPoissonRatioAxialTransverse = input;
  }

  /**
   * @brief Getter for default axial/transverse Shear modulus
   * @return The value of the default axial/transverse Shear modulus.
   */
  real64 getDefaultShearModulusAxialTransverse() const
  {
    return m_defaultShearModulusAxialTransverse;
  }

  /**
   * @brief Setter for the default axial/transverse Shear modulus.
   * @param[in] input New value for the default axial/transverse Shear modulus
   */
  void setDefaultShearModulusAxialTransverse( real64 const input )
  {
    m_defaultShearModulusAxialTransverse = input;
  }

  /**
   * @brief Accessor for effective bulk modulus
   * @return A const reference to arrayView1d<real64> containing the effective bulk
   *         modulus (at every element).
   */
  arrayView1d< real64 > const effectiveBulkModulus() { return m_effectiveBulkModulus; }

  /**
   * @brief Const accessor for effective bulk modulus
   * @return A const reference to arrayView1d<real64 const> containing the
   *         effective bulk modulus (at every element).
   */
  arrayView1d< real64 const > const effectiveBulkModulus() const { return m_effectiveBulkModulus; }

  /**
   * @brief Accessor for effective bulk modulus
   * @return A const reference to arrayView1d<real64> containing the effective bulk
   *         modulus (at every element).
   */
  arrayView1d< real64 > const effectiveShearModulus() { return m_effectiveShearModulus; }

  /**
   * @brief Const accessor for effective shear modulus
   * @return A const reference to arrayView1d<real64 const> containing the
   *         effective shear modulus (at every element).
   */
  arrayView1d< real64 const > const effectiveShearModulus() const { return m_effectiveShearModulus; }

  /**
   * @brief Getter for effective bulk modulus.
   * @return reference to mutable effective bulk modulus.
   */
  GEOS_HOST_DEVICE
  arrayView1d< real64 const > getEffectiveBulkModulus() const { return m_effectiveBulkModulus; }

  /**
   * @brief Getter for effective shear modulus.
   * @return reference to mutable effective shear modulus.
   */
  GEOS_HOST_DEVICE
  arrayView1d< real64 const > getEffectiveShearModulus() const { return m_effectiveShearModulus; }


protected:
  virtual void postInputInitialization() override;

  /// The default value of the transverse Young's modulus for any new
  /// allocations.
  real64 m_defaultYoungModulusTransverse;

  /// The default value of the axial Young's modulus for any new
  /// allocations.
  real64 m_defaultYoungModulusAxial;

  /// The default value of the transverse Poisson's ratio for any new
  /// allocations.
  real64 m_defaultPoissonRatioTransverse;

  /// The default value of the axial/transverse Poisson's ratio for any new
  /// allocations.
  real64 m_defaultPoissonRatioAxialTransverse;

  /// The default value of the axial/transverse Shear modulus for any new
  /// allocations.
  real64 m_defaultShearModulusAxialTransverse;

  /// The effective bulk modulus.
  array1d< real64 > m_effectiveBulkModulus;

  /// The effective shear modulus
  array1d< real64 > m_effectiveShearModulus;

  /// State variable: The material direction for each element/particle
  array3d< real64 > m_materialDirection;

  /// The default value of the transverse Young's modulus pressure derivative for new allocations.
  real64 m_defaultYoungModulusTransversePressureDerivative;

  /// The default value of the axial Young's modulus pressure derivative for new allocations.
  real64 m_defaultYoungModulusAxialPressureDerivative;

  /// The default value of the axial transverse Shear modulus pressure derivative for new allocations.
  real64 m_defaultShearModulusAxialTransversePressureDerivative;

  /// Pressure scale for saturating transverse Young's modulus pressure dependence.
  real64 m_defaultYoungModulusTransversePressureScale;

  /// Pressure scale for saturating axial Young's modulus pressure dependence.
  real64 m_defaultYoungModulusAxialPressureScale;

  /// Pressure scale for saturating axial-transverse shear modulus pressure dependence.
  real64 m_defaultShearModulusAxialTransversePressureScale;

  ///State variable: The velocity gradient for each element/particle
  array3d< real64 > m_velocityGradient;

  ///State variable: The plastic strain values for each quadrature point
  array3d< real64 > m_plasticStrain;

  /// State variable: The relaxation values for each quadrature point
  array2d< real64 > m_relaxation;

  // Crack-tip stress concentration variables:
  int m_enableCrackTipStressConcentration;
  array1d< real64 > m_crackTipStressConcentration;
  array1d< real64 > m_distanceToCrackTip;

  /// State variable: The accumulated plastic work values for each quadrature point
  array2d< real64 > m_basalPlanePlasticWork;
  array2d< real64 > m_plasticWork;

  /// Material parameter: lateral thermal expansion coefficient
  real64 m_alphaL;

  /// Material parameter: transverse thermal expansion coefficient
  real64 m_alphaT;

  /// State variable: Host-visible DFG damage values for each quadrature point.
  array2d< real64 > m_damage;

  /// State variable: basal-plane fracture damage for each quadrature point.
  array2d< real64 > m_basalPlaneDamage;

  /// State variable: comminution/powder damage for each quadrature point.
  array2d< real64 > m_comminutionDamage;

  /// Optional flag to activate stress-free closure of directional distension.
  int m_enableDistension;

  /// State variable: closeable distension normal to basal planes.
  array2d< real64 > m_basalNormalDistension;

  /// State variable: closeable distension in the basal/transverse plane.
  array2d< real64 > m_transverseDistension;

  /// State variable: distension-derived porosity for plotting/coupling.
  array2d< real64 > m_porosity;

  /// State variable: The temperature values for each element/particle
  array1d< real64 > m_temperature;

  /// State variable: The temperature rate values for each element/particle
  array1d< real64 > m_temperatureRate;


  /// State variable: The jacobian of the deformation
  array2d< real64 > m_jacobian;

  /// Discretization-sized variable: The length scale for each element/particle
  array1d< real64 > m_lengthScale;

  /// Discretization-sized variable: The strength scale for each element/particle
  array1d< real64 > m_strengthScale;

  /// Material parameter: The value of the failure strength
  real64 m_failureStrength;

  int m_maximumPrincipalStressDamage;
  real64 m_crackSpeed;

  /// Material parameter: The value of fracture energy release rate
  int m_scaleFractureEnergyReleaseRate;
  real64 m_fractureEnergyStrengthScaleExponent;
  real64 m_basalPlaneFractureEnergyReleaseRate;
  real64 m_totalFractureEnergyReleaseRate;
  real64 m_damageEvolutionExponent;

  /// Material parameter: The value of the damaged material frictional slope
  real64 m_damagedMaterialFrictionalSlope;

  /// Material parameters: The values controlling the failure envelope.
  /// The distortion branch is the negative signed-distortion branch and the
  /// optional positive branch defaults to these values when omitted.
  real64 m_distortionShearResponseX2;
  real64 m_distortionShearResponseY1;
  real64 m_distortionShearResponseY2;
  real64 m_distortionShearResponseM1;

  real64 m_positiveDistortionShearResponseX2;
  real64 m_positiveDistortionShearResponseY1;
  real64 m_positiveDistortionShearResponseY2;
  real64 m_positiveDistortionShearResponseM1;

  real64 m_inPlaneShearResponseX2;
  real64 m_inPlaneShearResponseY1;
  real64 m_inPlaneShearResponseY2;
  real64 m_inPlaneShearResponseM1;

  real64 m_coupledShearResponseX2;
  real64 m_coupledShearResponseY1;
  real64 m_coupledShearResponseY2;
  real64 m_coupledShearResponseM1;

  real64 m_distortionStrainHardeningC0;
  real64 m_inPlaneStrainHardeningC0;
  real64 m_coupledStrainHardeningC0;

  real64 m_maximumPlasticStrain;
};

} /* namespace constitutive */

} /* namespace geos */

#endif /* GEOSX_CONSTITUTIVE_SOLID_GRAPHITE_HPP_ */
