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
 *  @file JohnsonCook.hpp
 */


#ifndef GEOS_CONSTITUTIVE_SOLID_JOHNSONCOOK_HPP_
#define GEOS_CONSTITUTIVE_SOLID_JOHNSONCOOK_HPP_

#include "ElasticIsotropic.hpp"
#include "InvariantDecompositions.hpp"
#include "SolidModelDiscretizationOpsIsotropic.hpp"
#include "constitutive/ExponentialRelation.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geos
{
namespace constitutive
{

/**
 * @class JohnsonCookUpdates
 *
 * Class to provide elastic isotropic material updates that incorporate
 * Johnson–Cook plasticity with strain rate sensitivity and thermal effects.
 */
class JohnsonCookUpdates : public ElasticIsotropicUpdates
{
public:
  /**
   * @brief Constructor
   * @param[in] A                    ArrayView holding the initial yield strength (A parameter)
   * @param[in] B                    ArrayView holding the hardening modulus (B parameter)
   * @param[in] n                    ArrayView holding the strain hardening exponent (n parameter)
   * @param[in] C                    ArrayView holding the strain rate sensitivity coefficient (C parameter)
   * @param[in] m                    ArrayView holding the thermal softening exponent (m parameter)
   * @param[in] meltingTemperature   ArrayView holding the melting temperature (T_m)
   * @param[in] referenceTemperature ArrayView holding the reference temperature (T_0)
   * @param[in] deformationGradient  ArrayView holding the deformation gradient field
   * @param[in] velocityGradient     ArrayView holding the velocity gradient field
   * @param[in] plasticStrain        ArrayView holding the plastic strain field
   * @param[in] bulkModulus          ArrayView holding the bulk modulus data
   * @param[in] shearModulus         ArrayView holding the shear modulus data
   * @param[in] thermalExpansionCoefficient ArrayView holding the thermal expansion coefficient data
   * @param[in] newStress            ArrayView holding the new stress field
   * @param[in] oldStress            ArrayView holding the old stress field
   * @param[in] density              ArrayView holding the density data
   * @param[in] wavespeed            ArrayView holding the wave speed data
   * @param[in] temperature          ArrayView holding the temperature field
   */
  JohnsonCookUpdates( arrayView1d< real64 > const & A,
                      arrayView1d< real64 > const & B,
                      arrayView1d< real64 > const & n,
                      arrayView1d< real64 > const & C,
                      arrayView1d< real64 > const & m,
                      arrayView1d< real64 > const & meltingTemperature,
                      arrayView1d< real64 > const & referenceTemperature,
                      arrayView3d< real64 const > const & deformationGradient,
                      arrayView3d< real64 const > const & velocityGradient,
                      arrayView3d< real64 > const & plasticStrain,
                      arrayView1d< real64 const > const & bulkModulus,
                      arrayView1d< real64 const > const & shearModulus,
                      arrayView1d< real64 const > const & thermalExpansionCoefficient,
                      arrayView3d< real64, solid::STRESS_USD > const & newStress,
                      arrayView3d< real64, solid::STRESS_USD > const & oldStress,
                      arrayView2d< real64 > const & density,
                      arrayView2d< real64 > const & wavespeed
                      //arrayView1d< real64 > const & temperature,
                      //arrayView1d< real64 > const & yieldStrength,
                      //arrayView1d< real64 > const & defaultYieldStrength,   // new argument
                      //arrayView1d< real64 > const & referenceStrainRate  
                      ):
    ElasticIsotropicUpdates( bulkModulus,
                             shearModulus,
                             thermalExpansionCoefficient,
                             newStress,
                             oldStress,
                             density,
                             wavespeed,
                             /* disableInelasticity = */ false ),
    m_A( A ),
    m_B( B ),
    m_n( n ),
    m_C( C ),
    m_m( m ),
    m_meltingTemperature( meltingTemperature ),
    m_referenceTemperature( referenceTemperature ),
    m_deformationGradient( deformationGradient ),
    m_velocityGradient( velocityGradient ),
    m_plasticStrain( plasticStrain )
  
  {}

  // Deleted default constructor
  JohnsonCookUpdates() = delete;

  // Default copy and move constructors
  JohnsonCookUpdates( JohnsonCookUpdates const & ) = default;
  JohnsonCookUpdates( JohnsonCookUpdates && ) = default;

  // Deleted assignment operators
  JohnsonCookUpdates & operator=( JohnsonCookUpdates const & ) = delete;
  JohnsonCookUpdates & operator=( JohnsonCookUpdates && ) = delete;

  /// Use the "isotropic" form of inner product compression
  using DiscretizationOps = SolidModelDiscretizationOpsIsotropic;

  /// Use base version of saveConvergedState
  using ElasticIsotropicUpdates::saveConvergedState;

  GEOS_HOST_DEVICE
  virtual void smallStrainNoStateUpdate_StressOnly( localIndex const k,
                                                    localIndex const q,
                                                    real64 const ( &totalStrain )[6],
                                                    real64 ( &stress )[6] ) const override final;

  GEOS_HOST_DEVICE
  virtual void smallStrainNoStateUpdate( localIndex const k,
                                         localIndex const q,
                                         real64 const ( &totalStrain )[6],
                                         real64 ( &stress )[6],
                                         real64 ( &stiffness )[6][6] ) const override final;

  GEOS_HOST_DEVICE
  virtual void smallStrainNoStateUpdate( localIndex const k,
                                         localIndex const q,
                                         real64 const ( &totalStrain )[6],
                                         real64 ( &stress )[6],
                                         DiscretizationOps & stiffness ) const override final;

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
                                             real64 const ( & strainIncrement )[6],
                                             real64 ( & stress )[6] ) const override;

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
                                  DiscretizationOps & stiffness ) const override;

  GEOS_HOST_DEVICE
  virtual void getElasticStiffness( localIndex const k,
                                    localIndex const q,
                                    real64 ( &stiffness )[6][6] ) const override;

  GEOS_HOST_DEVICE
  virtual void getElasticStrain( localIndex const k,
                                 localIndex const q,
                                 real64 ( &elasticStrain )[6] ) const override final;

protected:
  /// Johnson–Cook parameters:
  arrayView1d< real64 > const m_A;                    ///< Initial yield strength parameter (A)
  arrayView1d< real64 > const m_B;                    ///< Hardening modulus (B)
  arrayView1d< real64 > const m_n;                    ///< Strain hardening exponent (n)
  arrayView1d< real64 > const m_C;                    ///< Strain rate sensitivity coefficient (C)
  arrayView1d< real64 > const m_m;                    ///< Thermal softening exponent (m)
  arrayView1d< real64 > const m_meltingTemperature;   ///< Melting temperature (T_m)
  arrayView1d< real64 > const m_referenceTemperature; ///< Reference temperature (T_0)

  /// Kinematic and state variables:
  arrayView3d< real64 const > const m_deformationGradient;
  arrayView3d< real64 const > const m_velocityGradient;
  arrayView3d< real64 > const m_plasticStrain;
};

//
// Inline definitions for JohnsonCookUpdates member functions
//

GEOS_HOST_DEVICE
inline
void JohnsonCookUpdates::getElasticStiffness( localIndex const k,
                                              localIndex const q,
                                              real64 ( & stiffness )[6][6] ) const
{
  GEOS_UNUSED_VAR( q );
  real64 const G = m_shearModulus[k];
  real64 const lambda = conversions::bulkModAndShearMod::toFirstLame( m_bulkModulus[k], G );

  LvArray::tensorOps::fill< 6, 6 >( stiffness, 0 );
  stiffness[0][0] = lambda + 2*G;
  stiffness[0][1] = lambda;
  stiffness[0][2] = lambda;
  stiffness[1][0] = lambda;
  stiffness[1][1] = lambda + 2*G;
  stiffness[1][2] = lambda;
  stiffness[2][0] = lambda;
  stiffness[2][1] = lambda;
  stiffness[2][2] = lambda + 2*G;
  stiffness[3][3] = G;
  stiffness[4][4] = G;
  stiffness[5][5] = G;
}

GEOS_HOST_DEVICE
inline
void JohnsonCookUpdates::getElasticStrain( localIndex const k,
                                           localIndex const q,
                                           real64 ( & elasticStrain)[6] ) const
{
    GEOS_UNUSED_VAR( k );
    GEOS_UNUSED_VAR( q );
    GEOS_UNUSED_VAR( elasticStrain );
    GEOS_ERROR( "getElasticStrain not implemented for JohnsonCook model" );
}

GEOS_HOST_DEVICE
inline
void JohnsonCookUpdates::smallStrainNoStateUpdate_StressOnly( localIndex const k,
                                                               localIndex const q,
                                                               real64 const ( &totalStrain )[6],
                                                               real64 ( & stress )[6] ) const
{
  GEOS_UNUSED_VAR( k );
  GEOS_UNUSED_VAR( q );
  GEOS_UNUSED_VAR( totalStrain );
  GEOS_UNUSED_VAR( stress );
  GEOS_ERROR( "smallStrainNoStateUpdate_StressOnly is not implemented for JohnsonCook model" );
}

GEOS_HOST_DEVICE
inline
void JohnsonCookUpdates::smallStrainNoStateUpdate( localIndex const k,
                                                   localIndex const q,
                                                   real64 const ( &totalStrain )[6],
                                                   real64 ( & stress )[6],
                                                   real64 ( & stiffness )[6][6] ) const
{
  smallStrainNoStateUpdate_StressOnly( k, q, totalStrain, stress );
  getElasticStiffness( k, q, stiffness );
}

GEOS_HOST_DEVICE
inline
void JohnsonCookUpdates::smallStrainNoStateUpdate( localIndex const k,
                                                   localIndex const q,
                                                   real64 const ( &totalStrain )[6],
                                                   real64 ( & stress )[6],
                                                   DiscretizationOps & stiffness ) const
{
  smallStrainNoStateUpdate_StressOnly( k, q, totalStrain, stress );
  stiffness.m_bulkModulus = m_bulkModulus[k];
  stiffness.m_shearModulus = m_shearModulus[k];
}

GEOS_HOST_DEVICE
inline
void JohnsonCookUpdates::smallStrainUpdate_StressOnly( localIndex const k,
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
  GEOS_ERROR( "smallStrainUpdate_StressOnly is not implemented for JohnsonCook model" );
}

GEOS_HOST_DEVICE
inline
void JohnsonCookUpdates::smallStrainUpdate_StressOnly( localIndex const k,
                                                       localIndex const q,
                                                       real64 const & timeIncrement,
                                                       real64 const ( & beginningRotation )[3][3],
                                                       real64 const ( & endRotation )[3][3],
                                                       real64 const ( & strainIncrement )[6],
                                                       real64 ( & stress )[6] ) const
{
  // Johnson–Cook model implementation:
  // 1. Compute trial stress via elastic predictor.
  real64 previousStress[6] = { 0 };
  LvArray::tensorOps::copy< 6 >( previousStress, m_oldStress[k][q] );
  
  real64 trialP;
  real64 trialQ;
  real64 oldDeviatoricStress[6] = { 0 };
  twoInvariant::stressDecomposition( previousStress,
                                      trialP,
                                      trialQ,
                                      oldDeviatoricStress );
  LvArray::tensorOps::scale< 6 >( oldDeviatoricStress, sqrt(2.0 / 3.0)*trialQ );
  
  // 2. Compute pressure based on the deformation gradient.
  real64 J = LvArray::tensorOps::determinant< 3 >( m_deformationGradient[k] );
  real64 pressure = -m_bulkModulus[k] * std::log( J );
  
  // 3. Compute the elastic deviatoric stress increment.
  real64 rotationTranspose[3][3] = { { 0 } };
  LvArray::tensorOps::transpose< 3, 3 >( rotationTranspose, beginningRotation );
  
  real64 unrotatedVelocityGradient[3][3]  = { { 0 } };
  LvArray::tensorOps::Rij_eq_AikBkj< 3, 3, 3 >( unrotatedVelocityGradient, rotationTranspose, m_velocityGradient[k] );
  
  real64 unrotatedVelocityGradientTranspose[3][3]  = { { 0 } };
  LvArray::tensorOps::transpose< 3, 3 >( unrotatedVelocityGradientTranspose, unrotatedVelocityGradient );
  
  real64 denseD[3][3] = { { 0 } };
  LvArray::tensorOps::copy< 3, 3 >( denseD, unrotatedVelocityGradient );
  LvArray::tensorOps::add< 3, 3 >( denseD, unrotatedVelocityGradientTranspose );
  LvArray::tensorOps::scale< 3, 3 >( denseD, 0.5 );
  
  real64 D[6] = { 0 };
  LvArray::tensorOps::denseToSymmetric< 3 >( D, denseD );
  
  real64 trD = LvArray::tensorOps::symTrace< 3 >( D ) / 3.0;
  D[0] -= trD;
  D[1] -= trD;
  D[2] -= trD;
  
  real64 deviatoricStressIncrement[6] = { 0 };
  LvArray::tensorOps::scaledCopy< 6 >( deviatoricStressIncrement, D, 2.0 * m_shearModulus[k] * timeIncrement );
  
  LvArray::tensorOps::copy< 6 >( stress, oldDeviatoricStress );
  LvArray::tensorOps::add< 6 >( stress, deviatoricStressIncrement );
  LvArray::tensorOps::symAddIdentity< 3 >( stress, -pressure );
  
  // 4. Johnson–Cook plasticity update:
  // Decompose the stress and compute the effective (trial) stress invariant.
  real64 deviator[6] = { 0 };
  twoInvariant::stressDecomposition( stress,
                                      trialP,
                                      trialQ,
                                      deviator );
  
  // Retrieve the current equivalent plastic strain (assumed stored in first component).
  real64 plasticStrain = m_plasticStrain[k][q][0];
  // Estimate a simple effective strain rate (this can be refined).
  real64 strainRate = sqrt(2.0/3.0) * plasticStrain;
  // Compute normalized temperature.
  real64 normalizedTemp = ( m_temperature[k] - m_referenceTemperature[k] ) /
                            ( m_meltingTemperature[k] - m_referenceTemperature[k] );
  normalizedTemp = std::min( 1.0, std::max( 0.0, normalizedTemp ) );
  
  // Compute current yield stress using the Johnson–Cook relation:
  // yieldStress = A + B*(plasticStrain)^n * (1 + C*log(strainRate)) * (1 - (normalizedTemp)^m)
  real64 currentYield = m_A[k] + m_B[k]*std::pow( plasticStrain, m_n[k] )*
                        ( 1.0 + m_C[k]*std::log( strainRate ) )*
                        ( 1.0 - std::pow( normalizedTemp, m_m[k] ) );
  
  if( trialQ > currentYield )
  {
    real64 oldStress[6] = { 0 };
    LvArray::tensorOps::copy< 6 >( oldStress, stress );
    
    twoInvariant::stressRecomposition( trialP,
                                       currentYield,
                                       deviator,
                                       stress );
    
    real64 stressIncrement[6] = {0};
    LvArray::tensorOps::copy< 6 >( stressIncrement, stress );
    LvArray::tensorOps::subtract< 6 >( stressIncrement, oldStress );
    
    // Compute plastic strain increment (simplified update).
    real64 plasticStrainIncrement[6] = {0};
    LvArray::tensorOps::copy< 6 >( plasticStrainIncrement, strainIncrement );
    LvArray::tensorOps::subtract< 6 >( plasticStrainIncrement, stressIncrement );
    
    // Update plastic strain
    real64 oldPlasticStrain[6] = { 0 };
    LvArray::tensorOps::copy< 6 >( oldPlasticStrain, m_plasticStrain[k][q] );
    LvArray::tensorOps::add< 6 >( oldPlasticStrain, plasticStrainIncrement );
    LvArray::tensorOps::copy< 6 >( m_plasticStrain[k][q], oldPlasticStrain );
  }
  
  saveStress( k, q, stress );
}

GEOS_HOST_DEVICE
inline
void JohnsonCookUpdates::smallStrainUpdate( localIndex const k,
                                            localIndex const q,
                                            real64 const & timeIncrement,
                                            real64 const ( &strainIncrement )[6],
                                            real64 ( & stress )[6],
                                            real64 ( & stiffness )[6][6] ) const
{
  smallStrainUpdate_StressOnly( k,
                                q,
                                timeIncrement,
                                strainIncrement,
                                stress );
  getElasticStiffness( k, q, stiffness );
}

GEOS_HOST_DEVICE
inline
void JohnsonCookUpdates::smallStrainUpdate( localIndex const k,
                                            localIndex const q,
                                            real64 const & timeIncrement,
                                            real64 const ( &strainIncrement )[6],
                                            real64 ( & stress )[6],
                                            DiscretizationOps & stiffness ) const
{
  smallStrainUpdate_StressOnly( k,
                                q,
                                timeIncrement,
                                strainIncrement,
                                stress );
  stiffness.m_bulkModulus = m_bulkModulus[k];
  stiffness.m_shearModulus = m_shearModulus[k];
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void JohnsonCookUpdates::computePlasticStrainIncrement( localIndex const k,
                                                        localIndex const q,
                                                        const real64 timeIncrement,
                                                        real64 const ( & strainIncrement )[6],
                                                        real64 const ( & stressIncrement )[6],
                                                        real64 ( & plasticStrainIncrement )[6] ) const
{ 
  GEOS_UNUSED_VAR( q );
  GEOS_UNUSED_VAR( timeIncrement );
  
  // For the Johnson–Cook model the elastic strain increment is computed 
  // via an isotropic–deviatoric decomposition of the stress increment.
  
  real64 trialP;
  real64 trialQ;
  real64 stressIncrementDeviator[6] = {0};
  twoInvariant::stressDecomposition( stressIncrement,
                                     trialP,
                                     trialQ,
                                     stressIncrementDeviator );
  
  real64 stressIncrementIsostatic[6] = {0};
  stressIncrementIsostatic[0] = trialP;
  stressIncrementIsostatic[1] = trialP;
  stressIncrementIsostatic[2] = trialP;
  
  // Compute the elastic strain increment based on the elastic moduli.
  real64 elasticStrainIncrement[6] = {0};
  for( int i = 0; i < 6; ++i )
  {
    if( m_bulkModulus[k] > 1.0e-12 )
    {
      // Off-diagonals get a factor of 2; here (1 + (i >= 3)) equals 2 for shear components.
      elasticStrainIncrement[i] += ( 1 + (i >= 3) ) * stressIncrementIsostatic[i] / ( 3.0 * m_bulkModulus[k] );
    }
    if( m_shearModulus[k] > 1.0e-12 )
    {
      elasticStrainIncrement[i] += ( 1 + (i >= 3) ) * sqrt(2.0/3.0) * trialQ * stressIncrementDeviator[i] / ( 2.0 * m_shearModulus[k] );
    }
  }
  
  // The plastic strain increment is the difference between the total strain increment 
  // and the computed elastic strain increment.
  LvArray::tensorOps::copy< 6 >( plasticStrainIncrement, strainIncrement );
  LvArray::tensorOps::subtract< 6 >( plasticStrainIncrement, elasticStrainIncrement );
}

/**
 * @class JohnsonCook
 *
 * Class to provide an elastic isotropic material response with Johnson–Cook plasticity.
 */
class JohnsonCook : public ElasticIsotropic
{
public:

  /// Alias for JohnsonCookUpdates
  using KernelWrapper = JohnsonCookUpdates;

  /**
   * Constructor
   * @param[in] name Name of the instance in the catalog
   * @param[in] parent The group which contains this instance
   */
  JohnsonCook( string const & name, Group * const parent );

  /**
   * Default Destructor
   */
  virtual ~JohnsonCook() override;

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  /**
   * @name Static Factory Catalog members and functions
   */
  ///@{

  /// String name to use for this class in the catalog
  static constexpr auto m_catalogNameString = "JohnsonCook";

  /**
   * @brief Static catalog string
   * @return A string that is used to register/lookup this class in the registry
   */
  static std::string catalogName() { return m_catalogNameString; }

  /**
   * @brief Get catalog name
   * @return Name string
   */
  virtual string getCatalogName() const override { return catalogName(); }

  ///@}

  /// Keys for data specified in this class.
  struct viewKeyStruct : public ElasticIsotropic::viewKeyStruct
  {
    /// string/key for initial yield strength (A parameter)
    static constexpr char const * initialYieldStrengthString() { return "initialYieldStrength"; }

    /// string/key for hardening modulus (B parameter)
    static constexpr char const * hardeningModulusString() { return "hardeningModulus"; }

    /// string/key for strain hardening exponent (n parameter)
    static constexpr char const * strainHardeningExponentString() { return "strainHardeningExponent"; }

    /// string/key for strain rate sensitivity coefficient (C parameter)
    static constexpr char const * strainRateCoeffString() { return "strainRateCoeff"; }

    /// string/key for thermal softening exponent (m parameter)
    static constexpr char const * thermalSofteningExponentString() { return "thermalSofteningExponent"; }

    /// string/key for melting temperature
    static constexpr char const * meltingTemperatureString() { return "meltingTemperature"; }

    /// string/key for reference temperature
    static constexpr char const * referenceTemperatureString() { return "referenceTemperature"; }

    /// string/key for deformation gradient
    static constexpr char const * deformationGradientString() { return "deformationGradient"; }

    /// string/key for velocity gradient
    static constexpr char const * velocityGradientString() { return "velocityGradient"; }

    /// string/key for plastic strain
    static constexpr char const * plasticStrainString() { return "plasticStrain"; }

    // Add the static function to return the string "yieldStrength"
    static constexpr char const * yieldStrengthString() { return "yieldStrength"; }

    // Add this method for reference strain rate
    static constexpr char const * referenceStrainRateString() { return "referenceStrainRate"; }

    // Add this method for reference strain rate
    static constexpr char const * plasticStrainRateString() { return "plasticStrainRate"; }

    static constexpr const char* defaultYieldStrengthString() {return "defaultYieldStrength";}  // Return the name of the field as a string

    


  };

  GEOS_HOST_DEVICE
  virtual arrayView1d< real64 const > getInitialYieldStrength() const final
  {
    return m_initialYieldStrength;
  }

  GEOS_HOST_DEVICE
  virtual arrayView1d< real64 const > getHardeningModulus() const final
  {
    return m_hardeningModulus;
  }

  GEOS_HOST_DEVICE
  virtual arrayView1d< real64 const > getStrainHardeningExponent() const final
  {
    return m_strainHardeningExponent;
  }

  GEOS_HOST_DEVICE
  virtual arrayView1d< real64 const > getStrainRateCoeff() const final
  {
    return m_strainRateCoeff;
  }

  GEOS_HOST_DEVICE
  virtual arrayView1d< real64 const > getThermalSofteningExponent() const final
  {
    return m_thermalSofteningExponent;
  }

  GEOS_HOST_DEVICE
  virtual arrayView1d< real64 const > getMeltingTemperature() const final
  {
    return m_meltingTemperature;
  }

  GEOS_HOST_DEVICE
  virtual arrayView1d< real64 const > getReferenceTemperature() const final
  {
    return m_referenceTemperature;
  }

  GEOS_HOST_DEVICE
  virtual arrayView3d< real64 const > getDeformationGradient() const final
  {
    return m_deformationGradient;
  }

  GEOS_HOST_DEVICE
  virtual arrayView3d< real64 const > getVelocityGradient() const final
  {
    return m_velocityGradient;
  }

  GEOS_HOST_DEVICE
  virtual arrayView3d< real64 const > getPlasticStrain() const final
  {
    return m_plasticStrain;
  }

  /**
   * @brief Create a instantiation of the JohnsonCookUpdates class
   *        that refers to the data in this.
   * @param includeState Flag whether to pass state arrays that may not be needed for "no-state" updates
   * @return An instantiation of JohnsonCookUpdates.
   */
  JohnsonCookUpdates createKernelUpdates( bool const includeState = true ) const
  {
    if( includeState )
    {
      return JohnsonCookUpdates( m_initialYieldStrength,
                                 m_hardeningModulus,
                                 m_strainHardeningExponent,
                                 m_strainRateCoeff,
                                 m_thermalSofteningExponent,
                                 m_meltingTemperature,
                                 m_referenceTemperature,
                                 m_deformationGradient,
                                 m_velocityGradient,
                                 m_plasticStrain,
                                 m_bulkModulus,
                                 m_shearModulus,
                                 m_thermalExpansionCoefficient,
                                 m_newStress,
                                 m_oldStress,
                                 m_density,
                                 m_wavespeed,
                                 /* disableInelasticity = */ false );
    }
    else // for "no state" updates, pass empty views to avoid transfer of stress data to device
    {
      return JohnsonCookUpdates( m_initialYieldStrength,
                                 m_hardeningModulus,
                                 m_strainHardeningExponent,
                                 m_strainRateCoeff,
                                 m_thermalSofteningExponent,
                                 m_meltingTemperature,
                                 m_referenceTemperature,
                                 m_deformationGradient,
                                 m_velocityGradient,
                                 m_plasticStrain,
                                 m_bulkModulus,
                                 m_shearModulus,
                                 m_thermalExpansionCoefficient,
                                 arrayView3d< real64, solid::STRESS_USD >(),
                                 arrayView3d< real64, solid::STRESS_USD >(),
                                 m_density,
                                 m_wavespeed,
                                 /* disableInelasticity = */ false );
    }
  }

  /**
   * @brief Construct an update kernel for a derived type.
   * @tparam UPDATE_KERNEL The type of update kernel from the derived type.
   * @tparam PARAMS The parameter pack to hold the constructor parameters for
   *   the derived update kernel.
   * @param constructorParams The constructor parameter for the derived type.
   * @return An @p UPDATE_KERNEL object.
   */
  template< typename UPDATE_KERNEL, typename ... PARAMS >
  UPDATE_KERNEL createDerivedKernelUpdates( PARAMS && ... constructorParams ) const
  {
    return UPDATE_KERNEL( std::forward< PARAMS >( constructorParams )...,
                          m_initialYieldStrength,
                          m_hardeningModulus,
                          m_strainHardeningExponent,
                          m_strainRateCoeff,
                          m_thermalSofteningExponent,
                          m_meltingTemperature,
                          m_referenceTemperature,
                          m_deformationGradient,
                          m_velocityGradient,
                          m_plasticStrain,
                          m_bulkModulus,
                          m_shearModulus,
                          m_thermalExpansionCoefficient,
                          m_newStress,
                          m_oldStress,
                          m_density,
                          m_wavespeed,
                          /* disableInelasticity = */ false );
  }

protected:

  /// Post-process XML data
  virtual void postInputInitialization() override;

  /// The default value of the initial yield strength for any new allocations
  real64 m_defaultInitialYieldStrength;

  /// The initial yield strength (A parameter) for each upper level dimension (i.e. cell) of *this
  array1d< real64 > m_initialYieldStrength;

  /// The hardening modulus (B parameter)
  array1d< real64 > m_hardeningModulus;

  /// The strain hardening exponent (n parameter)
  array1d< real64 > m_strainHardeningExponent;

  /// The strain rate sensitivity coefficient (C parameter)
  array1d< real64 > m_strainRateCoeff;

  /// The thermal softening exponent (m parameter)
  array1d< real64 > m_thermalSofteningExponent;

  /// The melting temperature (T_m) for each upper level dimension (i.e. cell) of *this
  array1d< real64 > m_meltingTemperature;

  /// The reference temperature (T_0)
  array1d< real64 > m_referenceTemperature;

  /// The deformation gradient for each upper level dimension (i.e. cell) of *this
  array3d< real64 > m_deformationGradient;

  /// The velocity gradient for each upper level dimension (i.e. cell) of *this
  array3d< real64 > m_velocityGradient;

  /// The plastic strain for each upper level dimension (i.e. cell) of *this
  array3d< real64 > m_plasticStrain;

};


} /* namespace constitutive */
} /* namespace geos */

#endif /* GEOS_CONSTITUTIVE_SOLID_JOHNSONCOOK_HPP_ */



