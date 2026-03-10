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
 * @file HyperViscoElastoPlastic.hpp
 * @brief A large strain hyperelastic visco plastic model (Nguyen et al. Int. J. Solids Struct. 2016)
 *
 */

#ifndef GEOS_CONSTITUTIVE_SOLID_HYPERVISCOELASTOPLASTIC_HPP_
#define GEOS_CONSTITUTIVE_SOLID_HYPERVISCOELASTOPLASTIC_HPP_

#include "SolidBase.hpp"
#include "InvariantDecompositions.hpp"
#include "PropertyConversions.hpp"
#include "SolidModelDiscretizationOpsFullyAnisotroipic.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geos
{

namespace constitutive
{

/**
 * @class HyperViscoElastoPlasticUpdates
 *
 * Class to provide material updates that may be
 * called from a kernel function.
 */
class HyperViscoElastoPlasticUpdates : public SolidBaseUpdates
{
public:

  /**
   * @brief Constructor
   * @param[in] deformationGradient The ArrayView holding the deformation gradient for each element/particle.
   * @param[in] plasticStrain The ArrayView holding the plastic strain for each quadrature point
   * @param[in] damage The ArrayView holding the damage for each quadrature point.
   * @param[in] jacobian The ArrayView holding the jacobian for each quadrature point.
   * @param[in] newStress The ArrayView holding the new stress data for each quadrature point.
   * @param[in] oldStress The ArrayView holding the old stress data for each quadrature point.
   */
  HyperViscoElastoPlasticUpdates( arrayView3d< real64 > const & deformationGradient,
                                  elasticStrain,
                                  plasticStrain,
                                  Ain,
                                  Bin,
                                  youngModuli,
                                  elasticPoissonRatio,
                                  relaxationTimes,
                                  arrayView1d< real64 const > const & thermalExpansionCoefficient,
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
    m_deformationGradient( deformationGradient ),
    m_plasticStrain( plasticStrain ),
    m_damage( damage ),
    m_jacobian( jacobian ),
    m_yieldStrength( yieldStrength ),
    m_strainHardeningSlope( strainHardeningSlope ),
    m_shearSofteningMagnitude( shearSofteningMagnitude ),
    m_shearSofteningShapeParameter1( shearSofteningShapeParameter1 ),
    m_shearSofteningShapeParameter2( shearSofteningShapeParameter2 ),
    m_maximumStretch( maximumStretch )
  {}

  /// Default copy constructor
  HyperViscoElastoPlasticUpdates( HyperViscoElastoPlasticUpdates const & ) = default;

  /// Default move constructor
  HyperViscoElastoPlasticUpdates( HyperViscoElastoPlasticUpdates && ) = default;

  /// Deleted default constructor
  HyperViscoElastoPlasticUpdates() = delete;

  /// Deleted copy assignment operator
  HyperViscoElastoPlasticUpdates & operator=( HyperViscoElastoPlasticUpdates const & ) = delete;

  /// Deleted move assignment operator
  HyperViscoElastoPlasticUpdates & operator=( HyperViscoElastoPlasticUpdates && ) =  delete;

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
                                             real64 const ( &beginningRotation )[3][3],
                                             real64 const ( &endRotation )[3][3],
                                             real64 const ( &strainIncrement )[6],
                                             real64 ( &stress )[6] ) const override;

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  virtual void saveConvergedState( localIndex const k,
                                   localIndex const q ) const override final
  {
    ElasticIsotropicUpdates::saveConvergedState( k, q );
  }

private:
  /// A reference to the ArrayView holding the deformation gradient for each element/particle.
  arrayView3d< real64 > const m_deformationGradient;

  /// A reference to the ArrayView holding the plastic strain for each quadrature point.
  arrayView3d< real64 > const m_plasticStrain;

  /// A reference to the ArrayView holding the damage for each quadrature point.
  arrayView2d< real64 > const m_damage;

  /// A reference to the ArrayView holding the jacobian for each quadrature point.
  arrayView2d< real64 > const m_jacobian;

  /// A reference to the ArrayView holding the yield strength for each element/particle
  arrayView1d< real64 > const m_yieldStrength;

  /// The strain hardening slope
  real64 const m_strainHardeningSlope;

  /// The shear softening magnitude
  real64 m_shearSofteningMagnitude;

  /// The shear softening shape parameter 1
  real64 m_shearSofteningShapeParameter1;

  /// The shear softening shape parameter 2
  real64 m_shearSofteningShapeParameter2;

  /// The compressive strength
  real64 const m_maximumStretch;

};


GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void HyperViscoElastoPlasticUpdates::smallStrainUpdate( localIndex const k,
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
  GEOS_ERROR( "smallStrainUpdate not implemented for StrainHardeningPolymer" );
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void HyperViscoElastoPlasticUpdates::smallStrainUpdate( localIndex const k,
                                                        localIndex const q,
                                                        real64 const & timeIncrement,
                                                        real64 const ( &strainIncrement )[6],
                                                        real64 ( & stress )[6],
                                                        DiscretizationOps & stiffness ) const
{
  GEOS_UNUSED_VAR( k );
  GEOS_UNUSED_VAR( q );
  GEOS_UNUSED_VAR( timeIncrement );
  GEOS_UNUSED_VAR( strainIncrement );
  GEOS_UNUSED_VAR( stress );
  GEOS_UNUSED_VAR( stiffness );
  GEOS_ERROR( "smallStrainUpdate not implemented for StrainHardeningPolymer" );
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void HyperViscoElastoPlasticUpdates::smallStrainUpdate_StressOnly( localIndex const k,
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
void HyperViscoElastoPlasticUpdates::smallStrainUpdate_StressOnly( localIndex const k,
                                                                   localIndex const q,
                                                                   real64 const & timeIncrement,
                                                                   real64 const ( &beginningRotation )[3][3],
                                                                   real64 const ( &endRotation )[3][3],
                                                                   real64 const ( &strainIncrement )[6],
                                                                   real64 ( & stress )[6] ) const
{
  // Elastic predictor (assume strainIncrement is all elastic)
  real64 dt = timeIncrement;
  localIndex numMaxwellElements = m_relaxationTimes[k].size();

  // Compute elastic predictor strain
  real64 rightCauchyDeformationTensor[3][3] = { { 0.0 } };
  LvArray::tensorOps::Rij_add_AikAjk< 3, 3 >( rightCauchyDeformationTensor, m_deformationGradient[k] );
  LvArray::tensorOps::transpose< 3 >( rightCauchyDeformationTensor );

  real64 rotatedOldInverseTransposePlasticDeformationGradient[3][3] = { { 0.0 } };
  LvArray::tensorOps::Rij_eq_AikBkj< 3, 3, 3 >( rotatedOldInverseTransposePlasticDeformationGradient, endRotation, m_plasticDeformationGradient[k] );
  LvArray::tensorOps::invert< 3 >( rotatedOldInverseTransposePlasticDeformationGradient );
  LvArray::tensorOps::transpose< 3 >( rotatedOldInverseTransposePlasticDeformationGradient );

  real64 elasticRigthCauchyDeformationTensorDense[3][3] = { { 0.0 } };
  LvArray::tensorOps::Rij_eq_AikSymBklAjl< 3 >( elasticRigthCauchyDeformationTensorDense, rotatedOldInverseTransposePlasticDeformationGradient, rightCauchyDeformationTensor );

  real64 elasticRightCauchyDeformationTensor[6] = { 0.0 };
  LvArray::tensorOps::denseToSymmetric< 3 >( elasticRightCauchyDeformationTensor, elasticRightCauchyDeformationTensorDense );

  // Compute the logarithmic strain measure
  // Compute eigen vectors
  real64 eigenValues[3] = { 0.0 };
  real64 eigenVectors[3][3] = { { 0.0 } }; // Are these normalized? TODO check!
  LvArray::tensorOps::symEigenvectors< 3 >( eigenValues, eigenVectors, elasticRightCauchyDeformationTensor );

  real64 currentTransform[3][3] = { {0.0} };
  tensorOps::Rij_eq_AikBkj< 3, 3, 3 >( currentTransform, eigenVectors, endRotation );

  real64 oldElasticStrain[6] = { 0.0 };
  LvArray::tensorOps::Rij_eq_AikSymBklAjl< 3 >( oldElasticStrain, currentTransform, m_elasticStrain[k] );

  // Check if unrotated deformation gradient is diagonal
  real64 elasticStrainPredictor[6] = { 0.0 }
  for( int i=0; i < 3; ++i )
  {
    elasticStrainPredictor[i] = 0.5*LvArray::math::log( eigenValues[i] );
  }

  real64 E_inf = m_youngModuli[0];
  real64 E_e = E_inf;

  for( localIndex i = 0; i < numMaxwellElements; ++i )
  {
    E_e += m_youngModuli[i+1] * LvArray::math::exp( -0.5*dt/m_relaxationTimes[i] );
  }

  real64 K_inf = E_inf/(3 ((1-2*m_elasticPoissonRatio)));
  real64 G_inf = E_inf/(2*(1+m_elasticPoissonRatio));
  real64 K_e = E_e/(3 ((1-2*m_elasticPoissonRatio)));
  real64 G_e = E_e/(2*(1+m_elasticPoissonRatio));

  real64 volumetricStress = K_e*LvArray::tensorOps::symTrace< 3 >( elasticStrainPredictor ) - (K_e - K_inf)*LvArray::tensorOps::symTrace< 3 >( oldElasticStrain );
  for( localIndex n = 0; n < numMaxwellElements; ++n )
  {
    volumetricStress += m_Bin[k][n]*LvArray::math::exp( -dt/m_relaxationTimes[n] );
  }

  localIndex voigtMap[3][3] = { {0, 5, 4}, {5, 1, 3}, {4, 3, 2} };
  real64 trialStress[3][3] = { { 0.0 } };
  for( localIndex i=0; i < 3; ++i )
  {
    for( localIndex j=0; j < 3; ++j )
    {
      // Check summation using voigt notation is correct for off axis elements
      localIndex voigtIndex = voigtMap[i][j];
      trialStress[voigtIndex] += 2*G_e*elasticStrainPredictor[voigtIndex] - 2*(G_e-G_inf)*oldElasticStrain[voigtIndex] +  i == j ? volumetricStress : 0.0;
      for( localIndex n=0; n < numMaxwellElements; ++n )
      {
        trialStress[voigtIndex] += m_Ain[k][n][voigtIndex]*LvArray::math::exp( -dt/m_relaxationTimes[n] );
      }
    }
  }

  real64 reverseTransform[3][3] = { { 0.0 } };
  LvArray::tensorOps::transpose< 3, 3 >( reverseTransform, currentTransform );

  real64 newElasticStrain[6] = { 0.0 };
  LvArray::tensorOps::Rij_eq_AikSymBklAjl< 3 >( newElasticStrain, reverseTransform, elasticStrainPredictor );

  LvArray::tensorOps::copy< 6 >( m_elasticStrain[k], newElasticStrain );

  saveStress( k, q, stress );

  if( m_disableInelasticity )
  {
    return;
  }


  // save new stress and return
  saveStress( k, q, stress );
  return;
}


/**
 * @class HyperViscoElastoPlastic
 *
 * Strain hardening polymer material model.
 */
class HyperViscoElastoPlastic : public ElasticIsotropic
{
public:

  /// @typedef Alias for HyperViscoElastoPlasticUpdates
  using KernelWrapper = HyperViscoElastoPlasticUpdates;

  /**
   * constructor
   * @param[in] name name of the instance in the catalog
   * @param[in] parent the group which contains this instance
   */
  HyperViscoElastoPlastic( string const & name, Group * const parent );

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
  static constexpr auto m_catalogNameString = "HyperViscoElastoPlastic";

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

    // string/key for quadrature point elastic strain value
    static constexpr char const * elasticStrainString() { return "elasticStrain"; }

    // string/key for quadrature point plastic strain value
    static constexpr char const * plasticStrainString() { return "plasticStrain"; }

    /// string/key for quadrature point deviatoric stress terms
    static constexpr char const * AinString() { return "A_in"; }

    /// string/key for quadrature point deviatoric stress terms
    static constexpr char const * BinString() { return "B_in"; }

    /// string/key for quadrature point damage value
    static constexpr char const * damageString() { return "damage"; }

    /// string/key for Young's moduli
    static constexpr char const * youngModuliString() { return "youngModuli"; }

    /// string/key for elastic Poisson's ratio
    static constexpr char const * elasticPoissonRatioString() { return "elasticPoissonRatio"; }

    /// string/key for relaxation times
    static constexpr char const * relaxationTimesString() { return "relaxationTimes"; }
  };

  /**
   * @brief Create a instantiation of the HyperViscoElastoPlasticUpdates class that refers to the data in this.
   * @return An instantiation of HyperViscoElastoPlasticUpdates.
   */
  HyperViscoElastoPlasticUpdates createKernelUpdates() const
  {
    return HyperViscoElastoPlastic( m_deformationGradient,
                                    m_elasticStrain,
                                    m_plasticStrain,
                                    m_Ain,
                                    m_Bin,
                                    m_youngModuli,
                                    m_elasticPoissonRatio,
                                    m_relaxationTimes,
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
                          m_elasticStrain,
                          m_plasticStrain,
                          m_Ain,
                          m_Bin,
                          m_youngModuli,
                          m_elasticPoissonRatio,
                          m_relaxationTimes,
                          m_thermalExpansionCoefficient,
                          m_newStress,
                          m_oldStress,
                          m_density,
                          m_wavespeed,
                          m_disableInelasticity );
  }


protected:
  virtual void postInputInitialization() override;

  /// State variable: The deformation gradient values for each element/particle.
  array3d< real64 > m_deformationGradient;

  /// State variable: The elastic strain from the previous time step for numerical integration stored in voigt notation
  array3d< real64 > m_elasticStrain;

  /// State variable: The plastic strain values for each quadrature point stored in voigt notation
  array3d< real64 > m_plasticStrain;

  /// State variable: The plastic deformation gradient for each quadrature point
  array2d< real64 > m_plasticDeformationGradient;

  /// State variable: The deviatoric stress terms, A_in, stored in voigt notation
  array3d< real64 > m_Ain;

  /// State variable: The volumetric stress terms, B_in, stored in voigt notation
  array2d< real64 > m_Bin;

  /// Material parameter: The Young's moduli
  array1d< real64 > m_youngModuli;

  /// Material parameter: The elastic Poisson's ratio
  real64 m_elasticPoissonRatio;

  /// Material parameter: The relaxation time constants for elastic properties
  array1d< real64 > m_relaxationTimes;
};

} /* namespace constitutive */

} /* namespace geos */

#endif /* GEOSX_CONSTITUTIVE_SOLID_HYPERVISCOELASTOPLASTIC_HPP_ */
