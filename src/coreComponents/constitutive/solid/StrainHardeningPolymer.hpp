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
 * @brief Simple flow stress model for strain hardening polymer
 *
 * This damage model is intended for use with damage-field partitioning (DFG) within the
 * MPM solver, but can also be used without DFG by any solver. It is only appropriate for
 * schemes implementing explicit time integration. CC: TODO add description of material model
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_KINEMATICDAMAGE_HPP
#define GEOSX_CONSTITUTIVE_SOLID_KINEMATICDAMAGE_HPP

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
   * @param[in] damage 
   * @param[in] jacobian The ArrayView holding the jacobian for each quardrature point.
   * @param[in] yieldStrength The yield strength
   * @param[in] maximumStretch The maximum stretch
   * @param[in] thermalSoftening not currently implemented (CC: TODO)
   * @param[in] bulkModulus The ArrayView holding the bulk modulus data for each element.
   * @param[in] shearModulus The ArrayView holding the shear modulus data for each element.
   * @param[in] thermalExpansionCoefficient The ArrayView holding the thermal expansion coefficient data for each element.
   * @param[in] newStress The ArrayView holding the new stress data for each quadrature point.
   * @param[in] oldStress The ArrayView holding the old stress data for each quadrature point.
   */
  StrainHardeningPolymerUpdates( arrayView2d< real64 > const & plasticStrain,
                                 arrayView2d< real64 > const & damage,
                                 arrayView2d< real64 > const & jacobian,
                                 arrayView2d< real64 > const & yieldStrength,
                                 arrayView2d< real64 > const & maximumStretch,
                                 // arrayView2d< real64 > const & thermalSoftening,
                                 arrayView1d< real64 const > const & bulkModulus,
                                 arrayView1d< real64 const > const & shearModulus,
                                 arrayView1d< real64 const > const & thermalExpansionCoefficient,
                                 arrayView3d< real64, solid::STRESS_USD > const & newStress,
                                 arrayView3d< real64, solid::STRESS_USD > const & oldStress,
                                 bool const & disableInelasticity ):
    ElasticIsotropicUpdates( bulkModulus, shearModulus, thermalExpansionCoefficient, newStress, oldStress, disableInelasticity ),
    m_plasticStrain( plasticStrain ),
    m_damage( damage ),
    m_jacobian( jacobian ),
    m_yieldStrength( yieldStrength ),
    m_maximumStretch( maximumStretch )
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
  void smallStrainUpdateHelper( localIndex const k,
                                localIndex const q,
                                real64 const dt,
                                real64 const ( &strainIncrement )[6],
                                real64 ( &stress )[6] ) const;

  GEOS_HOST_DEVICE
  void computePlasticStrainIncrement ( localIndex const k,
                                       localIndex const q,
                                       const realT dt,
                                       real64 const ( &strainIncrement )[6],
                                       real64 const ( &stressIncrement )[6],
                                       arrayView1d< real64 > plastic_strain_increment ) const;

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  virtual void saveConvergedState( localIndex const k,
                                   localIndex const q ) const override final
  {
    ElasticIsotropicUpdates::saveConvergedState( k, q );
  }

private:
  /// A reference to the ArrayView holding the plastic strain for each quadrature point.
  arrayView3d< real64 > const m_plasticStrain;

  /// A reference to the ArrayView holding the damage for each quadrature point.
  arrayView2d< real64 > const m_damage;

  /// A reference to the ArrayView holding the jacobian for each quadrature point.
  arrayView2d< real64 > const m_jacobian;

  /// The yield strength
  real64 const m_yieldStrength;

  /// The compressive strength
  real64 const m_maximumStretch;

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
  // elastic predictor (assume strainIncrement is all elastic)
  ElasticIsotropicUpdates::smallStrainUpdate( k, q, timeIncrement, strainIncrement, stress, stiffness );
  m_jacobian[k][q] *= exp( strainIncrement[0] + strainIncrement[1] + strainIncrement[2] );

  if( m_disableInelasticity )
  {
    return;
  }

  // call the constitutive model
  StrainHardeningPolymerUpdates::smallStrainUpdateHelper( k, q, timeIncrement, strainIncrement, stress );

  // It doesn't make sense to modify stiffness with this model

  // save new stress and return
  saveStress( k, q, stress );
  return;
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
  smallStrainUpdate( k, q, timeIncrement, strainIncrement, stress, stiffness.m_c );
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void StrainHardeningPolymerUpdates::smallStrainUpdate_StressOnly( localIndex const k,
                                                                  localIndex const q,
                                                                  real64 const & timeIncrement,
                                                                  real64 const ( &strainIncrement )[6],
                                                                  real64 ( & stress )[6] ) const
{
  // elastic predictor (assume strainIncrement is all elastic)
  ElasticIsotropicUpdates::smallStrainUpdate_StressOnly( k, q, timeIncrement, strainIncrement, stress );
  m_jacobian[k][q] *= exp( strainIncrement[0] + strainIncrement[1] + strainIncrement[2] );

  if( m_disableInelasticity )
  {
    return;
  }

  // call the constitutive model
  StrainHardeningPolymerUpdates::smallStrainUpdateHelper( k, q, timeIncrement, strainIncrement, stress );

  // save new stress and return
  saveStress( k, q, stress );
  return;
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void StrainHardeningPolymerUpdates::smallStrainUpdateHelper( localIndex const k,
                                                             localIndex const q,
                                                             real64 const dt,
                                                             real64 const ( &strainIncrement )[6],
                                                             real64 ( & stress )[6] ) const
{
  // Store old stress for plastic strain increment
  real64 oldStress[6];
  for(int i = 0; i < 6; i++)
  {
    oldStress[i] = stress[i];
  }

  // decompose into mean (P) and von Mises (Q) stress invariants
  real64 trialP;
  real64 trialQ;
  real64 deviator[6];
  twoInvariant::stressDecomposition( stress,
                                     trialP,
                                     trialQ,
                                     deviator );

  // check yield function
  real64 yield = trialQ / m_yieldStrength[k];
  if( yield < 1.0 ) // elasticity
  {
    return;
  } 
  else
  {
    //Compute stretches by finding eigenvalues
    // CC: TODO what function do I use to compute eigenvalues

    array1d< real64 > stretch( 3 );
    BlasLapackLA::matrixEigenvalues( H, stretch );

    // Find the largest eigenvalues
    real64 maximumStretch = 0.0;
    for( localIndex i = 0; i < 3; ++i )
    {
        maximumStretch = std::max( stretch[i], maximumStretch );
    }

    if(maximumStretch > m_maximumStretch[k])
    {
        m_damage[k][q] = 1.0;
    }

    // magnitude of plastic strain tensor
    real64 gamma_p = 0.0;
    for( int i = 0; i < 6; i++ )
    {
      gamma_p += ( 1 + (i >= 3) ) * m_plasticStrain[k][i] * m_plasticStrain[k][i];
    }
    gamma_p = sqrt( gamma_p );
    
    // This term starts at value r0 and decays with plastic shear strain to give plastic softening.
    // Put in a check to prevent roundoff error.
    realT gamma_by_r1_to_r2 = std::pow( gamma_p / m_shearSofteningShapeParameter1[k], m_shearSofteningShapeParameter2[k] );

    // Compute change in yield strength
    realT plasticSoftening = m_shearSofteningMagnitude[k] * std::exp( std::max( -1.0 * gamma_by_r1_to_r2, -16.0 ) );
    realT stretchHardening = m_strainHardeningSlope[k] * ( maximumStretch * maximumStretch - 1.0 / maximumStretch );
    m_yieldStrength[k] += plasticSoftening + stretchHardening;

    // // CC: need to add this later
    //real64 thermalStrengthReduction = 1.0;
    // if(m_thermalSoftening)
    // {
    //   thermalStrengthReduction = computeThermalStrengthReduction();
    // }
    // m_yieldStrength[k] *= thermalStrengthReduction;

    // CC: is there a minimum floating point value constant defined somewhere?
    m_yieldStrength[k] *= 0.0000001 + (1.0 - m_damage[k][q]);

    // re-construct stress = P*eye + sqrt(2/3)*Q*nhat
    twoInvariant::stressRecomposition( trialP,
                                       m_yieldStrength[k],
                                       deviator,
                                       stress );

    // Increment plastic strain
    real64 stressIncrement[6];
    for(int i = 0; i < 6; i++ )
    {
      stressIncrement[i] = oldStress[i]-stress[i];
    }

    // increment plastic strain
    array1d< real64 > plasticStrainIncrement;
    computePlasticStrainIncrement( k,
                                   q,
                                   dt,           
                                   strainIncrement,
                                   stressIncrement,
                                   plasticStrainIncrement );

    for(int i = 0; i < 3; i++)
    {
      m_plasticStrain[k][q][i] += plasticStrainIncrement[i];
    }
  }
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void StrainHardeningPolymerUpdates::computePlasticStrainIncrement ( localIndex const k,
                                                                    localIndex const q,
                                                                    const realT dt,
                                                                    real64 const ( &strainIncrement )[6],
                                                                    real64 const ( &stressIncrement )[6],
                                                                    arrayView1d< real64 > plastic_strain_increment ) const
{ // For hypo-elastic models we compute the increment in plastic strain assuming
  // for some increment in total strain and stress and elastic properties.

  // Isotroptic-deviatoric decomposition;
  real64 trialP;
  real64 trialQ;
  real64 stressIncrementDeviator[6];
  twoInvariant::stressDecomposition( stressIncrement,
                                     trialP,
                                     trialQ,
                                     stressIncrementDeviator );

  real64 stressIncrementIsostatic[6] = {0};
  stressIncrementIsostatic[0] = trialP;
  stressIncrementIsostatic[1] = trialP;
  stressIncrementIsostatic[2] = trialP;

  // For damage or softening it there may be cases where bulk or shear are approx 0, so we need
  // to be careful that we compute this
  real64 elasticStrainIncrement[6] = {0};
  for(int i = 0; i < 6; i++)
  {
    if (m_bulkModulus[k] > 1.0e-12)
    {
      stressIncrementIsostatic[i] *= 1.0/3.0/m_bulkModulus[k]; // this is now elastic strain
      elasticStrainIncrement[i] += stressIncrementIsostatic[i];
    }
    if (m_shearModulus[k] > 1.0e-12)
    {
      stressIncrementDeviator[i] *= 1.0/2.0/m_shearModulus[k]; // this is now elastic strain
      elasticStrainIncrement[i] += stressIncrementDeviator[i];
    }

    // add plastic strain to plastic strain tensor ep += D*dt-eedot*dt
    plasticStrainIncrement[i] = strainIncrement[i];
    plasticStrainIncrement[i] -= elasticStrainIncrement[i];
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
    // string/key for quadrature point plastic strain value
    static constexpr char const * plasticStrainString() { return "plasticStrain"; }

    /// string/key for quadrature point damage value
    static constexpr char const * damageString() { return "damage"; }

    /// string/key for quadrature point jacobian value
    static constexpr char const * jacobianString() { return "jacobian"; }

    /// string/key for strain hardening slope
    static constexpr char const * strainHardeningSlopeString() { return "strainHardeningSlope"; }

    /// string/key for shear softening magnitude
    static constexpr char const * shearSofteningMagnitudeString() { return "shearSofteningMagnitude"; }

    /// string/key for shear softening shape parameter 1
    static constexpr char const * shearSofteningShapeParameter1String() { return "shearSofteningShapeParameter1"; }

    /// string/key for shear softening shape parameter 2
    static constexpr char const * shearSofteningShapeParameter1String() { return "shearSofteningShapeParameter2"; }

    /// string/key for yield strength
    static constexpr char const * yieldStrengthString() { return "yieldStrength"; }

    /// string/key for maximum stretch
    static constexpr char const * maximumStretchString() { return "maximumStretch"; }
  };

  /**
   * @brief Create a instantiation of the StrainHardeningPolymerUpdates class that refers to the data in this.
   * @return An instantiation of StrainHardeningPolymerUpdates.
   */
  StrainHardeningPolymerUpdates createKernelUpdates() const
  {
    return StrainHardeningPolymerUpdates( m_plasticStrain,
                                          m_damage,
                                          m_jacobian,
                                          m_strainHardeningSlope,
                                          m_shearSofteningMagnitude,
                                          m_shearSofteningShapeParameter1,
                                          m_shearSofteningShapeParameter2,
                                          m_yieldStrength,
                                          m_maximumStretch,
                                          m_bulkModulus,
                                          m_shearModulus,
                                          m_thermalExpansionCoefficient,
                                          m_newStress,
                                          m_oldStress,
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
                          m_plasticStrain,
                          m_damage,
                          m_jacobian,
                          m_strainHardeningSlope,
                          m_shearSofteningMagnitude,
                          m_shearSofteningShapeParameter1,
                          m_shearSofteningShapeParameter2,
                          m_yieldStrength,
                          m_maximumStretch,
                          m_bulkModulus,
                          m_shearModulus,
                          m_thermalExpansionCoefficient,
                          m_newStress,
                          m_oldStress,
                          m_disableInelasticity );
  }


protected:
  virtual void postProcessInput() override;

  // State variable: The strain hardening slope values for each quadrature point
  array2d< real64 > m_strainHardeningSlope;

  // State variable: The shear softening magnitude values for each quadrature point
  array2d< real64 > m_shearSofteningMagnitude;

  // State variable: The shear softening shape parameter 1 values for each quadrature point
  array2d< real64 > m_shearSofteningShapeParameter1;

    // State variable: The shear softening shape parameter 2 values for each quadrature point
  array2d< real64 > m_shearSofteningShapeParameter2;

  /// State variable: The plastic strain values for each quadrature point
  array3d< real64 > m_plasticStrain;

  /// State variable: The damage values for each quadrature point
  array2d< real64 > m_damage;

  /// State variable: The jacobian of the deformation
  array2d< real64 > m_jacobian;

  /// Material parameter: The value of unconfined tensile strength
  real64 m_yieldStrength;

  /// Material parameter: The value of maximum theoretical strength
  real64 m_maximumStretch;
};

} /* namespace constitutive */

} /* namespace geos */

#endif /* GEOSX_CONSTITUTIVE_SOLID_KINEMATICDAMAGE_HPP_ */
