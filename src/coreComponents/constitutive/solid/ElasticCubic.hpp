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
 *  @file ElasticCubic.hpp
 */

#ifndef GEOS_CONSTITUTIVE_SOLID_ELASTICCUBIC_HPP_
#define GEOS_CONSTITUTIVE_SOLID_ELASTICCUBIC_HPP_

#include "SolidBase.hpp"
#include "constitutive/ExponentialRelation.hpp"
#include "LvArray/src/tensorOps.hpp"
#include "SolidModelDiscretizationOpsCubic.hpp"

namespace geos
{

namespace constitutive
{

/**
 * @class ElasticCubicUpdates
 *
 * Class to provide elastic transverse isotropic material updates that
 * may be called from a kernel function.
 */
class ElasticCubicUpdates : public SolidBaseUpdates
{
public:
  /**
   * @brief Constructor
   * @param[in] c11 The 11 component of the Voigt stiffness tensor.
   * @param[in] c12 The 12 component of the Voigt stiffness tensor.
   * @param[in] c44 The 44 component of the Voigt stiffness tensor.
   * @param[in] bulkModulus The effective bulk modulus for wavespeed calculations
   * @param[in] effectiveShearModulus The effective shear modulus for wavespeed calculations
   * @param[in] materailDirection The material direction for each point
   * @param[in] thermalExpansionCoefficient The ArrayView holding the thermal expansion coefficient data for each element.
   * @param[in] newStress The ArrayView holding the new stress data for each point.
   * @param[in] oldStress The ArrayView holding the old stress data for each point.
   * @param[in] disableInelasticity Flag to disable plastic response for inelastic models.
   */
  ElasticCubicUpdates( arrayView1d< real64 > const & c11,
                                     arrayView1d< real64 > const & c12,
                                     arrayView1d< real64 > const & c44,
                                     arrayView1d< real64 > const & bulkModulus,
                                     arrayView1d< real64 > const & effectiveShearModulus,
                                     arrayView3d< real64 > const & materialDirection,
                                     arrayView1d< real64 const > const & thermalExpansionCoefficient,
                                     arrayView3d< real64, solid::STRESS_USD > const & newStress,
                                     arrayView3d< real64, solid::STRESS_USD > const & oldStress,
                                     bool const & disableInelasticity ):
    SolidBaseUpdates( newStress, oldStress, thermalExpansionCoefficient, disableInelasticity ),
    m_c11( c11 ),
    m_c12( c12 ),
    m_c44( c44 ),
    m_bulkModulus( bulkModulus ),
    m_effectiveShearModulus( effectiveShearModulus ),
    m_materialDirection( materialDirection )
  {}

  /// Deleted default constructor
  ElasticCubicUpdates() = delete;

  /// Default copy constructor
  ElasticCubicUpdates( ElasticCubicUpdates const & ) = default;

  /// Default move constructor
  ElasticCubicUpdates( ElasticCubicUpdates && ) = default;

  /// Deleted copy assignment operator
  ElasticCubicUpdates & operator=( ElasticCubicUpdates const & ) = delete;

  /// Deleted move assignment operator
  ElasticCubicUpdates & operator=( ElasticCubicUpdates && ) =  delete;

  // Use transverse isotropic form of inner product compression
  using DiscretizationOps = SolidModelDiscretizationOpsCubic;

  /// Use base version of saveConvergedState
  using SolidBaseUpdates::saveConvergedState;

  // total strain interfaces

  GEOS_HOST_DEVICE
  virtual void smallStrainNoStateUpdate_StressOnly( localIndex const k,
                                                    localIndex const q,
                                                    real64 const ( &totalStrain )[6],
                                                    real64 ( &stress )[6] ) const override;

  GEOS_HOST_DEVICE
  virtual void smallStrainNoStateUpdate( localIndex const k,
                                         localIndex const q,
                                         real64 const ( &totalStrain )[6],
                                         real64 ( &stress )[6],
                                         real64 ( &stiffness )[6][6] ) const override;

  GEOS_HOST_DEVICE
  virtual void smallStrainNoStateUpdate( localIndex const k,
                                         localIndex const q,
                                         real64 const ( &totalStrain )[6],
                                         real64 ( &stress )[6],
                                         DiscretizationOps & stiffness ) const;

  // incremental strain interfaces

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
                                             real64 const ( & strainIncrement)[6],
                                             real64 ( & stress )[6]) const override;

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

  // miscellaneous getters

  GEOS_HOST_DEVICE
  virtual void getElasticStiffness( localIndex const k, localIndex const q, real64 ( &stiffness )[6][6] ) const override;

protected:

  /// A reference to the ArrayView holding c11 for each element.
  arrayView1d< real64 > const m_c11;

  /// A reference to the ArrayView holding c12 for each element.
  arrayView1d< real64 > const m_c12;

  /// A reference to the ArrayView holding c44 for each element.
  arrayView1d< real64 > const m_c44;

  /// A reference to the ArrayView holding the bulkModulus for each element
  arrayView1d< real64 > const m_bulkModulus;

  /// A reference to the ArrayView holding the effective shear modulus for each element
  arrayView1d< real64 > const m_effectiveShearModulus;

  // A reference to the ArrayView holding the material direction for each element/particle.
  arrayView3d< real64 > const m_materialDirection;
};

inline
GEOS_HOST_DEVICE
void ElasticCubicUpdates::getElasticStiffness( localIndex const k,
                                                localIndex const q,
                                                real64 ( & stiffness )[6][6] ) const
{
  GEOS_UNUSED_VAR( q );
  LvArray::tensorOps::fill< 6, 6 >( stiffness, 0 );

  stiffness[0][0] = m_c11[k];
  stiffness[0][1] = m_c12[k];
  stiffness[0][2] = m_c12[k];

  stiffness[1][0] = m_c12[k];
  stiffness[1][1] = m_c11[k];
  stiffness[1][2] = m_c12[k];

  stiffness[2][0] = m_c12[k];
  stiffness[2][1] = m_c12[k];
  stiffness[2][2] = m_c11[k];

  stiffness[3][3] = m_c44[k];
  stiffness[4][4] = m_c44[k];
  stiffness[5][5] = m_c44[k];
}

inline
GEOS_HOST_DEVICE
void ElasticCubicUpdates::smallStrainNoStateUpdate_StressOnly( localIndex const k,
                                                               localIndex const q,
                                                               real64 const ( &totalStrain )[6],
                                                               real64 ( & stress )[6] ) const
{
  GEOS_UNUSED_VAR( q );
  stress[0] = m_c11[k] * totalStrain[0] + m_c12[k] * totalStrain[1] + m_c12[k]*totalStrain[2];
  stress[1] = m_c12[k] * totalStrain[0] + m_c11[k] * totalStrain[1] + m_c12[k]*totalStrain[2];
  stress[2] = m_c12[k] * totalStrain[0] + m_c12[k] * totalStrain[1] + m_c11[k]*totalStrain[2];

  stress[3] = m_c44[k]*totalStrain[3];
  stress[4] = m_c44[k]*totalStrain[4];
  stress[5] = m_c44[k]*totalStrain[5];
}

inline
GEOS_HOST_DEVICE
void ElasticCubicUpdates::smallStrainNoStateUpdate( localIndex const k,
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
void ElasticCubicUpdates::smallStrainNoStateUpdate( localIndex const k,
                                                    localIndex const q,
                                                    real64 const ( &totalStrain )[6],
                                                    real64 ( & stress )[6],
                                                    DiscretizationOps & stiffness ) const
{
  smallStrainNoStateUpdate_StressOnly( k, q, totalStrain, stress );
  stiffness.m_c11 = m_c11[k];
  stiffness.m_c12 = m_c12[k];
  stiffness.m_c44 = m_c44[k];
}

inline
GEOS_HOST_DEVICE
void ElasticCubicUpdates::smallStrainUpdate_StressOnly( localIndex const k,
                                                        localIndex const q,
                                                        real64 const & timeIncrement,
                                                        real64 const ( &strainIncrement )[6],
                                                        real64 ( & stress )[6] ) const
{
  GEOS_UNUSED_VAR( timeIncrement );
  smallStrainNoStateUpdate_StressOnly( k, q, strainIncrement, stress ); // stress =  incrementalStress
  LvArray::tensorOps::add< 6 >( stress, m_oldStress[k][q] );            // stress += m_oldStress
  saveStress( k, q, stress );                                           // m_newStress = stress
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void ElasticCubicUpdates::smallStrainUpdate_StressOnly( localIndex const k,
                                                        localIndex const q,
                                                        real64 const & timeIncrement,
                                                        real64 const ( & beginningRotation )[3][3],
                                                        real64 const ( & endRotation )[3][3],
                                                        real64 const ( & strainIncrement )[6],
                                                        real64 ( & stress )[6] ) const
{
  GEOS_UNUSED_VAR( timeIncrement );
  GEOS_UNUSED_VAR( endRotation );

  // Unrotate the material direction first so transforming the stiffness matrix only occurs once

  // Get inverse of rotation matrix using transpose
  real64 beginningRotationTranspose[3][3] = { { 0 } };
  LvArray::tensorOps::transpose< 3, 3 >( beginningRotationTranspose, beginningRotation );
  
  // Normalize the material direction for safety
  real64 materialDirection[3][3] ={ { 0 } };
  LvArray::tensorOps::copy< 3, 3 >( materialDirection, m_materialDirection[k] ); // currently assume directions are already normalized
  // LvArray::tensorOps::normalize< 3 >( materialDirection );

  // Use beginning rotation to unrotate material direction
  real64 unrotatedMaterialDirection[3][3] = { { 0 } };
  LvArray::tensorOps::Rij_eq_AikBjk< 3, 3, 3 >( unrotatedMaterialDirection, beginningRotationTranspose, materialDirection );

  // // Compute the rotational axis between the z direction of the stiffness tensor and the material direction
  // real64 axis[3] = {0};
  // axis[2] = 1;
  // real64 rotationAxis[3] = { 0 };
  // LvArray::tensorOps::crossProduct( rotationAxis, axis, unrotatedMaterialDirection );

  // // Compute the rotation matrix to transform the stiffness tensor
  // real64 v[3][3] = { { 0 } };
  // v[0][1] = -rotationAxis[2];
  // v[0][2] = rotationAxis[1];
  // v[1][2] = -rotationAxis[0];
  // v[1][0] = rotationAxis[2];
  // v[2][0] = -rotationAxis[1];
  // v[2][1] = rotationAxis[0];

  // real64 c = LvArray::tensorOps::AiBi< 3 >( axis, unrotatedMaterialDirection );
  // real64 s = LvArray::tensorOps::l2Norm< 3 >( rotationAxis );

  real64 R[3][3] = { { 0 } };
  LvArray::tensorOps::copy< 3, 3 >( R, unrotatedMaterialDirection );
  // LvArray::tensorOps::addIdentity< 3 >( R, 1);
  // LvArray::tensorOps::add< 3, 3 >( R, v);

  // real64 temp[3][3] = { {0} };
  // LvArray::tensorOps::Rij_eq_AikBkj< 3, 3, 3 >( temp, v, v );
  // LvArray::tensorOps::scale< 3, 3 >( temp, ( 1 - c ) / ( s * s ));
  // LvArray::tensorOps::add< 3, 3 >( R, temp );

  // Get M to transform stiffness matrix in Voigt notation
  real64 M[6][6] = { {0} };
  M[0][0] = R[0][0] * R[0][0];
  M[0][1] = R[0][1] * R[0][1];
  M[0][2] = R[0][2] * R[0][2];
  M[0][3] = 2 * R[0][1] * R[0][2];
  M[0][4] = 2 * R[0][0] * R[0][2];
  M[0][5] = 2 * R[0][0] * R[0][1]; 

  M[1][0] = R[1][0] * R[1][0];
  M[1][1] = R[1][1] * R[1][1];
  M[1][2] = R[1][2] * R[1][2];
  M[1][3] = 2 * R[1][1] * R[1][2];
  M[1][4] = 2 * R[1][0] * R[1][2];
  M[1][5] = 2 * R[1][0] * R[1][1];

  M[2][0] = R[2][0] * R[2][0];
  M[2][1] = R[2][1] * R[2][1];
  M[2][2] = R[2][2] * R[2][2];
  M[2][3] = 2 * R[2][1] * R[2][2];
  M[2][4] = 2 * R[2][0] * R[2][2];
  M[2][5] = 2 * R[2][0] * R[2][1];

  M[3][0] = R[1][0] * R[2][0];
  M[3][1] = R[1][1] * R[2][1];
  M[3][2] = R[1][2] * R[2][2];
  M[3][3] = R[1][1] * R[2][2] + R[1][2] * R[2][1];
  M[3][4] = R[1][0] * R[2][2] + R[1][2] * R[2][0];
  M[3][5] = R[1][0] * R[2][1] + R[2][2] * R[2][0];

  M[4][0] = R[0][0] * R[2][0];
  M[4][1] = R[0][1] * R[2][1];
  M[4][2] = R[0][2] * R[2][2];
  M[4][3] = R[0][1] * R[2][2] + R[0][2] * R[2][1];
  M[4][4] = R[0][0] * R[2][2] + R[0][2] * R[2][0];
  M[4][5] = R[0][0] * R[2][1] + R[0][1] * R[2][0];

  M[5][0] = R[0][0] * R[1][0];
  M[5][1] = R[0][1] * R[1][1];
  M[5][2] = R[0][2] * R[1][2];
  M[5][3] = R[0][1] * R[1][2] + R[0][2] * R[1][1];
  M[5][4] = R[0][0] * R[1][2] + R[0][2] * R[1][0];
  M[5][5] = R[0][0] * R[1][1] + R[0][1] * R[1][0];

  // Build Voigt stiffness matrix
  real64 const c11 = m_c11[k];
  real64 const c12 = m_c12[k];
  real64 const c44 = m_c44[k];

  real64 stiffnessMatrix[6][6] = { { 0 } };
  stiffnessMatrix[0][0] = c11;
  stiffnessMatrix[0][1] = c12;
  stiffnessMatrix[0][2] = c12;

  stiffnessMatrix[1][0] = c12;
  stiffnessMatrix[1][1] = c11;
  stiffnessMatrix[1][2] = c12;

  stiffnessMatrix[2][0] = c12;
  stiffnessMatrix[2][1] = c12;
  stiffnessMatrix[2][2] = c11;

  stiffnessMatrix[3][3] = c44;
  stiffnessMatrix[4][4] = c44;
  stiffnessMatrix[5][5] = c44;

  real64 temp2[6][6] = { { 0 } };
  LvArray::tensorOps::Rij_eq_AikBkj< 6, 6, 6 >( temp2, M, stiffnessMatrix ); // M * S
  real64 rotatedStiffnessMatrix[6][6] = { { 0 } };
  LvArray::tensorOps::Rij_eq_AikBjk< 6, 6, 6 >( rotatedStiffnessMatrix, temp2, M ); // ( M * S ) * M^T

  LvArray::tensorOps::Ri_eq_AijBj< 6, 6 >( stress, rotatedStiffnessMatrix, strainIncrement );
  LvArray::tensorOps::add< 6 >( stress, m_oldStress[k][q] );
  saveStress( k, q, stress );  
}

inline
GEOS_HOST_DEVICE
void ElasticCubicUpdates::smallStrainUpdate( localIndex const k,
                                                           localIndex const q,
                                                           real64 const & timeIncrement,
                                                           real64 const ( &strainIncrement )[6],
                                                           real64 ( & stress )[6],
                                                           real64 ( & stiffness )[6][6] ) const
{
  smallStrainUpdate_StressOnly( k, q, timeIncrement, strainIncrement, stress );
  getElasticStiffness( k, q, stiffness );
}

GEOS_HOST_DEVICE
inline
void ElasticCubicUpdates::smallStrainUpdate( localIndex const k,
                                                           localIndex const q,
                                                           real64 const & timeIncrement,
                                                           real64 const ( &strainIncrement )[6],
                                                           real64 ( & stress )[6],
                                                           DiscretizationOps & stiffness ) const
{
  smallStrainUpdate_StressOnly( k, q, timeIncrement, strainIncrement, stress );
  stiffness.m_c11 = m_c11[k];
  stiffness.m_c12 = m_c12[k];
  stiffness.m_c44 = m_c44[k];
}

/**
 * @class ElasticCubic
 *
 * Class to provide a  elastic transverse isotropic material response.
 */
class ElasticCubic : public SolidBase
{
public:

  /// @typedef Alias for ElasticCubicUpdates
  using KernelWrapper = ElasticCubicUpdates;

  /**
   * @brief constructor
   * @param[in]name name of the instance in the catalog
   * @param[in]parent the group which contains this instance
   */
  ElasticCubic( string const & name, Group * const parent );

  /**
   * Destructor
   */
  virtual ~ElasticCubic() override;

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  /**
   * @name Static Factory Catalog members and functions
   */
  ///@{

  /// string name to use for this class in the catalog
  static constexpr auto m_catalogNameString = "ElasticCubic";

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
    /// string/key for default c11 component of Voigt stiffness tensor
    static constexpr char const * defaultC11String() { return "defaultC11"; };

    /// string/key for default c12 component of Voigt stiffness tensor
    static constexpr char const * defaultC12String() { return "defaultC12"; };

    /// string/key for default c44 component of Voigt stiffness tensor
    static constexpr char const * defaultC44String() { return "defaultC44"; };

    /// string/key for c11 component of Voigt stiffness tensor
    static constexpr char const * c11String() { return "c11"; }

    /// string/key for c12 component of Voigt stiffness tensor
    static constexpr char const * c12String() { return "c12"; }

    /// string/key for c44 component of Voigt stiffness tensor
    static constexpr char const * c44String() { return "c44"; }

    /// string/key for effective bulk modulus 
    static constexpr char const * bulkModulusString() { return "bulkModulus"; }

    /// string/key for effective shear modulus 
    static constexpr char const * effectiveShearModulusString() { return "effectiveShearModulus"; }

    /// string/key for material direction value
    static constexpr char const * materialDirectionString() { return "materialDirection"; }
  };

  /**
   * @brief Const-Getter for 11 component of Voigt stiffness tensor.
   * @return reference to immutable 11 component of Voigt stiffness tensor.
   */
  arrayView1d< real64 const > getC11() const { return m_c11; }

  /**
   * @brief Getter for 11 component of Voigt stiffness tensor.
   * @return reference to mutable 11 component of Voigt stiffness tensor.
   */
  arrayView1d< real64 > getC11() { return m_c11; }

  /**
   * @brief Const-Getter for 12 component of Voigt stiffness tensor.
   * @return reference to immutable 12 component of Voigt stiffness tensor.
   */
  arrayView1d< real64 const > getC12() const { return m_c12; }

  /**
   * @brief Getter for 12 component of Voigt stiffness tensor.
   * @return reference to mutable 12 component of Voigt stiffness tensor.
   */
  arrayView1d< real64 > getC12() { return m_c12; }

  /**
   * @brief Const-Getter for 44 component of Voigt stiffness tensor.
   * @return reference to immutable 44 component of Voigt stiffness tensor.
   */
  arrayView1d< real64 const > getC44() const { return m_c44; }

  /**
   * @brief Getter for 44 component of Voigt stiffness tensor.
   * @return reference to mutable 44 component of Voigt stiffness tensor.
   */
  arrayView1d< real64 > getC44() { return m_c44; }

  /**
   * @brief Accessor for bulk modulus
   * @return A const reference to arrayView1d<real64> containing the effective bulk
   *         modulus (at every element).
   */
  arrayView1d< real64 > const bulkModulus() { return m_bulkModulus; }

  /**
   * @brief Const accessor for effective bulk modulus
   * @return A const reference to arrayView1d<real64 const> containing the
   *         effective bulk modulus (at every element).
   */
  arrayView1d< real64 const > const bulkModulus() const { return m_bulkModulus; }

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
  arrayView1d< real64 const > getBulkModulus() const override final { return m_bulkModulus; }

  /**
  * @brief Getter for effective shear modulus.
  * @return reference to mutable effective shear modulus.
  */
  GEOS_HOST_DEVICE
  arrayView1d< real64 const > getEffectiveShearModulus() const { return m_effectiveShearModulus; }


  /**
   * @brief Create a instantiation of the
   *        ElasticCubicUpdates class that refers to the
   *        data in this.
   * @return An instantiation of ElasticCubicUpdates.
   */
  ElasticCubicUpdates createKernelUpdates() const
  {
    return ElasticCubicUpdates( m_c11,
                                m_c12,
                                m_c44,
                                m_bulkModulus,
                                m_effectiveShearModulus,
                                m_materialDirection,
                                m_thermalExpansionCoefficient,
                                m_newStress,
                                m_oldStress,
                                m_disableInelasticity );
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
  UPDATE_KERNEL createDerivedKernelUpdates( PARAMS && ... constructorParams )
  {
    return UPDATE_KERNEL( std::forward< PARAMS >( constructorParams )...,
                          m_c11,
                          m_c12,
                          m_c44,
                          m_bulkModulus,
                          m_effectiveShearModulus,
                          m_materialDirection,
                          m_thermalExpansionCoefficient,
                          m_newStress,
                          m_oldStress,
                          m_disableInelasticity );
  }

protected:
  virtual void postProcessInput() override;

  /// The default value of the C11 for any new allocations.
  real64 m_defaultC11;

  /// The default value of the C12 for any new allocations.
  real64 m_defaultC12;

  /// The default value of the C44 for any new allocations.
  real64 m_defaultC44;

  /// The 11 component of the Voigt stiffness tensor.
  array1d< real64 > m_c11;

  /// The 12 component of the Voigt stiffness tensor.
  array1d< real64 > m_c12;

  /// The 44 component of the Voigt stiffness tensor.
  array1d< real64 > m_c44;

  /// The bulk modulus.
  array1d< real64 > m_bulkModulus;

  /// The effective shear modulus
  array1d< real64 > m_effectiveShearModulus;
  
  /// State variable: The material direction for each element/particle
  array3d< real64 > m_materialDirection;
};

} /* namespace constitutive */

} /* namespace geos */

#endif /* GEOS_CONSTITUTIVE_SOLID_ELASTICCUBIC_HPP_ */
