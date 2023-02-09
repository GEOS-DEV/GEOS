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
 * @file CeramicDamage.hpp
 * @brief Simple damage model for modeling material failure in brittle materials.
 * 
 * This damage model is intended for use with damage-field partitioning (DFG) within the
 * MPM solver. We don't use the damage to modify the stress. Rather, it flags the
 * creation of fracture surfaces which can ceramically separate due to DFG.
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_KINEMATICDAMAGE_HPP
#define GEOSX_CONSTITUTIVE_SOLID_KINEMATICDAMAGE_HPP

#include "ElasticIsotropic.hpp"
#include "InvariantDecompositions.hpp"
#include "PropertyConversions.hpp"
#include "SolidModelDiscretizationOpsFullyAnisotroipic.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geosx
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
   * @param[in] bulkModulus The ArrayView holding the bulk modulus data for each element.
   * @param[in] shearModulus The ArrayView holding the shear modulus data for each element.
   * @param[in] newStress The ArrayView holding the new stress data for each quadrature point.
   * @param[in] oldStress The ArrayView holding the old stress data for each quadrature point.
   */
  CeramicDamageUpdates( arrayView2d< real64 > const & damage,
                          real64 const & defaultFailureStress,
                          arrayView1d< real64 const > const & bulkModulus,
                          arrayView1d< real64 const > const & shearModulus,
                          arrayView3d< real64, solid::STRESS_USD > const & newStress,
                          arrayView3d< real64, solid::STRESS_USD > const & oldStress,
                          bool const & disableInelasticity ):
    ElasticIsotropicUpdates( bulkModulus, shearModulus, newStress, oldStress, disableInelasticity ),
    m_damage( damage ),
    m_defaultFailureStress( defaultFailureStress )
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

  GEOSX_HOST_DEVICE
  virtual void smallStrainUpdate( localIndex const k,
                                  localIndex const q,
                                  real64 const ( &strainIncrement )[6],
                                  real64 ( &stress )[6],
                                  real64 ( &stiffness )[6][6] ) const override final;

  GEOSX_HOST_DEVICE
  virtual void smallStrainUpdate( localIndex const k,
                                  localIndex const q,
                                  real64 const ( &strainIncrement )[6],
                                  real64 ( &stress )[6],
                                  DiscretizationOps & stiffness ) const final;

  GEOSX_HOST_DEVICE
  virtual void smallStrainUpdate_StressOnly( localIndex const k,
                                             localIndex const q,
                                             real64 const ( &strainIncrement )[6],
                                             real64 ( &stress )[6] ) const override final;

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  virtual void saveConvergedState( localIndex const k,
                                   localIndex const q ) const override final
  {
    ElasticIsotropicUpdates::saveConvergedState( k, q );
  }

private:
  /// A reference to the ArrayView holding the damage for each quadrature point.
  arrayView2d< real64 > const m_damage;

  /// The failure stress
  real64 const m_defaultFailureStress;
};


GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void CeramicDamageUpdates::smallStrainUpdate( localIndex const k,
                                                localIndex const q,
                                                real64 const ( &strainIncrement )[6],
                                                real64 ( & stress )[6],
                                                real64 ( & stiffness )[6][6] ) const
{
  // elastic predictor (assume strainIncrement is all elastic)
  ElasticIsotropicUpdates::smallStrainUpdate( k, q, strainIncrement, stress, stiffness );

  if( m_disableInelasticity )
  {
    return;
  }

  // get principal stresses
  real64 eigenValues[3] = {};
  real64 eigenVectors[3][3] = {};
  LvArray::tensorOps::symEigenvectors< 3 >( eigenValues, eigenVectors, stress );

  // check for damage
  real64 maxPrincipalStress = fmax( eigenValues[0], fmax( eigenValues[1], eigenValues[2] ) ); // TODO: Check if eigenvalues are returned in order
  if( maxPrincipalStress > m_defaultFailureStress )
  {
    m_damage[k][q] = 1.0;
  }

  // It doesn't make sense to modify stiffness with this model

  // save new stress and return
  saveStress( k, q, stress );
  return;
}


GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void CeramicDamageUpdates::smallStrainUpdate( localIndex const k,
                                                localIndex const q,
                                                real64 const ( &strainIncrement )[6],
                                                real64 ( & stress )[6],
                                                DiscretizationOps & stiffness ) const
{
  smallStrainUpdate( k, q, strainIncrement, stress, stiffness.m_c );
}


GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE // TODO: Is there a way to not have to re-write the constitutive model twice for regular and StressOnly?
void CeramicDamageUpdates::smallStrainUpdate_StressOnly( localIndex const k,
                                                           localIndex const q,
                                                           real64 const ( &strainIncrement )[6],
                                                           real64 ( & stress )[6] ) const
{
  // elastic predictor (assume strainIncrement is all elastic)
  ElasticIsotropicUpdates::smallStrainUpdate_StressOnly( k, q, strainIncrement, stress );

  if( m_disableInelasticity )
  {
    return;
  }

  // get principal stresses
  real64 eigenValues[3] = {};
  real64 eigenVectors[3][3] = {};
  LvArray::tensorOps::symEigenvectors< 3 >( eigenValues, eigenVectors, stress );

  // check for damage
  real64 maxPrincipalStress = fmax( eigenValues[0], fmax( eigenValues[1], eigenValues[2] ) ); // TODO: Check if eigenvalues are returned in order
  if( maxPrincipalStress > m_defaultFailureStress )
  {
    m_damage[k][q] = 1.0;
  }

  // TODO: construct consistent tangent stiffness? Not possible for spectral damage.
  // We don't modify the stiffness for now...

  // save new stress and return
  saveStress( k, q, stress );
  return;
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

  arrayView2d< real64 const > getDamage() const { return m_damage; }

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
    /// string/key for failure stress
    static constexpr char const * defaultFailureStressString() { return "defaultFailureStress"; }

    /// string/key for quadrature point damage value
    static constexpr char const * damageString() { return "damage"; }
  };

  /**
   * @brief Create a instantiation of the CeramicDamageUpdate class that refers to the data in this.
   * @return An instantiation of CeramicDamageUpdate.
   */
  CeramicDamageUpdates createKernelUpdates() const
  {
    return CeramicDamageUpdates( m_damage,
                                   m_defaultFailureStress,
                                   m_bulkModulus,
                                   m_shearModulus,
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
                          m_damage,
                          m_defaultFailureStress,
                          m_bulkModulus,
                          m_shearModulus,
                          m_newStress,
                          m_oldStress,
                          m_disableInelasticity );
  }


protected:
  virtual void postProcessInput() override;

  /// State variable: The damage values for each quadrature point
  array2d< real64 > m_damage;

  /// Material parameter: The default value of failure stress
  real64 m_defaultFailureStress;
};

} /* namespace constitutive */

} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_SOLID_KINEMATICDAMAGE_HPP_ */
