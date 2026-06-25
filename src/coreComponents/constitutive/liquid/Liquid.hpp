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
 *  @file Liquid.hpp
 */

#ifndef GEOS_CONSTITUTIVE_LIQUID_LIQUID_HPP_
#define GEOS_CONSTITUTIVE_LIQUID_LIQUID_HPP_

#include "constitutive/ContinuumBase.hpp"
// #include "SolidModelDiscretizationOpsIsotropic.hpp"
// #include "constitutive/ExponentialRelation.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geos
{

namespace constitutive
{

/**
 * @class LiquidUpdates
 *
 * Class to provide elastic isotropic material updates that may be
 * called from a kernel function.
 */
class LiquidUpdates : public ContinuumBaseUpdates
{
public:
  /**
   * @brief Constructor
   * @param[in] bulkModulus          The ArrayView holding the bulkd modulus data for each element.
   * @param[in] viscosity            The ArrayView holding the viscosity data for each element.
   * @param[in] jacobian             The ArrayView holding the jacobian data for each element.
   * @param[in] newStress            The ArrayView holding the new stress data for each quadrature point.
   * @param[in] oldStress            The ArrayView holding the old stress data for each quadrature point.
   */
  LiquidUpdates( arrayView1d< real64 const > const & bulkModulus,
                 arrayView1d< real64 const > const & viscosity,
                 arrayView2d< real64 const > const & jacobian,
                 arrayView3d< real64 const > const & velocityGradient,
                 arrayView3d< real64, solid::STRESS_USD > const & newStress,
                 arrayView3d< real64, solid::STRESS_USD > const & oldStress,
                 arrayView2d< real64 > const & density,
                 arrayView2d< real64 > const & wavespeed ):
    ContinuumBaseUpdates( newStress,
                          oldStress,
                          density,
                          wavespeed ),
    m_bulkModulus( bulkModulus ),
    m_viscosity( viscosity ),
    m_jacobian( jacobian ),
    m_velocityGradient( velocityGradient )
  {}

  /// Deleted default constructor
  LiquidUpdates() = delete;

  /// Default copy constructor
  LiquidUpdates( LiquidUpdates const & ) = default;

  /// Default move constructor
  LiquidUpdates( LiquidUpdates && ) = default;

  /// Deleted copy assignment operator
  LiquidUpdates & operator=( LiquidUpdates const & ) = delete;

  /// Deleted move assignment operator
  LiquidUpdates & operator=( LiquidUpdates && ) =  delete;

  // /// Use the "isotropic" form of inner product compression
  // using DiscretizationOps = SolidModelDiscretizationOpsIsotropic;

  /// Use base version of saveConvergedState
  using ContinuumBaseUpdates::saveConvergedState;

  GEOS_HOST_DEVICE
  virtual void smallStrainNoStateUpdate_StressOnly( localIndex const k,
                                                    localIndex const q,
                                                    real64 const ( &totalStrain )[6],
                                                    real64 ( &stress )[6] ) const override;

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
  void smallStrainUpdate( localIndex const k,
                          localIndex const q,
                          real64 const & timeIncrement,
                          real64 const ( &strainIncrement )[6],
                          real64 ( &stress )[6],
                          real64 ( &stiffness )[6][6] ) const;

  // TODO: confirm hyper stress/strain measures before activatiing
  // CC: need hyperelastic model for hyperelasticMMS constitutive model in MPM
  GEOS_HOST_DEVICE
  virtual void hyperUpdate( localIndex const k,
                            localIndex const q,
                            real64 const ( &FminusI )[3][3],
                            real64 ( &stress )[6] ) const override final;

  GEOS_HOST_DEVICE
  virtual void hyperUpdate( localIndex const k,
                            localIndex const q,
                            real64 const ( &FminusI )[3][3],
                            real64 ( &stress )[6],
                            real64 ( &stiffness )[6][6] ) const override final;

protected:

  /// A reference to the ArrayView holding the bulk modulus for each element.
  arrayView1d< real64 const > const m_bulkModulus;

  /// A reference to the ArrayView holding the viscosity for each element.
  arrayView1d< real64 const > const m_viscosity;

  /// A reference to the ArrayView holding the jacobian for each element.
  arrayView2d< real64 const > const m_jacobian;

  /// A reference to the ArrayView holding the velocity gradient for each element.
  arrayView3d< real64 const > const m_velocityGradient;

};


GEOS_HOST_DEVICE
inline
void LiquidUpdates::smallStrainNoStateUpdate_StressOnly( localIndex const k,
                                                      localIndex const q,
                                                      real64 const ( &totalStrain )[6],
                                                      real64 ( & stress )[6] ) const
{
  GEOS_UNUSED_VAR( k );
  GEOS_UNUSED_VAR( q );
  GEOS_UNUSED_VAR( totalStrain );
  GEOS_UNUSED_VAR( stress );
  GEOS_ERROR( "smallStrainNoStateUpdate_StressOnly overload not implemented for Liquid" );
}


GEOS_HOST_DEVICE
inline
void LiquidUpdates::smallStrainUpdate_StressOnly( localIndex const k,
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
  GEOS_ERROR( "smallStrainUpdate_StressOnly overload not implemented for Liquid" );
}


GEOS_HOST_DEVICE
inline
void LiquidUpdates::smallStrainUpdate_StressOnly( localIndex const k,
                                                  localIndex const q,
                                                  real64 const & timeIncrement,
                                                  real64 const ( &beginningRotation )[3][3],
                                                  real64 const ( &endRotation )[3][3],
                                                  real64 const ( &strainIncrement )[6],
                                                  real64 ( & stress )[6] ) const
{
  GEOS_UNUSED_VAR( q );
  GEOS_UNUSED_VAR( endRotation );
  GEOS_UNUSED_VAR( timeIncrement );
  GEOS_UNUSED_VAR( strainIncrement );

  real64 beginningRotationTranspose[3][3] = { };
  LvArray::tensorOps::transpose< 3, 3 >( beginningRotationTranspose, beginningRotation );

  //  temp matrix for computing rotations.
  real64 tempMat[ 3 ][ 3 ]= { };
  
  // Unrotate velocity gradient:
  real64 unrotatedVelocityGradient[3][3] = { };
  LvArray::tensorOps::Rij_eq_AkiBkj< 3, 3, 3 >( tempMat, beginningRotation, m_velocityGradient[k] );
  LvArray::tensorOps::Rij_eq_AikBkj< 3, 3, 3 >( unrotatedVelocityGradient, tempMat, beginningRotation );
  real64 unrotatedVelocityGradientTranspose[3][3] = { };
  LvArray::tensorOps::transpose< 3, 3 >( unrotatedVelocityGradientTranspose, unrotatedVelocityGradient );

  // Symmetric part of unrotated velocity gradient, the strain "rate" for the step
  real64 denseD[3][3] = { };
  LvArray::tensorOps::copy< 3, 3 >( denseD, unrotatedVelocityGradient );
  LvArray::tensorOps::add< 3, 3 >( denseD, unrotatedVelocityGradientTranspose );
  LvArray::tensorOps::scale< 3, 3 >( denseD, 0.5 );
  real64 D[6] = { };
  LvArray::tensorOps::denseToSymmetric<3>( D, denseD );


  real64 devD[6] = { };
  LvArray::tensorOps::copy< 6 >( devD, D );
  real64 trD = LvArray::tensorOps::symTrace< 3 >( D ) / 3.0 ;
  LvArray::tensorOps::symAddIdentity< 3 >( devD, -trD );

  real64 K = m_bulkModulus[k];
  real64 mu = m_viscosity[k];

  real64 pressure =  LvArray::math::min( m_bulkModulus[k] * LvArray::math::log( m_jacobian[k][0] ), 0.0);

  stress[0] = -pressure + mu*devD[0];
  stress[1] = -pressure + mu*devD[1];
  stress[2] = -pressure + mu*devD[2];

  stress[3] = mu*devD[3];
  stress[4] = mu*devD[4];
  stress[5] = mu*devD[5];

  m_wavespeed[k][0] = sqrt( K / m_density[k][0] );
}

// CC: place filler for MPM hyperelasticMMS model
// There is a commented out implementation below from another person
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void LiquidUpdates::hyperUpdate( localIndex const k,
                              localIndex const q,
                              real64 const ( &FminusI )[3][3],
                              real64 ( & stress )[6] ) const
{
  GEOS_UNUSED_VAR( k );
  GEOS_UNUSED_VAR( q );
  GEOS_UNUSED_VAR( FminusI );
  GEOS_UNUSED_VAR( stress );
  GEOS_ERROR( "hyperUpdate() not implemented for this model" );
}

// CC: place filler for MPM hyperelasticMMS model
// There is a commented out implementation below from another person
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void LiquidUpdates::hyperUpdate( localIndex const k,
                              localIndex const q,
                              real64 const ( &FminusI )[3][3],
                              real64 ( & stress )[6],
                              real64 ( & stiffness )[6][6] ) const
{
  GEOS_UNUSED_VAR( k );
  GEOS_UNUSED_VAR( q );
  GEOS_UNUSED_VAR( FminusI );
  GEOS_UNUSED_VAR( stress );
  GEOS_UNUSED_VAR( stiffness );
  GEOS_ERROR( "hyperUpdate() not implemented for this model" );
}

/**
 * @class Liquid
 *
 * Class to provide an elastic isotropic material response.
 */
class Liquid : public ContinuumBase
{
public:

  /// Alias for LiquidUpdates
  using KernelWrapper = LiquidUpdates;

  /**
   * constructor
   * @param[in] name name of the instance in the catalog
   * @param[in] parent the group which contains this instance
   */
  Liquid( string const & name,
       Group * const parent );

  /**
   * Default Destructor
   */
  virtual ~Liquid() override;

  /**
   * @name Static Factory Catalog members and functions
   */
  ///@{

  /// string name to use for this class in the catalog
  static constexpr auto m_catalogNameString = "Liquid";

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
  struct viewKeyStruct : public ContinuumBase::viewKeyStruct
  {
    /// string/key for default bulk modulus
    static constexpr char const * defaultBulkModulusString() { return "defaultBulkModulus"; }

    /// string/key for bulk modulus
    static constexpr char const * bulkModulusString() { return "bulkModulus"; }

    /// string/key for default viscosity
    static constexpr char const * defaultViscosityString() { return "defaultViscosity"; }

    /// string/key for viscosity
    static constexpr char const * viscosityString() { return "viscosity"; }

    /// string/key for jacobian
    static constexpr char const * jacobianString() { return "jacobian"; }

    /// string/key for velocity gradient
    static constexpr char const * velocityGradientString() { return "velocityGradient"; }
  };

  /**
   * @brief Accessor for bulk modulus
   * @return A const reference to arrayView1d<real64> containing the bulk modulus (at every element).
   */
  arrayView1d< real64 > const bulkModulus() { return m_bulkModulus; }

  /**
   * @brief Const accessor for bulk modulus
   * @return A const reference to arrayView1d<real64 const> containing bulk modulus (at every element).
   */
  arrayView1d< real64 const > const bulkModulus() const { return m_bulkModulus; }

  /**
   * @brief Accessor for viscosity
   * @return A const reference to arrayView1d<real64> containing the viscosity (at every element).
   */
  arrayView1d< real64 > const viscosity() { return m_viscosity; }

  /**
   * @brief Const accessor for viscosity
   * @return A const reference to arrayView1d<real64 const> containing viscosity (at every element).
   */
  arrayView1d< real64 const > const viscosity() const { return m_viscosity; }

  /**
   * @brief Accessor for jacobian
   * @return A const reference to arrayView1d<real64> containing the jacobian (at every element).
   */
  arrayView2d< real64 > const jacobian() { return m_jacobian; }

  /**
   * @brief Const accessor for jacobian
   * @return A const reference to arrayView1d<real64 const> containing jacobian (at every element).
   */
  arrayView2d< real64 const > const jacobian() const { return m_jacobian; }

  /**
   * @brief Accessor for velocity gradient
   * @return A const reference to arrayView1d<real64> containing the velocity gradient (at every element).
   */
  arrayView3d< real64 > const velocityGradient() { return m_velocityGradient; }

  /**
   * @brief Const accessor for velocity gradient
   * @return A const reference to arrayView1d<real64 const> containing velocity gradient (at every element).
   */
  arrayView3d< real64 const > const velocityGradient() const { return m_velocityGradient; }

  GEOS_HOST_DEVICE
  virtual arrayView1d< real64 const > getBulkModulus() const final
  {
    return m_bulkModulus;
  }

  GEOS_HOST_DEVICE
  virtual arrayView1d< real64 const > getViscosity() const final
  {
    return m_viscosity;
  }

  GEOS_HOST_DEVICE
  virtual arrayView2d< real64 const > getJacobian() const final
  {
    return m_jacobian;
  }

  GEOS_HOST_DEVICE
  virtual arrayView3d< real64 const > getVelocityGradient() const final
  {
    return m_velocityGradient;
  }

  /**
   * @brief Create a instantiation of the LiquidUpdate class
   *        that refers to the data in this.
   * @param includeState Flag whether to pass state arrays that may not be needed for "no-state" updates
   * @return An instantiation of LiquidUpdate.
   */
  LiquidUpdates createKernelUpdates( bool const includeState = true ) const
  {
    if( includeState )
    {
      return LiquidUpdates( m_bulkModulus,
                            m_viscosity,
                            m_jacobian,
                            m_velocityGradient,
                            m_newStress,
                            m_oldStress,
                            m_density,
                            m_wavespeed );
    }
    else // for "no state" updates, pass empty views to avoid transfer of stress data to device
    {
      return LiquidUpdates( m_bulkModulus,
                            m_viscosity,
                            m_jacobian,
                            m_velocityGradient,
                            arrayView3d< real64, solid::STRESS_USD >(),
                            arrayView3d< real64, solid::STRESS_USD >(),
                            m_density,
                            m_wavespeed );
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
                          m_bulkModulus,
                          m_viscosity,
                          m_jacobian,
                          m_velocityGradient,
                          m_newStress,
                          m_oldStress,
                          m_wavespeed );
  }

  /**
   * @brief Allocate constitutive arrays
   * @param parent Object's parent group (element subregion)
   * @param numConstitutivePointsPerParentIndex Number of quadrature points per element
   */
  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

protected:

  /// Post-process XML data
  virtual void postInputInitialization() override;

  /// The default value of the bulk modulus
  real64 m_defaultBulkModulus;

  /// The bulk modulus for each upper level dimension (i.g. cell) of *this
  array1d< real64 > m_bulkModulus;

  /// The default value of the viscosity
  real64 m_defaultViscosity;

  /// The viscosity for each upper level dimension (i.e. cell) of *this
  array1d< real64 > m_viscosity;

  /// The jacobian for each upper level dimension (i.e. cell) of *this
  array2d< real64 > m_jacobian;

  // The velocity gradient for each upper level dimension (i.e. cell) of *this
  array3d< real64 > m_velocityGradient;

};

} /* namespace constitutive */

} /* namespace geos */

#endif /* GEOS_CONSTITUTIVE_LIQUID_LIQUID_HPP_ */
