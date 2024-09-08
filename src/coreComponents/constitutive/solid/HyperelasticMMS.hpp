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
 *  @file HyperelasticMMS.hpp
 */

#ifndef GEOS_CONSTITUTIVE_SOLID_HYPERELASTICMMS_HPP_
#define GEOS_CONSTITUTIVE_SOLID_HYPERELASTICMMS_HPP_

#include "SolidBase.hpp"
#include "PropertyConversions.hpp"
#include "SolidModelDiscretizationOpsIsotropic.hpp"
#include "constitutive/ExponentialRelation.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geos
{

namespace constitutive
{

/**
 * @class HyperelasticMMSUpdates
 *
 * Class to provide hyperelastic material updates that may be
 * called from a kernel function for method of manufactured solutions.
 */
class HyperelasticMMSUpdates : public SolidBaseUpdates
{
public:
  /**
   * @brief Constructor
   * @param[in] lambda  The ArrayView holding the first Lame constant data for each element.
   * @param[in] shearModulus The ArrayView holding the second Lame constant data for each element.
   * @param[in] thermalExpansionCoefficient The ArrayView holding the thermalExpansionCoefficient.
   * @param[in] newStress    The ArrayView holding the new stress data for each quadrature point.
   * @param[in] oldStress    The ArrayView holding the old stress data for each quadrature point.
   * @param[in] disableInelasticity Flag to disable plasticity for inelastic models
   */
  HyperelasticMMSUpdates( arrayView1d< real64 const > const & lambda,
                          arrayView1d< real64 const > const & shearModulus,
                          arrayView1d< real64 const > const & thermalExpansionCoefficient,
                          arrayView3d< real64, solid::STRESS_USD > const & newStress,
                          arrayView3d< real64, solid::STRESS_USD > const & oldStress,
                          arrayView2d< real64 > const & density,
                          arrayView2d< real64 > const & wavespeed,
                          const bool & disableInelasticity ):
    SolidBaseUpdates( newStress,
                      oldStress,
                      density,
                      wavespeed,
                      thermalExpansionCoefficient,
                      disableInelasticity ),
    m_lambda( lambda ),
    m_shearModulus( shearModulus )
  {}

  /// Deleted default constructor
  HyperelasticMMSUpdates() = delete;

  /// Default copy constructor
  HyperelasticMMSUpdates( HyperelasticMMSUpdates const & ) = default;

  /// Default move constructor
  HyperelasticMMSUpdates( HyperelasticMMSUpdates && ) = default;

  /// Deleted copy assignment operator
  HyperelasticMMSUpdates & operator=( HyperelasticMMSUpdates const & ) = delete;

  /// Deleted move assignment operator
  HyperelasticMMSUpdates & operator=( HyperelasticMMSUpdates && ) =  delete;

  /// Use the "isotropic" form of inner product compression
  using DiscretizationOps = SolidModelDiscretizationOpsIsotropic;

  /// Use base version of saveConvergedState
  using SolidBaseUpdates::saveConvergedState;

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
  void smallStrainNoStateUpdate( localIndex const k,
                                         localIndex const q,
                                         real64 const ( &totalStrain )[6],
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
                                             real64 const ( & beginningRotation )[3][3],
                                             real64 const ( & endRotation )[3][3],
                                             real64 const ( &strainIncrement )[6],
                                             real64 ( &stress )[6] ) const override;

  GEOS_HOST_DEVICE
  virtual void smallStrainUpdate( localIndex const k,
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
  virtual void getElasticStiffness( localIndex const k,
                                    localIndex const q,
                                    real64 ( &stiffness )[6][6] ) const override;

  GEOS_HOST_DEVICE
  virtual void getElasticStrain( localIndex const k,
                                 localIndex const q,
                                 real64 ( &elasticStrain )[6] ) const override;

  GEOS_HOST_DEVICE
  virtual real64 getLambda( localIndex const k ) const final
  {
    return m_lambda[k];
  }

  GEOS_HOST_DEVICE
  virtual real64 getShearModulus( localIndex const k ) const override final
  {
    return m_shearModulus[k];
  }

  GEOS_HOST_DEVICE
  virtual void viscousStateUpdate( localIndex const k,
                                   localIndex const q,
                                   real64 beta ) const override;

  GEOS_HOST_DEVICE
  virtual void hyperUpdate( localIndex const k,
                            localIndex const q,
                            real64 const ( & FminusI )[3][3],
                            real64 ( & stress )[6] ) const override;

  GEOS_HOST_DEVICE
  virtual void hyperUpdate( localIndex const k,
                            localIndex const q,
                            real64 const ( & FminusI )[3][3],
                            real64 ( & stress )[6],
                            real64 ( & stiffness )[6][6] ) const override;

protected:

  /// A reference to the ArrayView holding the first Lame constant for each element.
  arrayView1d< real64 const > const m_lambda;

  /// A reference to the ArrayView holding the shear modulus for each element.
  arrayView1d< real64 const > const m_shearModulus;

};


GEOS_HOST_DEVICE
inline
void HyperelasticMMSUpdates::getElasticStiffness( localIndex const k,
                                                   localIndex const q,
                                                   real64 ( & stiffness )[6][6] ) const
{
  GEOS_UNUSED_VAR( q );
  real64 const G = m_shearModulus[k];
  real64 const lambda = m_lambda[k];

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
void HyperelasticMMSUpdates::getElasticStrain( localIndex const k,
                                               localIndex const q,
                                               real64 ( & elasticStrain )[6] ) const
{
  real64 const E = conversions::lameConstants::toYoungMod( m_lambda[k], m_shearModulus[k] );
  real64 const nu = conversions::lameConstants::toPoissonRatio( m_lambda[k], m_shearModulus[k] );

  elasticStrain[0] = (    m_newStress[k][q][0] - nu*m_newStress[k][q][1] - nu*m_newStress[k][q][2])/E;
  elasticStrain[1] = (-nu*m_newStress[k][q][0] +    m_newStress[k][q][1] - nu*m_newStress[k][q][2])/E;
  elasticStrain[2] = (-nu*m_newStress[k][q][0] - nu*m_newStress[k][q][1] +    m_newStress[k][q][2])/E;

  elasticStrain[3] = m_newStress[k][q][3] / m_shearModulus[k];
  elasticStrain[4] = m_newStress[k][q][4] / m_shearModulus[k];
  elasticStrain[5] = m_newStress[k][q][5] / m_shearModulus[k];
}


GEOS_HOST_DEVICE
inline
void HyperelasticMMSUpdates::smallStrainNoStateUpdate_StressOnly( localIndex const k,
                                                                   localIndex const q,
                                                                   real64 const ( &totalStrain )[6],
                                                                   real64 ( & stress )[6] ) const
{
  GEOS_UNUSED_VAR( k );
  GEOS_UNUSED_VAR( q );
  GEOS_UNUSED_VAR( totalStrain );
  GEOS_UNUSED_VAR( stress );
}


GEOS_HOST_DEVICE
inline
void HyperelasticMMSUpdates::smallStrainNoStateUpdate( localIndex const k,
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
void HyperelasticMMSUpdates::smallStrainNoStateUpdate( localIndex const k,
                                                        localIndex const q,
                                                        real64 const ( &totalStrain )[6],
                                                        real64 ( & stress )[6],
                                                        DiscretizationOps & stiffness ) const
{
  smallStrainNoStateUpdate_StressOnly( k, q, totalStrain, stress );

  stiffness.m_bulkModulus = conversions::lameConstants::toBulkMod(m_lambda[k], m_shearModulus[k]);
  stiffness.m_shearModulus = m_shearModulus[k];
}


GEOS_HOST_DEVICE
inline
void HyperelasticMMSUpdates::smallStrainUpdate_StressOnly( localIndex const k,
                                                            localIndex const q,
                                                            real64 const & timeIncrement,
                                                            real64 const ( &strainIncrement )[6],
                                                            real64 ( & stress )[6] ) const
{
  GEOS_UNUSED_VAR( timeIncrement );
  smallStrainNoStateUpdate_StressOnly( k, q, strainIncrement, stress ); // stress  = incrementalStress
  LvArray::tensorOps::add< 6 >( stress, m_oldStress[k][q] );            // stress += m_oldStress
  saveStress( k, q, stress );                                           // m_newStress = stress
}

GEOS_HOST_DEVICE
inline 
void HyperelasticMMSUpdates::smallStrainUpdate_StressOnly( localIndex const k,
                                                           localIndex const q,
                                                           real64 const & timeIncrement,
                                                           real64 const ( & beginningRotation )[3][3],
                                                           real64 const ( & endRotation )[3][3],
                                                           real64 const ( & strainIncrement )[6],
                                                           real64 ( &stress )[6] ) const
{
  // CC: confirm this is correct
  GEOS_UNUSED_VAR( beginningRotation );
  GEOS_UNUSED_VAR( endRotation );
  smallStrainUpdate_StressOnly( k,
                                q,
                                timeIncrement,
                                strainIncrement,
                                stress );
}


GEOS_HOST_DEVICE
inline
void HyperelasticMMSUpdates::smallStrainUpdate( localIndex const k,
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
void HyperelasticMMSUpdates::smallStrainUpdate( localIndex const k,
                                                 localIndex const q,
                                                 real64 const & timeIncrement,
                                                 real64 const ( &strainIncrement )[6],
                                                 real64 ( & stress )[6],
                                                 DiscretizationOps & stiffness ) const
{
  smallStrainUpdate_StressOnly( k, q, timeIncrement, strainIncrement, stress );
  stiffness.m_bulkModulus = conversions::lameConstants::toBulkMod(m_lambda[k], m_shearModulus[k]);
  stiffness.m_shearModulus = m_shearModulus[k];
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void HyperelasticMMSUpdates::viscousStateUpdate( localIndex const k,
                                                  localIndex const q,
                                                  real64 beta ) const
{
  GEOS_UNUSED_VAR( k );
  GEOS_UNUSED_VAR( q );
  GEOS_UNUSED_VAR( beta );
}

// CC: hyperelastic update for model
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void HyperelasticMMSUpdates::hyperUpdate( localIndex const k,
                                          localIndex const q,
                                          real64 const ( & FminusI )[3][3],
                                          real64 ( & stress )[6] ) const
{
  GEOS_UNUSED_VAR( q );

  real64 lambda = m_lambda[k];
  real64 G = m_shearModulus[k];

  real64 F [3][3];
  LvArray::tensorOps::copy< 3, 3 >( F, FminusI );
  F[0][0] += 1;
  F[1][1] += 1;
  F[2][2] += 1;

  real64 C[3][3] = {{0}};
  for(int i = 0; i < 3; i++)
  {
    for(int j = 0; j < 3; j++)
    {
      for(int kk = 0; kk < 3; kk++)
      {
        C[i][j] += F[i][kk] * F[j][kk];
      }
    }
  }

  real64 J = F[0][0] * ( F[1][1] * F[2][2] - F[1][2] * F[2][1] ) -
             F[0][1] * ( F[1][0] * F[2][2] - F[1][2] * F[2][0] ) +
             F[0][2] * ( F[1][0] * F[2][1] - F[1][1] * F[2][0] );

  real64 const x1 = lambda * std::log(J) / J;
  real64 const x2 = G / J;

  stress[0] = x1 +  x2 * ( C[0][0] - 1 );
  stress[1] = x1 +  x2 * ( C[1][1] - 1 );
  stress[2] = x1 +  x2 * ( C[2][2] - 1 );

  stress[3] = x2 * C[1][2];
  stress[4] = x2 * C[0][2];
  stress[5] = x2 * C[0][1];

  m_wavespeed[k][0] = sqrt( (conversions::lameConstants::toBulkMod( m_lambda[k], m_shearModulus[k] ) + (4.0/3.0) * m_shearModulus[k] ) / m_density[k][0] );
}

// CC: hyperelastic update for model
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void HyperelasticMMSUpdates::hyperUpdate( localIndex const k,
                                          localIndex const q,
                                          real64 const ( & FminusI )[3][3],
                                          real64 ( & stress )[6],
                                          real64 ( & stiffness )[6][6] ) const
{
  hyperUpdate(k, q, FminusI, stress);
  getElasticStiffness( k, q, stiffness );
}


/**
 * @class HyperelasticMMS
 *
 * Class to provide an elastic isotropic material response.
 */
class HyperelasticMMS : public SolidBase
{
public:

  /// Alias for HyperelasticMMSUpdates
  using KernelWrapper = HyperelasticMMSUpdates;

  /**
   * constructor
   * @param[in] name name of the instance in the catalog
   * @param[in] parent the group which contains this instance
   */
  HyperelasticMMS( string const & name, Group * const parent );

  /**
   * Default Destructor
   */
  virtual ~HyperelasticMMS() override;

  /**
   * @name Static Factory Catalog members and functions
   */
  ///@{

  /// string name to use for this class in the catalog
  static constexpr auto m_catalogNameString = "HyperelasticMMS";

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
  struct viewKeyStruct : public SolidBase::viewKeyStruct
  {
    /// string/key for default bulk modulus
    static constexpr char const * defaultBulkModulusString() { return "defaultBulkModulus"; }

    /// string/key for default poisson ratio
    static constexpr char const * defaultPoissonRatioString() { return "defaultPoissonRatio"; }

    /// string/key for default shear modulus
    static constexpr char const * defaultShearModulusString() { return "defaultShearModulus"; }

    /// string/key for default Young's modulus
    static constexpr char const * defaultYoungModulusString() { return "defaultYoungModulus"; }

    /// string/key for default first Lame constant
    static constexpr char const * defaultLambdaString() { return "defaultLambda"; }

    /// string/key for first Lame constant
    static constexpr char const * lambdaString() { return "lambda"; }

    /// string/key for shear modulus
    static constexpr char const * shearModulusString() { return "shearModulus"; }

  };

  /**
   * @brief Accessor for first Lame constant
   * @return A const reference to arrayView1d<real64> containing the first Lame constant
   *         (at every element).
   */
  arrayView1d< real64 > const lambda() { return m_lambda; }

  /**
   * @brief Const accessor for first Lame constant
   * @return A const reference to arrayView1d<real64 const> containing the first Lame constant
   *         (at every element).
   */
  arrayView1d< real64 const > const lambda() const { return m_lambda; }

  /**
   * @brief Accessor for shear modulus
   * @return A const reference to arrayView1d<real64> containing the shear modulus (at every element).
   */
  arrayView1d< real64 > const shearModulus() { return m_shearModulus; }

  /**
   * @brief Const accessor for shear modulus
   * @return A const reference to arrayView1d<real64 const> containing the
   *         shear modulus (at every element).
   */
  arrayView1d< real64 const > const shearModulus() const { return m_shearModulus; }

  GEOS_HOST_DEVICE
  virtual arrayView1d< real64 const > getLambda() const final
  {
    return m_lambda;
  }
  GEOS_HOST_DEVICE
  virtual arrayView1d< real64 const > getShearModulus() const override final
  {
    return m_shearModulus;
  }

  /**
   * @brief Create a instantiation of the HyperelasticMMSUpdate class
   *        that refers to the data in this.
   * @param includeState Flag whether to pass state arrays that may not be needed for "no-state" updates
   * @return An instantiation of HyperelasticMMSUpdate.
   */
  HyperelasticMMSUpdates createKernelUpdates( bool const includeState = true ) const
  {
    if( includeState )
    {
      return HyperelasticMMSUpdates( m_lambda,
                                     m_shearModulus,
                                     m_thermalExpansionCoefficient,
                                     m_newStress,
                                     m_oldStress,
                                     m_density,
                                     m_wavespeed,
                                     m_disableInelasticity );
    }
    else // for "no state" updates, pass empty views to avoid transfer of stress data to device
    {
      return HyperelasticMMSUpdates( m_lambda,
                                     m_shearModulus,
                                     m_thermalExpansionCoefficient,
                                     arrayView3d< real64, solid::STRESS_USD >(),
                                     arrayView3d< real64, solid::STRESS_USD >(),
                                     m_density,
                                     m_wavespeed,
                                     m_disableInelasticity );
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
                          m_lambda,
                          m_shearModulus,
                          m_newStress,
                          m_oldStress,
                          m_density,
                          m_wavespeed,
                          m_disableInelasticity );
  }

protected:

  /// Post-process XML data
  virtual void postInputInitialization() override;

  /// The default value of the first Lame constant for any new allocations.
  real64 m_defaultLambda;

  /// The default value of the second Lame constant for any new allocations.
  real64 m_defaultShearModulus;

  /// The first lame constant for each upper level dimension (i.e. cell) of *this
  array1d< real64 > m_lambda;

  /// The shearModulus for each upper level dimension (i.e. cell) of *this
  array1d< real64 > m_shearModulus;

};

} /* namespace constitutive */

} /* namespace geos */

#endif /* GEOS_CONSTITUTIVE_SOLID_HYPERELASTICMMS_HPP_ */
