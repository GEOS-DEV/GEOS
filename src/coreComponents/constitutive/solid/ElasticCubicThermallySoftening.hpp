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
 *  @file ElasticCubicThermallySoftening.hpp
 */

#ifndef GEOS_CONSTITUTIVE_SOLID_ELASTICCUBICTHERMALLYSOFTENING_HPP_
#define GEOS_CONSTITUTIVE_SOLID_ELASTICCUBICTHERMALLYSOFTENING_HPP_

#include "ElasticCubic.hpp"
#include "constitutive/ExponentialRelation.hpp"
#include "LvArray/src/tensorOps.hpp"
#include "SolidModelDiscretizationOpsCubic.hpp"

namespace geos
{

namespace constitutive
{

/**
 * @class ElasticCubicThermallySofteningUpdates
 *
 * Class to provide elastic transverse isotropic material updates that
 * may be called from a kernel function.
 */
class ElasticCubicThermallySofteningUpdates : public ElasticCubicUpdates
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
  ElasticCubicThermallySofteningUpdates( real64 const & referenceC11,
                                         real64 const & referenceC12,
                                         real64 const & referenceC44,
                                         real64 const & firstOrderC11ThermalCoefficient,
                                         real64 const & firstOrderC12ThermalCoefficient,
                                         real64 const & firstOrderC44ThermalCoefficient,
                                         real64 const & secondOrderC11ThermalCoefficient,
                                         real64 const & secondOrderC12ThermalCoefficient,
                                         real64 const & secondOrderC44ThermalCoefficient,
                                         real64 const & referenceTemperature,
                                         arrayView1d< real64 > const & temperature,   
                                         arrayView1d< real64 > const & c11,
                                         arrayView1d< real64 > const & c12,
                                         arrayView1d< real64 > const & c44,
                                         arrayView1d< real64 > const & bulkModulus,
                                         arrayView1d< real64 > const & effectiveShearModulus,
                                         arrayView3d< real64 > const & materialDirection,
                                         arrayView1d< real64 const > const & thermalExpansionCoefficient,
                                         arrayView3d< real64, solid::STRESS_USD > const & newStress,
                                         arrayView3d< real64, solid::STRESS_USD > const & oldStress,
                                         bool const & disableInelasticity ):
    ElasticCubicUpdates( c11, c12, c44, bulkModulus, effectiveShearModulus, materialDirection, thermalExpansionCoefficient, newStress, oldStress, disableInelasticity ),
    m_referenceC11( referenceC11 ),
    m_referenceC12( referenceC12 ),
    m_referenceC44( referenceC44 ),
    m_firstOrderC11ThermalCoefficient( firstOrderC11ThermalCoefficient ),
    m_firstOrderC12ThermalCoefficient( firstOrderC12ThermalCoefficient ),
    m_firstOrderC44ThermalCoefficient( firstOrderC44ThermalCoefficient ),
    m_secondOrderC11ThermalCoefficient( secondOrderC11ThermalCoefficient ),
    m_secondOrderC12ThermalCoefficient( secondOrderC12ThermalCoefficient ),
    m_secondOrderC44ThermalCoefficient( secondOrderC44ThermalCoefficient ),
    m_referenceTemperature( referenceTemperature ),
    m_temperature( temperature )
  {}

  /// Deleted default constructor
  ElasticCubicThermallySofteningUpdates() = delete;

  /// Default copy constructor
  ElasticCubicThermallySofteningUpdates( ElasticCubicThermallySofteningUpdates const & ) = default;

  /// Default move constructor
  ElasticCubicThermallySofteningUpdates( ElasticCubicThermallySofteningUpdates && ) = default;

  /// Deleted copy assignment operator
  ElasticCubicThermallySofteningUpdates & operator=( ElasticCubicThermallySofteningUpdates const & ) = delete;

  /// Deleted move assignment operator
  ElasticCubicThermallySofteningUpdates & operator=( ElasticCubicThermallySofteningUpdates && ) =  delete;

  // Use transverse isotropic form of inner product compression
  using DiscretizationOps = SolidModelDiscretizationOpsCubic;

  /// Use base version of saveConvergedState
  using ElasticCubicUpdates::saveConvergedState;

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
                                         DiscretizationOps & stiffness ) const override;

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
                                  DiscretizationOps & stiffness ) const override;

  // miscellaneous getters

  GEOS_HOST_DEVICE
  virtual void getElasticStiffness( localIndex const k, localIndex const q, real64 ( &stiffness )[6][6] ) const override;

protected:


  /// A reference value for the c11.
  real64 const m_referenceC11;

  /// A reference value for the c12.
  real64 const m_referenceC12;

  /// A reference value for the c44.
  real64 const m_referenceC44;

  /// A value for the first order thermal coefficient for c11.
  real64 const m_firstOrderC11ThermalCoefficient;

  /// A value for the first order thermal coefficient for c12.
  real64 const m_firstOrderC12ThermalCoefficient;
  
  /// A value for the first order thermal coefficient for c44.
  real64 const m_firstOrderC44ThermalCoefficient;

  /// A value for the second order thermal coefficient for c11.
  real64 const m_secondOrderC11ThermalCoefficient;

  /// A value for the second order thermal coefficient for c12.
  real64 const m_secondOrderC12ThermalCoefficient;
  
  /// A value for the second order thermal coefficient for c44.
  real64 const m_secondOrderC44ThermalCoefficient;

  /// A value for the reference temperature.
  real64 const m_referenceTemperature;

  /// A reference to the ArrayView holding the temperature for each element
  arrayView1d< real64 > const m_temperature;
};

inline
GEOS_HOST_DEVICE
void ElasticCubicThermallySofteningUpdates::getElasticStiffness( localIndex const k,
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
void ElasticCubicThermallySofteningUpdates::smallStrainNoStateUpdate_StressOnly( localIndex const k,
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
void ElasticCubicThermallySofteningUpdates::smallStrainNoStateUpdate( localIndex const k,
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
void ElasticCubicThermallySofteningUpdates::smallStrainNoStateUpdate( localIndex const k,
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
void ElasticCubicThermallySofteningUpdates::smallStrainUpdate_StressOnly( localIndex const k,
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
void ElasticCubicThermallySofteningUpdates::smallStrainUpdate_StressOnly( localIndex const k,
                                                        localIndex const q,
                                                        real64 const & timeIncrement,
                                                        real64 const ( & beginningRotation )[3][3],
                                                        real64 const ( & endRotation )[3][3],
                                                        real64 const ( & strainIncrement )[6],
                                                        real64 ( & stress )[6] ) const
{
  real64 dT = m_temperature[k] - m_referenceTemperature;

  m_c11[k] = m_referenceC11 * ( 1.0 + m_firstOrderC11ThermalCoefficient * dT  + m_secondOrderC11ThermalCoefficient * dT * dT );
  m_c12[k] = m_referenceC12 * ( 1.0 + m_firstOrderC12ThermalCoefficient * dT  + m_secondOrderC12ThermalCoefficient * dT * dT );
  m_c44[k] = m_referenceC44 * ( 1.0 + m_firstOrderC44ThermalCoefficient * dT  + m_secondOrderC44ThermalCoefficient * dT * dT );

  ElasticCubicUpdates::smallStrainUpdate_StressOnly( k, 
                                                     q, 
                                                     timeIncrement, 
                                                     beginningRotation, 
                                                     endRotation, 
                                                     strainIncrement,
                                                     stress );
}

inline
GEOS_HOST_DEVICE
void ElasticCubicThermallySofteningUpdates::smallStrainUpdate( localIndex const k,
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
void ElasticCubicThermallySofteningUpdates::smallStrainUpdate( localIndex const k,
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
 * @class ElasticCubicThermallySoftening
 *
 * Class to provide a  elastic transverse isotropic material response.
 */
class ElasticCubicThermallySoftening : public ElasticCubic
{
public:

  /// @typedef Alias for ElasticCubicThermallySofteningUpdates
  using KernelWrapper = ElasticCubicThermallySofteningUpdates;

  /**
   * @brief constructor
   * @param[in]name name of the instance in the catalog
   * @param[in]parent the group which contains this instance
   */
  ElasticCubicThermallySoftening( string const & name, Group * const parent );

  /**
   * Destructor
   */
  virtual ~ElasticCubicThermallySoftening() override;

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  /**
   * @name Static Factory Catalog members and functions
   */
  ///@{

  /// string name to use for this class in the catalog
  static constexpr auto m_catalogNameString = "ElasticCubicThermallySoftening";

  /**
   * @return A string that is used to register/lookup this class in the registry
   */
  static string catalogName() { return m_catalogNameString; }

  virtual string getCatalogName() const override { return catalogName(); }
  ///@}

  /**
   * Keys for data specified in this class.
   */
  struct viewKeyStruct : public ElasticCubic::viewKeyStruct
  {

    /// string/key for the first order thermal coefficient for the c11 component of Voigt stiffness tensor
    static constexpr char const * firstOrderC11ThermalCoefficientString() { return "firstOrderC11ThermalCoefficient"; };

    /// string/key for the first order thermal coefficient for the c12 component of Voigt stiffness tensor
    static constexpr char const * firstOrderC12ThermalCoefficientString() { return "firstOrderC12ThermalCoefficient"; };

    /// string/key for the first order thermal coefficient for the c44 component of Voigt stiffness tensor
    static constexpr char const * firstOrderC44ThermalCoefficientString() { return "firstOrderC44ThermalCoefficient"; };

    /// string/key for the second order thermal coefficient for the c11 component of Voigt stiffness tensor
    static constexpr char const * secondOrderC11ThermalCoefficientString() { return "secondOrderC11ThermalCoefficient"; };

    /// string/key for the second order thermal coefficient for the c12 component of Voigt stiffness tensor
    static constexpr char const * secondOrderC12ThermalCoefficientString() { return "secondOrderC12ThermalCoefficient"; };

    /// string/key for the second order thermal coefficient for the c44 component of Voigt stiffness tensor
    static constexpr char const * secondOrderC44ThermalCoefficientString() { return "secondOrderC44ThermalCoefficient"; };
  
    /// string/key for the reference temperature
    static constexpr char const * referenceTemperatureString() { return "referenceTemperature"; };

    /// string/key for the temperature
    static constexpr char const * temperatureString() { return "temperature"; };

  };

  /**
   * @brief Const-Getter for element temperature
   * @return reference to immutable temperature.
   */
  arrayView1d< real64 const > getTemperature() const { return m_temperature; }

  /**
   * @brief Getter for temperature.
   * @return reference to mutable temperature.
   */
  arrayView1d< real64 > getTemperature() { return m_temperature; }


  /**
   * @brief Create a instantiation of the
   *        ElasticCubicThermallySofteningUpdates class that refers to the
   *        data in this.
   * @return An instantiation of ElasticCubicThermallySofteningUpdates.
   */
  ElasticCubicThermallySofteningUpdates createKernelUpdates() const
  {
    return ElasticCubicThermallySofteningUpdates( m_defaultC11,
                                                  m_defaultC12,
                                                  m_defaultC44,
                                                  m_firstOrderC11ThermalCoefficient,
                                                  m_firstOrderC12ThermalCoefficient,
                                                  m_firstOrderC44ThermalCoefficient,
                                                  m_secondOrderC11ThermalCoefficient,
                                                  m_secondOrderC12ThermalCoefficient,
                                                  m_secondOrderC44ThermalCoefficient,
                                                  m_referenceTemperature,
                                                  m_temperature,
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
                          m_defaultC11,
                          m_defaultC12,
                          m_defaultC44,
                          m_firstOrderC11ThermalCoefficient,
                          m_firstOrderC12ThermalCoefficient,
                          m_firstOrderC44ThermalCoefficient,
                          m_secondOrderC11ThermalCoefficient,
                          m_secondOrderC12ThermalCoefficient,
                          m_secondOrderC44ThermalCoefficient,
                          m_referenceTemperature,
                          m_temperature,
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

  /// The default value of the first order thermal C11 coefficient for any new allocations.
  real64 m_firstOrderC11ThermalCoefficient;

  /// The default value of the first order thermal C12 coefficient for any new allocations.
  real64 m_firstOrderC12ThermalCoefficient;

  /// The default value of the first order thermal C44 coefficient for any new allocations.
  real64 m_firstOrderC44ThermalCoefficient;

  /// The default value of the second order thermal C11 coefficient for any new allocations.
  real64 m_secondOrderC11ThermalCoefficient;

  /// The default value of the second order thermal C12 coefficient for any new allocations.
  real64 m_secondOrderC12ThermalCoefficient;

  /// The default value of the second order thermal C44 coefficient for any new allocations.
  real64 m_secondOrderC44ThermalCoefficient;

  /// The default value of the reference temperature
  real64 m_referenceTemperature;

  // Array storing element temperature
  array1d< real64 > m_temperature;
};

} /* namespace constitutive */

} /* namespace geos */

#endif /* GEOS_CONSTITUTIVE_SOLID_ELASTICCUBIC_HPP_ */
