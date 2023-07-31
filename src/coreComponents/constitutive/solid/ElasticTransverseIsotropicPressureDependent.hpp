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
 *  @file ElasticTransverseIsotropicPressureDependent.hpp
 */

#ifndef GEOS_CONSTITUTIVE_SOLID_ElasticTransverseIsotropicPressureDependent_HPP_
#define GEOS_CONSTITUTIVE_SOLID_ElasticTransverseIsotropicPressureDependent_HPP_

#include "ElasticTransverseIsotropic.hpp"
#include "constitutive/ExponentialRelation.hpp"
#include "LvArray/src/tensorOps.hpp"
#include "SolidModelDiscretizationOpsTransverseIsotropic.hpp"

namespace geos
{

namespace constitutive
{

/**
 * @class ElasticTransverseIsotropicPressureDependentUpdates
 *
 * Class to provide elastic transverse isotropic material updates that
 * may be called from a kernel function.
 *
 * @note The "transverse" directions are 1 and 2 (or 0 and 1 in C-index)
 */
class ElasticTransverseIsotropicPressureDependentUpdates : public ElasticTransverseIsotropicUpdates
{
public:
  /**
   * @brief Constructor
   * @param[in] c11 The 11 component of the Voigt stiffness tensor.
   * @param[in] c13 The 13 component of the Voigt stiffness tensor.
   * @param[in] c33 The 33 component of the Voigt stiffness tensor.
   * @param[in] c44 The 44 component of the Voigt stiffness tensor.
   * @param[in] c66 The 66 component of the Voigt stiffness tensor.
   * @param[in] dc11dp The pressure derivative of the 11 component of the Voigt stiffness tensor.
   * @param[in] dc13dp The pressure derivative of the 13 component of the Voigt stiffness tensor.
   * @param[in] dc33dp The pressure derivative of the 33 component of the Voigt stiffness tensor.
   * @param[in] dc44dp The pressure derivative of the 44 component of the Voigt stiffness tensor.
   * @param[in] dc66dp The pressure derivative of the 66 component of the Voigt stiffness tensor.
   * @param[in] thermalExpansionCoefficient The ArrayView holding the thermal expansion coefficient data for each element.
   * @param[in] newStress The ArrayView holding the new stress data for each point.
   * @param[in] oldStress The ArrayView holding the old stress data for each point.
   * @param[in] disableInelasticity Flag to disable plastic response for inelastic models.
   */
  ElasticTransverseIsotropicPressureDependentUpdates( real64 const & dc11dp,
                                                      real64 const & dc13dp,
                                                      real64 const & dc33dp,
                                                      real64 const & dc44dp,
                                                      real64 const & dc66dp,
                                                      arrayView1d< real64 const > const & c11,
                                                      arrayView1d< real64 const > const & c13,
                                                      arrayView1d< real64 const > const & c33,
                                                      arrayView1d< real64 const > const & c44,
                                                      arrayView1d< real64 const > const & c66,
                                                      arrayView1d< real64 const > const & thermalExpansionCoefficient,
                                                      arrayView3d< real64, solid::STRESS_USD > const & newStress,
                                                      arrayView3d< real64, solid::STRESS_USD > const & oldStress,
                                                      bool const & disableInelasticity ):
    ElasticTransverseIsotropicUpdates( c11, c13, c33, c44, c66, thermalExpansionCoefficient, newStress, oldStress, disableInelasticity ),
    m_dc11dp( dc11dp ),
    m_dc13dp( dc13dp ),
    m_dc33dp( dc33dp ),
    m_dc44dp( dc44dp ),
    m_dc66dp( dc66dp )
  {}

  /// Deleted default constructor
  ElasticTransverseIsotropicPressureDependentUpdates() = delete;

  /// Default copy constructor
  ElasticTransverseIsotropicPressureDependentUpdates( ElasticTransverseIsotropicPressureDependentUpdates const & ) = default;

  /// Default move constructor
  ElasticTransverseIsotropicPressureDependentUpdates( ElasticTransverseIsotropicPressureDependentUpdates && ) = default;

  /// Deleted copy assignment operator
  ElasticTransverseIsotropicPressureDependentUpdates & operator=( ElasticTransverseIsotropicPressureDependentUpdates const & ) = delete;

  /// Deleted move assignment operator
  ElasticTransverseIsotropicPressureDependentUpdates & operator=( ElasticTransverseIsotropicPressureDependentUpdates && ) =  delete;

  // Use transverse isotropic form of inner product compression
  using DiscretizationOps = SolidModelDiscretizationOpsTransverseIsotropic;

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
                                         DiscretizationOps & stiffness ) const override;

  // incremental strain interfaces

  GEOS_HOST_DEVICE
  virtual void smallStrainUpdate_StressOnly( localIndex const k,
                                             localIndex const q,
                                             real64 const & timeIncrement,
                                             real64 const ( &strainIncrement )[6],
                                             real64 ( &stress )[6] ) const override;

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
  virtual void getElasticStiffness( localIndex const k, localIndex const q, real64 ( &stiffness )[6][6] ) const override;

protected:

  /// The pressure derivative of 11 component of the Voigt stiffness tensor.
  real64 const m_dc11dp;

  /// The pressure derivative of 13 component of the Voigt stiffness tensor.
  real64 const m_dc13dp;

  /// The pressure derivative of 33 component of the Voigt stiffness tensor.
  real64 const m_dc33dp;

  /// The pressure derivative of 44 component of the Voigt stiffness tensor.
  real64 const m_dc44dp;

  /// The pressure derivative of 66 component of the Voigt stiffness tensor.
  real64 const m_dc66dp;

};

inline
GEOS_HOST_DEVICE
void ElasticTransverseIsotropicPressureDependentUpdates::getElasticStiffness( localIndex const k,
                                                                              localIndex const q,
                                                                              real64 ( & stiffness )[6][6] ) const
{
  // GEOS_UNUSED_VAR( q );
  LvArray::tensorOps::fill< 6, 6 >( stiffness, 0 );

  real64 pressure = (-1.0/3.0)*(m_newStress[k][q][0] + m_newStress[k][q][1] + m_newStress[k][q][2]);

  real64 c11 = m_c11[k] + m_dc11dp * pressure;
  real64 c13 = m_c13[k] + m_dc13dp * pressure;
  real64 c33 = m_c33[k] + m_dc33dp * pressure;
  real64 c44 = m_c44[k] + m_dc44dp * pressure;
  real64 c66 = m_c66[k] + m_dc66dp * pressure;

  stiffness[0][0] = c11;
  stiffness[0][1] = c11 - 2 * c66;
  stiffness[0][2] = c13;

  stiffness[1][0] = stiffness[0][1];
  stiffness[1][1] = c11;
  stiffness[1][2] = c13;

  stiffness[2][0] = stiffness[0][2];
  stiffness[2][1] = stiffness[1][2];
  stiffness[2][2] = c33;

  stiffness[3][3] = c44;
  stiffness[4][4] = c44;
  stiffness[5][5] = c66;
}

inline
GEOS_HOST_DEVICE
void ElasticTransverseIsotropicPressureDependentUpdates::smallStrainNoStateUpdate_StressOnly( localIndex const k,
                                                                                              localIndex const q,
                                                                                              real64 const ( &totalStrain )[6],
                                                                                              real64 ( & stress )[6] ) const
{
  GEOS_UNUSED_VAR( q );

  real64 pressure = (-1.0/3.0)*(stress[0] + stress[1] + stress[2]);

  real64 c11 = m_c11[k] + m_dc11dp * pressure;
  real64 c13 = m_c13[k] + m_dc13dp * pressure;
  real64 c33 = m_c33[k] + m_dc33dp * pressure;
  real64 c44 = m_c44[k] + m_dc44dp * pressure;
  real64 c66 = m_c66[k] + m_dc66dp * pressure;

  real64 const c12temp = ( c11 - 2.0 * c66 );
  stress[0] = c11 * totalStrain[0] +  c12temp * totalStrain[1] + c13 * totalStrain[2];
  stress[1] = c12temp * totalStrain[0] + c11 * totalStrain[1] + c13 * totalStrain[2];
  stress[2] = c13 * totalStrain[0] + c13 * totalStrain[1] + c33 * totalStrain[2];

  stress[3] = c44 * totalStrain[3];
  stress[4] = c44 * totalStrain[4];
  stress[5] = c66 * totalStrain[5];
}

inline
GEOS_HOST_DEVICE
void ElasticTransverseIsotropicPressureDependentUpdates::smallStrainNoStateUpdate( localIndex const k,
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
void ElasticTransverseIsotropicPressureDependentUpdates::smallStrainNoStateUpdate( localIndex const k,
                                                                                   localIndex const q,
                                                                                   real64 const ( &totalStrain )[6],
                                                                                   real64 ( & stress )[6],
                                                                                   DiscretizationOps & stiffness ) const
{
  smallStrainNoStateUpdate_StressOnly( k, q, totalStrain, stress );

  real64 pressure = (-1.0/3.0)*(stress[0] + stress[1] + stress[2]);

  stiffness.m_c11 = m_c11[k] + m_dc11dp * pressure;
  stiffness.m_c13 = m_c13[k] + m_dc13dp * pressure;
  stiffness.m_c33 = m_c33[k] + m_dc33dp * pressure;
  stiffness.m_c44 = m_c44[k] + m_dc44dp * pressure;
  stiffness.m_c66 = m_c66[k] + m_dc66dp * pressure;
}

inline
GEOS_HOST_DEVICE
void ElasticTransverseIsotropicPressureDependentUpdates::smallStrainUpdate_StressOnly( localIndex const k,
                                                                      localIndex const q,
                                                                      real64 const & timeIncrement,
                                                                      real64 const ( &strainIncrement )[6],
                                                                      real64 ( & stress )[6] ) const
{
  GEOS_UNUSED_VAR( timeIncrement );
  smallStrainNoStateUpdate_StressOnly( k, q, strainIncrement, stress ); // stress = incrementalStress
  LvArray::tensorOps::add< 6 >( stress, m_oldStress[k][q] );            // stress += m_oldStress
  saveStress( k, q, stress );                                           // m_newStress = stress
}

inline
GEOS_HOST_DEVICE
void ElasticTransverseIsotropicPressureDependentUpdates::smallStrainUpdate( localIndex const k,
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
void ElasticTransverseIsotropicPressureDependentUpdates::smallStrainUpdate( localIndex const k,
                                                                            localIndex const q,
                                                                            real64 const & timeIncrement,
                                                                            real64 const ( &strainIncrement )[6],
                                                                            real64 ( & stress )[6],
                                                                            DiscretizationOps & stiffness ) const
{
  smallStrainUpdate_StressOnly( k, q, timeIncrement, strainIncrement, stress );

  real64 pressure = (-1.0/3.0)*(stress[0] + stress[1] + stress[2]);

  stiffness.m_c11 = m_c11[k] + m_dc11dp * pressure;
  stiffness.m_c13 = m_c13[k] + m_dc13dp * pressure;
  stiffness.m_c33 = m_c33[k] + m_dc33dp * pressure;
  stiffness.m_c44 = m_c44[k] + m_dc44dp * pressure;
  stiffness.m_c66 = m_c66[k] + m_dc66dp * pressure;
}

/**
 * @class ElasticTransverseIsotropicPressureDependent
 *
 * Class to provide a  elastic transverse isotropic material response.
 */
class ElasticTransverseIsotropicPressureDependent : public ElasticTransverseIsotropic
{
public:

  /// @typedef Alias for ElasticTransverseIsotropicPressureDependentUpdates
  using KernelWrapper = ElasticTransverseIsotropicPressureDependentUpdates;

  /**
   * @brief constructor
   * @param[in]name name of the instance in the catalog
   * @param[in]parent the group which contains this instance
   */
  ElasticTransverseIsotropicPressureDependent( string const & name, Group * const parent );

  /**
   * Destructor
   */
  virtual ~ElasticTransverseIsotropicPressureDependent() override;

  /**
   * @name Static Factory Catalog members and functions
   */
  ///@{

  /// string name to use for this class in the catalog
  static constexpr auto m_catalogNameString = "ElasticTransverseIsotropicPressureDependent";

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
    /// string/key for c11 component of Voigt stiffness tensor
    static constexpr char const * dc11dpString() { return "dc11dp"; }

    /// string/key for c13 component of Voigt stiffness tensor
    static constexpr char const * dc13dpString() { return "dc13dp"; }

    /// string/key for c33 component of Voigt stiffness tensor
    static constexpr char const * dc33dpString() { return "dc33dp"; }

    /// string/key for c44 component of Voigt stiffness tensor
    static constexpr char const * dc44dpString() { return "dc44dp"; }

    /// string/key for c66 component of Voigt stiffness tensor
    static constexpr char const * dc66dpString() { return "dc66dp"; }
  };

  /**
   * @brief Const-Getter for 11 component of Voigt stiffness tensor.
   * @return reference to immutable 11 component of Voigt stiffness tensor.
   */
  real64 const & getdC11dp() const { return m_dc11dp; }

  /**
   * @brief Getter for 11 component of Voigt stiffness tensor.
   * @return reference to mutable 11 component of Voigt stiffness tensor.
   */
  real64 & getdC11dp() { return m_dc11dp; }

  /**
   * @brief Const-Getter for 13 component of Voigt stiffness tensor.
   * @return reference to immutable 13 component of Voigt stiffness tensor.
   */
  real64 const & getdC13dp() const { return m_dc13dp; }

  /**
   * @brief Getter for 13 component of Voigt stiffness tensor.
   * @return reference to mutable 13 component of Voigt stiffness tensor.
   */
  real64 & getdC13dp() { return m_dc13dp; }

  /**
   * @brief Const-Getter for 33 component of Voigt stiffness tensor.
   * @return reference to immutable 33 component of Voigt stiffness tensor.
   */
  real64 const & getdC33dp() const { return m_dc33dp; }

  /**
   * @brief Getter for 33 component of Voigt stiffness tensor.
   * @return reference to mutable 33 component of Voigt stiffness tensor.
   */
  real64 & getdC33dp() { return m_dc33dp; }

  /**
   * @brief Const-Getter for 44 component of Voigt stiffness tensor.
   * @return reference to immutable 44 component of Voigt stiffness tensor.
   */
  real64 const & getdC44dp() const { return m_dc44dp; }

  /**
   * @brief Getter for 44 component of Voigt stiffness tensor.
   * @return reference to mutable 44 component of Voigt stiffness tensor.
   */
  real64 & getdC44dp() { return m_dc44dp; }

  /**
   * @brief Const-Getter for 66 component of Voigt stiffness tensor.
   * @return reference to immutable 66 component of Voigt stiffness tensor.
   */
  real64 const & getdC66dp() const { return m_dc66dp; }

  /**
   * @brief Getter for 66 component of Voigt stiffness tensor.
   * @return reference to mutable 66 component of Voigt stiffness tensor.
   */
  real64 & getdC66dp() { return m_dc66dp; }

  /**
   * @brief Create a instantiation of the
   *        ElasticTransverseIsotropicPressureDependentUpdates class that refers to the
   *        data in this.
   * @return An instantiation of ElasticTransverseIsotropicPressureDependentUpdates.
   */
  ElasticTransverseIsotropicPressureDependentUpdates createKernelUpdates() const
  {
    return ElasticTransverseIsotropicPressureDependentUpdates( m_dc11dp,
                                                               m_dc13dp,
                                                               m_dc33dp,
                                                               m_dc44dp,
                                                               m_dc66dp,
                                                               m_c11,
                                                               m_c13,
                                                               m_c33,
                                                               m_c44,
                                                               m_c66,
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
                          m_dc11dp,
                          m_dc13dp,
                          m_dc33dp,
                          m_dc44dp,
                          m_dc66dp,
                          m_c11,
                          m_c13,
                          m_c33,
                          m_c44,
                          m_c66,
                          m_thermalExpansionCoefficient,
                          m_newStress,
                          m_oldStress,
                          m_disableInelasticity );
  }

protected:
  virtual void postProcessInput() override;

  /// The pressure derivative of the 11 component of the Voigt stiffness tensor.
  real64 m_dc11dp;

  /// The pressure derivative of the 13 component of the Voigt stiffness tensor.
  real64 m_dc13dp;

  /// The pressure derivative of the 33 component of the Voigt stiffness tensor.
  real64 m_dc33dp;

  /// The pressure derivative of the 44 component of the Voigt stiffness tensor.
  real64 m_dc44dp;

  /// The pressure derivative of the 66 component of the Voigt stiffness tensor.
  real64 m_dc66dp;
};

} /* namespace constitutive */

} /* namespace geos */

#endif /* GEOS_CONSTITUTIVE_SOLID_ElasticTransverseIsotropicPressureDependent_HPP_ */
