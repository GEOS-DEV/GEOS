/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/*
 * @file WellInjectionConstraint.hpp
 */


#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLINJECTIONCONSTRAINT_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLINJECTIONCONSTRAINT_HPP

#include "common/format/EnumStrings.hpp"
#include "dataRepository/Group.hpp"
#include "functions/TableFunction.hpp"

namespace geos
{

using namespace dataRepository;
/**
 * @class InjectionConstraint
 * @brief This class describes   constraint used to control a injection well.
 */

template< typename ConstraintType >
class InjectionConstraint : public ConstraintType
{
public:
  typedef InjectionConstraint< ConstraintType > classtype;
  /**
   * @name Constructor / Destructor
   */
  ///@{

  /**
   * @brief Constructor for WellControls Objects.
   * @param[in] name the name of this instantiation of WellControls in the repository
   * @param[in] parent the parent group of this instantiation of WellControls
   */
  explicit InjectionConstraint( string const & name, dataRepository::Group * const parent );

  /**
   * @brief Default destructor.
   */
  ~InjectionConstraint() override;

  /**
   * @brief Deleted default constructor.
   */
  InjectionConstraint() = delete;

  /**
   * @brief Deleted copy constructor.
   */
  InjectionConstraint( InjectionConstraint const & ) = delete;

  /**
   * @brief Deleted move constructor.
   */
  InjectionConstraint( InjectionConstraint && ) = delete;

  /**
   * @brief Deleted assignment operator.
   * @return a reference to a constraint object
   */
  InjectionConstraint & operator=( InjectionConstraint const & ) = delete;

  /**
   * @brief Deleted move operator.
   * @return a reference to a constraint object
   */
  InjectionConstraint & operator=( InjectionConstraint && ) = delete;

  ///@}

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new Constraint object through the object catalog.
   */
  static string catalogName()
  {
    return "Injection"+ConstraintType::catalogName();
  }

  virtual string getCatalogName() const override { return catalogName(); }

  struct injectionStreamKey
  {
    /// String key for the well injection stream
    static constexpr char const * injectionStreamString() { return "injectionStream"; }
    /// String key for the well injection temperature
    static constexpr char const * injectionTemperatureString() { return "injectionTemperature"; }
  };

  /**
   * @brief Const accessor for the composition of the injection stream
   * @return a global component fraction vector
   */
  arrayView1d< real64 const > getInjectionStream() const { return m_injectionStream; }

  /**
   * @brief Const accessor for the temperature of the injection stream
   * @return the temperature of the injection stream
   */
  real64 getInjectionTemperature() const { return m_injectionTemperature; }

  /**
   * @brief Set composition of the injection stream
   * @param[in] injectionStream a global component fraction vector
   */
  void setInjectionStream( arrayView1d< real64 const > const & injectionStream ) { m_injectionStream = injectionStream; }

  /**
   * @brief Set temperature of the injection stream
   * @param[in] temperature the temperature of the injection stream
   */
  void setInjectionTemperature( real64 temperature ) { m_injectionTemperature = temperature; }

protected:

  virtual void postInputInitialization() override;
  static bool isViolated( const real64 & currentValue, const real64 & constraintValue )
  { return currentValue > constraintValue; }

  void validateInjectionStream();
private:

  /// Vector with global component fractions at the injector
  array1d< real64 > m_injectionStream;

  /// Temperature at the injector
  real64 m_injectionTemperature;

};


} //namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLINJECTIONCONSTRAINT_HPP
