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
 * @file WellMassRateConstraints.hpp
 */


#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLMASSRATECONSTRAINTS_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLMASSRATECONSTRAINTS_HPP

#include "common/format/EnumStrings.hpp"
#include "dataRepository/Group.hpp"
#include "functions/TableFunction.hpp"
#include "WellConstraintsBase.hpp"
namespace geos
{
namespace dataRepository
{
namespace keys
{
static constexpr auto MassProductionConstraint = "MassProductionConstraint";
static constexpr auto MassInjectionConstraint = "MassInjectionConstraint";
}
}

namespace constraintViewStruct
{

struct MassConstraintKey
{
  // String key for the well target mass rate
  static constexpr char const * targetMassRateString() { return "targetMassRate"; }
};

}

/**
 * @class MassConstraint
 * @brief This class describes a mass rate constraint used to control a  well.
 */

class MassConstraint : public WellConstraintBase
{
public:

  /**
   * @name Constructor / Destructor
   */
  ///@{

  /**
   * @brief Constructor for WellControls Objects.
   * @param[in] name the name of this instantiation of WellControls in the repository
   * @param[in] parent the parent group of this instantiation of WellControls
   */
  explicit MassConstraint( string const & name, dataRepository::Group * const parent );


  /**
   * @brief Default destructor.
   */
  ~MassConstraint() override;

  /**
   * @brief Deleted default constructor.
   */
  MassConstraint() = delete;

  /**
   * @brief Deleted copy constructor.
   */
  MassConstraint( MassConstraint const & ) = delete;

  /**
   * @brief Deleted move constructor.
   */
  MassConstraint( MassConstraint && ) = delete;

  /**
   * @brief Deleted assignment operator.
   * @return a reference to a constraint object
   */
  MassConstraint & operator=( MassConstraint const & ) = delete;

  /**
   * @brief Deleted move operator.
   * @return a reference to a constraint object
   */
  MassConstraint & operator=( MassConstraint && ) = delete;

  ///@}

  /**
   * @name Getters / Setters
   */

  // Temp interface - tjb
  virtual ConstraintTypeId getControl() const override { return ConstraintTypeId::MASSRATE; };
  ///@}

  // Mass constraint defintion keys
  constraintViewStruct::MassConstraintKey viewKeysMassConstraint;


protected:

  virtual void postInputInitialization() override;



};


/**
 * @class MassProductionConstraint
 * @brief This class describes a mass rate constraint used to control a production well.
 */

class MassProductionConstraint : public MassConstraint
{
public:

  /**
   * @name Constructor / Destructor
   */
  ///@{

  /**
   * @brief Constructor for WellControls Objects.
   * @param[in] name the name of this instantiation of WellControls in the repository
   * @param[in] parent the parent group of this instantiation of WellControls
   */
  explicit MassProductionConstraint( string const & name, dataRepository::Group * const parent );


  /**
   * @brief Default destructor.
   */
  ~MassProductionConstraint() override;

  /**
   * @brief Deleted default constructor.
   */
  MassProductionConstraint() = delete;

  /**
   * @brief Deleted copy constructor.
   */
  MassProductionConstraint( MassProductionConstraint const & ) = delete;

  /**
   * @brief Deleted move constructor.
   */
  MassProductionConstraint( MassProductionConstraint && ) = delete;

  /**
   * @brief Deleted assignment operator.
   * @return a reference to a constraint object
   */
  MassProductionConstraint & operator=( MassProductionConstraint const & ) = delete;

  /**
   * @brief Deleted move operator.
   * @return a reference to a constraint object
   */
  MassProductionConstraint & operator=( MassProductionConstraint && ) = delete;

  ///@}

  /**
   * @name Getters / Setters
   */

  /**
   * @brief Get name of constraint
   * @return constraint key
   */
  virtual std::string  getConstraintKey( ) const override { return dataRepository::keys::MassProductionConstraint; };
  ///@}

  // Mass constraint defintion keys
  constraintViewStruct::MassConstraintKey viewKeysMassConstraint;



  virtual bool checkViolation( WellConstraintBase const & currentConstraint, real64 const & currentTime ) const override;

protected:

  virtual void postInputInitialization() override;



};

/**
 * @class MassInjectionConstraint
 * @brief This class describes a Mass rate constraint used to control a injection well.
 */

class MassInjectionConstraint : public MassConstraint
{
public:


  /**
   * @name Constructor / Destructor
   */
  ///@{

  /**
   * @brief Constructor for MassInjectionConstraint Objects.
   * @param[in] name the name of this instantiation of WellControls in the repository
   * @param[in] parent the parent group of this instantiation of WellControls
   */
  explicit MassInjectionConstraint( string const & name, dataRepository::Group * const parent );


  /**
   * @brief Default destructor.
   */
  ~MassInjectionConstraint() override;

  /**
   * @brief Deleted default constructor.
   */
  MassInjectionConstraint() = delete;

  /**
   * @brief Deleted copy constructor.
   */
  MassInjectionConstraint( MassInjectionConstraint const & ) = delete;

  /**
   * @brief Deleted move constructor.
   */
  MassInjectionConstraint( MassInjectionConstraint && ) = delete;

  /**
   * @brief Deleted assignment operator.
   * @return a reference to a constraint object
   */
  MassInjectionConstraint & operator=( MassInjectionConstraint const & ) = delete;

  /**
   * @brief Deleted move operator.
   * @return a reference to a constraint object
   */
  MassInjectionConstraint & operator=( MassInjectionConstraint && ) = delete;

  ///@}

  /**
   * @name Getters / Setters
   */
  ///@{
  /**
   * @brief Get name of constraint
   * @return constraint key
   */
  virtual std::string getConstraintKey( ) const override { return dataRepository::keys::MassInjectionConstraint; };
  ///@}

  // Mass constraint defintion keys
  constraintViewStruct::MassConstraintKey viewKeysMassConstraint;



  // Injection stream definition keys
  constraintViewStruct::injectionStreamKey viewKeysInjectionStream;

  virtual bool checkViolation( WellConstraintBase const & currentConstraint, real64 const & currentTime ) const override;

protected:

  virtual void postInputInitialization() override;

private:

  /// Vector with global component fractions at the injector
  array1d< real64 > m_injectionStream;

  /// Temperature at the injector
  real64 m_injectionTemperature;

};



} //namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLCONSTRAINT_HPP
