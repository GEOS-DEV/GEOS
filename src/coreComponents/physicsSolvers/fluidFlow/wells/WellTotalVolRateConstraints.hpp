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
 * @file WellVolumeRateConstraints.hpp
 */


#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLTOTALVOLRATECONSTRAINTS_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLTOTALVOLRATECONSTRAINTS_HPP

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
static constexpr auto totalVolProductionConstraint = "TotalVolProductionConstraint";
static constexpr auto totalVolInjectionConstraint = "TotalVolInjectionConstraint";
}
}


/**
 * @class TotalVolConstraint
 * @brief This class describes a volume rate constraint used to control a well.
 */

class TotalVolConstraint : public WellConstraintBase
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
  explicit TotalVolConstraint( string const & name, dataRepository::Group * const parent );


  /**
   * @brief Default destructor.
   */
  ~TotalVolConstraint() override;

  /**
   * @brief Deleted default constructor.
   */
  TotalVolConstraint() = delete;

  /**
   * @brief Deleted copy constructor.
   */
  TotalVolConstraint( TotalVolConstraint const & ) = delete;

  /**
   * @brief Deleted move constructor.
   */
  TotalVolConstraint( TotalVolConstraint && ) = delete;

  /**
   * @brief Deleted assignment operator.
   * @return a reference to a constraint object
   */
  TotalVolConstraint & operator=( TotalVolConstraint const & ) = delete;

  /**
   * @brief Deleted move operator.
   * @return a reference to a constraint object
   */
  TotalVolConstraint & operator=( TotalVolConstraint && ) = delete;

  ///@}
  /**
   * @brief Struct to serve as a container for variable strings and keys.
   * @struct viewKeyStruct
   */
  struct viewKeyStruct
  {
    /// String key for the volume rate
    static constexpr char const * volumeRateString() { return "volumeRate"; }
  };
  /**
   * @name Getters / Setters
   */
  ///@{
  /**
   * @brief Get name of constraint
   * @return constraint key
   */
  virtual std::string getConstraintKey( ) const override { return "TotalVolInjectionConstraint"; };
  ///@}

  // Temp interface - tjb
  virtual ConstraintTypeId getControl() const override { return ConstraintTypeId::TOTALVOLRATE; };

protected:

  virtual void postInputInitialization() override;

};


/**
 * @class TotalVolProductionConstraint
 * @brief This class describes a volume rate constraint used to control a production well.
 */

class TotalVolProductionConstraint : public TotalVolConstraint
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
  explicit TotalVolProductionConstraint( string const & name, dataRepository::Group * const parent );


  /**
   * @brief Default destructor.
   */
  ~TotalVolProductionConstraint() override;

  /**
   * @brief Deleted default constructor.
   */
  TotalVolProductionConstraint() = delete;

  /**
   * @brief Deleted copy constructor.
   */
  TotalVolProductionConstraint( TotalVolProductionConstraint const & ) = delete;

  /**
   * @brief Deleted move constructor.
   */
  TotalVolProductionConstraint( TotalVolProductionConstraint && ) = delete;

  /**
   * @brief Deleted assignment operator.
   * @return a reference to a constraint object
   */
  TotalVolProductionConstraint & operator=( TotalVolProductionConstraint const & ) = delete;

  /**
   * @brief Deleted move operator.
   * @return a reference to a constraint object
   */
  TotalVolProductionConstraint & operator=( TotalVolProductionConstraint && ) = delete;

  ///@}

  /**
   * @name Getters / Setters
   */
  ///@{
  /**
   * @brief Get name of constraint
   * @return constraint key
   */
  virtual std::string getConstraintKey( ) const override { return "TotalVolProductionConstraint"; };
  ///@}

  virtual bool checkViolation( WellConstraintBase const & currentConstraint, real64 const & currentTime ) const override;
protected:

  virtual void postInputInitialization() override;

};

/**
 * @class TotalVolInjectionConstraint
 * @brief This class describes a volume rate constraint used to control a injection well.
 */

class TotalVolInjectionConstraint : public TotalVolConstraint
{
public:


  /**
   * @name Constructor / Destructor
   */
  ///@{

  /**
   * @brief Constructor for TotalVolInjectionConstraint Objects.
   * @param[in] name the name of this instantiation of WellControls in the repository
   * @param[in] parent the parent group of this instantiation of WellControls
   */
  explicit TotalVolInjectionConstraint( string const & name, dataRepository::Group * const parent );


  /**
   * @brief Default destructor.
   */
  ~TotalVolInjectionConstraint() override;

  /**
   * @brief Deleted default constructor.
   */
  TotalVolInjectionConstraint() = delete;

  /**
   * @brief Deleted copy constructor.
   */
  TotalVolInjectionConstraint( TotalVolInjectionConstraint const & ) = delete;

  /**
   * @brief Deleted move constructor.
   */
  TotalVolInjectionConstraint( TotalVolInjectionConstraint && ) = delete;

  /**
   * @brief Deleted assignment operator.
   * @return a reference to a constraint object
   */
  TotalVolInjectionConstraint & operator=( TotalVolInjectionConstraint const & ) = delete;

  /**
   * @brief Deleted move operator.
   * @return a reference to a constraint object
   */
  TotalVolInjectionConstraint & operator=( TotalVolInjectionConstraint && ) = delete;

  ///@}

  /**
   * @name Getters / Setters
   */
  ///@{


  ///@}

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
