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


#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLVOLUMERATECONSTRAINTS_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLVOLUMERATECONSTRAINTS_HPP

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
static constexpr auto volumeProductionConstraint = "VolumeProductionConstraint";
static constexpr auto volumeInjectionConstraint = "VolumeInjectionConstraint";
}
}


/**
 * @class VolumeConstraint
 * @brief This class describes a volume rate constraint used to control a well.
 */

class VolumeConstraint : public WellConstraintBase
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
  explicit VolumeConstraint( string const & name, dataRepository::Group * const parent );


  /**
   * @brief Default destructor.
   */
  ~VolumeConstraint() override;

  /**
   * @brief Deleted default constructor.
   */
  VolumeConstraint() = delete;

  /**
   * @brief Deleted copy constructor.
   */
  VolumeConstraint( VolumeConstraint const & ) = delete;

  /**
   * @brief Deleted move constructor.
   */
  VolumeConstraint( VolumeConstraint && ) = delete;

  /**
   * @brief Deleted assignment operator.
   * @return a reference to a constraint object
   */
  VolumeConstraint & operator=( VolumeConstraint const & ) = delete;

  /**
   * @brief Deleted move operator.
   * @return a reference to a constraint object
   */
  VolumeConstraint & operator=( VolumeConstraint && ) = delete;

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
  virtual std::string getConstraintKey( ) const override { return "VolumeInjectionConstraint"; };
  ///@}

  // Temp interface - tjb
  virtual ConstraintTypeId getControl() const override { return ConstraintTypeId::TOTALVOLRATE; };

protected:

  virtual void postInputInitialization() override;

};


/**
 * @class VolumeProductionConstraint
 * @brief This class describes a volume rate constraint used to control a production well.
 */

class VolumeProductionConstraint : public VolumeConstraint
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
  explicit VolumeProductionConstraint( string const & name, dataRepository::Group * const parent );


  /**
   * @brief Default destructor.
   */
  ~VolumeProductionConstraint() override;

  /**
   * @brief Deleted default constructor.
   */
  VolumeProductionConstraint() = delete;

  /**
   * @brief Deleted copy constructor.
   */
  VolumeProductionConstraint( VolumeProductionConstraint const & ) = delete;

  /**
   * @brief Deleted move constructor.
   */
  VolumeProductionConstraint( VolumeProductionConstraint && ) = delete;

  /**
   * @brief Deleted assignment operator.
   * @return a reference to a constraint object
   */
  VolumeProductionConstraint & operator=( VolumeProductionConstraint const & ) = delete;

  /**
   * @brief Deleted move operator.
   * @return a reference to a constraint object
   */
  VolumeProductionConstraint & operator=( VolumeProductionConstraint && ) = delete;

  ///@}

  /**
   * @name Getters / Setters
   */
  ///@{
  /**
   * @brief Get name of constraint
   * @return constraint key
   */
  virtual std::string getConstraintKey( ) const override { return "VolumeProductionConstraint"; };
  ///@}

  virtual bool checkViolation( WellConstraintBase const & currentConstraint, real64 const & currentTime ) const override;
protected:

  virtual void postInputInitialization() override;

};

/**
 * @class VolumeInjectionConstraint
 * @brief This class describes a volume rate constraint used to control a injection well.
 */

class VolumeInjectionConstraint : public VolumeConstraint
{
public:


  /**
   * @name Constructor / Destructor
   */
  ///@{

  /**
   * @brief Constructor for VolumeInjectionConstraint Objects.
   * @param[in] name the name of this instantiation of WellControls in the repository
   * @param[in] parent the parent group of this instantiation of WellControls
   */
  explicit VolumeInjectionConstraint( string const & name, dataRepository::Group * const parent );


  /**
   * @brief Default destructor.
   */
  ~VolumeInjectionConstraint() override;

  /**
   * @brief Deleted default constructor.
   */
  VolumeInjectionConstraint() = delete;

  /**
   * @brief Deleted copy constructor.
   */
  VolumeInjectionConstraint( VolumeInjectionConstraint const & ) = delete;

  /**
   * @brief Deleted move constructor.
   */
  VolumeInjectionConstraint( VolumeInjectionConstraint && ) = delete;

  /**
   * @brief Deleted assignment operator.
   * @return a reference to a constraint object
   */
  VolumeInjectionConstraint & operator=( VolumeInjectionConstraint const & ) = delete;

  /**
   * @brief Deleted move operator.
   * @return a reference to a constraint object
   */
  VolumeInjectionConstraint & operator=( VolumeInjectionConstraint && ) = delete;

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
