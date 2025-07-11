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
 * @file WellControls.hpp
 */


#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLPRESSURECONSTRAINTS_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLPRESSURECONSTRAINTS_HPP

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
static constexpr auto minimumBHPConstraint = "MinimumBHPConstraint";
static constexpr auto maximumBHPConstraint = "MaximumBHPConstraint";
}
}

/**
 * @class MinimumBHPConstraint
 * @brief This class describes a minimum pressure constraint used to control a injection well.
 */
class MinimumBHPConstraint : public WellConstraintBase
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
  explicit MinimumBHPConstraint( string const & name, dataRepository::Group * const parent );


  /**
   * @brief Default destructor.
   */
  ~MinimumBHPConstraint() override;

  /**
   * @brief Deleted default constructor.
   */
  MinimumBHPConstraint() = delete;

  /**
   * @brief Deleted copy constructor.
   */
  MinimumBHPConstraint( MinimumBHPConstraint const & ) = delete;

  /**
   * @brief Deleted move constructor.
   */
  MinimumBHPConstraint( MinimumBHPConstraint && ) = delete;

  /**
   * @brief Deleted assignment operator.
   * @return a reference to a constraint object
   */
  MinimumBHPConstraint & operator=( MinimumBHPConstraint const & ) = delete;

  /**
   * @brief Deleted move operator.
   * @return a reference to a constraint object
   */
  MinimumBHPConstraint & operator=( MinimumBHPConstraint && ) = delete;

  ///@}


  /**
   * @name Getters / Setters
   */
  ///@{

  // Temp interface - tjb
  virtual WellControls::Control getControl() const override { return WellControls::Control::BHP; };
  /**
   * @brief Get name of constraint
   * @return constraint key
   */
  virtual std::string  getConstraintKey( ) const override { return dataRepository::keys::minimumBHPConstraint; };
  ///@}
  /**
   * @brief Struct to serve as a container for variable strings and keys.
   * @struct viewKeyStruct
   */
  struct viewKeyStruct
  {
    /// String key for the well target BHP
    static constexpr char const * targetBHPString() { return "targetBHP"; }
  }
  viewKeysWellBHPConstraint;

protected:

  virtual void postInputInitialization() override;


};

/**
 * @class WellMinimumBHPConstraint
 * @brief This class describes a maximum pressure constraint used to control a injection well.
 */
class MaximumBHPConstraint : public WellConstraintBase
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
  explicit MaximumBHPConstraint( string const & name, dataRepository::Group * const parent );


  /**
   * @brief Default destructor.
   */
  ~MaximumBHPConstraint() override;

  /**
   * @brief Deleted default constructor.
   */
  MaximumBHPConstraint() = delete;

  /**
   * @brief Deleted copy constructor.
   */
  MaximumBHPConstraint( MaximumBHPConstraint const & ) = delete;

  /**
   * @brief Deleted move constructor.
   */
  MaximumBHPConstraint( MaximumBHPConstraint && ) = delete;

  /**
   * @brief Deleted assignment operator.
   * @return a reference to a constraint object
   */
  MaximumBHPConstraint & operator=( MaximumBHPConstraint const & ) = delete;

  /**
   * @brief Deleted move operator.
   * @return a reference to a constraint object
   */
  MaximumBHPConstraint & operator=( MaximumBHPConstraint && ) = delete;

  ///@}


  /**
   * @name Getters / Setters
   */
  ///@{
  // Temp interface - tjb
  virtual WellControls::Control getControl() const override { return WellControls::Control::BHP; };
  /**
   * @brief Get name of constraint
   * @return constraint key
   */
  virtual std::string getConstraintKey( ) const override { return dataRepository::keys::maximumBHPConstraint; };

  ///@}
  /**
   * @brief Struct to serve as a container for variable strings and keys.
   * @struct viewKeyStruct
   */
  struct viewKeyStruct
  {
    /// String key for the well target BHP
    static constexpr char const * targetBHPString() { return "targetBHP"; }
  }
  viewKeysWellBHPConstraint;

protected:

  virtual void postInputInitialization() override;

private:


};


} //namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLPRESSURECONSTRAINTS_HPP
