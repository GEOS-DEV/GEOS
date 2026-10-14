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


#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLWHPCONSTRAINTS_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLWHPCONSTRAINTS_HPP

#include "common/format/EnumStrings.hpp"
#include "dataRepository/Group.hpp"
#include "functions/TableFunction.hpp"
#include "WellConstraintsBase.hpp"
namespace geos
{

/**
 * @class WHPConstraint
 * @brief This class describes a minimum pressure constraint used to control a injection well.
 */
class WHPConstraint : public WellConstraintBase
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
  explicit WHPConstraint( string const & name, dataRepository::Group * const parent );


  /**
   * @brief Default destructor.
   */
  ~WHPConstraint() override;

  /**
   * @brief Deleted default constructor.
   */
  WHPConstraint() = delete;

  /**
   * @brief Deleted copy constructor.
   */
  WHPConstraint( WHPConstraint const & ) = delete;

  /**
   * @brief Deleted move constructor.
   */
  WHPConstraint( WHPConstraint && ) = delete;

  /**
   * @brief Deleted assignment operator.
   * @return a reference to a constraint object
   */
  WHPConstraint & operator=( WHPConstraint const & ) = delete;

  /**
   * @brief Deleted move operator.
   * @return a reference to a constraint object
   */
  WHPConstraint & operator=( WHPConstraint && ) = delete;

  ///@}


  /**
   * @name Getters / Setters
   */
  ///@{

  // Temp interface - tjb
  virtual ConstraintTypeId getControl() const override { return ConstraintTypeId::WHP; };

  /**
   * @name Getters / Setters
   */
  ///@{

  ///@}

  ///@}
  /**
   * @brief Struct to serve as a container for variable strings and keys.
   * @struct viewKeyStruct
   */
  struct viewKeyStruct
  {
    /// String key for the well target WHP
    static constexpr char const * targetWHPString() { return "targetWHP"; }
    /// String key for the well reference elevation (for WHP control)
    static constexpr char const * refElevString() { return "referenceElevation"; }
    /// string key for constraint values entered table name
    static constexpr char const * flowTableNameString() { return "flowTableName"; }
  }
  viewKeysWellWHPConstraint;

  //virtual bool checkViolation( WellConstraintBase const & currentConstraint, real64 const & currentTime ) const override;

  /**
   * @brief Getter for the reference elevation where the WHP control is enforced
   * @return the reference elevation
   */
  real64 getReferenceElevation() const { return m_refElevation; }

  /**
   * @brief Getter for the reference gravity coefficient
   * @return the reference gravity coefficient
   */
  real64 getReferenceGravityCoef() const { return m_refGravCoef; }

  /**
   * @brief Getter for the flow table name
   * @return the name of the flow table
   */
  std::string  getFlowTableName() const { return m_flowTableName; }

  /**
   * @brief Setter for the reference gravity
   */
  void setReferenceGravityCoef( real64 const & refGravCoef ) { m_refGravCoef = refGravCoef; }



protected:

  virtual void postInputInitialization() override;

  /// Reference elevation
  real64 m_refElevation;

  /// Gravity coefficient of the reference elevation
  real64 m_refGravCoef;

  /// Flow table name
  string m_flowTableName;

};

/**
 * @class MinimumWHPConstraint
 * @brief This class describes a minimum pressure constraint used to control a injection well.
 */
class MinimumWHPConstraint : public WHPConstraint
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
  explicit MinimumWHPConstraint( string const & name, dataRepository::Group * const parent );


  /**
   * @brief Default destructor.
   */
  ~MinimumWHPConstraint() override;

  /**
   * @brief Deleted default constructor.
   */
  MinimumWHPConstraint() = delete;

  /**
   * @brief Deleted copy constructor.
   */
  MinimumWHPConstraint( MinimumWHPConstraint const & ) = delete;

  /**
   * @brief Deleted move constructor.
   */
  MinimumWHPConstraint( MinimumWHPConstraint && ) = delete;

  /**
   * @brief Deleted assignment operator.
   * @return a reference to a constraint object
   */
  MinimumWHPConstraint & operator=( MinimumWHPConstraint const & ) = delete;

  /**
   * @brief Deleted move operator.
   * @return a reference to a constraint object
   */
  MinimumWHPConstraint & operator=( MinimumWHPConstraint && ) = delete;

  ///@}


  /**
   * @name Getters / Setters
   */
  ///@{

  // Temp interface - tjb
  virtual ConstraintTypeId getControl() const override { return ConstraintTypeId::WHP; };

  /**
   * @brief Struct to serve as a container for variable strings and keys.
   * @struct viewKeyStruct
   */
  struct viewKeyStruct
  {
    /// String key for the well target WHP
    static constexpr char const * targetWHPString() { return "targetWHP"; }
  }
  viewKeysWellWHPConstraint;
  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new Constraint object through the object catalog.
   */
  static string catalogName()
  {
    return "MinimumWHPConstraint";
  }
  virtual string getCatalogName() const override { return catalogName(); }
  virtual bool checkViolation( WellConstraintBase const & currentConstraint, real64 const & currentTime ) const override;



protected:

  virtual void postInputInitialization() override;


};

/**
 * @class WellMinimumWHPConstraint
 * @brief This class describes a maximum pressure constraint used to control a injection well.
 */
class MaximumWHPConstraint : public WHPConstraint
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
  explicit MaximumWHPConstraint( string const & name, dataRepository::Group * const parent );


  /**
   * @brief Default destructor.
   */
  ~MaximumWHPConstraint() override;

  /**
   * @brief Deleted default constructor.
   */
  MaximumWHPConstraint() = delete;

  /**
   * @brief Deleted copy constructor.
   */
  MaximumWHPConstraint( MaximumWHPConstraint const & ) = delete;

  /**
   * @brief Deleted move constructor.
   */
  MaximumWHPConstraint( MaximumWHPConstraint && ) = delete;

  /**
   * @brief Deleted assignment operator.
   * @return a reference to a constraint object
   */
  MaximumWHPConstraint & operator=( MaximumWHPConstraint const & ) = delete;

  /**
   * @brief Deleted move operator.
   * @return a reference to a constraint object
   */
  MaximumWHPConstraint & operator=( MaximumWHPConstraint && ) = delete;

  ///@}


  /**
   * @name Getters / Setters
   */
  ///@{
  // Temp interface - tjb
  virtual ConstraintTypeId getControl() const override { return ConstraintTypeId::WHP; };


  ///@}
  /**
   * @brief Struct to serve as a container for variable strings and keys.
   * @struct viewKeyStruct
   */
  struct viewKeyStruct
  {
    /// String key for the well target WHP
    static constexpr char const * targetWHPString() { return "targetWHP"; }
  } viewKeysWellWHPConstraint;
  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new Constraint object through the object catalog.
   */
  static string catalogName()
  {
    return "MaximumWHPConstraint";
  }
  virtual string getCatalogName() const override { return catalogName(); }
  virtual bool checkViolation( WellConstraintBase const & currentConstraint, real64 const & currentTime ) const override;
protected:

  virtual void postInputInitialization() override;



private:


};


} //namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLWHPCONSTRAINTS_HPP
