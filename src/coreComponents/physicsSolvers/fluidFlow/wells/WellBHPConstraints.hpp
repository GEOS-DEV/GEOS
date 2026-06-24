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
 * @file WellBHPConstraints.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLBHPCONSTRAINTS_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLBHPCONSTRAINTS_HPP

#include "common/format/EnumStrings.hpp"
#include "dataRepository/Group.hpp"
#include "functions/TableFunction.hpp"
#include "WellConstraintsBase.hpp"
namespace geos
{

enum class BHPConstraintTypeId : integer
{
  MIN,    /**< The well operates at a specified minimum bottom hole pressure (BHP) */
  MAX,    /**< The well operates at a specified maximum bottom hole pressure (BHP) */
  UNINITIALIZED,   /**< This is the current well control before postInputInitialization (needed to restart from file properly) */
};


/**
 * @class BHPConstraint
 * @brief This class describes a minimum pressure constraint used to control a injection well.
 */
template< BHPConstraintTypeId T >
class BHPConstraint : public WellConstraintBase
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
  explicit BHPConstraint( string const & name, dataRepository::Group * const parent );


  /**
   * @brief Default destructor.
   */
  ~BHPConstraint() override;

  /**
   * @brief Deleted default constructor.
   */
  BHPConstraint() = delete;

  /**
   * @brief Deleted copy constructor.
   */
  BHPConstraint( BHPConstraint const & ) = delete;

  /**
   * @brief Deleted move constructor.
   */
  BHPConstraint( BHPConstraint && ) = delete;

  /**
   * @brief Deleted assignment operator.
   * @return a reference to a constraint object
   */
  BHPConstraint & operator=( BHPConstraint const & ) = delete;

  /**
   * @brief Deleted move operator.
   * @return a reference to a constraint object
   */
  BHPConstraint & operator=( BHPConstraint && ) = delete;

  ///@}

  /**
   * @name Getters / Setters
   */
  ///@{

  // Temp interface - tjb
  virtual ConstraintTypeId getControl() const override { return ConstraintTypeId::BHP; };

  ///@}
  /**
   * @brief Struct to serve as a container for variable strings and keys.
   * @struct viewKeyStruct
   */
  struct viewKeyStruct
  {
    /// String key for the well target BHP
    static constexpr char const * targetBHPString() { return "targetBHP"; }
    /// String key for the well reference elevation (for BHP control)
    static constexpr char const * refElevString() { return "referenceElevation"; }
  }
  viewKeysWellBHPConstraint;

  virtual bool checkViolation( WellConstraintBase const & currentConstraint, real64 const & currentTime ) const override;

  /**
   * @brief Getter for the reference elevation where the BHP control is enforced
   * @return the reference elevation
   */
  real64 getReferenceElevation() const { return m_refElevation; }

  /**
   * @brief Set the reference elevation where the BHP control is enforced
   * @return the reference elevation
   */
  void setReferenceElevation( real64 const & refElevation ) { m_refElevation=refElevation; }

  /**
   * @brief Getter for the reference gravity coefficient
   * @return the reference gravity coefficient
   */
  real64 getReferenceGravityCoef() const { return m_refGravCoef; }

  /**
   * @brief Setter for the reference gravity
   */
  void setReferenceGravityCoef( real64 const & refGravCoef ) { m_refGravCoef = refGravCoef; }


  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new Constraint object through the object catalog.
   */
  static string catalogName()
  {
    if constexpr (T == BHPConstraintTypeId::MAX)
    {
      return "MaximumBHPConstraint";
    }
    else
    {
      return "MinimumBHPConstraint";
    }

  }


  virtual string getCatalogName() const override { return catalogName(); }
protected:

  virtual void postInputInitialization() override;

  /// Reference elevation
  real64 m_refElevation;

  /// Gravity coefficient of the reference elevation
  real64 m_refGravCoef;

};

using MinimumBHPConstraint = BHPConstraint< BHPConstraintTypeId::MIN >;
using MaximumBHPConstraint = BHPConstraint< BHPConstraintTypeId::MAX >;

} //namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLBHPCONSTRAINTS_HPP
