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
 * @file WellProductionConstraints.hpp
 */


#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLPRODUCTIONCONSTRAINT_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLPRODUCTIONCONSTRAINT_HPP

#include "common/format/EnumStrings.hpp"
#include "dataRepository/Group.hpp"
#include "functions/TableFunction.hpp"

namespace geos
{
using namespace dataRepository;
/**
 * @class ProductionConstraint
 * @brief This class describes constraint used to control a production well.
 */

template< typename ConstraintType >
class ProductionConstraint : public ConstraintType
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
  explicit ProductionConstraint( string const & name, dataRepository::Group * const parent );

  /**
   * @brief Default destructor.
   */
  ~ProductionConstraint() override;

  /**
   * @brief Deleted default constructor.
   */
  ProductionConstraint() = delete;

  /**
   * @brief Deleted copy constructor.
   */
  ProductionConstraint( ProductionConstraint const & ) = delete;

  /**
   * @brief Deleted move constructor.
   */
  ProductionConstraint( ProductionConstraint && ) = delete;

  /**
   * @brief Deleted assignment operator.
   * @return a reference to a constraint object
   */
  ProductionConstraint & operator=( ProductionConstraint const & ) = delete;

  /**
   * @brief Deleted move operator.
   * @return a reference to a constraint object
   */
  ProductionConstraint & operator=( ProductionConstraint && ) = delete;

  ///@}

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new Constraint object through the object catalog.
   */
  static string catalogName()
  {
    return "Production"+ConstraintType::catalogName();
  }
  virtual string getCatalogName() const override { return catalogName(); }
protected:

  virtual void postInputInitialization() override;

  static bool isViolated( const real64 & currentValue, const real64 & constraintValue )
  { return currentValue < constraintValue; }
};


} //namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLPRODUCTIONCONSTRAINT_HPP
