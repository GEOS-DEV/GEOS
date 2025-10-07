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

/**
 * @class MassConstraint
 * @brief This class describes a mass rate constraint used to control a  well.
 */
template< typename WellConstraintType >
class MassConstraint : public WellConstraintType
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

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new Constraint object through the object catalog.
   */
  static string catalogName()
  {
    if constexpr ( std::is_same_v< WellConstraintType, InjectionConstraint > )    // special case
    {
      return "MassInjectionConstraint";
    }
    else   // default
    {
      return "MassProductionConstraint";
    }
  }
  ///@}

  /**
   * @name Getters / Setters
   */

  // Temp interface - tjb
  virtual ConstraintTypeId getControl() const override { return ConstraintTypeId::MASSRATE; };
  ///@}

  virtual bool checkViolation( WellConstraintBase const & currentConstraint, real64 const & currentTime ) const override;

protected:

  virtual void postInputInitialization() override;

};


} //namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLCONSTRAINT_HPP
