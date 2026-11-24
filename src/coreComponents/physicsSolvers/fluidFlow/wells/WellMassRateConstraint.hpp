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
 * @file WellMassRateConstraint.hpp
 */


#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLMASSRATECONSTRAINT_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLMASSRATECONSTRAINT_HPP

#include "common/format/EnumStrings.hpp"
#include "dataRepository/Group.hpp"
#include "functions/TableFunction.hpp"
#include "WellConstraintsBase.hpp"
namespace geos
{

/**
 * @class MassRateConstraint
 * @brief This class describes a mass rate constraint used to control a  well.
 */

class MassRateConstraint : public WellConstraintBase
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
  explicit MassRateConstraint( string const & name, dataRepository::Group * const parent );


  /**
   * @brief Default destructor.
   */
  ~MassRateConstraint() override;

  /**
   * @brief Deleted default constructor.
   */
  MassRateConstraint() = delete;

  /**
   * @brief Deleted copy constructor.
   */
  MassRateConstraint( MassRateConstraint const & ) = delete;

  /**
   * @brief Deleted move constructor.
   */
  MassRateConstraint( MassRateConstraint && ) = delete;

  /**
   * @brief Deleted assignment operator.
   * @return a reference to a constraint object
   */
  MassRateConstraint & operator=( MassRateConstraint const & ) = delete;

  /**
   * @brief Deleted move operator.
   * @return a reference to a constraint object
   */
  MassRateConstraint & operator=( MassRateConstraint && ) = delete;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new Constraint object through the object catalog.
   */
  static string catalogName()
  {
    return "MassRateConstraint";
  }
  ///@}

  struct viewKeyStruct
  {
    /// String key for the well target  rate
    static constexpr char const * massRateString() { return "massRate"; }
  };

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

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLMASSRATECONSTRAINT_HPP
