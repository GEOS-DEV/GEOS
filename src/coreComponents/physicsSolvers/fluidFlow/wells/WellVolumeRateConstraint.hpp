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



/**
 * @class VolumeRateConstraint
 * @brief This class describes a volume rate constraint used to control a well.
 */

class VolumeRateConstraint : public WellConstraintBase
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
  explicit VolumeRateConstraint( string const & name, dataRepository::Group * const parent );


  /**
   * @brief Default destructor.
   */
  ~VolumeRateConstraint() override;

  /**
   * @brief Deleted default constructor.
   */
  VolumeRateConstraint() = delete;

  /**
   * @brief Deleted copy constructor.
   */
  VolumeRateConstraint( VolumeRateConstraint const & ) = delete;

  /**
   * @brief Deleted move constructor.
   */
  VolumeRateConstraint( VolumeRateConstraint && ) = delete;

  /**
   * @brief Deleted assignment operator.
   * @return a reference to a constraint object
   */
  VolumeRateConstraint & operator=( VolumeRateConstraint const & ) = delete;

  /**
   * @brief Deleted move operator.
   * @return a reference to a constraint object
   */
  VolumeRateConstraint & operator=( VolumeRateConstraint && ) = delete;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new Constraint object through the object catalog.
   */
  static string catalogName()
  {
    return "VolumeRateConstraint";
  }
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

  // Temp interface - tjb
  virtual ConstraintTypeId getControl() const override { return ConstraintTypeId::TOTALVOLRATE; };
  ///@}

  virtual bool checkViolation( WellConstraintBase const & currentConstraint, real64 const & currentTime ) const override;
protected:

  virtual void postInputInitialization() override;

};


} //namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLCONSTRAINT_HPP
