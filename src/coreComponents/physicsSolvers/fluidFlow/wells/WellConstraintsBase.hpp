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
 * @file WellConstraintBase.hpp
 */


#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLCONSTRAINTBASE_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLCONSTRAINTBASE_HPP

#include "common/format/EnumStrings.hpp"

#include "functions/TableFunction.hpp"
#include "dataRepository/Group.hpp"
namespace geos
{



enum class ConstraintTypeId : integer
{
  BHP,    /**< The well operates at a specified bottom hole pressure (BHP) */
  PHASEVOLRATE,   /**< The well operates at a specified phase volumetric flow rate */
  TOTALVOLRATE,   /**< The well operates at a specified total volumetric flow rate */
  MASSRATE,   /**<The well operates at a specified mass rate */
  WHP,  /**< The well operates at a specified wellhead   pressure (WHP) */
  LIQUIDRATE,  /**< The well operates at a specified liquid rate */
  UNINITIALIZED,   /**< This is the current well control before postInputInitialization (needed to restart from file properly) */
};

/**
 * @class WellContraint
 * @brief This class describes a constraint used to control a well.
 */
class WellConstraintBase : public dataRepository::Group
{
public:
  friend class WellControls;

  /// Catalog interface specific to WellConstraintBase
  using CatalogInterface = dataRepository::CatalogInterface< WellConstraintBase, string const &, dataRepository::Group * const >;

  /// Get the singleton catalog for WellConstraintBase
  static CatalogInterface::CatalogType & getCatalog();

  /**
   * @brief function to return the catalog name of the derived class
   * @return a string that contains the catalog name of the derived class
   */
  virtual string getCatalogName() const = 0;

  /**
   * @name Getters / Setters
   */
  ///@{

  /**
   * @brief Constructor for WellControls Objects.
   * @param[in] name the name of this instantiation of WellControls in the repository
   * @param[in] parent the parent group of this instantiation of WellControls
   */
  explicit WellConstraintBase( string const & name, dataRepository::Group * const parent );


  /**
   * @brief Default destructor.
   */
  ~WellConstraintBase() override;

  /**
   * @brief Deleted default constructor.
   */
  WellConstraintBase() = delete;

  /**
   * @brief Deleted copy constructor.
   */
  WellConstraintBase( WellConstraintBase const & ) = delete;

  /**
   * @brief Deleted move constructor.
   */
  WellConstraintBase( WellConstraintBase && ) = delete;

  /**
   * @brief Deleted assignment operator.
   * @return a reference to a constraint object
   */
  WellConstraintBase & operator=( WellConstraintBase const & ) = delete;

  /**
   * @brief Deleted move operator.
   * @return a reference to a constraint object
   */
  WellConstraintBase & operator=( WellConstraintBase && ) = delete;

  ///@}


  /**
   * @name Getters / Setters
   */
  ///@{

  // Temp interface - tjb
  virtual ConstraintTypeId getControl() const = 0;

  /**
   * @brief Defines whether the constraint should be evaluated or not
   * @brief Some workflows require the well model to define a constraint
   * @brief of similar type to user defined constraints. For example,
   * @brief rate constraints to evaluated WHP constraints.
   * @return true if the constraint is active, false otherwise
   */
  bool isConstraintActive( ) const { return m_isConstraintActive; }

  /**
   * @brief Sets constraint active status
   * @param[in] constraintActive true if the constraint is active, false otherwise
   */
  void setConstraintActive( bool const & constraintActive ) { m_isConstraintActive=constraintActive; }

  /**
   * @brief Sets constraint value
   * @param[in] constraint value
   */
  void setConstraintValue( real64 const & constraintValue )
  {
    m_constraintValue = constraintValue;
  }

  /**
   * @brief Get the target bottom hole pressure value.
   * @return a value for the target bottom hole pressure
   */
  real64 getConstraintValue( real64 const & currentTime ) const
  {
    if( m_constraintScheduleTableName.empty() )
    {
      return m_rateSign*m_constraintValue;
    }

    return m_rateSign*m_constraintScheduleTable->evaluate( &currentTime );
  }

  ///@}

  /**
   * @brief Struct to serve as a container for variable strings and keys.
   * @struct viewKeyStruct
   */
  struct viewKeyStruct
  {
    /// string key for schedule table name
    static constexpr char const * constraintScheduleTableNameString() { return "constraintScheduleTableName"; }

    /// String key for the well constraint value
    static constexpr char const * constraintValueString() { return "value"; }

    /// String key for the well constraint active flag
    static constexpr char const * constraintActiveString() { return "constraintActive"; }

  }
  /// ViewKey struct for the WellControls class
  viewKeysWellConstraint;

  // Quantities computed from well constraint solve with this boundary condition
  // Until we have a more general interface for constraints to return these quantities,
  // we will store them in the constraint object itself.
  // This is not ideal but it is a temporary solution to avoid having to solve the constraint multiple times in the well solver and in the
  // test.
  void setBHP( real64 bhp ){ m_BHP=bhp;};
  void setPhaseVolumeRates( array1d< real64 > const & phaseVolumeRates ) { m_phaseVolumeRates = phaseVolumeRates; };
  void setTotalVolumeRate( real64 totalVolumeRate ){ m_totalVolumeRate = totalVolumeRate; };
  void setMassRate( real64 massRate ){ m_massRate = massRate; };

  /**
   * @brief Getter for the bottom hole pressure
   * @return bottom hole pressure
   */
  real64 bottomHolePressure() const { return m_BHP; }

  /**
   * @brief Getter for the phase volume rates
   * @return an arrayView1d storing the phase volume rates
   */
  arrayView1d< real64 const > phaseVolumeRates() const { return m_phaseVolumeRates; }

  /**
   * @brief Getter for the total volume rate
   * @return mass rate
   */
  real64 totalVolumeRate() const { return m_totalVolumeRate; }

  /**
   * @brief Getter for the liquid rate
   * @return liquid rate
   */
  real64 liquidRate() const { return m_liquidRate; }

  /**
   * @brief Getter for the mass rate
   * @return mass rate
   */
  real64 massRate() const { return m_massRate; }

  // endof This needs to be somewhere else tjb
  /**
   * @brief Check if this constraint is violated
   * @return true if limiting constraint, false otherwise
   */
  virtual bool checkViolation( WellConstraintBase const & currentConstraint, real64 const & currentTime ) const = 0;

protected:

  virtual void postInputInitialization() override;

  /**
   * @brief set next time step based on tables intervals
   * @param[in] currentTime the current time
   * @param[inout] nextDt the time step
   */
  void setNextDtFromTables( real64 const currentTime, real64 & nextDt );


protected:

  /// Constraint status
  integer m_isConstraintActive;

  /// Flag to indicate whether a schedule table should be generated for constraint value;
  bool m_useScheduleTable;

  /// Constraint value
  real64 m_constraintValue;

  void setNextDtFromTable( TableFunction const * table, real64 const currentTime, real64 & nextDt );

  /// Constraint schedule table name
  string m_constraintScheduleTableName;

  /// Constraint values versus time
  TableFunction const * m_constraintScheduleTable;

  // Quantities computed from well constraint solve with this boundary condition

  // botton hole pressure
  real64 m_BHP;

  // phase rates
  array1d< real64 >  m_phaseVolumeRates;

  // liquid rate
  real64 m_liquidRate;

  // total volume rate
  real64 m_totalVolumeRate;

  // mass rate
  real64 m_massRate;

  /// Rate sign. +1 for injector, -1 for producer
  real64 m_rateSign;
};


} //namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLCONSTRAINTBASE_HPP
