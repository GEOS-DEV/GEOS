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


#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLCONSTRAINT_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLCONSTRAINT_HPP

#include "common/format/EnumStrings.hpp"
#include "dataRepository/Group.hpp"
#include "functions/TableFunction.hpp"
#include "WellControls.hpp"
namespace geos
{
namespace dataRepository
{
namespace keys
{
static constexpr auto wellConstraint = "WellConstraint";
}
}

/*
   class WellContraint : public dataRepository::Group
   {
   public:

   virtual bool estimateWellSolution()=0;
   virtual void assembleConstraintEquation()=0;
   virtual void getBHP()=0;
   virtual void getPhaseRate(integer phase)=0;
   }
 */


/**
 * @class WellContraint
 * @brief This class describes a constraint used to control a well.
 */
class WellConstraint : public dataRepository::Group
{
public:

  /** Types of well controls
   * Used to specifiy a well's operating conditions
   */
  enum class Control : integer
  {
    BHP,  /**< The well operates at a specified bottom hole pressure (BHP) */
    PHASEVOLRATE, /**< The well operates at a specified phase volumetric flow rate */
    TOTALVOLRATE, /**< The well operates at a specified total volumetric flow rate */
    MASSRATE, /**<The well operates at a specified mass rate */
    UNINITIALIZED, /**< This is the current well control before postInputInitialization (needed to restart from file properly) */
  };


  /**
   * @name Constructor / Destructor
   */
  ///@{

  /**
   * @brief Constructor for WellControls Objects.
   * @param[in] name the name of this instantiation of WellControls in the repository
   * @param[in] parent the parent group of this instantiation of WellControls
   */
  explicit WellConstraint( string const & name, dataRepository::Group * const parent );


  /**
   * @brief Default destructor.
   */
  ~WellConstraint() override;

  /**
   * @brief Deleted default constructor.
   */
  WellConstraint() = delete;

  /**
   * @brief Deleted copy constructor.
   */
  WellConstraint( WellConstraint const & ) = delete;

  /**
   * @brief Deleted move constructor.
   */
  WellConstraint( WellConstraint && ) = delete;

  /**
   * @brief Deleted assignment operator.
   * @return a reference to a constraint object
   */
  WellConstraint & operator=( WellConstraint const & ) = delete;

  /**
   * @brief Deleted move operator.
   * @return a reference to a constraint object
   */
  WellConstraint & operator=( WellConstraint && ) = delete;

  ///@}

  /**
   * @name Getters / Setters
   */
  ///@{

  void setWellType( WellControls::Type wellType ) { m_wellType=wellType;}
  /**
   * @brief Get the target bottom hole pressure value.
   * @return a value for the target bottom hole pressure
   */
  real64 getConstraintValue( real64 const & currentTime ) const
  {
    return m_constraintScheduleTable->evaluate( &currentTime );
  }

  ///@}

  /**
   * @brief Is the constraint for an injector?
   * @return a boolean
   */
  bool isInjectionConstraint() const { return ( m_wellType == WellControls::Type::INJECTOR ); }

  /**
   * @brief Is the constraint for a producer?
   * @return a boolean
   */
  bool isProductionConstraint() const { return ( m_wellType == WellControls::Type::PRODUCER ); }

  /**
   * @brief Struct to serve as a container for variable strings and keys.
   * @struct viewKeyStruct
   */
  struct viewKeyStruct
  {
    /// string key for schedule table name
    static constexpr char const * constraintScheduleTableNameString() { return "constraintScheduleTableName"; }

  }
  /// ViewKey struct for the WellControls class
  viewKeysWellConstraint;

protected:

  virtual void postInputInitialization() override;

  /**
   * @brief set next time step based on tables intervals
   * @param[in] currentTime the current time
   * @param[inout] nextDt the time step
   */
  void setNextDtFromTables( real64 const currentTime, real64 & nextDt );

protected:


  void setNextDtFromTable( TableFunction const * table, real64 const currentTime, real64 & nextDt );

  /// Total rate table name
  string m_constraintScheduleTableName;

  /// Constraint values versus time
  TableFunction const * m_constraintScheduleTable;

private:
  /// Well type (as Type enum)
  WellControls::Type m_wellType;
};

/**
 * @class WellBHPContraint
 * @brief This class describes a pressure constraint used to control a well.
 */

class WellBHPConstraint : public WellConstraint
{
public:

  /** Types of well controls
   * Used to specifiy a well's operating conditions
   */
  enum class Control : integer
  {
    BHP,  /**< The well operates at a specified bottom hole pressure (BHP) */
    PHASEVOLRATE, /**< The well operates at a specified phase volumetric flow rate */
    TOTALVOLRATE, /**< The well operates at a specified total volumetric flow rate */
    MASSRATE, /**<The well operates at a specified mass rate */
    UNINITIALIZED, /**< This is the current well control before postInputInitialization (needed to restart from file properly) */
  };


  /**
   * @name Constructor / Destructor
   */
  ///@{

  /**
   * @brief Constructor for WellControls Objects.
   * @param[in] name the name of this instantiation of WellControls in the repository
   * @param[in] parent the parent group of this instantiation of WellControls
   */
  explicit WellBHPConstraint( string const & name, dataRepository::Group * const parent );


  /**
   * @brief Default destructor.
   */
  ~WellBHPConstraint() override;

  /**
   * @brief Deleted default constructor.
   */
  WellBHPConstraint() = delete;

  /**
   * @brief Deleted copy constructor.
   */
  WellBHPConstraint( WellBHPConstraint const & ) = delete;

  /**
   * @brief Deleted move constructor.
   */
  WellBHPConstraint( WellBHPConstraint && ) = delete;

  /**
   * @brief Deleted assignment operator.
   * @return a reference to a constraint object
   */
  WellBHPConstraint & operator=( WellBHPConstraint const & ) = delete;

  /**
   * @brief Deleted move operator.
   * @return a reference to a constraint object
   */
  WellBHPConstraint & operator=( WellBHPConstraint && ) = delete;

  ///@}

  /**
   * @name Getters / Setters
   */
  ///@{

  ///@}

  /**
   * @brief Struct to serve as a container for variable strings and keys.
   * @struct viewKeyStruct
   */
  struct viewKeyStruct
  {
    /// String key for the well reference elevation (for BHP control)
    static constexpr char const * refElevString() { return "referenceElevation"; }
    /// String key for the well type
    static constexpr char const * typeString() { return "type"; }
    /// String key for the well input control
    static constexpr char const * inputControlString() { return "control"; }
    /// String key for the well current control
    static constexpr char const * currentControlString() { return "currentControl"; }
    /// String key for the well target BHP
    static constexpr char const * targetBHPString() { return "targetBHP"; }
    /// String key for the well target rate
    static constexpr char const * targetTotalRateString() { return "targetTotalRate"; }
    /// String key for the well target phase rate
    static constexpr char const * targetPhaseRateString() { return "targetPhaseRate"; }
    /// String key for the well target phase name
    static constexpr char const * targetPhaseNameString() { return "targetPhaseName"; }
    /// String key for the well target phase name
    static constexpr char const * targetMassRateString() { return "targetMassRate"; }
    /// String key for the well injection stream
    static constexpr char const * injectionStreamString() { return "injectionStream"; }
    /// String key for the well injection temperature
    static constexpr char const * injectionTemperatureString() { return "injectionTemperature"; }
    /// String key for checking the rates at surface conditions
    static constexpr char const * useSurfaceConditionsString() { return "useSurfaceConditions"; }
    /// String key for the surface pressure
    static constexpr char const * surfacePressureString() { return "surfacePressure"; }
    /// String key for the surface temperature
    static constexpr char const * surfaceTemperatureString() { return "surfaceTemperature"; }
    /// string key for total rate table name
    static constexpr char const * targetTotalRateTableNameString() { return "targetTotalRateTableName"; }
    /// string key for phase rate table name
    static constexpr char const * targetPhaseRateTableNameString() { return "targetPhaseRateTableName"; }
    /// string key for mass rate table name
    static constexpr char const * targetMassRateTableNameString() { return "targetMassRateTableName"; }
    /// string key for BHP table name
    static constexpr char const * targetBHPTableNameString() { return "targetBHPTableName"; }


  }
  /// ViewKey struct for the WellControls class
  viewKeysWellConstraint;

protected:

  virtual void postInputInitialization() override;

private:

  /// Target bottom hole pressure value
  real64 m_targetBHP;

};

/**
 * @class WellBHPContraint
 * @brief This class describes a pressure constraint used to control a well.
 */

class WellPhaseRateConstraint : public WellConstraint
{
public:

  /** Types of well controls
   * Used to specifiy a well's operating conditions
   */
  enum class Control : integer
  {
    BHP,  /**< The well operates at a specified bottom hole pressure (BHP) */
    PHASEVOLRATE, /**< The well operates at a specified phase volumetric flow rate */
    TOTALVOLRATE, /**< The well operates at a specified total volumetric flow rate */
    MASSRATE, /**<The well operates at a specified mass rate */
    UNINITIALIZED, /**< This is the current well control before postInputInitialization (needed to restart from file properly) */
  };


  /**
   * @name Constructor / Destructor
   */
  ///@{

  /**
   * @brief Constructor for WellControls Objects.
   * @param[in] name the name of this instantiation of WellControls in the repository
   * @param[in] parent the parent group of this instantiation of WellControls
   */
  explicit WellPhaseRateConstraint( string const & name, dataRepository::Group * const parent );


  /**
   * @brief Default destructor.
   */
  ~WellPhaseRateConstraint() override;

  /**
   * @brief Deleted default constructor.
   */
  WellPhaseRateConstraint() = delete;

  /**
   * @brief Deleted copy constructor.
   */
  WellPhaseRateConstraint( WellPhaseRateConstraint const & ) = delete;

  /**
   * @brief Deleted move constructor.
   */
  WellPhaseRateConstraint( WellPhaseRateConstraint && ) = delete;

  /**
   * @brief Deleted assignment operator.
   * @return a reference to a constraint object
   */
  WellBHPConstraint & operator=( WellBHPConstraint const & ) = delete;

  /**
   * @brief Deleted move operator.
   * @return a reference to a constraint object
   */
  WellPhaseRateConstraint & operator=( WellPhaseRateConstraint && ) = delete;

  ///@}

  /**
   * @name Getters / Setters
   */
  ///@{

  ///@}

  /**
   * @brief Struct to serve as a container for variable strings and keys.
   * @struct viewKeyStruct
   */
  struct viewKeyStruct
  {
    /// String key for the well reference elevation (for BHP control)
    static constexpr char const * refElevString() { return "referenceElevation"; }
    /// String key for the well type
    static constexpr char const * typeString() { return "type"; }
    /// String key for the well input control
    static constexpr char const * inputControlString() { return "control"; }
    /// String key for the well current control
    static constexpr char const * currentControlString() { return "currentControl"; }
    /// String key for the well target BHP
    static constexpr char const * targetBHPString() { return "targetBHP"; }
    /// String key for the well target rate
    static constexpr char const * targetTotalRateString() { return "targetTotalRate"; }
    /// String key for the well target phase rate
    static constexpr char const * targetPhaseRateString() { return "targetPhaseRate"; }
    /// String key for the well target phase name
    static constexpr char const * targetPhaseNameString() { return "targetPhaseName"; }
    /// String key for the well target phase name
    static constexpr char const * targetMassRateString() { return "targetMassRate"; }
    /// String key for the well injection stream
    static constexpr char const * injectionStreamString() { return "injectionStream"; }
    /// String key for the well injection temperature
    static constexpr char const * injectionTemperatureString() { return "injectionTemperature"; }
    /// String key for checking the rates at surface conditions
    static constexpr char const * useSurfaceConditionsString() { return "useSurfaceConditions"; }
    /// String key for the surface pressure
    static constexpr char const * surfacePressureString() { return "surfacePressure"; }
    /// String key for the surface temperature
    static constexpr char const * surfaceTemperatureString() { return "surfaceTemperature"; }
    /// string key for total rate table name
    static constexpr char const * targetTotalRateTableNameString() { return "targetTotalRateTableName"; }
    /// string key for phase rate table name
    static constexpr char const * targetPhaseRateTableNameString() { return "targetPhaseRateTableName"; }
    /// string key for mass rate table name
    static constexpr char const * targetMassRateTableNameString() { return "targetMassRateTableName"; }
    /// string key for BHP table name
    static constexpr char const * targetBHPTableNameString() { return "targetBHPTableName"; }


  }
  /// ViewKey struct for the WellControls class
  viewKeysWellConstraint;

protected:

  virtual void postInputInitialization() override;

private:
  /// Name of the targeted phase
  string m_targetPhaseName;

  /// Target phase rate
  real64 m_targetPhaseRate;

  /// Flag to decide whether rates are controlled at rates or surface conditions
  integer m_useSurfaceConditions;
};

} //namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLCONSTRAINT_HPP
