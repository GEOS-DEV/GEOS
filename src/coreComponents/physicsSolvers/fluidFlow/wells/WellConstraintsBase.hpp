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
static constexpr auto wellConstraintBase = "WellConstraintBase";
}
}

namespace constraintViewStruct
{

struct constraintValueKey
{
  /// String key for the well constraint value
  static constexpr char const * constraintValueString() { return "constraintValue"; }
  /// string key for constraint values entered table name
  static constexpr char const * constraintScheduleTableNameString() { return "constraintScheduleTableName"; }
};

struct surfaceConditionsKey
{
  /// String key for checking the rates at surface conditions
  static constexpr char const * useSurfaceConditionsString() { return "useSurfaceConditions"; }
  /// String key for the surface pressure
  static constexpr char const * surfacePressureString() { return "surfacePressure"; }
  /// String key for the surface temperature
  static constexpr char const * surfaceTemperatureString() { return "surfaceTemperature"; }
};

struct injectionStreamKey
{
  /// String key for the well injection stream
  static constexpr char const * injectionStreamString() { return "injectionStream"; }
  /// String key for the well injection temperature
  static constexpr char const * injectionTemperatureString() { return "injectionTemperature"; }
};

}

/**
 * @brief Validate the injection stream and temperature.
 * @param[in] injectionStream the injection stream vector
 * @param[in] injectionTemperature the injection temperature
 * @param[in] dataContext context for error messages
 */
template< typename T >
void validateInjectionStream( array1d< real64 > const & injectionStream,
                              real64 const & injectionTemperature,
                              std::string const & className,
                              T const & context,
                              dataRepository::DataContext const & dataContext )
{

  GEOS_THROW_IF( (injectionStream.empty()  && injectionTemperature >= 0) ||
                 (!injectionStream.empty() && injectionTemperature < 0),
                 className << " "  << dataContext << ": Both "
                           << constraintViewStruct::injectionStreamKey::injectionStreamString() << " and " << constraintViewStruct::injectionStreamKey::injectionTemperatureString()
                           << " must be specified for multiphase simulations",
                 InputError );

  if( !injectionStream.empty())
  {
    real64 sum = 0.0;
    for( localIndex ic = 0; ic < injectionStream.size(); ++ic )
    {
      GEOS_ERROR_IF( injectionStream[ic] < 0.0 || injectionStream[ic] > 1.0,
                     context.getWrapperDataContext( constraintViewStruct::injectionStreamKey::injectionStreamString() ) << ": Invalid injection stream" );
      sum += injectionStream[ic];
    }
    GEOS_THROW_IF( LvArray::math::abs( 1.0 - sum ) > std::numeric_limits< real64 >::epsilon(),
                   context.getWrapperDataContext( constraintViewStruct::injectionStreamKey::injectionStreamString() ) << ": Invalid injection stream",
                   InputError );
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
class WellConstraintBase : public dataRepository::Group
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


  /**
   * @brief Get name of constraint
   * @return constraint key
   */
  virtual std::string getConstraintKey( ) const = 0;


  /**
   * @brief Get the target bottom hole pressure value.
   * @return a value for the target bottom hole pressure
   */
  real64 getConstraintValue( real64 const & currentTime ) const
  {
    return m_constraintScheduleTable->evaluate( &currentTime );
  }

  ///@}

  // Phase constraint defintion keys
  constraintViewStruct::constraintValueKey viewKeysConstraintValue;

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


  /// Constraint value
  real64 m_constraintValue;

  void setNextDtFromTable( TableFunction const * table, real64 const currentTime, real64 & nextDt );

  /// Total rate table name
  string m_constraintScheduleTableName;

  /// Constraint values versus time
  TableFunction const * m_constraintScheduleTable;

};



} //namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLCONSTRAINT_HPP
