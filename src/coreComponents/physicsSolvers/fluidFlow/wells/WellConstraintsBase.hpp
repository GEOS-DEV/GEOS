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
 * @brief Register fields required to define surface conditions for constraint
 * @param[in] useSurfaceConditions 0 - use reservoir conditions, 1 use specified P &  T
 * @param[in] surfacePres surface pressure
 * @param[in] surfaceTemp surface pressure
 */
template< typename T >
void registerSurfaceConditions( integer & useSurfaceConditions,
                                real64 & surfacePres,
                                real64 & surfaceTemp,
                                T & context )
{
  context.registerWrapper( constraintViewStruct::surfaceConditionsKey::useSurfaceConditionsString(), &useSurfaceConditions ).
    setDefaultValue( 0 ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "Flag to specify whether rates are checked at surface or reservoir conditions.\n"
                    "Equal to 1 for surface conditions, and to 0 for reservoir conditions" );

  context.registerWrapper( constraintViewStruct::surfaceConditionsKey::surfacePressureString(), &surfacePres ).
    setDefaultValue( 0 ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "Surface pressure used to compute volumetric rates when surface conditions are used [Pa]" );

  context.registerWrapper( constraintViewStruct::surfaceConditionsKey::surfaceTemperatureString(), &surfaceTemp ).
    setDefaultValue( 0 ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "Surface temperature used to compute volumetric rates when surface conditions are used [K]" );
}

/**
 * @brief Register fields required to define an injection stream.
 * @param[in] injectionStream the injection stream vector
 * @param[in] injectionTemperature the injection temperature
 * @param[in] context class needing fields
 */
template< typename T >
void registerInjectionStream( array1d< real64 > & injectionStream,
                              real64 & injectionTemperature,
                              T & context )
{
  context.registerWrapper( constraintViewStruct::injectionStreamKey::injectionStreamString(), &injectionStream ).
    setDefaultValue( -1 ).
    setSizedFromParent( 0 ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "Global component densities of the injection stream [moles/m^3 or kg/m^3]" );

  context.registerWrapper( constraintViewStruct::injectionStreamKey::injectionTemperatureString(), &injectionTemperature ).
    setDefaultValue( -1 ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "Temperature of the injection stream [K]" );
}
/**
 * @brief Validate the surface conditions
 * @param[in] useSurfaceConditions 0 - use reservoir conditions, 1 use specified P &  T
 * @param[in] className owner of fields
 * @param[in] dataContext context for error messages
 */
template< typename T >
void validateSurfaceConditions( integer useSurfaceConditions,
                                std::string const & className,
                                T const & context )
{
  GEOS_THROW_IF( useSurfaceConditions != 0 &&  useSurfaceConditions != 1,
                 className << " "  << context.getDataContext()  << ": The flag to select surface/reservoir conditions must be equal to 0 or 1",
                 InputError );
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
                              T const & context )
{
  GEOS_THROW_IF( (injectionStream.empty()  && injectionTemperature >= 0) ||
                 (!injectionStream.empty() && injectionTemperature < 0),
                 className << " "  << context.getDataContext() << ": Both "
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

  // Temp interface - tjb
  virtual WellControls::Control getControl() const = 0;

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

  //
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
   * @brief Getter for the mass rate
   * @return mass rate
   */
  real64 massRate() const { return m_massRate; }

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


  /// Constraint value
  real64 m_constraintValue;

  void setNextDtFromTable( TableFunction const * table, real64 const currentTime, real64 & nextDt );

  /// Total rate table name
  string m_constraintScheduleTableName;

  /// Constraint values versus time
  TableFunction const * m_constraintScheduleTable;

  // Quantities computed from well constraint solve with this boundary condition

  // botton hole pressure
  real64 m_BHP;

  // phase rates
  array1d< real64 >  m_phaseVolumeRates;

  // total volume rate
  real64 m_totalVolumeRate;

  // mass rate
  real64 m_massRate;

};



} //namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLCONSTRAINT_HPP
