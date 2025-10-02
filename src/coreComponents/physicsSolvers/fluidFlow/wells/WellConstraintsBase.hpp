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
  virtual ConstraintTypeId getControl() const = 0;

  /**
   * @brief Get name of constraint
   * @return constraint key
   */
  virtual std::string getConstraintKey( ) const = 0;

  /**
   * @brief Defines whether the constraint should be evaluated or not
   * @brief Some workflows require the well model to define a constraint
   * @brief of similar type to user defined constraints. For example,
   * @brief rate constraints to evaluated WHP constraints.
   * @return true if the constraint is active, false otherwise
   */
  bool isConstraintActive( ) { return m_isConstraintActive; }

  /**
   * @brief Sets constraint active status
   * @param[in] constraintActive true if the constraint is active, false otherwise
   */
  bool setConstraintActive( bool const & constraintActive ) { return m_isConstraintActive=constraintActive; }

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

  // Quantities computed from well constraint solve with this boundary condition
  // This needs to be somewhere else tjb
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
  bool m_isConstraintActive;

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

/**
 * @class ProductionConstraint
 * @brief This class describes constraint used to control a production well.
 */

class ProductionConstraint : public WellConstraintBase
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

protected:

  virtual void postInputInitialization() override;


};

/**
 * @class InjectionConstraint
 * @brief This class describes   constraint used to control a injection well.
 */

class InjectionConstraint : public WellConstraintBase
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
  explicit InjectionConstraint( string const & name, dataRepository::Group * const parent );

  /**
   * @brief Default destructor.
   */
  ~InjectionConstraint() override;

  /**
   * @brief Deleted default constructor.
   */
  InjectionConstraint() = delete;

  /**
   * @brief Deleted copy constructor.
   */
  InjectionConstraint( InjectionConstraint const & ) = delete;

  /**
   * @brief Deleted move constructor.
   */
  InjectionConstraint( InjectionConstraint && ) = delete;

  /**
   * @brief Deleted assignment operator.
   * @return a reference to a constraint object
   */
  InjectionConstraint & operator=( InjectionConstraint const & ) = delete;

  /**
   * @brief Deleted move operator.
   * @return a reference to a constraint object
   */
  InjectionConstraint & operator=( InjectionConstraint && ) = delete;

  ///@}

  struct injectionStreamKey
  {
    /// String key for the well injection stream
    static constexpr char const * injectionStreamString() { return "injectionStream"; }
    /// String key for the well injection temperature
    static constexpr char const * injectionTemperatureString() { return "injectionTemperature"; }
  };

  /**
   * @brief Const accessor for the composition of the injection stream
   * @return a global component fraction vector
   */
  arrayView1d< real64 const > getInjectionStream() const { return m_injectionStream; }

  /**
   * @brief Const accessor for the temperature of the injection stream
   * @return the temperature of the injection stream
   */
  real64 getInjectionTemperature() const { return m_injectionTemperature; }

protected:

  virtual void postInputInitialization() override;

  void validateInjectionStream();
private:

  /// Vector with global component fractions at the injector
  array1d< real64 > m_injectionStream;

  /// Temperature at the injector
  real64 m_injectionTemperature;

};


} //namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLCONSTRAINT_HPP
