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


#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLCONTROLS_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLCONTROLS_HPP

#include "common/format/EnumStrings.hpp"
#include "dataRepository/Group.hpp"
#include "functions/TableFunction.hpp"
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"

#include "physicsSolvers/fluidFlow/wells/WellBHPConstraints.hpp"
#include "physicsSolvers/fluidFlow/wells/WellTotalVolRateConstraints.hpp"
#include "physicsSolvers/fluidFlow/wells/WellPhaseRateConstraints.hpp"
#include "physicsSolvers/fluidFlow/wells/WellMassRateConstraints.hpp"
#include "physicsSolvers/fluidFlow/wells/WellLiquidRateConstraints.hpp"
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidBase.hpp"
namespace geos
{
namespace dataRepository
{
namespace keys
{
static constexpr auto wellControls = "WellControls";
}
}


/**
 * @class WellControls
 * @brief This class describes the controls used to operate a well.
 */
class WellControls : public dataRepository::Group
{
public:

  /** Type of wells
   * Either producer or injector.
   */
  enum class Type : integer
  {
    PRODUCER,  /**< A production well */
    INJECTOR   /**< An injection well */
  };

  /** Status of wells
   * Either open or closed
   */
  enum class Status : integer
  {
    OPEN,  /**< flowing well */
    CLOSED   /**< shutin well */
  };

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
  explicit WellControls( string const & name, dataRepository::Group * const parent );


  /**
   * @brief Default destructor.
   */
  ~WellControls() override;

  /**
   * @brief Deleted default constructor.
   */
  WellControls() = delete;

  /**
   * @brief Deleted copy constructor.
   */
  WellControls( WellControls const & ) = delete;

  /**
   * @brief Deleted move constructor.
   */
  WellControls( WellControls && ) = delete;

  /**
   * @brief Deleted assignment operator.
   * @return a reference to a perforation object
   */
  WellControls & operator=( WellControls const & ) = delete;

  /**
   * @brief Deleted move operator.
   * @return a reference to a perforation object
   */
  WellControls & operator=( WellControls && ) = delete;

  ///@}

  /**
   * @brief Create a new geometric object (box, plane, etc) as a child of this group.
   * @param childKey the catalog key of the new geometric object to create
   * @param childName the name of the new geometric object in the repository
   * @return the group child
   */
  virtual Group * createChild( string const & childKey, string const & childName ) override;
  /// Expand catalog for schema generation

  virtual void expandObjectCatalogs() override;


  /**
   * @name Getters / Setters
   */
  ///@{


  /**
   * @brief Get the control type for the well.
   * @return the Control enum enforced at the well
   */
  Control getControl() const { return m_currentControl; }

  /**
   * @brief Set the control type for the well.
   * @param[in] newControl type
   */
  void setControl( Control const & newControl )  {  m_currentControl = newControl; }

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
   * @brief Returns the target bottom hole pressure value.
   * @param[in] targetTime time at which to evaluate the constraint
   * @return the injector maximum bottom hole pressure or producer minimum bottom hole pressure
   */
  real64 getTargetBHP( real64 const & targetTime ) const
  {
    if( isProducer())
    {
      return m_minBHPConstraint->getConstraintValue( targetTime );
    }
    return m_maxBHPConstraint->getConstraintValue( targetTime );
  }

  /**
   * @brief Const accessor for the temperature of the injection stream
   * @return the temperature of the injection stream
   */
  real64 getInjectionTemperature() const
  {
    real64 injectionTemperature = 0.0;
    for( auto c : m_injectionRateConstraintList )
      if( c->isConstraintActive())
      {
        injectionTemperature =  dynamic_cast< InjectionConstraint * >(c)->getInjectionTemperature();
        break;
      }
    return injectionTemperature;
  }
  /**
   * @brief Const accessor for the  injection stream
   * @return the injection stream
   */
  arrayView1d< real64 const > getInjectionStream() const
  {
    arrayView1d< real64 const > injectionStream;
    for( auto c : m_injectionRateConstraintList )
      if( c->isConstraintActive())
      {
        injectionStream =  dynamic_cast< InjectionConstraint * >(c)->getInjectionStream();
        break;
      }
    return injectionStream;
  }
  /**
   * @brief Return the reference elvation where pressure constraint is measured
   * @return  vertical location of constraint
   */
  real64 getReferenceElevation() const
  {
    if( isProducer () )
    {
      return getMinBHPConstraint()->getReferenceElevation();
    }
    return getMaxBHPConstraint()->getReferenceElevation();
  }
  /**
   * @brief Getter for the flag specifying whether we check rates at surface or reservoir conditions
   * @return 1 if we use surface conditions, and 0 otherwise
   */
  integer useSurfaceConditions() const { return m_useSurfaceConditions; }

  /**
   * @brief Getter for the reservoir region associated with reservoir volume constraint
   * @return name of reservoir region
   */
  string referenceReservoirRegion() const { return m_referenceReservoirRegion; }

  /**
   * @brief Getter for the surface pressure when m_useSurfaceConditions == 1
   * @return the surface pressure
   */
  const real64 & getSurfacePressure() const { return m_surfacePres; }

  /**
   * @brief Getter for the surface temperature when m_useSurfaceConditions == 1
   * @return the surface temperature
   */
  const real64 & getSurfaceTemperature() const { return m_surfaceTemp; }

  /**
   * @brief Is the well an injector?
   * @return a boolean
   */
  bool isInjector() const { return ( m_type == Type::INJECTOR ); }

  /**
   * @brief Is the well a producer?
   * @return a boolean
   */
  bool isProducer() const { return ( m_type == Type::PRODUCER ); }

  /**
   * @brief Is the well open (or shut) at currentTime, status initalized in WellSolverBase::implicitStepSetup
   * @return a boolean
   */
  bool isWellOpen() const;

  void setWellState( bool open );
  bool getWellState() const;


  void setConstraintSwitch( bool constraintSwitch );
  bool getConstraintSwitch() const;

  void setCurrentConstraint( WellConstraintBase * currentConstraint ) { m_currentConstraint = currentConstraint;}
  WellConstraintBase *  getCurrentConstraint() { return m_currentConstraint; }
  WellConstraintBase const *  getCurrentConstraint() const { return m_currentConstraint; }


  /**
   * @brief Getter for the flag to enable crossflow
   * @return the flag deciding whether crossflow is allowed or not
   */
  bool isCrossflowEnabled() const { return m_isCrossflowEnabled; }

  /**
   * @brief Getter for the initial pressure coefficient
   * @return the initial pressure coefficient
   */
  real64 getInitialPressureCoefficient() const { return m_initialPressureCoefficient; }

  /**
   * @brief set next time step based on tables intervals
   * @param[in] currentTime the current time
   * @param[inout] nextDt the time step
   */
  void setNextDtFromTables( real64 const & currentTime, real64 & nextDt );


  /**
   * @brief setter for multi fluid separator
   * @param[in] fluidSeparatorPtr single or multiphase separator
   */
  void setFluidSeparator( std::unique_ptr< constitutive::ConstitutiveBase > fluidSeparatorPtr )  {  m_fluidSeparatorPtr = std::move( fluidSeparatorPtr );}
  /**
   * @brief Getter for multi fluid separator
   * @return reference to separator
   */
  constitutive::MultiFluidBase & getMultiFluidSeparator()  { return dynamicCast< constitutive::MultiFluidBase & >( *m_fluidSeparatorPtr ); }

  /**
   * @brief Getter for single fluid separator
   * @return reference to separator
   */
  constitutive::SingleFluidBase & getSingleFluidSeparator()  { return dynamicCast< constitutive::SingleFluidBase & >( *m_fluidSeparatorPtr ); }


  /**
   * @brief Getter for the reservoir average pressure when m_useSurfaceConditions == 0
   * @return the pressure
   */
  real64 getRegionAveragePressure() const { return m_regionAveragePressure; }

  /**
   * @brief Set the reservoir average pressure when m_useSurfaceConditions == 0
   * @param[in] regionAveragePressure value for pressure
   */
  void setRegionAveragePressure( real64 regionAveragePressure ) { m_regionAveragePressure = regionAveragePressure; }

  /**
   * @brief Getter for the reservoir average temperature when m_useSurfaceConditions == 0
   * @return the temperature
   */
  real64 getRegionAverageTemperature() const { return m_regionAverageTemperature; }

  /**
   * @brief Set the reservoir average temperature when m_useSurfaceConditions == 0
   * @param[in] regionAverageTemperature value for temperature
   */
  void setRegionAverageTemperature( real64 regionAverageTemperature ) { m_regionAverageTemperature = regionAverageTemperature; }

  /**
   * @brief Set well status from time and internal action, eg. all perfs closed
   * @param[in] currentTime the current time
   * @param[in] status
   */
  void setWellStatus ( real64 const & currentTime, WellControls::Status status );

  /**
   * @brief Is the well open (or shut) based on internal action
   * @return a Status
   */
  WellControls::Status getWellStatus () const { return m_wellStatus; }
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

    /// String key for the well current control
    static constexpr char const * currentControlString() { return "currentControl"; }

    /// String key for checking the rates at surface conditions
    static constexpr char const * useSurfaceConditionsString() { return "useSurfaceConditions"; }
    /// String key for reference reservoir region
    static constexpr char const * referenceReservoirRegionString() { return "referenceReservoirRegion"; }
    /// String key for the surface pressure
    static constexpr char const * surfacePressureString() { return "surfacePressure"; }
    /// String key for the surface temperature
    static constexpr char const * surfaceTemperatureString() { return "surfaceTemperature"; }

    /// string key for status table name
    static constexpr char const * statusTableNameString() { return "statusTableName"; }
    /// string key for perforation status table name
    static constexpr char const * perfStatusTableNameString() { return "perfStatusTableName"; }
    /// string key for the crossflow flag
    static constexpr char const * enableCrossflowString() { return "enableCrossflow"; }
    /// string key for the initial pressure coefficient
    static constexpr char const * initialPressureCoefficientString() { return "initialPressureCoefficient"; }
    /// string key for the esitmate well solution flag
    static constexpr char const * estimateWellSolutionString() { return "estimateWellSolution"; }

    /// string key for the minimum BHP presssure for a producer
    static constexpr char const * minimumBHPConstraintString() { return "MinimumBHPConstraint"; }
    /// string key for the maximum BHP presssure for a injection
    static constexpr char const * maximumBHPConstraintString() { return "MaximumBHPConstraint"; }
    /// string key for the maximum phase rate for a producer
    static constexpr char const * phaseProductionConstraintString() { return "PhaseProductionConstraint"; }
    /// string key for the maximum phase rate for a injection
    static constexpr char const * phaseInjectionConstraintString() { return "PhaseInjectionConstraint"; }
    /// string key for the maximum volume rate for a producer
    static constexpr char const * totalVolProductionConstraintString() { return "TotalVolProductionConstraint"; }
    /// string key for the maximum volume rate for a injector
    static constexpr char const * totalVolInjectionConstraintString() { return "TotalVolInjectionConstraint"; }
    /// string key for the maximum mass rate for a producer
    static constexpr char const * massProductionConstraintString() { return "massProductionConstraint"; }
    /// string key for the maximum mass rate for a injector
    static constexpr char const * massInjectionConstraintString() { return "massInjectionConstraint"; }
    /// string key for the liquid rate for a producer
    static constexpr char const * liquidProductionConstraintString() { return "liquidProductionConstraint"; }
  }
  /// ViewKey struct for the WellControls class
  viewKeysWellControls;

  static void setNextDtFromTable( TableFunction const * table, real64 const currentTime, real64 & nextDt );

  /**
   * @brief Create a constraint
   * @tparam ConstraintType the type of constraint to create
   * @param[in] constraintName name to assign to the constraint
   */
  template< typename ConstraintType > void createConstraint ( string const & constraintName );


  /**
   * @brief Getters for constraints
   */
  MinimumBHPConstraint * getMinBHPConstraint() { return m_minBHPConstraint; };
  MinimumBHPConstraint * getMinBHPConstraint() const { return m_minBHPConstraint; };
  MaximumBHPConstraint * getMaxBHPConstraint() { return m_maxBHPConstraint; };
  MaximumBHPConstraint * getMaxBHPConstraint() const { return m_maxBHPConstraint; };
  // Lists of rate constraints
  std::vector< WellConstraintBase * >  getProdRateConstraints() { return m_productionRateConstraintList; };
  std::vector< WellConstraintBase * >  getProdRateConstraints() const { return m_productionRateConstraintList; };
  std::vector< WellConstraintBase * >  getInjRateConstraints() { return m_injectionRateConstraintList; }
  std::vector< WellConstraintBase * >  getInjRateConstraints() const { return m_injectionRateConstraintList; }
protected:

  virtual void postInputInitialization() override;



private:

  /// Well type (as Type enum)
  Type m_type;

  /// Reference elevation
  real64 m_refElevation;

  /// Gravity coefficient of the reference elevation
  real64 m_refGravCoef;

  /// Well controls as a Control enum
  Control m_currentControl;

  /// Flag to decide whether rates are controlled at rates or surface conditions
  integer m_useSurfaceConditions;

  // Fuild model to compute properties for constraint equation user specified conditions
  std::unique_ptr< constitutive::ConstitutiveBase >  m_fluidSeparatorPtr;

  /// Reservoir region associated with reservoir volume constraint
  string m_referenceReservoirRegion;

  /// Surface pressure
  real64 m_surfacePres;

  /// Surface temperature
  real64 m_surfaceTemp;


  /// Well status table name
  string m_statusTableName;

  /// Perforation status table name
  string m_perfStatusTableName;

  /// Flag to enable crossflow
  integer m_isCrossflowEnabled;

  /// Tuning coefficient for the initial well pressure
  real64 m_initialPressureCoefficient;

  /// Rate sign. +1 for injector, -1 for producer
  real64 m_rateSign;


  /// Status table
  TableFunction const * m_statusTable;

  bool m_wellOpen;



  /// List of constraints
  //constraint_array m_ConstraintList;
  // Bool to trigger old/new constraint switch logic
  bool m_constraintSwitch;

  // Current constraint
  WellConstraintBase * m_currentConstraint;
  // Minimum and maximum BHP and WHP constraints
  MinimumBHPConstraint *  m_minBHPConstraint;
  MaximumBHPConstraint * m_maxBHPConstraint;


  // Lists of rate constraints
  std::vector< WellConstraintBase * > m_productionRateConstraintList;
  std::vector< WellConstraintBase * > m_injectionRateConstraintList;


  /// Well status
  WellControls::Status m_wellStatus;


  /// Region average pressure used in volume rate constraint calculations
  real64 m_regionAveragePressure;

  /// Region average temperature used in volume rate constraint calculations
  real64 m_regionAverageTemperature;

};

ENUM_STRINGS( WellControls::Type,
              "producer",
              "injector" );

ENUM_STRINGS( WellControls::Control,
              "BHP",
              "phaseVolRate",
              "totalVolRate",
              "massRate",
              "uninitialized" );


} //namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLCONTROLS_HPP
