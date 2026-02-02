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
#include "physicsSolvers/PhysicsSolverBase.hpp"
#include "common/format/EnumStrings.hpp"
#include "dataRepository/Group.hpp"
#include "functions/TableFunction.hpp"
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"

#include "physicsSolvers/fluidFlow/wells/WellInjectionConstraint.hpp"
#include "physicsSolvers/fluidFlow/wells/WellProductionConstraint.hpp"
#include "physicsSolvers/fluidFlow/wells/WellBHPConstraints.hpp"
#include "physicsSolvers/fluidFlow/wells/WellVolumeRateConstraint.hpp"
#include "physicsSolvers/fluidFlow/wells/WellPhaseVolumeRateConstraint.hpp"
#include "physicsSolvers/fluidFlow/wells/WellMassRateConstraint.hpp"
#include "physicsSolvers/fluidFlow/wells/WellLiquidRateConstraint.hpp"
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidBase.hpp"
#include "physicsSolvers/fluidFlow/wells/WellNewtonSolver.hpp"

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
class WellControls :  public dataRepository::Group //public PhysicsSolverBase
{
public:

  /** Type of wells
   * Either producer or injector
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


  /// String used to form the solverName used to register single-physics solvers in CoupledSolver
  static string coupledSolverAttributePrefix() { return "well"; }
  /**
   * @brief Create a new geometric object (box, plane, etc) as a child of this group.
   * @param childKey the catalog key of the new geometric object to create
   * @param childName the name of the new geometric object in the repository
   * @return the group child
   */
  virtual Group * createChild( string const & childKey, string const & childName ) override;
  /// Expand catalog for schema generation

  virtual void expandObjectCatalogs() override;

  virtual void registerWellDataOnMesh( WellElementSubRegion & subRegion ) = 0;
  virtual void setConstitutiveNames( ElementSubRegionBase & subRegion ) const = 0;
  /**
   * @brief Create well separator
   */
  virtual void  createSeparator( WellElementSubRegion & subRegion ) = 0;

  /**
   * @defgroup WellManager Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   * The "Well" versions apply to individual well subRegions, whereas the others apply to all wells
   */
  /**@{*/

  virtual void validateWellConstraints( real64 const & time_n,
                                        real64 const & dt,
                                        WellElementSubRegion const & subRegion ) = 0;

  virtual bool isCompositional() const = 0;
  /**
   *   * @brief Initialize well for the beginning of a simulation or restart
   *   @param domain the domain
   *   @param mesh the mesh level
   *   @param subRegion the well subRegion
   *  @param time_n the current time
   */
  virtual void initializeWell( DomainPartition & domain, MeshLevel & mesh, WellElementSubRegion & subRegion, real64 const & time_n ) = 0;
  /**
   * @brief function to set the next time step size
   * @param[in] currentTime the current time
   * @param[in] currentDt the current time step size
   * @param[in] domain the domain object
   * @return the prescribed time step size
   */
  real64 setNextDt( real64 const & currentTime,
                    real64 const & currentDt,
                    WellElementSubRegion & subRegion );
  // Bring the base class implicitStepSetup into scope to avoid hiding the overloaded virtual function


  virtual void implicitStepSetup( real64 const & time_n,
                                  real64 const & GEOS_UNUSED_PARAM( dt ),
                                  ElementRegionManager & elemManager,
                                  WellElementSubRegion & subRegion ) = 0;

  virtual void
  implicitStepComplete( real64 const & time,
                        real64 const & dt,
                        WellElementSubRegion const & subRegion ) = 0;
  virtual real64 updateSubRegionState( ElementRegionManager const & elemManager, WellElementSubRegion & subRegion ) = 0;

  /**
   * @brief Function to evaluate well constraints after applying the solution update
   * @param time_n the time at the beginning of the time step
   * @param subRegion the well subRegion
   * @return true if all constraints are satisfied, false otherwise
   */
  bool evaluateConstraints( real64 const & time_n,
                            WellElementSubRegion & subRegion );

  /**
   * @brief Function to evaluate well constraints after applying the solution update
   * @param time_n the time at the beginning of the time step
   * @param subRegion the well subRegion
   * @return true if all constraints are satisfied, false otherwise
   */
  bool evaluateConstraints( real64 const & time_n,
                            real64 const & dt,
                            integer const cycleNumber,
                            integer const coupledIterationNumber,
                            DomainPartition & domain,
                            MeshLevel & mesh,
                            ElementRegionManager & elemManager,
                            WellElementSubRegion & subRegion,
                            DofManager const & dofManager );

  void solveConstraint( WellConstraintBase *constraint,
                        real64 const & time_n,
                        real64 const & dt,
                        integer const cycleNumber,
                        integer const coupledIterationNumber,
                        DomainPartition & domain,
                        MeshLevel & mesh,
                        ElementRegionManager & elemManager,
                        WellElementSubRegion & subRegion,
                        DofManager const & dofManager );
  /**
   * @brief assembles the accumulation term for an individual well
   * @param time_n time at the beginning of the time step
   * @param dt the time step size
   * @param subRegion the well subregion containing all the primary and dependent fields
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   */
  virtual void assembleWellAccumulationTerms( real64 const & time,
                                              real64 const & dt,
                                              WellElementSubRegion & subRegion,
                                              DofManager const & dofManager,
                                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                              arrayView1d< real64 > const & localRhs ) = 0;
  /**
   * @brief assembles the well momentum terms for an individual well
   * @param time_n time at the beginning of the time step
   * @param dt the time step size
   * @param subRegion the well subregion containing all the primary and dependent fields
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   */
  virtual void assembleWellPressureRelations( real64 const & time_n,
                                              real64 const & dt,
                                              WellElementSubRegion const & subRegion,
                                              DofManager const & dofManager,
                                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                              arrayView1d< real64 > const & localRhs ) = 0;
  /**
   * @brief assembles the well constraint terms for an individual well
   * @param time_n time at the beginning of the time step
   * @param dt the time step size
   * @param subRegion the well subregion containing all the primary and dependent fields
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   */
  virtual void assembleWellConstraintTerms( real64 const & time_n,
                                            real64 const & dt,
                                            WellElementSubRegion const & subRegion,
                                            DofManager const & dofManager,
                                            CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                            arrayView1d< real64 > const & localRhs ) = 0;

  /**
   * @brief Recompute the perforation rates for all the wells
   * @param time_n the time at the beginning of the time step
   * @param dt the time step size
   * @param elemManager the element region manager
   * @param subRegion the well subregion containing all the primary and dependent fields
   */
  virtual void computeWellPerforationRates( real64 const & time_n,
                                            real64 const & GEOS_UNUSED_PARAM( dt ),
                                            ElementRegionManager & elemManager,
                                            WellElementSubRegion & subRegion ) = 0;

  /**
   * @brief assembles the flux terms for individual well for all connections between well elements
   * @param time_n previous time value
   * @param dt time step
   * @param subRegion the well subregion containing all the primary and dependent fields
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   */
  virtual void assembleWellFluxTerms( real64 const & time,
                                      real64 const & dt,
                                      WellElementSubRegion const & subRegion,
                                      DofManager const & dofManager,
                                      CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                      arrayView1d< real64 > const & localRhs ) = 0;
  virtual real64
  calculateWellResidualNorm( real64 const & time_n,
                             real64 const & dt,
                             NonlinearSolverParameters const & nonlinearSolverParameters,
                             WellElementSubRegion const & subRegion,
                             DofManager const & dofManager,
                             arrayView1d< real64 const > const & localRhs ) = 0;

  virtual array1d< real64 >
  calculateLocalWellResidualNorm( real64 const & time_n,
                                  real64 const & dt,
                                  NonlinearSolverParameters const & nonlinearSolverParameters,
                                  WellElementSubRegion const & subRegion,
                                  DofManager const & dofManager,
                                  arrayView1d< real64 const > const & localRhs ) = 0;

  virtual real64
  scalingForWellSystemSolution( WellElementSubRegion & subRegion,
                                DofManager const & dofManager,
                                arrayView1d< real64 const > const & localSolution ) = 0;

  virtual bool
  checkWellSystemSolution( WellElementSubRegion & subRegion,
                           DofManager const & dofManager,
                           arrayView1d< real64 const > const & localSolution,
                           real64 const scalingFactor ) = 0;

  virtual void
  applyWellSystemSolution( DofManager const & dofManager,
                           arrayView1d< real64 const > const & localSolution,
                           real64 const scalingFactor,
                           real64 const dt,
                           DomainPartition & domain,
                           MeshLevel & mesh,
                           WellElementSubRegion & subRegion ) = 0;

  virtual void applyWellBoundaryConditions( real64 const time_n,
                                            real64 const dt,
                                            ElementRegionManager & elemManager,
                                            WellElementSubRegion & subRegion,
                                            DofManager const & dofManager,
                                            arrayView1d< real64 > const & localRhs,
                                            CRSMatrixView< real64, globalIndex const > const & localMatrix ) = 0;
  /**
   * @brief Recompute all dependent quantities from primary variables (including constitutive models)
   * @param subRegion the well subregion containing all the primary and dependent fields
   */
  virtual real64 updateWellState( ElementRegionManager const & elemManager, WellElementSubRegion & subRegion ) = 0;

  /**
   * @brief Reset the well state to the beginning of the time step
   * @param subRegion the well subregion containing all the primary and dependent fields
   */
  virtual void resetStateToBeginningOfStep( ElementRegionManager const & elemManager, WellElementSubRegion & subRegion ) = 0;

  virtual void postInputInitialization() override;

  virtual void initializeWellPostInitialConditionsPreSubGroups( WellElementSubRegion & subRegion ) = 0;
  virtual void printRates( real64 const & time_n,
                           real64 const & dt,
                           WellElementSubRegion const & subRegion ) = 0;
  /**@}*/
  /**
   * @brief Apply a given functor to a container if the container can be
   *        cast to one of the specified types.
   * @tparam CASTTYPE      the first type that will be used in the attempted casting of container
   * @tparam CASTTYPES     a variadic list of types that will be used in the attempted casting of container
   * @tparam CONTAINERTYPE the type of container
   * @tparam LAMBDA        the type of lambda function to call in the function
   * @param[in] container  a pointer to the container which will be passed to the lambda function
   * @param[in] lambda     the lambda function to call in the function
   * @return               a boolean to indicate whether the lambda was successfully applied to the container.
   */
  template< typename T0, typename T1, typename ... CASTTYPES, typename CONTAINERTYPE, typename LAMBDA >
  static bool applyLambdaToContainer( CONTAINERTYPE container, LAMBDA && lambda )
  {
    using Pointee = std::remove_pointer_t< std::remove_reference_t< CONTAINERTYPE > >;
    using T = std::conditional_t< std::is_const< Pointee >::value, T0 const, T0 >;
    T * const castedContainer = dynamic_cast< T * >( container );

    if( castedContainer != nullptr )
    {
      lambda( *castedContainer );
      return true;
    }

    return applyLambdaToContainer< T1, CASTTYPES... >( container, std::forward< LAMBDA >( lambda ) );
  }

  // Base case: no more types to try
  template< typename CONTAINERTYPE, typename LAMBDA >
  static bool applyLambdaToContainer( CONTAINERTYPE /*container*/, LAMBDA && /*lambda*/ )
  {
    return false;
  }

  // Single-type overload: try only T0 and stop
  template< typename T0, typename CONTAINERTYPE, typename LAMBDA >
  static bool applyLambdaToContainer( CONTAINERTYPE container, LAMBDA && lambda )
  {
    using Pointee = std::remove_pointer_t< std::remove_reference_t< CONTAINERTYPE > >;
    using T = std::conditional_t< std::is_const< Pointee >::value, T0 const, T0 >;
    T * const castedContainer = dynamic_cast< T * >( container );

    if( castedContainer != nullptr )
    {
      lambda( *castedContainer );
      return true;
    }

    return false;
  }


  /**
   * @copydoc forInjectionConstraints(LAMBDA &&)
   */
  template< typename GROUPTYPE = Group, typename ... GROUPTYPES, typename LAMBDA >
  void forInjectionConstraints( LAMBDA && lambda ) const
  {
    for( auto const * constraintIter : m_injectionRateConstraintList )
    {
      applyLambdaToContainer< GROUPTYPE, GROUPTYPES... >( constraintIter, [&]( auto const & castedSubGroup )
      {
        lambda( castedSubGroup );
      } );
    }
  }

  // non-const overload
  template< typename GROUPTYPE = Group, typename ... GROUPTYPES, typename LAMBDA >
  void forInjectionConstraints( LAMBDA && lambda )
  {
    for( auto * constraintIter : m_injectionRateConstraintList )
    {
      applyLambdaToContainer< GROUPTYPE, GROUPTYPES... >( constraintIter, [&]( auto & castedSubGroup )
      {
        lambda( castedSubGroup );
      } );
    }
  }
  /**
   * @copydoc forProductionConstraints(LAMBDA &&)
   */
  template< typename GROUPTYPE = Group, typename ... GROUPTYPES, typename LAMBDA >
  void forProductionConstraints( LAMBDA && lambda ) const
  {
    for( auto const * constraintIter : m_productionRateConstraintList )
    {
      applyLambdaToContainer< GROUPTYPE, GROUPTYPES... >( constraintIter, [&]( auto const & castedSubGroup )
      {
        lambda( castedSubGroup );
      } );
    }
  }

  // non-const overload
  template< typename GROUPTYPE = Group, typename ... GROUPTYPES, typename LAMBDA >
  void forProductionConstraints( LAMBDA && lambda )
  {
    for( auto * constraintIter : m_productionRateConstraintList )
    {
      applyLambdaToContainer< GROUPTYPE, GROUPTYPES... >( constraintIter, [&]( auto & castedSubGroup )
      {
        lambda( castedSubGroup );
      } );
    }
  }
  /**
   * @name Getters / Setters
   */
  ///@{

  /**
   * @brief Get the Constitutive Name object
   *
   * @tparam CONSTITUTIVE_BASE_TYPE the base type of the constitutive model.
   * @param subRegion the element subregion on which the constitutive model is registered
   * @return the name name of the constitutive model of type CONSTITUTIVE_BASE_TYPE registered on the subregion.
   */
  template< typename CONSTITUTIVE_BASE_TYPE >
  static string getConstitutiveName( ElementSubRegionBase const & subRegion );

  /**
   * @brief Register wrapper with given name and store constitutive model name on the subregion
   *
   * @tparam CONSTITUTIVE the base type of the constitutive model.
   * @param subRegion the subregion on which the constitutive model is registered.
   * @param wrapperName the wrapper name to register.
   * @param constitutiveType the type description of the constitutive model.
   */
  template< typename CONSTITUTIVE >
  void setConstitutiveName( ElementSubRegionBase & subRegion, string const & wrapperName, string const & constitutiveType ) const;

  /**
   * @brief return the list of target regions
   * @return the array of region names
   */
  string_array const & getTargetRegionNames() const {return m_targetRegionNames;}
  /**
   * @brief Get the control type for the well.
   * @return the Control enum enforced at the well
   */
  std::string getFlowSolverName() const { return m_flowSolverName; }

  /**
   * @brief Set the control type for the well.
   * @param[in] flowSolverName  the name of the flow solver
   */
  void setFlowSolverName( const std::string & flowSolverName )   { m_flowSolverName = flowSolverName;    }

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
   * @brief Get the input control type for the well.
   * @return the Control enum enforced at the well
   */
  Control getInputControl() const { return m_inputControl; }

  /**
   * @brief getter for esitmator switch
   * @return True if estimate well solution
   */
  integer estimateSolution() const { return m_estimateSolution; }

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
  real64 getTargetBHP( real64 const & targetTime ) const;


  /**
   * @brief Const accessor for the temperature of the injection stream
   * @return the temperature of the injection stream
   */
  real64 getInjectionTemperature() const;

  /**
   * @brief Const accessor for the  injection stream
   * @return the injection stream
   */
  arrayView1d< real64 const > getInjectionStream() const;

  /**
   * @brief Const accessor for the phase constraint index
   * @return phase index associated with phase constraint
   */
  integer getConstraintPhaseIndex() const;

  /**
   * @brief Return the reference elvation where pressure constraint is measured
   * @return  vertical location of constraint
   */
  real64 getReferenceElevation() const;


  /**
   * @brief Getter for the flag specifying whether we check rates at surface or reservoir conditions
   * @return 1 if we use surface conditions, and 0 otherwise
   */
  integer useSurfaceConditions() const { return m_useSurfaceConditions; }

  /**
   * @brief Getter for the reservoir region associated with reservoir volume constraint
   * @return name of reservoir region
   */
  string getReferenceReservoirRegion() const { return m_referenceReservoirRegion; }

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
   * @brief getter for iso/thermal switch
   * @return True if thermal
   */
  integer isThermal() const { return m_isThermal; }

  /**
   * @brief setter for iso/thermal switch
   * @param[in] isThermal
   */

  void setThermal( bool isThermal )   {  m_isThermal=isThermal; }
  /**
   * @brief Is the well open (or shut) at currentTime, status initalized in WellSolverBase::implicitStepSetup
   * @return a boolean
   */
  bool isWellOpen() const;

  /**
   * @brief Set the well state
   * @param[in] open boolean
   */
  void setWellState( bool open );
  /**
   * @brief Get the well state
   * @return a boolean
   */
  bool getWellState() const;


  /**
   * @brief Set the current consrtaint
   * @param[in] currentConstraint pointer to constraint
   */
  void setCurrentConstraint( WellConstraintBase * currentConstraint ) { m_currentConstraint = currentConstraint;}
  /**
   * @brief Get the current consrtaint
   * @return pointer to constraint
   */
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
   * @brief Utility function to keep the well variables during a time step (used in
   * poromechanics simulations)
   * @param[in] keepVariablesConstantDuringInitStep flag to tell the solver to freeze its
   * primary variables during a time step
   * @detail This function is meant to be called by a specific task before/after the
   * initialization step
   */
  void setKeepVariablesConstantDuringInitStep( bool const keepVariablesConstantDuringInitStep )
  { m_keepVariablesConstantDuringInitStep = keepVariablesConstantDuringInitStep; }


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


  virtual string wellElementDofName() const = 0;

  virtual string resElementDofName() const = 0;

  virtual localIndex numFluidComponents() const = 0;

  virtual localIndex numFluidPhases() const = 0;
  /**
   * @brief Struct to serve as a container for variable strings and keys.
   * @struct viewKeyStruct
   */
  struct viewKeyStruct
  {
    /// String key for the fluid model names
    static constexpr char const * fluidNamesString() { return "fluidNames"; }
    ///   String key for the write CSV flag
    static constexpr char const * writeCSVFlagString() { return "writeCSV"; }
    static constexpr char const * timeStepFromTablesFlagString() { return "timeStepFromTables"; }
/// @return string for the   targetRegions wrapper
    static constexpr char const * targetRegionsString() { return "targetRegions"; }

    /// String key for the well reference elevation (for BHP control)
    static constexpr char const * refElevString() { return "referenceElevation"; }
    /// String key for the well type
    static constexpr char const * typeString() { return "type"; }
    /// String key for the well current control
    static constexpr char const * currentControlString() { return "currentControl"; }
    /// String key for the well input control
    static constexpr char const * inputControlString() { return "control"; }
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

    /// string key for the minimum BHP presssure for a producer
    static constexpr char const * minimumBHPConstraintString() { return "MinimumBHPConstraint"; }
    /// string key for the maximum BHP presssure for a injection
    static constexpr char const * maximumBHPConstraintString() { return "MaximumBHPConstraint"; }
    /// string key for the maximum phase rate for a producer
    static constexpr char const * productionPhaseVolumeRateConstraintString() { return "ProductionPhaseVolumeRateConstraint"; }
    /// string key for the maximum phase rate for a injection
    static constexpr char const * injectionPhaseVolumeRateConstraint() { return "InjectionPhaseVolumeRateConstraint"; }
    /// string key for the maximum volume rate for a producer
    static constexpr char const * productionVolumeRateConstraint() { return "ProductionVolumeRateConstraint"; }
    /// string key for the maximum volume rate for a injector
    static constexpr char const * injectionVolumeRateConstraint() { return "InjectionVolumeRateConstraint"; }
    /// string key for the maximum mass rate for a producer
    static constexpr char const * productionMassRateConstraint() { return "ProductionMassRateConstraint"; }
    /// string key for the maximum mass rate for a injector
    static constexpr char const * injectionMassRateConstraint() { return "InjectionMassRateConstraint"; }
    /// string key for the liquid rate for a producer
    static constexpr char const * productionLiquidRateConstraint() { return "ProductionLiquidRateConstraint"; }
    /// string key for the minimum BHP presssure for a producer
    static constexpr char const * wellNewtonSolverString() { return "WellNewtonSolver"; }

    /// string key for the esitmate well solution flag
    static constexpr char const * estimateWellSolutionString() { return "estimateWellSolution"; }
    /// string key for the enable iso thermal estimator flag
    static constexpr char const * enableIsoThermalEstimatorString() { return "enableIsoThermalEstimator"; }

    // control data (not registered on the mesh)

    static constexpr char const * massDensityString() { return "massDensity";}

    static constexpr char const * currentBHPString() { return "currentBHP"; }

    static constexpr char const * currentPhaseVolRateString() { return "currentPhaseVolumetricRate"; }
    static constexpr char const * currentVolRateString() { return "currentVolRate"; }

    static constexpr char const * currentTotalVolRateString() { return "currentTotalVolumetricRate"; }

    static constexpr char const * currentMassRateString() { return "currentMassRate"; }

  }
  /// ViewKey struct for the WellControls class
  viewKeysWellControls;
  void setPerforationStatus( real64 const & time_n, WellElementSubRegion & subRegion );
  void setGravCoef( WellElementSubRegion & subRegion, R1Tensor const & gravVector );
  /**
   * @brief Set next time step based on a table function
   * @param[in] table the table function
   * @param[in] currentTime the current time
   * @param[inout] nextDt the time step
   */

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

  /**
   * @brief Getters for constraint lists
   */
  std::vector< WellConstraintBase * >  getProdRateConstraints() { return m_productionRateConstraintList; };
  std::vector< WellConstraintBase * >  getProdRateConstraints() const { return m_productionRateConstraintList; };
  std::vector< WellConstraintBase * >  getInjRateConstraints() { return m_injectionRateConstraintList; }
  std::vector< WellConstraintBase * >  getInjRateConstraints() const { return m_injectionRateConstraintList; }

  /**
   * @brief Set thermal effects enable
   * @param[in] true/false
   */
  void enableThermalEffects ( bool enable ) { m_thermalEffectsEnabled = enable; };

  /**
   * @brief Are thermal effects enabled
   * @return true if thermal effects are enabled, false otherwise
   */
  bool thermalEffectsEnabled() const { return m_thermalEffectsEnabled; }

  /**
   * @brief Is isoThermalEstimator  enabled
   * @return true if isoThermalEstimator is enabled, false otherwise
   */
  bool isoThermalEstimatorEnabled() const { return m_enableIsoThermalEstimator; }

  WellNewtonSolver & getWellNewtonSolver() { return m_wellNewtonSolver; }

  void selectWellConstraint( real64 const & time_n,
                             real64 const & dt,
                             integer const cycleNumber,
                             integer const coupledIterationNumber,
                             DomainPartition & domain,
                             MeshLevel & mesh,
                             ElementRegionManager & elemManager,
                             WellElementSubRegion & subRegion,
                             DofManager const & dofManager );


protected:

  virtual void postRestartInitialization( )override;
private:
  /// List of names of regions the solver will be applied to
  string_array m_targetRegionNames;

protected:

  /// Well type (as Type enum)
  Type m_type;

  /// Name of the flow solver managing this well
  std::string m_flowSolverName;

  /// the max number of fluid phases
  integer m_numPhases;

  /// the number of fluid components
  integer m_numComponents;

  /// the number of Degrees of Freedom per well element
  integer m_numDofPerWellElement;

  /// the number of Degrees of Freedom per reservoir element
  integer m_numDofPerResElement;

  /// flag indicating whether thermal formulation is used
  integer m_isThermal;
  /// flag to freeze the initial state during initialization in coupled problems
  bool m_keepVariablesConstantDuringInitStep;
  /// rates output
  integer m_writeCSV;
  string const m_ratesOutputDir;

  // flag to enable time step selection base on rates/bhp tables coordinates
  integer m_timeStepFromTables;
  /// Reference elevation
  real64 m_refElevation;

  /// Gravity coefficient of the reference elevation
  real64 m_refGravCoef;

  /// Input well controls as a Control enum
  Control m_inputControl;

  /// Well controls as a Control enum
  Control m_currentControl;

  /// Flag to decide whether rates are controlled at rates or surface conditions
  integer m_useSurfaceConditions;

  // Fuild model to compute properties for constraint equation user specified conditions
  std::unique_ptr< constitutive::ConstitutiveBase >  m_fluidSeparatorPtr;

  /// name of the fluid constitutive model used as a reference for component/phase description on subregion
  string m_referenceFluidModelName;

  /// Reservoir region associated with reservoir volume constraint
  string m_referenceReservoirRegion;

  /// Surface pressure
  real64 m_surfacePres;

  /// Surface temperature
  real64 m_surfaceTemp;

  /// Perforation status table name
  string m_perfStatusTableName;

  /// Flag to enable crossflow
  integer m_isCrossflowEnabled;

  /// Tuning coefficient for the initial well pressure
  real64 m_initialPressureCoefficient;

  // Current constrint
  WellConstraintBase * m_currentConstraint;

  // Minimum and maximum BHP and WHP constraints
  MinimumBHPConstraint *  m_minBHPConstraint;
  MaximumBHPConstraint * m_maxBHPConstraint;

  // Lists of rate constraints
  std::vector< WellConstraintBase * > m_productionRateConstraintList;
  std::vector< WellConstraintBase * > m_injectionRateConstraintList;

  /// Well status
  WellControls::Status m_wellStatus;

  /// Well open flag
  bool m_wellOpen;

  /// Well status table name
  string m_statusTableName;

  /// Status table
  TableFunction const * m_statusTable;

  /// Region average pressure used in volume rate constraint calculations
  real64 m_regionAveragePressure;

  /// Region average temperature used in volume rate constraint calculations
  real64 m_regionAverageTemperature;

  integer m_estimateSolution;
  integer m_enableIsoThermalEstimator;
  bool m_thermalEffectsEnabled;

  WellNewtonSolver m_wellNewtonSolver;
};


// Use local aliases to avoid accidental macro expansion of the tokens 'Type' or 'Control'
using WellControls_Type = WellControls::Type;
ENUM_STRINGS( WellControls_Type,
              "producer",
              "injector" );

using WellControls_Control = WellControls::Control;
ENUM_STRINGS( WellControls_Control,
              "BHP",
              "phaseVolRate",
              "totalVolRate",
              "massRate",
              "uninitialized" );


template< typename CONSTITUTIVE >
void WellControls::setConstitutiveName( ElementSubRegionBase & subRegion, string const & wrapperName, string const & constitutiveType ) const
{
  subRegion.registerWrapper< string >( wrapperName ).
    setPlotLevel( dataRepository::PlotLevel::NOPLOT ).
    setRestartFlags( dataRepository::RestartFlags::NO_WRITE ).
    setSizedFromParent( 0 );

  string & constitutiveName = subRegion.getReference< string >( wrapperName );
  constitutiveName = getConstitutiveName< CONSTITUTIVE >( subRegion );
  GEOS_ERROR_IF( constitutiveName.empty(), GEOS_FMT( "{}: {} constitutive model not found on subregion {}",
                                                     getDataContext(), constitutiveType, subRegion.getName() ) );
}
template< typename CONSTITUTIVE_BASE_TYPE >
string WellControls::getConstitutiveName( ElementSubRegionBase const & subRegion )
{
  string validName;
  dataRepository::Group const & constitutiveModels = subRegion.getConstitutiveModels();

  constitutiveModels.forSubGroups< CONSTITUTIVE_BASE_TYPE >( [&]( dataRepository::Group const & model )
  {
    GEOS_ERROR_IF( !validName.empty(), "A valid constitutive model was already found." );
    validName = model.getName();
  } );

  return validName;
}

/**
 * @brief Get the Constitutive Model object
 * @tparam BASETYPE the base type of the constitutive model.
 * @tparam LOOKUP_TYPE the type of the key used to look up the constitutive model.
 * @param dataGroup the data group containing the constitutive models.
 * @param key the key used to look up the constitutive model.
 * @return the constitutive model of type @p BASETYPE registered on the @p dataGroup with the key @p key.
 */
template< typename BASETYPE = constitutive::ConstitutiveBase, typename LOOKUP_TYPE >
static BASETYPE const & getConstitutiveModel( dataRepository::Group const & dataGroup, LOOKUP_TYPE const & key )
{
  dataRepository::Group const & constitutiveModels = dataGroup.getGroup( ElementSubRegionBase::groupKeyStruct::constitutiveModelsString() );
  return constitutiveModels.getGroup< BASETYPE >( key );
}
/**
 * @brief Get the Constitutive Model object
 * @tparam BASETYPE the base type of the constitutive model.
 * @tparam LOOKUP_TYPE the type of the key used to look up the constitutive model.
 * @param dataGroup the data group containing the constitutive models.
 * @param key the key used to look up the constitutive model.
 * @return the constitutive model of type @p BASETYPE registered on the @p dataGroup with the key @p key.
 */
template< typename BASETYPE = constitutive::ConstitutiveBase, typename LOOKUP_TYPE >
static BASETYPE & getConstitutiveModel( dataRepository::Group & dataGroup, LOOKUP_TYPE const & key )
{
  dataRepository::Group & constitutiveModels = dataGroup.getGroup( ElementSubRegionBase::groupKeyStruct::constitutiveModelsString() );
  return constitutiveModels.getGroup< BASETYPE >( key );
}
} //namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLCONTROLS_HPP
