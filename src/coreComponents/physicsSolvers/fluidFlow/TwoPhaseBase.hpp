/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file TwoPhaseBase.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_TWOPHASEBASE_HPP_
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_TWOPHASEBASE_HPP_

#include "physicsSolvers/fluidFlow/FlowSolverBase.hpp"

namespace geosx
{


/**
 * @class TwoPhaseBase
 *
 * class to perform an immiscible two-phase finite volume solve.
 */
class TwoPhaseBase : public FlowSolverBase
{
public:

  static constexpr localIndex NUM_PHASES = 2;
  static constexpr localIndex NUM_DOF    = 2;
  
  // define the column offset of the derivatives
  struct ColOffset
  {
    static constexpr integer DPRES = 0;
    static constexpr integer DSAT  = 1;
  };

  // define the row offset for the mass conservation eqns 
  struct RowOffset
  {
    static constexpr integer NONWETTING = 0;
    static constexpr integer WETTING    = 1;
  };

  
  /**
   * @brief main constructor for Group Objects
   * @param name the name of this instantiation of Group in the repository
   * @param parent the parent group of this instantiation of Group
   */
  TwoPhaseBase( const std::string & name,
                    Group * const parent );


  /// deleted default constructor
  TwoPhaseBase() = delete;

  /// deleted copy constructor
  TwoPhaseBase( TwoPhaseBase const & ) = delete;

  /// default move constructor
  TwoPhaseBase( TwoPhaseBase && ) = default;

  /// deleted assignment operator
  TwoPhaseBase & operator=( TwoPhaseBase const & ) = delete;

  /// deleted move operator
  TwoPhaseBase & operator=( TwoPhaseBase && ) = delete;

  /**
   * @brief default destructor
   */
  virtual ~TwoPhaseBase() override = default;

  virtual void RegisterDataOnMesh( Group * const MeshBodies ) override;

  virtual real64
  SolverStep( real64 const & time_n,
              real64 const & dt,
              integer const cycleNumber,
              DomainPartition * domain ) override;

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/

  virtual void
  ImplicitStepSetup( real64 const & time_n,
                     real64 const & dt,
                     DomainPartition * const domain,
                     DofManager & dofManager,
                     ParallelMatrix & matrix,
                     ParallelVector & rhs,
                     ParallelVector & solution ) override;

  virtual void
  AssembleSystem( real64 const time_n,
                  real64 const dt,
                  DomainPartition * const domain,
                  DofManager const & dofManager,
                  ParallelMatrix & matrix,
                  ParallelVector & rhs ) override;

  virtual bool
  CheckSystemSolution( DomainPartition const * const domain,
                       DofManager const & dofManager,
                       ParallelVector const & solution,
                       real64 const scalingFactor ) override;

  virtual void
  SolveSystem( DofManager const & dofManager,
               ParallelMatrix & matrix,
               ParallelVector & rhs,
               ParallelVector & solution ) override;

  virtual void
  ResetStateToBeginningOfStep( DomainPartition * const domain ) override;

  virtual void
  ImplicitStepComplete( real64 const & time,
                        real64 const & dt,
                        DomainPartition * const domain ) override;

  /**
   * @brief assembles the accumulation terms for all cells
   * @param time_n previous time value
   * @param dt time step
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   */
  void AssembleAccumulationTerms( DomainPartition const * const domain,
                                  DofManager const * const dofManager,
                                  ParallelMatrix * const matrix,
                                  ParallelVector * const rhs );

  /**
   * @brief assembles the flux terms for all cells
   * @param time_n previous time value
   * @param dt time step
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   */
  virtual void
  AssembleFluxTerms( real64 const time_n,
                     real64 const dt,
                     DomainPartition const * const domain,
                     DofManager const * const dofManager,
                     ParallelMatrix * const matrix,
                     ParallelVector * const rhs ) = 0;



  /**@}*/

  /**
   * @brief Function to update all constitutive state and dependent variables
   * @param dataGroup group that contains the fields
   */
  void UpdateState( Group * dataGroup ) const;

  /**
   * @brief Function to update all constitutive models
   * @param dataGroup group that contains the fields
   */
  void UpdateFluidModel( Group * const dataGroup ) const;


  struct viewKeyStruct : FlowSolverBase::viewKeyStruct
  {
    static constexpr auto elemDofFieldString = "elemVariables";

    // inputs
    static constexpr auto relPermNameString  = "relPermName";
    static constexpr auto relPermIndexString = "relPermIndex";
    
    // primary variables    
    static constexpr auto wettingPhaseSatString = "wettingPhaseSaturation";
    static constexpr auto deltaWettingPhaseSatString = "deltaWettingPhaseSaturation";
    static constexpr auto phaseSatString = "phaseSaturation"; // currently only used to update relperms
    
    // intermediate fields
    static constexpr auto phaseMobilityString = "phaseMobility";
    static constexpr auto dPhaseMobility_dPressureString = "dPhaseMobility_dPressure";
    static constexpr auto dPhaseMobility_dSaturationString = "dPhaseMobility_dSaturation";

    //backup fields
    static constexpr auto porosityOldString = "porosityOld";
    static constexpr auto phaseDensityOldString = "phaseDensityOld";

  } viewKeysTwoPhaseBase;

  viewKeyStruct & viewKeys()
  { return viewKeysTwoPhaseBase; }

  viewKeyStruct const & viewKeys() const
  { return viewKeysTwoPhaseBase; }

  struct groupKeyStruct : SolverBase::groupKeyStruct
  {
  } groupKeysTwoPhaseBase;

  groupKeyStruct & groupKeys()
  { return groupKeysTwoPhaseBase; }

  groupKeyStruct const & groupKeys() const
  { return groupKeysTwoPhaseBase; }

protected:

  virtual void PostProcessInput() override;

  virtual void InitializePreSubGroups( Group * const rootGroup ) override;
 
  virtual void InitializePostInitialConditions_PreSubGroups( dataRepository::Group * const rootGroup ) override;

protected:

  /**
   * @brief Resize the allocated multidimensional fields
   * @param domain the domain containing the mesh and fields
   *
   * Resize fields along dimensions 1 and 2 (0 is the size of containing object, i.e. element subregion)
   * once the number of phases/components is known (e.g. component fractions)
   */
  void ResizeFields( MeshLevel * const meshLevel );

  /**
   * @brief Function to update all constitutive models
   * @param dataGroup group that contains the fields
   */
  void UpdateSolidModel( Group * const dataGroup ) const;

  /**
   * @brief Function to update fluid mobility
   * @param dataGroup group that contains the fields
   */
  void UpdatePhaseMobility( Group * const dataGroup ) const;

  /**
   * @brief Update all relevant fluid models using current values of pressure and composition
   * @param dataGroup the group storing the required fields
   */
  void UpdateRelPermModel( Group * const dataGroup ) const; 
  
  /**
   * @brief Backup current values of all constitutive fields that participate in the accumulation term
   * @param domain the domain containing the mesh and fields
   */
  void BackupFields( DomainPartition * const domain );

  /**
   * @brief Setup stored views into domain data for the current step
   */
  void ResetViews( DomainPartition * const domain ) override;

  /// views into primary variable fields

  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> m_pressure;
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> m_deltaPressure;

  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> m_wettingPhaseSat;
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> m_deltaWettingPhaseSat;

  /// views into auxiliary data
  
  ElementRegionManager::ElementViewAccessor<arrayView2d<real64>> m_phaseMob;
  ElementRegionManager::ElementViewAccessor<arrayView2d<real64>> m_dPhaseMob_dPres;
  ElementRegionManager::ElementViewAccessor<arrayView2d<real64>> m_dPhaseMob_dSat;

  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> m_pvMult;
  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> m_dPvMult_dPres;
  
  ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> m_phaseDens;
  ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> m_dPhaseDens_dPres;

  /// views into backup fields

  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> m_porosityOld;
  ElementRegionManager::ElementViewAccessor<arrayView2d<real64>> m_phaseDensOld;

  /// name of the rel perm constitutive model
  string m_relPermName;

  /// index of the rel perm constitutive model
  localIndex m_relPermIndex;

  /// map from the phase indices to the row indices
  array1d<localIndex> m_phaseToRow;
  
  /// index of the wetting phase in the MaterialViewAccessors
  localIndex m_ipw;

  /// index of the non-wetting phase in the MaterialViewAccessors
  localIndex m_ipnw;

};


} /* namespace geosx */

#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_TWOPHASEBASE_HPP_
