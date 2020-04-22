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
   * @param targetIndex index for the constitutive models
   */
  void UpdateState( Group & dataGroup, localIndex const targetIndex );

  /**
   * @brief Function to update the fluid constitutive model
   * @param dataGroup group that contains the fields
   * @param targetIndex index for the fluid constitutive model
   */
  void UpdateFluidModel( Group & dataGroup, localIndex const targetIndex );


  struct viewKeyStruct : FlowSolverBase::viewKeyStruct
  {
    static constexpr auto elemDofFieldString = "elemVariables";

    // inputs
    static constexpr auto relPermNamesString  = "relPermNames";

    // primary variables
    static constexpr auto phaseSatString = "phaseSaturation";
    static constexpr auto deltaPhaseSatString = "deltaPhaseSaturation";
    static constexpr auto newPhaseSatString = "newPhaseSaturation";

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
  {} groupKeysTwoPhaseBase;

  groupKeyStruct & groupKeys()
  { return groupKeysTwoPhaseBase; }

  groupKeyStruct const & groupKeys() const
  { return groupKeysTwoPhaseBase; }

protected:

  virtual void PostProcessInput() override;

  virtual void InitializePreSubGroups( Group * const rootGroup ) override;

  virtual void InitializePostInitialConditions_PreSubGroups( Group * const rootGroup ) override;

  /**
   * @brief Resize the allocated multidimensional fields
   * @param domain the domain containing the mesh and fields
   *
   * Resize fields along dimension 1 (0 is the size of containing object, i.e. element subregion)
   */
  virtual void ResizeFields( MeshLevel & meshLevel );

  /**
   * @brief Checks constitutive models for consistency
   * @param cm        reference to the global constitutive model manager
   */
  void ValidateConstitutiveModels( constitutive::ConstitutiveManager const & cm ) const;

  /**
   * @brief Function to update the solid model
   * @param dataGroup group that contains the fields
   * @param targetIndex index for the solid constitutive model
   */
  void UpdateSolidModel( Group & dataGroup, localIndex const targetIndex );

  /**
   * @brief Function to update fluid mobility
   * @param dataGroup group that contains the fields
   * @param targetIndex index for the fluid and relperm constitutive models
   */
  void UpdatePhaseMobility( Group & dataGroup, localIndex const targetIndex ) const;

  /**
   * @brief Function to update relative permeability
   * @param dataGroup the group storing the required fields
   * @param targetIndex index for the relperm constitutive model
   */
  void UpdateRelPermModel( Group & dataGroup, localIndex const targetIndex );

  /**
   * @brief Backup current values of all constitutive fields that participate in the accumulation term
   * @param domain the domain containing the mesh and fields
   */
  void BackupFields( DomainPartition * const domain );

  /**
   * @brief Setup stored views into domain data for the current step
   */
  void ResetViews( DomainPartition * const domain ) override;

  /// name of the rel perm constitutive model
  array1d< string > m_relPermModelNames;

  /// map from the phase indices to the row indices
  array1d< localIndex > m_phaseToRow;

  /// views into primary variable fields

  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > m_pressure;
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > m_deltaPressure;

  // when we move to the old saturation / new saturation format, I will be able to remove one of these three vectors
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > m_phaseSat;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > m_deltaPhaseSat;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > m_newPhaseSat;

  /// views into auxiliary data

  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > m_phaseMob;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > m_dPhaseMob_dPres;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > m_dPhaseMob_dSat;

  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > m_pvMult;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > m_dPvMult_dPres;

  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 > > m_phaseDens;
  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 > > m_dPhaseDens_dPres;

  /// views into backup fields

  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > m_porosityOld;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > m_phaseDensOld;

};


} /* namespace geosx */

#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_TWOPHASEBASE_HPP_
