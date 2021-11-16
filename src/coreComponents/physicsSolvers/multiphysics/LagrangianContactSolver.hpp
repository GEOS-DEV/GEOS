/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file LagrangianContactSolver.hpp
 *
 */

#ifndef GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_LAGRANGIANCONTACTSOLVER_HPP_
#define GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_LAGRANGIANCONTACTSOLVER_HPP_

#include "physicsSolvers/SolverBase.hpp"

namespace geosx
{

class SolidMechanicsLagrangianFEM;

class LagrangianContactSolver : public SolverBase
{
public:

  LagrangianContactSolver( const string & name,
                           Group * const parent );

  ~LagrangianContactSolver() override;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string catalogName()
  {
    return "LagrangianContact";
  }

  virtual void initializePreSubGroups() override;

  virtual void registerDataOnMesh( Group & MeshBodies ) override final;

  virtual void
  setupDofs( DomainPartition const & domain,
             DofManager & dofManager ) const override;

  virtual void
  setupSystem( DomainPartition & domain,
               DofManager & dofManager,
               CRSMatrix< real64, globalIndex > & localMatrix,
               array1d< real64 > & localRhs,
               array1d< real64 > & localSolution,
               bool const setSparsity = true ) override;

  virtual void
  implicitStepSetup( real64 const & time_n,
                     real64 const & dt,
                     DomainPartition & domain ) override final;

  virtual void
  implicitStepComplete( real64 const & time_n,
                        real64 const & dt,
                        DomainPartition & domain ) override final;

  virtual void
  assembleSystem( real64 const time,
                  real64 const dt,
                  DomainPartition & domain,
                  DofManager const & dofManager,
                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                  arrayView1d< real64 > const & localRhs ) override;

  virtual void
  applyBoundaryConditions( real64 const time,
                           real64 const dt,
                           DomainPartition & domain,
                           DofManager const & dofManager,
                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                           arrayView1d< real64 > const & localRhs ) override;

  virtual real64
  calculateResidualNorm( DomainPartition const & domain,
                         DofManager const & dofManager,
                         arrayView1d< real64 const > const & localRhs ) override;

  virtual void
  solveSystem( DofManager const & dofManager,
               ParallelMatrix & matrix,
               ParallelVector & rhs,
               ParallelVector & solution ) override;

  virtual void
  applySystemSolution( DofManager const & dofManager,
                       arrayView1d< real64 const > const & localSolution,
                       real64 const scalingFactor,
                       DomainPartition & domain ) override;

  virtual void
  resetStateToBeginningOfStep( DomainPartition & domain ) override;

  virtual real64
  solverStep( real64 const & time_n,
              real64 const & dt,
              int const cycleNumber,
              DomainPartition & domain ) override;

  virtual void
  setNextDt( real64 const & currentDt,
             real64 & nextDt ) override;


  virtual real64
  explicitStep( real64 const & time_n,
                real64 const & dt,
                integer const cycleNumber,
                DomainPartition & domain ) override;

  virtual real64
  nonlinearImplicitStep( real64 const & time_n,
                         real64 const & dt,
                         integer const cycleNumber,
                         DomainPartition & domain ) override;

  virtual bool
  lineSearch( real64 const & time_n,
              real64 const & dt,
              integer const cycleNumber,
              DomainPartition & domain,
              DofManager const & dofManager,
              CRSMatrixView< real64, globalIndex const > const & localMatrix,
              arrayView1d< real64 > const & localRhs,
              arrayView1d< real64 const > const & localSolution,
              real64 const scaleFactor,
              real64 & lastResidual ) override;

  void computeFaceDisplacementJump( DomainPartition & domain );

  void assembleForceResidualDerivativeWrtTraction( DomainPartition & domain,
                                                   DofManager const & dofManager,
                                                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                   arrayView1d< real64 > const & localRhs );

  void assembleTractionResidualDerivativeWrtDisplacementAndTraction( DomainPartition const & domain,
                                                                     DofManager const & dofManager,
                                                                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                     arrayView1d< real64 > const & localRhs );

  void assembleStabilization( DomainPartition const & domain,
                              DofManager const & dofManager,
                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                              arrayView1d< real64 > const & localRhs );

  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    constexpr static char const * solidSolverNameString() { return "solidSolverName"; }
    constexpr static char const * stabilizationNameString() { return "stabilizationName"; }
    constexpr static char const * contactRelationNameString() { return "contactRelationName"; }
    constexpr static char const * activeSetMaxIterString() { return "activeSetMaxIter"; }

    constexpr static char const * rotationMatrixString() { return "rotationMatrix"; }

    constexpr static char const * tractionString() { return "traction"; }
    constexpr static char const * deltaTractionString() { return "deltaTraction"; }
    constexpr static char const * fractureStateString() { return "fractureState"; }
    constexpr static char const * previousFractureStateString() { return "previousFractureState"; }
    constexpr static char const * dispJumpString() { return "displacementJump"; }
    constexpr static char const * previousDispJumpString() { return "previousLocalJump"; }

    constexpr static char const * slidingCheckToleranceString() { return "slidingCheckTolerance"; }
    constexpr static char const * normalDisplacementToleranceString() { return "normalDisplacementTolerance"; }
    constexpr static char const * normalTractionToleranceString() { return "normalTractionTolerance"; }
    constexpr static char const * slidingToleranceString() { return "slidingTolerance"; }

  };

  string const & getContactRelationName() const { return m_contactRelationName; }

  SolidMechanicsLagrangianFEM const * getSolidSolver() const { return m_solidSolver; }

  SolidMechanicsLagrangianFEM * getSolidSolver() { return m_solidSolver; }

  integer const & getActiveSetMaxIter() const { return m_activeSetMaxIter; }

protected:
  virtual void postProcessInput() override final;

  virtual void
  initializePostInitialConditionsPreSubGroups() override final;

private:

  string m_solidSolverName;
  SolidMechanicsLagrangianFEM * m_solidSolver;

  string m_stabilizationName;

  string m_contactRelationName;
  localIndex m_contactRelationFullIndex;

  integer m_activeSetMaxIter;

  integer m_activeSetIter = 0;

  real64 const m_slidingCheckTolerance = 0.05;

  string const m_tractionKey = viewKeyStruct::tractionString();

  real64 m_initialResidual[3] = {0.0, 0.0, 0.0};

  /**
   * @struct FractureState
   *
   * A struct for the fracture states
   */
  struct FractureState
  {
    static constexpr integer STICK = 0;    ///< element is closed: no jump across the discontinuity
    static constexpr integer SLIP = 1;     ///< element is sliding: no normal jump across the discontinuity, but sliding is allowed for
    static constexpr integer NEW_SLIP = 2; ///< element just starts sliding: no normal jump across the discontinuity, but sliding is allowed
                                           ///< for
    static constexpr integer OPEN = 3;     ///< element is open: no constraints are imposed
  };

  string fractureStateToString( integer const & state ) const
  {
    string stringState;
    switch( state )
    {
      case FractureState::STICK:
      {
        stringState = "stick";
        break;
      }
      case FractureState::SLIP:
      {
        stringState = "slip";
        break;
      }
      case FractureState::NEW_SLIP:
      {
        stringState = "new_slip";
        break;
      }
      case FractureState::OPEN:
      {
        stringState = "open";
        break;
      }
    }
    return stringState;
  }

  void createPreconditioner( DomainPartition const & domain );

public:

  void initializeFractureState( MeshLevel & mesh,
                                string const & fieldName ) const;

  void setFractureStateForElasticStep( DomainPartition & domain ) const;

  bool updateFractureState( DomainPartition & domain ) const;

  void synchronizeFractureState( DomainPartition & domain ) const;

  bool isFractureAllInStickCondition( DomainPartition const & domain ) const;

  void computeFractureStateStatistics( DomainPartition const & domain,
                                       globalIndex & numStick,
                                       globalIndex & numSlip,
                                       globalIndex & numOpen,
                                       bool printAll = false ) const;

  bool isElementInOpenState( FaceElementSubRegion const & subRegion,
                                                 localIndex const kfe ) const;

  void computeRotationMatrices( DomainPartition & domain ) const;

  void computeTolerances( DomainPartition & domain ) const;

  void computeFaceNodalArea( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition,
                             ArrayOfArraysView< localIndex const > const & faceToNodeMap,
                             localIndex const kf0,
                             array1d< real64 > & nodalArea ) const;

  real64 const machinePrecision = std::numeric_limits< real64 >::epsilon();

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static bool compareFractureStates( integer const state0,
                                     integer const state1 )
  {
    return state0 == state1
           || ( state0 == FractureState::NEW_SLIP && state1 == FractureState::SLIP )
           || ( state0 == FractureState::SLIP && state1 == FractureState::NEW_SLIP );
  }
};

} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_LAGRANGIANCONTACTSOLVER_HPP_ */
