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
  LagrangianContactSolver( const std::string & name,
                           Group * const parent );

  ~LagrangianContactSolver() override;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string CatalogName()
  {
    return "LagrangianContact";
  }

  virtual void InitializePreSubGroups( Group * const rootGroup ) override;

  virtual void RegisterDataOnMesh( dataRepository::Group * const MeshBodies ) override final;

  virtual void SetupDofs( DomainPartition const * const domain,
                          DofManager & dofManager ) const override;

  virtual void SetupSystem( DomainPartition * const domain,
                            DofManager & dofManager,
                            ParallelMatrix & matrix,
                            ParallelVector & rhs,
                            ParallelVector & solution,
                            bool const setSparsity = true ) override;

  virtual void
  ImplicitStepSetup( real64 const & time_n,
                     real64 const & dt,
                     DomainPartition * const domain,
                     DofManager & dofManager,
                     ParallelMatrix & matrix,
                     ParallelVector & rhs,
                     ParallelVector & solution ) override final;

  virtual void ImplicitStepComplete( real64 const & time_n,
                                     real64 const & dt,
                                     DomainPartition * const domain ) override final;

  virtual void AssembleSystem( real64 const time,
                               real64 const dt,
                               DomainPartition * const domain,
                               DofManager const & dofManager,
                               ParallelMatrix & matrix,
                               ParallelVector & rhs ) override;

  virtual void ApplyBoundaryConditions( real64 const time,
                                        real64 const dt,
                                        DomainPartition * const domain,
                                        DofManager const & dofManager,
                                        ParallelMatrix & matrix,
                                        ParallelVector & rhs ) override;

  virtual real64
  CalculateResidualNorm( DomainPartition const * const domain,
                         DofManager const & dofManager,
                         ParallelVector const & rhs ) override;

  virtual void SolveSystem( DofManager const & dofManager,
                            ParallelMatrix & matrix,
                            ParallelVector & rhs,
                            ParallelVector & solution ) override;

  virtual void
  ApplySystemSolution( DofManager const & dofManager,
                       ParallelVector const & solution,
                       real64 const scalingFactor,
                       DomainPartition * const domain ) override;

  virtual void ResetStateToBeginningOfStep( DomainPartition * const domain ) override;

  virtual real64 SolverStep( real64 const & time_n,
                             real64 const & dt,
                             int const cycleNumber,
                             DomainPartition * const domain ) override;

  virtual void SetNextDt( real64 const & currentDt,
                          real64 & nextDt ) override;


  virtual real64 ExplicitStep( real64 const & time_n,
                               real64 const & dt,
                               integer const cycleNumber,
                               DomainPartition * const domain ) override;

  virtual real64 NonlinearImplicitStep( real64 const & time_n,
                                        real64 const & dt,
                                        integer const cycleNumber,
                                        DomainPartition * const domain,
                                        DofManager const & dofManager,
                                        ParallelMatrix & matrix,
                                        ParallelVector & rhs,
                                        ParallelVector & solution ) override;

  virtual bool LineSearch( real64 const & time_n,
                           real64 const & dt,
                           integer const GEOSX_UNUSED_PARAM( cycleNumber ),
                           DomainPartition * const domain,
                           DofManager const & dofManager,
                           ParallelMatrix & matrix,
                           ParallelVector & rhs,
                           ParallelVector const & solution,
                           real64 const scaleFactor,
                           real64 & lastResidual ) override;

  void UpdateDeformationForCoupling( DomainPartition * const domain );

  void AssembleForceResidualDerivativeWrtTraction( DomainPartition * const domain,
                                                   DofManager const & dofManager,
                                                   ParallelMatrix * const matrix,
                                                   ParallelVector * const rhs );

  void AssembleTractionResidualDerivativeWrtDisplacementAndTraction( DomainPartition const * const domain,
                                                                     DofManager const & dofManager,
                                                                     ParallelMatrix * const matrix,
                                                                     ParallelVector * const rhs );

  void AssembleStabilization( DomainPartition const * const domain,
                              DofManager const & dofManager,
                              ParallelMatrix * const matrix,
                              ParallelVector * const rhs );

  real64 SplitOperatorStep( real64 const & time_n,
                            real64 const & dt,
                            integer const cycleNumber,
                            DomainPartition * const domain );

  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    constexpr static auto solidSolverNameString = "solidSolverName";
    constexpr static auto stabilizationNameString = "stabilizationName";
    constexpr static auto contactRelationNameString = "contactRelationName";
    constexpr static auto activeSetMaxIterString = "activeSetMaxIter";

    constexpr static auto tractionString = "traction";
    constexpr static auto deltaTractionString = "deltaTraction";
    constexpr static auto fractureStateString = "fractureState";
    constexpr static auto previousFractureStateString = "previousFractureState";
    constexpr static auto localJumpString = "localJump";
    constexpr static auto previousLocalJumpString = "previousLocalJump";

    constexpr static auto slidingCheckToleranceString = "slidingCheckTolerance";
    constexpr static auto normalDisplacementToleranceString = "normalDisplacementTolerance";
    constexpr static auto normalTractionToleranceString = "normalTractionTolerance";
    constexpr static auto slidingToleranceString = "slidingTolerance";

  } LagrangianContactSolverViewKeys;

  string const & getContactRelationName() const { return m_contactRelationName; }

protected:
  virtual void PostProcessInput() override final;

  virtual void
  InitializePostInitialConditions_PreSubGroups( dataRepository::Group * const problemManager ) override final;

private:

  string m_solidSolverName;
  SolidMechanicsLagrangianFEM * m_solidSolver;

  string m_stabilizationName;

  string m_contactRelationName;
  localIndex m_contactRelationFullIndex;

  integer m_activeSetMaxIter;

  integer m_activeSetIter = 0;

  real64 const m_slidingCheckTolerance = 0.05;
  //real64 m_normalDisplacementTolerance = 1.e-7;
  //real64 m_normalTractionTolerance = 1.e-4;
  //real64 m_slidingTolerance = 1.e-7;

  string const m_tractionKey = viewKeyStruct::tractionString;

  real64 m_initialResidual[3] = {0.0, 0.0, 0.0};

  /**
   * @struct FractureState
   *
   * A struct for the fracture states
   */
  struct FractureState
  {
    static constexpr integer STICK = 0;    ///< element is closed: no jump across the discontinuity
    static constexpr integer SLIP = 1;     ///< element is sliding: no normal jump across the discontinuity, but sliding
                                           ///< is allowed for
    static constexpr integer NEW_SLIP = 2; ///< element just starts sliding: no normal jump across the discontinuity,
                                           ///< but sliding is allowed for
    static constexpr integer OPEN = 3;     ///< element is open: no constraints are imposed
  };

  string FractureStateToString( integer const & state ) const
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

  bool CompareFractureStates( integer const & state0, integer const & state1 ) const
  {
    if( state0 == state1 )
    {
      return true;
    }
    else if( state0 == FractureState::NEW_SLIP && state1 == FractureState::SLIP )
    {
      return true;
    }
    else if( state0 == FractureState::SLIP && state1 == FractureState::NEW_SLIP )
    {
      return true;
    }
    return false;
  }

  void ComputeTolerances( DomainPartition * const domain ) const;

public:

  void InitializeFractureState( MeshLevel * const mesh,
                                string const fieldName ) const;

  void SetFractureStateForElasticStep( DomainPartition * const domain ) const;

  bool UpdateFractureState( DomainPartition * const domain ) const;

  void SynchronizeFractureState( DomainPartition * const domain ) const;

  bool IsFractureAllInStickCondition( DomainPartition const * const domain ) const;

  void ComputeFractureStateStatistics( DomainPartition const * const domain,
                                       globalIndex & numStick,
                                       globalIndex & numSlip,
                                       globalIndex & numOpen,
                                       bool printAll = false ) const;
};

} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_LAGRANGIANCONTACTSOLVER_HPP_ */
