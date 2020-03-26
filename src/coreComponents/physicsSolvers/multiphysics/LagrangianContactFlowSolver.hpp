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
 * @file LagrangianContactFlowSolver.hpp
 *
 */

#ifndef GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_LAGRANGIANCONTACTFLOWSOLVER_HPP_
#define GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_LAGRANGIANCONTACTFLOWSOLVER_HPP_

#include "physicsSolvers/SolverBase.hpp"
#include "physicsSolvers/multiphysics/LagrangianContactSolver.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBase.hpp"

namespace geosx
{

class LagrangianContactSolver;
class FlowSolverBase;

class LagrangianContactFlowSolver : public SolverBase
{
public:
  LagrangianContactFlowSolver( const std::string & name,
                               Group * const parent );

  ~LagrangianContactFlowSolver() override;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string CatalogName()
  {
    return "LagrangianContactWithFlow";
  }

  virtual void InitializePreSubGroups( Group * const rootGroup ) override;

  virtual void SetupDofs( DomainPartition const * const domain,
                          DofManager & dofManager ) const override;

  virtual void SetupSystem( DomainPartition * const domain,
                            DofManager & dofManager,
                            ParallelMatrix & matrix,
                            ParallelVector & rhs,
                            ParallelVector & solution ) override;

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

  void AssembleForceResidualDerivativeWrtPressure( DomainPartition * const domain,
                                                   DofManager const & dofManager,
                                                   ParallelMatrix * const matrix,
                                                   ParallelVector * const rhs );

  void AssembleFluidMassResidualDerivativeWrtDisplacement( real64 const dt,
                                                           DomainPartition const * const domain,
                                                           DofManager const & dofManager,
                                                           ParallelMatrix * const matrix,
                                                           ParallelVector * const rhs );

  void AssembleStabiliziation( real64 const dt,
                               DomainPartition const * const domain,
                               DofManager const & dofManager,
                               ParallelMatrix * const matrix,
                               ParallelVector * const rhs );

  real64 SplitOperatorStep( real64 const & time_n,
                            real64 const & dt,
                            integer const cycleNumber,
                            DomainPartition * const domain );

  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    constexpr static auto contactSolverNameString = "contactSolverName";
    constexpr static auto flowSolverNameString = "flowSolverName";
    constexpr static auto stabilizationNameString = "stabilizationName";

//    constexpr static auto stabilizationNameString = "stabilizationName";
//    constexpr static auto contactRelationNameString = "contactRelationName";
//    constexpr static auto activeSetMaxIterString = "activeSetMaxIter";

//    constexpr static auto tractionString = "traction";
//    constexpr static auto deltaTractionString = "deltaTraction";
//    constexpr static auto fractureStateString = "fractureState";
//    constexpr static auto integerFractureStateString = "integerFractureState";
//    constexpr static auto previousFractureStateString = "previousFractureState";
//    constexpr static auto localJumpString = "localJump";
//    constexpr static auto previousLocalJumpString = "previousLocalJump";

    constexpr static auto defaultConductivityString = "defaultConductivity";

//    constexpr static auto slidingCheckToleranceString = "slidingCheckTolerance";
//    constexpr static auto normalDisplacementToleranceString = "normalDisplacementTolerance";
//    constexpr static auto normalTractionToleranceString = "normalTractionTolerance";
//    constexpr static auto slidingToleranceString = "slidingTolerance";

  } LagrangianContactFlowSolverViewKeys;

  //string const & getContactRelationName() const { return m_contactRelationName; }

protected:
  virtual void PostProcessInput() override final;

  virtual void
  InitializePostInitialConditions_PreSubGroups( dataRepository::Group * const problemManager ) override final;

private:

  string m_contactSolverName;
  string m_flowSolverName;

  LagrangianContactSolver * m_contactSolver;
  FlowSolverBase * m_flowSolver;

  string m_stabilizationName;

  real64 m_defaultConductivity;

  integer m_activeSetIter = 0;

  string const m_tractionKey = LagrangianContactSolver::viewKeyStruct::tractionString;
  string const m_fractureStateKey = LagrangianContactSolver::viewKeyStruct::fractureStateString;
  string const m_localJumpKey = LagrangianContactSolver::viewKeyStruct::localJumpString;
  string const m_pressureKey = FlowSolverBase::viewKeyStruct::pressureString;
  string const m_deltaPressureKey = FlowSolverBase::viewKeyStruct::deltaPressureString;

  real64 m_initialResidual[4] = {0.0, 0.0, 0.0, 0.0};

};

} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_LAGRANGIANCONTACTFLOWSOLVER_HPP_ */
