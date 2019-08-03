/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file HydrofractureSolver.hpp
 *
 */

#ifndef HydrofractureSolver_HPP_
#define HydrofractureSolver_HPP_

#include "../SolverBase.hpp"

namespace geosx
{

class FlowSolverBase;
class SolidMechanicsLagrangianFEM;

class HydrofractureSolver : public SolverBase
{
public:
  HydrofractureSolver( const std::string& name,
                       ManagedGroup * const parent );

  ~HydrofractureSolver() override;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string CatalogName()
  {
    return "Hydrofracture";
  }

  virtual void RegisterDataOnMesh( dataRepository::ManagedGroup * const MeshBodies ) override final;

  void SetupSystem ( DomainPartition * const domain,
                     systemSolverInterface::EpetraBlockSystem * const blockSystem );

  void SetSparsityPattern( DomainPartition const * const domain,
                           Epetra_FECrsGraph * const sparsity );

  void SetNumRowsAndTrilinosIndices( MeshLevel * const meshLevel,
                                     localIndex & numLocalRows,
                                     globalIndex & numGlobalRows,
                                     localIndex offset );

  virtual void ImplicitStepSetup( real64 const& time_n,
                                  real64 const& dt,
                                  DomainPartition * const domain,
                                  systemSolverInterface::EpetraBlockSystem * const blockSystem ) override final;

  virtual void ImplicitStepComplete( real64 const& time_n,
                                     real64 const& dt,
                                     DomainPartition * const domain ) override final;

  virtual void AssembleSystem( DomainPartition * const domain,
                               systemSolverInterface::EpetraBlockSystem * const blockSystem,
                               real64 const time,
                               real64 const dt ) override;

  virtual void ApplyBoundaryConditions( DomainPartition * const domain,
                                        systemSolverInterface::EpetraBlockSystem * const blockSystem,
                                        real64 const time_n,
                                        real64 const dt ) override;

  virtual real64
  CalculateResidualNorm( systemSolverInterface::EpetraBlockSystem const *const blockSystem,
                         DomainPartition *const domain) override;

  virtual void SolveSystem( systemSolverInterface::EpetraBlockSystem * const blockSystem,
                            SystemSolverParameters const * const params ) override;

  virtual void
  ApplySystemSolution( systemSolverInterface::EpetraBlockSystem const * const blockSystem,
                       real64 const scalingFactor,
                       DomainPartition * const domain ) override;

  virtual void ResetStateToBeginningOfStep( DomainPartition * const domain ) override;

  virtual real64 SolverStep( real64 const & time_n,
                             real64 const & dt,
                             int const cycleNumber,
                             DomainPartition * domain ) override;

  virtual real64 ExplicitStep( real64 const & time_n,
                               real64 const & dt,
                               integer const cycleNumber,
                               DomainPartition * const domain ) override;

  void UpdateDeformationForCoupling( DomainPartition * const domain );

  void ApplyFractureFluidCoupling( DomainPartition * const domain,
                                   systemSolverInterface::EpetraBlockSystem & blockSystem );

  void AssembleForceResidualDerivativeWrtPressure( DomainPartition * const domain,
                                   systemSolverInterface::EpetraBlockSystem & blockSystem );

  void AssembleFluidMassResidualDerivativeWrtDisplacement( DomainPartition * const domain,
                                                           systemSolverInterface::EpetraBlockSystem & blockSystem );


  real64 SplitOperatorStep( real64 const& time_n,
                            real64 const& dt,
                            integer const cycleNumber,
                            DomainPartition * const domain );

  enum class couplingTypeOption : int
  {
    FixedStress,
    TightlyCoupled
  };

  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    constexpr static auto couplingTypeOptionString = "couplingTypeOptionEnum";
    constexpr static auto couplingTypeOptionStringString = "couplingTypeOption";

    constexpr static auto totalMeanStressString = "totalMeanStress";
    constexpr static auto oldTotalMeanStressString = "oldTotalMeanStress";

    constexpr static auto solidSolverNameString = "solidSolverName";
    constexpr static auto fluidSolverNameString = "fluidSolverName";
  } HydrofractureSolverViewKeys;

protected:
  virtual void PostProcessInput() override final;

  virtual void
  InitializePostInitialConditions_PreSubGroups( dataRepository::ManagedGroup * const problemManager ) override final;

private:

  string m_solidSolverName;
  string m_flowSolverName;
  string m_couplingTypeOptionString;
  couplingTypeOption m_couplingTypeOption;

  SolidMechanicsLagrangianFEM * m_solidSolver;
  FlowSolverBase * m_flowSolver;


};

} /* namespace geosx */

#endif /* HydrofractureSolver_HPP_ */
