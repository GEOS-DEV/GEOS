// Copyright (c) 2018, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-746361. All Rights
// reserved. See file COPYRIGHT for details.
//
// This file is part of the GEOSX Simulation Framework.

//
// GEOSX is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.

/**
 * @file SinglePhaseFlow_TPFA.hpp
 */

#ifndef SINGLE_PHASE_FLOW_TPFA_HPP_
#define SINGLE_PHASE_FLOW_TPFA_HPP_

#include "physicsSolvers/SolverBase.hpp"
#include "systemSolverInterface/LinearSolverWrapper.hpp"


struct stabledt
{
  double m_maxdt;
};

namespace ML_Epetra
{ class MultiLevelPreconditioner; }

namespace geosx
{
namespace dataRepository
{
class ManagedGroup;
}
class BoundaryConditionBase;
class FiniteElementBase;
class DomainPartition;

/**
 * @class SinglePhaseFlow_TPFA
 *
 * class to perform a single phase, two point flux approximation finite volume solve.
 */
class SinglePhaseFlow_TPFA : public SolverBase
{
public:
  SinglePhaseFlow_TPFA( const std::string& name,
                                ManagedGroup * const parent );


  virtual ~SinglePhaseFlow_TPFA();

  static string CatalogName() { return "SinglePhaseFlow_TPFA"; }

  virtual void FillDocumentationNode() override final;

  virtual void FillOtherDocumentationNodes( dataRepository::ManagedGroup * const group ) override final;

  virtual void InitializePreSubGroups( dataRepository::ManagedGroup * const problemManager ) override final;

//  virtual void ReadXML_PostProcess() override final;

  virtual void TimeStep( real64 const& time_n,
                         real64 const& dt,
                         integer const cycleNumber,
                         dataRepository::ManagedGroup * domain ) override;

  void TimeStepExplicit( real64 const& time_n,
                         real64 const& dt,
                         integer const cycleNumber,
                         DomainPartition * domain );

  void TimeStepQuasiStatic( real64 const& time_n,
                            real64 const& dt,
                            integer const cycleNumber,
                            DomainPartition& domain );

  real64 TimeStepImplicit( real64 const & time_n,
                           real64 const & dt,
                           integer const cycleNumber,
                           DomainPartition * const domain );

  void TimeStepImplicitSetup( real64 const& time_n,
                              real64 const& dt,
                              DomainPartition * const domain );

  void TimeStepImplicitComplete( real64 const & time,
                                 real64 const & dt,
                                 DomainPartition * const domain );

  void SetupSystem ( DomainPartition * const domain,
                     systemSolverInterface::EpetraBlockSystem * const blockSystem );

  void SetSparsityPattern( DomainPartition const * const domain,
                           Epetra_FECrsGraph * const sparsity );

  void SetNumRowsAndTrilinosIndices( MeshLevel * const meshLevel,
                                     localIndex & numLocalRows,
                                     localIndex & numGlobalRows,
                                     localIndex_array& localIndices,
                                     localIndex offset );

  void SetupMLPreconditioner( DomainPartition const & domain,
                              ML_Epetra::MultiLevelPreconditioner* MLPrec );


  real64 Assemble ( DomainPartition * const domain,
                    systemSolverInterface::EpetraBlockSystem * const blockSystem,
                    real64 const time,
                    real64 const dt );

  realT CalculateElementResidualAndDerivative( real64 const density,
                                               FiniteElementBase const * const fe,
                                               const Array2dT<R1Tensor>& dNdX,
                                               const realT* const detJ,
                                               R2SymTensor const * const refStress,
                                               array<R1Tensor> const & u,
                                               array<R1Tensor> const & uhat,
                                               array<R1Tensor> const & uhattilde,
                                               array<R1Tensor> const & vtilde,
                                               realT const dt,
                                               Epetra_SerialDenseMatrix& dRdU,
                                               Epetra_SerialDenseVector& R,
                                               real64 c[6][6] );

  void ApplyDirichletBC_implicit( ManagedGroup * object,
                                  BoundaryConditionBase const * const bc,
                                  lSet const & set,
                                  real64 const time_n,
                                  systemSolverInterface::EpetraBlockSystem & blockSystem );

  void ApplySystemSolution( systemSolverInterface::EpetraBlockSystem const * const blockSystem,
                            real64 const scalingFactor,
                            localIndex const dofOffset,
                            dataRepository::ManagedGroup * const nodeManager );


  void MakeGeometryParameters( DomainPartition * const  domain );


  enum class timeIntegrationOption
  {
    SteadyState,
    ImplicitTransient,
    ExplicitTransient
  };

  struct viewKeyStruct : ObjectManagerBase::viewKeyStruct
  {
    constexpr static auto trilinosIndexString = "trilinosIndex_SinglePhaseFlow_TPFA";
    constexpr static auto fluidPressureString = "fluidPressure";
    constexpr static auto deltaFluidPressureString = "deltaFluidPressure";
    constexpr static auto volumeString = "volume";
    constexpr static auto deltaVolumeString = "deltaVolume";
    constexpr static auto porosityString = "Porosity";
    constexpr static auto deltaPorosityString = "deltaPorosity";
    constexpr static auto permeabilityString = "permeablity";
    constexpr static auto faceAreaString = "faceArea";
    constexpr static auto faceCenterString = "faceCenter";

    dataRepository::ViewKey trilinosIndex = { trilinosIndexString };
    dataRepository::ViewKey timeIntegrationOption = { "timeIntegrationOption" };
    dataRepository::ViewKey fieldVarName = { "fieldName" };
    dataRepository::ViewKey functionalSpace = { "functionalSpace" };
  } viewKeys;

  struct groupKeyStruct
  {
    dataRepository::GroupKey systemSolverParameters = { "SystemSolverParameters" };
  } groupKeys;


  systemSolverInterface::LinearSolverWrapper m_linearSolverWrapper;
  systemSolverInterface::EpetraBlockSystem m_linearSystem;

  SystemSolverParameters * getSystemSolverParameters() {return this->GetGroup<SystemSolverParameters>(groupKeys.systemSolverParameters); }

private:

  stabledt m_stabledt;
  timeIntegrationOption m_timeIntegrationOption;
  SinglePhaseFlow_TPFA();
  Array2dT<real64> m_faceToElemLOverA;

  constexpr static int m_dim = 1;
  localIndex_array m_faceConnectors;

};


} /* namespace geosx */

#endif /*  */
