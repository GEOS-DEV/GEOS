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
 * @file SolidMechanicsLagrangianFEM.hpp
 */

#ifndef SOLID_MECHANICS_LAGRANGIAN_FEM_HPP_
#define SOLID_MECHANICS_LAGRANGIAN_FEM_HPP_

#include "physicsSolvers/SolverBase.hpp"
#include "systemSolverInterface/LinearSolverWrapper.hpp"
#include "common/TimingMacros.hpp"
#include "MPI_Communications/CommunicationTools.hpp"

#include "mesh/MeshForLoopInterface.hpp"
//#include "rajaInterface/GEOSX_RAJA_Interface.hpp"

namespace ML_Epetra
{ class MultiLevelPreconditioner; }

namespace geosx
{
namespace dataRepository
{
class ManagedGroup;
}
class FieldSpecificationBase;
class FiniteElementBase;
class DomainPartition;

class SolidMechanicsLagrangianFEM : public SolverBase
{
public:
  SolidMechanicsLagrangianFEM( const std::string& name,
                                ManagedGroup * const parent );


  SolidMechanicsLagrangianFEM( SolidMechanicsLagrangianFEM const & ) = delete;
  SolidMechanicsLagrangianFEM( SolidMechanicsLagrangianFEM && ) = default ;

  SolidMechanicsLagrangianFEM& operator=( SolidMechanicsLagrangianFEM const & ) = delete;
  SolidMechanicsLagrangianFEM& operator=( SolidMechanicsLagrangianFEM && ) = delete ;

  virtual ~SolidMechanicsLagrangianFEM() override;

  static string CatalogName() { return "SolidMechanics_LagrangianFEM"; }

  virtual void InitializePreSubGroups(ManagedGroup * const rootGroup) override;

  virtual void RegisterDataOnMesh( ManagedGroup * const MeshBody ) override final;

  virtual real64 SolverStep( real64 const& time_n,
                         real64 const& dt,
                         integer const cycleNumber,
                         DomainPartition * domain ) override;


  void TimeStepQuasiStatic( real64 const& time_n,
                            real64 const& dt,
                            integer const cycleNumber,
                            DomainPartition& domain );

  real64 TimeStepImplicit( real64 const & time_n,
                           real64 const & dt,
                           integer const cycleNumber,
                           DomainPartition * const domain );


  void SetupSystem ( DomainPartition * const domain,
                     systemSolverInterface::EpetraBlockSystem * const blockSystem );

  void SetSparsityPattern( DomainPartition const * const domain,
                           Epetra_FECrsGraph * const sparsity );

  void SetNumRowsAndTrilinosIndices( ManagedGroup * const domain,
                                     localIndex & numLocalRows,
                                     globalIndex & numGlobalRows,
                                     localIndex_array& localIndices,
                                     localIndex offset );

  void SetupMLPreconditioner( DomainPartition const & domain,
                              ML_Epetra::MultiLevelPreconditioner* MLPrec );

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/
  virtual real64 ExplicitStep( real64 const& time_n,
                                   real64 const& dt,
                                   integer const cycleNumber,
                                   DomainPartition * domain ) override;

  virtual void ImplicitStepSetup( real64 const& time_n,
                              real64 const& dt,
                              DomainPartition * const domain,
                              systemSolverInterface::EpetraBlockSystem * const blockSystem ) override;

  virtual void AssembleSystem ( DomainPartition * const domain,
                                  systemSolverInterface::EpetraBlockSystem * const blockSystem,
                                  real64 const time,
                                  real64 const dt ) override;

  virtual void SolveSystem( systemSolverInterface::EpetraBlockSystem * const blockSystem,
                            SystemSolverParameters const * const params ) override;

  virtual void ApplySystemSolution( systemSolverInterface::EpetraBlockSystem const * const blockSystem,
                            real64 const scalingFactor,
                            DomainPartition * const domain  ) override;

  virtual void ApplyBoundaryConditions( DomainPartition * const domain,
                                        systemSolverInterface::EpetraBlockSystem * const blockSystem,
                                        real64 const time,
                                        real64 const dt ) override;

  virtual real64
  CalculateResidualNorm(systemSolverInterface::EpetraBlockSystem const *const blockSystem, DomainPartition *const domain) override;

  virtual void ResetStateToBeginningOfStep( DomainPartition * const domain ) override;

  virtual void ImplicitStepComplete( real64 const & time,
                                 real64 const & dt,
                                 DomainPartition * const domain ) override;

  /**@}*/
  virtual real64
  ExplicitElementKernelLaunchSelector( localIndex const er,
                                       localIndex const esr,
                                       set<localIndex> const & elementList,
                                       arrayView2d<localIndex> const & elemsToNodes,
                                       arrayView3d< R1Tensor > const & dNdX,
                                       arrayView2d<real64> const & detJ,
                                       arrayView1d<R1Tensor> const & u,
                                       arrayView1d<R1Tensor> const & uhat,
                                       arrayView1d<R1Tensor> const & vel,
                                       arrayView1d<R1Tensor> & acc,
                                       ElementRegionManager::ConstitutiveRelationAccessor<constitutive::ConstitutiveBase>& constitutiveRelations,
                                       ElementRegionManager::MaterialViewAccessor< arrayView2d<real64> > const & meanStress,
                                       ElementRegionManager::MaterialViewAccessor< arrayView2d<R2SymTensor> > const & devStress,
                                       real64 const dt,
                                       localIndex NUM_NODES_PER_ELEM,
                                       localIndex NUM_QUADRATURE_POINTS );

  template< typename LEAFCLASS >
  real64
  ExplicitElementKernelLaunchSelectorT( localIndex const er,
                                        localIndex const esr,
                                        set<localIndex> const & elementList,
                                        arrayView2d<localIndex> const & elemsToNodes,
                                        arrayView3d< R1Tensor > const & dNdX,
                                        arrayView2d<real64> const & detJ,
                                        arrayView1d<R1Tensor> const & u,
                                        arrayView1d<R1Tensor> const & uhat,
                                        arrayView1d<R1Tensor> const & vel,
                                        arrayView1d<R1Tensor> & acc,
                                        ElementRegionManager::ConstitutiveRelationAccessor<constitutive::ConstitutiveBase>& constitutiveRelations,
                                        ElementRegionManager::MaterialViewAccessor< arrayView2d<real64> > const & meanStress,
                                        ElementRegionManager::MaterialViewAccessor< arrayView2d<R2SymTensor> > const & devStress,
                                        real64 const dt,
                                        localIndex NUM_NODES_PER_ELEM,
                                        localIndex NUM_QUADRATURE_POINTS );

  template< localIndex NUM_NODES_PER_ELEM, localIndex NUM_QUADRATURE_POINTS >
  static real64
  ExplicitElementKernelLaunch( localIndex const er,
                               localIndex const esr,
                               set<localIndex> const & elementList,
                               arrayView2d<localIndex> const & elemsToNodes,
                               arrayView3d< R1Tensor > const & dNdX,
                               arrayView2d<real64> const & detJ,
                               arrayView1d<R1Tensor> const & u,
                               arrayView1d<R1Tensor> const & uhat,
                               arrayView1d<R1Tensor> const & vel,
                               arrayView1d<R1Tensor> & acc,
                               ElementRegionManager::ConstitutiveRelationAccessor<constitutive::ConstitutiveBase> constitutiveRelations,
                               ElementRegionManager::MaterialViewAccessor< arrayView2d<real64> > const & meanStress,
                               ElementRegionManager::MaterialViewAccessor< arrayView2d<R2SymTensor> > const & devStress,
                               real64 const dt );

  realT ElementResidualAndJacobianKernel( real64 const density,
                                          FiniteElementBase const * const fe,
                                          arraySlice2d<R1Tensor const> const& dNdX,
                                          arraySlice1d<realT const> const& detJ,
                                          R2SymTensor const * const refStress,
                                          r1_array const& u,
                                          r1_array const& uhat,
                                          r1_array const& uhattilde,
                                          r1_array const& vtilde,
                                          realT const dt,
                                          Epetra_SerialDenseMatrix& dRdU,
                                          Epetra_SerialDenseVector& R,
                                          real64 c[6][6] );








  void ApplyDisplacementBC_implicit( real64 const time,
                                     DomainPartition & domain,
                                     systemSolverInterface::EpetraBlockSystem & blockSystem  );

  void ForceBC( dataRepository::ManagedGroup * const object,
                FieldSpecificationBase const* const bc,
                set<localIndex> const & set,
                real64 time,
                systemSolverInterface::EpetraBlockSystem & blockSystem );


  void ApplyTractionBC( DomainPartition * const domain,
                        real64 const time,
                        systemSolverInterface::EpetraBlockSystem & blockSystem );


  enum class timeIntegrationOption : int
  {
    QuasiStatic,
    ImplicitDynamic,
    ExplicitDynamic
  };

  void SetTimeIntegrationOption( string const & stringVal )
  {
    if( stringVal == "ExplicitDynamic" )
    {
      this->m_timeIntegrationOption = timeIntegrationOption::ExplicitDynamic;
    }
    else if( stringVal == "ImplicitDynamic" )
    {
      this->m_timeIntegrationOption = timeIntegrationOption::ImplicitDynamic;
    }
    else if ( stringVal == "QuasiStatic" )
    {
      this->m_timeIntegrationOption = timeIntegrationOption::QuasiStatic;
    }
    else
    {
      GEOS_ERROR("Invalid time integration option: " << stringVal);
    }
  }

  struct viewKeyStruct
  {
    static constexpr auto vTildeString = "velocityTilde";
    static constexpr auto uhatTildeString = "uhatTilde";
    static constexpr auto cflFactorString = "cflFactor";
    static constexpr auto newmarkGammaString = "newmarkGamma";
    static constexpr auto newmarkBetaString = "newmarkBeta";
    static constexpr auto massDampingString = "massDamping";
    static constexpr auto stiffnessDampingString = "stiffnessDamping";
    static constexpr auto useVelocityEstimateForQSString = "useVelocityForQS";
    static constexpr auto trilinosIndexString = "trilinosIndex";
    static constexpr auto timeIntegrationOptionStringString = "timeIntegrationOption";
    static constexpr auto timeIntegrationOptionString = "timeIntegrationOptionEnum";
    static constexpr auto maxNumResolvesString = "maxNumResolves";
    static constexpr auto strainTheoryString = "strainTheory";


    dataRepository::ViewKey vTilde = { vTildeString };
    dataRepository::ViewKey uhatTilde = { uhatTildeString };
    dataRepository::ViewKey newmarkGamma = { newmarkGammaString };
    dataRepository::ViewKey newmarkBeta = { newmarkBetaString };
    dataRepository::ViewKey massDamping = { massDampingString };
    dataRepository::ViewKey stiffnessDamping = { stiffnessDampingString };
    dataRepository::ViewKey useVelocityEstimateForQS = { useVelocityEstimateForQSString };
    dataRepository::ViewKey trilinosIndex = { trilinosIndexString };
    dataRepository::ViewKey timeIntegrationOption = { timeIntegrationOptionString };
  } solidMechanicsViewKeys;

  struct groupKeyStruct
  {
    dataRepository::GroupKey systemSolverParameters = { "SystemSolverParameters" };
  } solidMechanicsGroupKeys;

protected:
  virtual void PostProcessInput() override final;

  virtual void InitializePostInitialConditions_PreSubGroups( dataRepository::ManagedGroup * const problemManager ) override final;


protected:

  real64 m_newmarkGamma;
  real64 m_newmarkBeta;
  real64 m_massDamping;
  real64 m_stiffnessDamping;
  string m_timeIntegrationOptionString;
  timeIntegrationOption m_timeIntegrationOption;
  integer m_useVelocityEstimateForQS;
  real64 m_maxForce = 0.0;
  integer m_maxNumResolves;
  integer m_strainTheory;

  array1d< array1d < set<localIndex> > > m_elemsAttachedToSendOrReceiveNodes;
  array1d< array1d < set<localIndex> > > m_elemsNotAttachedToSendOrReceiveNodes;
  set<localIndex> m_sendOrRecieveNodes;
  set<localIndex> m_nonSendOrRecieveNodes;
  MPI_iCommData m_icomm;

  SolidMechanicsLagrangianFEM();

};



template< int N >
void Integrate( const R2SymTensor& fieldvar,
                arraySlice1d<R1Tensor const> const & dNdX,
                real64 const& detJ,
                real64 const& detF,
                const R2Tensor& Finv,
                arraySlice1d<R1Tensor> & result)
{
  real64 const integrationFactor = detJ * detF;

  R2Tensor P;
  P.AijBkj( fieldvar,Finv);
  P *= integrationFactor;

  for( int a=0 ; a<N ; ++a )  // loop through all shape functions in element
  {
    result[a].minusAijBj( P, dNdX[a] );
  }
}


template< typename LEAFCLASS >
real64 SolidMechanicsLagrangianFEM::
ExplicitElementKernelLaunchSelectorT( localIndex const er,
                                      localIndex const esr,
                                      set<localIndex> const & elementList,
                                      arrayView2d<localIndex> const & elemsToNodes,
                                      arrayView3d< R1Tensor > const & dNdX,
                                      arrayView2d<real64> const & detJ,
                                      arrayView1d<R1Tensor> const & u,
                                      arrayView1d<R1Tensor> const & uhat,
                                      arrayView1d<R1Tensor> const & vel,
                                      arrayView1d<R1Tensor> & acc,
                                      ElementRegionManager::ConstitutiveRelationAccessor<constitutive::ConstitutiveBase>& constitutiveRelations,
                                      ElementRegionManager::MaterialViewAccessor< arrayView2d<real64> > const & meanStress,
                                      ElementRegionManager::MaterialViewAccessor< arrayView2d<R2SymTensor> > const & devStress,
                                      real64 const dt,
                                      localIndex NUM_NODES_PER_ELEM,
                                      localIndex NUM_QUADRATURE_POINTS )
{
  real64 rval = 0;

  if( NUM_NODES_PER_ELEM==8 && NUM_QUADRATURE_POINTS==8 )
  {
    rval = LEAFCLASS::template ExplicitElementKernelLaunch<8,8>( er,
                                                                 esr,
                                                                 elementList,
                                                                 elemsToNodes,
                                                                 dNdX,
                                                                 detJ,
                                                                 u,
                                                                 uhat,
                                                                 vel,
                                                                 acc,
                                                                 constitutiveRelations,
                                                                 meanStress,
                                                                 devStress,
                                                                 dt );
  }

  return rval;
}

} /* namespace geosx */

#endif /* SOLID_MECHANICS_LAGRANGIAN_FEM_HPP_ */
