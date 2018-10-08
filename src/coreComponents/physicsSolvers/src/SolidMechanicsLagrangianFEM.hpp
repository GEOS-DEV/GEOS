/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
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

#include "mesh/MeshForLoopInterface.hpp"
//#include "rajaInterface/GEOSX_RAJA_Interface.hpp"


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

class SolidMechanics_LagrangianFEM : public SolverBase
{
public:
  SolidMechanics_LagrangianFEM( const std::string& name,
                                ManagedGroup * const parent );


  SolidMechanics_LagrangianFEM( SolidMechanics_LagrangianFEM const & ) = delete;
  SolidMechanics_LagrangianFEM( SolidMechanics_LagrangianFEM && ) = default ;

  SolidMechanics_LagrangianFEM& operator=( SolidMechanics_LagrangianFEM const & ) = delete;
  SolidMechanics_LagrangianFEM& operator=( SolidMechanics_LagrangianFEM && ) = delete ;

  virtual ~SolidMechanics_LagrangianFEM() override;

  static string CatalogName() { return "SolidMechanics_LagrangianFEM"; }

  virtual void FillDocumentationNode() override final;
  
  virtual void FillOtherDocumentationNodes( dataRepository::ManagedGroup * const group ) override final;

  virtual void FinalInitializationPreSubGroups( dataRepository::ManagedGroup * const problemManager ) override final;

  virtual void ReadXML_PostProcess() override final;

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

  template< localIndex NUM_NODES_PER_ELEM, localIndex NUM_QUADRATURE_POINTS >
  real64 ExplicitElementKernel( localIndex const er,
                                localIndex const esr,
                                set<localIndex> const & elementList,
                                array2d<localIndex> const & elemsToNodes,
                                array3d< R1Tensor > const & dNdX,
                                array2d<real64> const & detJ,
                                r1_array const & u,
                                r1_array const & uhat,
                                r1_array const & acc,
                                ElementRegionManager::ConstitutiveRelationAccessor<constitutive::ConstitutiveBase> constitutiveRelations,
                                ElementRegionManager::MaterialViewAccessor< array2d<real64> > meanStress,
                                ElementRegionManager::MaterialViewAccessor< array2d<R2SymTensor> > devStress,
                                real64 const dt );


  realT CalculateElementResidualAndDerivative( real64 const density,
                                               FiniteElementBase const * const fe,
                                               const array_view<R1Tensor,2>& dNdX,
                                               const realT* const detJ,
                                               R2SymTensor const * const refStress,
                                               array1d<R1Tensor> const & u,
                                               array1d<R1Tensor> const & uhat,
                                               array1d<R1Tensor> const & uhattilde,
                                               array1d<R1Tensor> const & vtilde,
                                               realT const dt,
                                               Epetra_SerialDenseMatrix& dRdU,
                                               Epetra_SerialDenseVector& R,
                                               real64 c[6][6] );

  void ApplyDisplacementBC_implicit( real64 const time,
                                     DomainPartition & domain,
                                     systemSolverInterface::EpetraBlockSystem & blockSystem  );

  void ForceBC( dataRepository::ManagedGroup * const object,
                BoundaryConditionBase const* const bc,
                set<localIndex> const & set,
                real64 time,
                systemSolverInterface::EpetraBlockSystem & blockSystem );

  void TractionBC( dataRepository::ManagedGroup * const object,
                   BoundaryConditionBase const* const bc,
                   set<localIndex> const & set,
                   real64 time,
                   systemSolverInterface::EpetraBlockSystem & blockSystem );


  enum class timeIntegrationOption
  {
    QuasiStatic,
    ImplicitDynamic,
    ExplicitDynamic
  };

  struct viewKeyStruct
  {
    dataRepository::ViewKey vTilde = { "velocityTilde" };
    dataRepository::ViewKey uhatTilde = { "uhatTilde" };
    dataRepository::ViewKey newmarkGamma = { "newmarkGamma" };
    dataRepository::ViewKey newmarkBeta = { "newmarkBeta" };
    dataRepository::ViewKey massDamping = { "massDamping" };
    dataRepository::ViewKey stiffnessDamping = { "stiffnessDamping" };
    dataRepository::ViewKey useVelocityEstimateForQS = { "useVelocityEstimateForQuasiStatic" };
    dataRepository::ViewKey trilinosIndex = { "trilinosIndex" };
    dataRepository::ViewKey ghostRank = { "ghostRank" };
    dataRepository::ViewKey timeIntegrationOption = { "timeIntegrationOption" };


  } solidMechanicsViewKeys;

  struct groupKeyStruct
  {
    dataRepository::GroupKey systemSolverParameters = { "SystemSolverParameters" };
  } solidMechanicsGroupKeys;

private:

  real64 m_maxForce;
  stabledt m_stabledt;
  timeIntegrationOption m_timeIntegrationOption;

  array1d< array1d < set<localIndex> > > m_elemsAttachedToSendOrReceiveNodes;
  array1d< array1d < set<localIndex> > > m_elemsNotAttachedToSendOrReceiveNodes;

  SolidMechanics_LagrangianFEM();

};

namespace Integration
{

template<typename T>
#ifdef USE_CONTAINERARRAY_RETURN_PTR
void OnePoint( T const * const __restrict__ dydx,
               T * const __restrict__ dy,
               T * const __restrict__ y,
               real64 const dx,
               localIndex const length )
#else
void OnePoint( T const & dydx,
               T & dy,
               T & y,
               real64 const dx,
               localIndex const length )
#endif
{

  forall_in_range(0, length, GEOSX_LAMBDA (localIndex a){
      
    dy[a][0] = dydx[a][0] * dx;
    dy[a][1] = dydx[a][1] * dx;
    dy[a][2] = dydx[a][2] * dx;
    
    y[a][0] += dy[a][0];
    y[a][1] += dy[a][1];
    y[a][2] += dy[a][2];
    
  });

}

template<typename T, typename U>
void OnePoint (U dydx,
               T dy_1, T dy_2, T dy_3,
               T y_1, T y_2, T y_3,
               real64 const dx,
               localIndex const length)
{
  
  forall_in_range(0, length, GEOSX_LAMBDA (localIndex a) {
      
    dy_1[a] = dydx[a][0] * dx;
    dy_2[a] = dydx[a][1] * dx;
    dy_3[a] = dydx[a][2] * dx;
    
    y_1[a] += dy_1[a];
    y_2[a] += dy_2[a];
    y_3[a] += dy_3[a];
  });
  
}



template<typename T>
#ifdef USE_CONTAINERARRAY_RETURN_PTR
void OnePoint( T const * const __restrict__ dydx,
               T * const __restrict__ y,
               real64 const dx,
               localIndex const length )
#else
void OnePoint( T const &  dydx,
               T & y,
               real64 const dx,
               localIndex const length )
#endif
{

  forall_in_range(0, length, GEOSX_LAMBDA (localIndex a) {
      
    y[a].plus_cA( dx, dydx[a] );
    
    });

}

//Three arrays present
template<typename T, typename U>
void OnePoint( T const dydx_0,
               T const dydx_1,
               T const dydx_2,
               U &  y,
               real64 const dx,
               localIndex const length )
{
  
  forall_in_range(0, length, GEOSX_LAMBDA (localIndex a) {
    //y[a].plus_cA( dx, dydx[a] );
    y[a][0] += dx*dydx_0[a];
    y[a][1] += dx*dydx_1[a];
    y[a][2] += dx*dydx_2[a];
  });

}




} // namespace integration

} /* namespace geosx */

#endif /* SOLID_MECHANICS_LAGRANGIAN_FEM_HPP_ */
