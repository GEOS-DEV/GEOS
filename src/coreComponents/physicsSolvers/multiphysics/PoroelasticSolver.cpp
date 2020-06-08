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
 * @file PoroelasticSolver.cpp
 *
 */


#include "PoroelasticSolver.hpp"

#include "common/DataLayouts.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/solid/PoroElastic.hpp"
#include "constitutive/fluid/SingleFluidBase.hpp"
#include "managers/NumericalMethodsManager.hpp"
#include "finiteElement/Kinematics.h"
#include "linearAlgebra/solvers/BlockPreconditioner.hpp"
#include "linearAlgebra/solvers/SeparateComponentPreconditioner.hpp"
#include "managers/DomainPartition.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

PoroelasticSolver::PoroelasticSolver( const std::string & name,
                                      Group * const parent ):
  SolverBase( name, parent ),
  m_solidSolverName(),
  m_flowSolverName(),
  m_couplingTypeOptionString(),
  m_couplingTypeOption()

{
  registerWrapper( viewKeyStruct::solidSolverNameString, &m_solidSolverName )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Name of the solid mechanics solver to use in the poroelastic solver" );

  registerWrapper( viewKeyStruct::fluidSolverNameString, &m_flowSolverName )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Name of the fluid mechanics solver to use in the poroelastic solver" );

  registerWrapper( viewKeyStruct::couplingTypeOptionStringString, &m_couplingTypeOptionString )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Coupling option: (FIM, SIM_FixedStress)" );
}

void PoroelasticSolver::RegisterDataOnMesh( dataRepository::Group * const MeshBodies )
{
  for( auto & mesh : MeshBodies->GetSubGroups() )
  {
    ElementRegionManager * const elemManager = mesh.second->group_cast< MeshBody * >()->getMeshLevel( 0 )->getElemManager();

    elemManager->forElementSubRegions< CellElementSubRegion, FaceElementSubRegion >( [&]( ElementSubRegionBase & elementSubRegion )
    {
      elementSubRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::totalMeanStressString )->
        setDescription( "Total Mean Stress" );
      elementSubRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::oldTotalMeanStressString )->
        setDescription( "Total Mean Stress" );
    } );
  }
}

void PoroelasticSolver::SetupDofs( DomainPartition const * const domain,
                                   DofManager & dofManager ) const
{
  GEOSX_MARK_FUNCTION;
  m_solidSolver->SetupDofs( domain, dofManager );
  m_flowSolver->SetupDofs( domain, dofManager );

  dofManager.addCoupling( keys::TotalDisplacement,
                          FlowSolverBase::viewKeyStruct::pressureString,
                          DofManager::Connector::Elem );
}

void PoroelasticSolver::SetupSystem( DomainPartition * const domain,
                                     DofManager & dofManager,
                                     ParallelMatrix & matrix,
                                     ParallelVector & rhs,
                                     ParallelVector & solution )
{
  // setup monolithic coupled system
  dofManager.setMesh( domain, 0, 0 );
  SetupDofs( domain, dofManager );
  dofManager.reorderByRank();

  if( m_precond )
  {
    m_precond->clear();
  }

  localIndex const numLocalDof = dofManager.numLocalDofs();
  matrix.createWithLocalSize( numLocalDof, numLocalDof, 30, MPI_COMM_GEOSX );
  rhs.createWithLocalSize( numLocalDof, MPI_COMM_GEOSX );
  solution.createWithLocalSize( numLocalDof, MPI_COMM_GEOSX );

  dofManager.setSparsityPattern( m_matrix, true );

  if( !m_precond && m_linearSolverParameters.get().solverType != "direct" )
  {
    CreatePreconditioner();
  }
}

void PoroelasticSolver::ImplicitStepSetup( real64 const & time_n,
                                           real64 const & dt,
                                           DomainPartition * const domain,
                                           DofManager & dofManager,
                                           ParallelMatrix & matrix,
                                           ParallelVector & rhs,
                                           ParallelVector & solution )
{
  if( m_couplingTypeOption == couplingTypeOption::SIM_FixedStress )
  {
    MeshLevel * const mesh = domain->getMeshBodies()->GetGroup< MeshBody >( 0 )->getMeshLevel( 0 );
    ElementRegionManager * const elemManager = mesh->getElemManager();

    ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > const totalMeanStress =
      elemManager->ConstructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( viewKeyStruct::totalMeanStressString );

    ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > oldTotalMeanStress =
      elemManager->ConstructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( viewKeyStruct::oldTotalMeanStressString );

    //***** loop over all elements and initialize the derivative arrays *****
    forAllElemsInMesh( mesh, [&]( localIndex const er,
                                  localIndex const esr,
                                  localIndex const k )->void
    {
      oldTotalMeanStress[er][esr][k] = totalMeanStress[er][esr][k];
    } );
  }
  else if( m_couplingTypeOption == couplingTypeOption::FIM )
  {
    m_flowSolver->ImplicitStepSetup( time_n, dt, domain,
                                     dofManager,
                                     matrix,
                                     rhs,
                                     solution );
    m_solidSolver->ImplicitStepSetup( time_n, dt, domain,
                                      dofManager,
                                      matrix,
                                      rhs,
                                      solution );
  }
}

void PoroelasticSolver::ImplicitStepComplete( real64 const & time_n,
                                              real64 const & dt,
                                              DomainPartition * const domain )
{
  if( m_couplingTypeOption == couplingTypeOption::FIM )
  {
    m_solidSolver->updateStress( domain ); // TODO: to be moved in m_solidSolver->ImplicitStepComplete
    m_solidSolver->ImplicitStepComplete( time_n, dt, domain );
    m_flowSolver->ImplicitStepComplete( time_n, dt, domain );
  }
}

void PoroelasticSolver::PostProcessInput()
{
  SolverBase::PostProcessInput();

  m_flowSolver  = this->getParent()->GetGroup< SinglePhaseBase >( m_flowSolverName );
  m_solidSolver = this->getParent()->GetGroup< SolidMechanicsLagrangianFEM >( m_solidSolverName );

  GEOSX_ERROR_IF( m_flowSolver == nullptr, "Flow solver not found or invalid type: " << m_flowSolverName );
  GEOSX_ERROR_IF( m_solidSolver == nullptr, "Solid solver not found or invalid type: " << m_solidSolverName );

  string ctOption = this->getReference< string >( viewKeyStruct::couplingTypeOptionStringString );

  if( ctOption == "SIM_FixedStress" )
  {
    m_couplingTypeOption = couplingTypeOption::SIM_FixedStress;

    // For this coupled solver the minimum number of Newton Iter should be 0 for both flow and solid solver,
    // otherwise it will never converge.
    m_flowSolver->getNonlinearSolverParameters().m_minIterNewton = 0;
    m_solidSolver->getNonlinearSolverParameters().m_minIterNewton = 0;
  }
  else if( ctOption == "FIM" )
  {
    m_couplingTypeOption = couplingTypeOption::FIM;
  }
  else
  {
    GEOSX_ERROR( "invalid coupling type option: " + ctOption );
  }

}

void PoroelasticSolver::InitializePostInitialConditions_PreSubGroups( Group * const problemManager )
{
  if( m_couplingTypeOption == couplingTypeOption::SIM_FixedStress )
  {
    this->getParent()->GetGroup( m_flowSolverName )->group_cast< SinglePhaseBase * >()->setPoroElasticCoupling();
    // Calculate initial total mean stress
    this->UpdateDeformationForCoupling( problemManager->GetGroup< DomainPartition >( keys::domain ));
  }
}

PoroelasticSolver::~PoroelasticSolver()
{
  // TODO Auto-generated destructor stub
}

void PoroelasticSolver::ResetStateToBeginningOfStep( DomainPartition * const domain )
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup< MeshBody >( 0 )->getMeshLevel( 0 );
  ElementRegionManager * const elemManager = mesh->getElemManager();

  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > const totalMeanStress =
    elemManager->ConstructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( viewKeyStruct::totalMeanStressString );

  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > oldTotalMeanStress =
    elemManager->ConstructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( viewKeyStruct::oldTotalMeanStressString );

  //***** loop over all elements and initialize the derivative arrays *****
  forAllElemsInMesh( mesh, [&]( localIndex const er,
                                localIndex const esr,
                                localIndex const k )->void
  {
    totalMeanStress[er][esr][k] = oldTotalMeanStress[er][esr][k];
  } );
}

real64 PoroelasticSolver::SolverStep( real64 const & time_n,
                                      real64 const & dt,
                                      int const cycleNumber,
                                      DomainPartition * const domain )
{
  real64 dt_return = dt;
  if( m_couplingTypeOption == couplingTypeOption::SIM_FixedStress )
  {
    dt_return = SplitOperatorStep( time_n, dt, cycleNumber, domain->group_cast< DomainPartition * >() );
  }
  else if( m_couplingTypeOption == couplingTypeOption::FIM )
  {
    SetupSystem( domain,
                 m_dofManager,
                 m_matrix,
                 m_rhs,
                 m_solution );

    ImplicitStepSetup( time_n,
                       dt,
                       domain,
                       m_dofManager,
                       m_matrix,
                       m_rhs,
                       m_solution );

    dt_return = this->NonlinearImplicitStep( time_n, dt, cycleNumber, domain,
                                             m_dofManager,
                                             m_matrix,
                                             m_rhs,
                                             m_solution );

    ImplicitStepComplete( time_n,
                          dt_return,
                          domain );
  }
  return dt_return;
}

void PoroelasticSolver::UpdateDeformationForCoupling( DomainPartition * const domain )
{

  MeshLevel & mesh = *domain->getMeshBody( 0 )->getMeshLevel( 0 );
  NodeManager & nodeManager = *mesh.getNodeManager();

  NumericalMethodsManager const & numericalMethodManager = domain->getNumericalMethodManager();

  FiniteElementDiscretizationManager const &
  feDiscretizationManager = numericalMethodManager.getFiniteElementDiscretizationManager();

  FiniteElementDiscretization const &
  feDiscretization = *feDiscretizationManager.GetGroup< FiniteElementDiscretization >( m_discretizationName );

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X = nodeManager.referencePosition();
  arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const & u = nodeManager.totalDisplacement();

  forTargetSubRegionsComplete< CellElementSubRegion >( mesh, [&]( localIndex const,
                                                                  localIndex const,
                                                                  localIndex const,
                                                                  ElementRegionBase & elemRegion,
                                                                  CellElementSubRegion & elementSubRegion )
  {
    string const & solidName = m_solidSolver->solidMaterialNames()[m_solidSolver->targetRegionIndex( elemRegion.getName() )];
    SolidBase const & solid = GetConstitutiveModel< SolidBase >( elementSubRegion, solidName );

    arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes = elementSubRegion.nodeList();

    arrayView1d< real64 > const &
    totalMeanStress = elementSubRegion.getReference< array1d< real64 > >( viewKeyStruct::totalMeanStressString );

    arrayView1d< real64 > const &
    oldTotalMeanStress = elementSubRegion.getReference< array1d< real64 > >( viewKeyStruct::oldTotalMeanStressString );

    arrayView1d< real64 const > const &
    pres = elementSubRegion.getReference< array1d< real64 > >( FlowSolverBase::viewKeyStruct::pressureString );

    arrayView1d< real64 const > const &
    dPres = elementSubRegion.getReference< array1d< real64 > >( FlowSolverBase::viewKeyStruct::deltaPressureString );

    arrayView1d< real64 > const &
    poro = elementSubRegion.getReference< array1d< real64 > >( SinglePhaseBase::viewKeyStruct::porosityString );

    arrayView1d< real64 const > const &
    poroOld = elementSubRegion.getReference< array1d< real64 > >( SinglePhaseBase::viewKeyStruct::porosityOldString );

    arrayView1d< real64 const > const &
    volume = elementSubRegion.getReference< array1d< real64 > >( CellBlock::viewKeyStruct::elementVolumeString );

    arrayView1d< real64 > const &
    dVol = elementSubRegion.getReference< array1d< real64 > >( SinglePhaseBase::viewKeyStruct::deltaVolumeString );

    arrayView1d< real64 const > const & bulkModulus = solid.getReference< array1d< real64 > >( "BulkModulus" );

    real64 const biotCoefficient = solid.getReference< real64 >( "BiotCoefficient" );

    arrayView3d< real64 const, solid::STRESS_USD > const & stress = solid.getStress();


    localIndex const numNodesPerElement = elemsToNodes.size( 1 );
    localIndex const numQuadraturePoints = feDiscretization.m_finiteElement->n_quadrature_points();

    forAll< parallelHostPolicy >( elementSubRegion.size(), [=] ( localIndex const ei )
    {

      R1Tensor u_local[10];

      for( localIndex i = 0; i < numNodesPerElement; ++i )
      {
        localIndex const nodeIndex = elemsToNodes( ei, i );
        u_local[ i ] = u[ nodeIndex ];
      }

      real64 effectiveMeanStress = 0.0;
      for( localIndex q=0; q<numQuadraturePoints; ++q )
      {
        effectiveMeanStress += ( stress( ei, q, 0 ) + stress( ei, q, 1 ) + stress( ei, q, 2 ) );
      }
      effectiveMeanStress /= ( 3 * numQuadraturePoints );

      totalMeanStress[ei] = effectiveMeanStress - biotCoefficient * (pres[ei] + dPres[ei]);

      poro[ei] = poroOld[ei] + (biotCoefficient - poroOld[ei]) / bulkModulus[ei]
                 * (totalMeanStress[ei] - oldTotalMeanStress[ei] + dPres[ei]);

      // update element volume
      R1Tensor Xlocal[ElementRegionManager::maxNumNodesPerElem];
      for( localIndex a = 0; a < elemsToNodes.size( 1 ); ++a )
      {
        Xlocal[a] = X[elemsToNodes[ei][a]];
        Xlocal[a] += u[elemsToNodes[ei][a]];
      }

      dVol[ei] = computationalGeometry::HexVolume( Xlocal ) - volume[ei];
    } );
  } );
}

void PoroelasticSolver::AssembleSystem( real64 const time_n,
                                        real64 const dt,
                                        DomainPartition * const domain,
                                        DofManager const & dofManager,
                                        ParallelMatrix & matrix,
                                        ParallelVector & rhs )
{

  // assemble J_SS
  m_solidSolver->AssembleSystem( time_n, dt, domain,
                                 dofManager,
                                 matrix,
                                 rhs );

  // assemble J_FF
  m_flowSolver->AssembleSystem( time_n, dt, domain,
                                dofManager,
                                matrix,
                                rhs );

  // assemble J_SF
  AssembleCouplingTerms( domain,
                         dofManager,
                         &matrix,
                         &rhs );

}

void PoroelasticSolver::AssembleCouplingTerms( DomainPartition * const domain,
                                               DofManager const & dofManager,
                                               ParallelMatrix * const matrix,
                                               ParallelVector * const rhs )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const & mesh = *domain->getMeshBodies()->GetGroup< MeshBody >( 0 )->getMeshLevel( 0 );
  NodeManager const * const nodeManager = mesh.getNodeManager();

  NumericalMethodsManager const & numericalMethodManager = domain->getNumericalMethodManager();

  FiniteElementDiscretizationManager const &
  feDiscretizationManager = numericalMethodManager.getFiniteElementDiscretizationManager();

  string const uDofKey = dofManager.getKey( keys::TotalDisplacement );
  arrayView1d< globalIndex const > const & uDofNumber = nodeManager->getReference< globalIndex_array >( uDofKey );

  arrayView2d< real64 const, nodes::INCR_DISPLACEMENT_USD > const & incr_disp = nodeManager->incrementalDisplacement();

  string const pDofKey = dofManager.getKey( FlowSolverBase::viewKeyStruct::pressureString );

  FiniteElementDiscretization const *
    feDiscretization = feDiscretizationManager.GetGroup< FiniteElementDiscretization >( m_discretizationName );

  matrix->open();
  rhs->open();

  // begin subregion loop
  forTargetSubRegionsComplete< CellElementSubRegion >( mesh, [&]( localIndex const,
                                                                  localIndex const,
                                                                  localIndex const,
                                                                  ElementRegionBase const & region,
                                                                  CellElementSubRegion const & elementSubRegion )
  {
    string const & fluidName = m_flowSolver->fluidModelNames()[m_flowSolver->targetRegionIndex( region.getName() )];
    SingleFluidBase const & fluid = GetConstitutiveModel< SingleFluidBase >( elementSubRegion, fluidName );

    string const & solidName = m_solidSolver->solidMaterialNames()[m_solidSolver->targetRegionIndex( region.getName() )];
    SolidBase const & solid = GetConstitutiveModel< SolidBase >( elementSubRegion, solidName );

    arrayView4d< real64 const > const & dNdX = elementSubRegion.dNdX();

    arrayView1d< integer const > const & ghostRank = elementSubRegion.ghostRank();

    arrayView2d< real64 const > const & detJ = elementSubRegion.detJ();

    arrayView1d< globalIndex const > const & pDofNumber = elementSubRegion.getReference< globalIndex_array >( pDofKey );

    arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes = elementSubRegion.nodeList();
    localIndex const numNodesPerElement = elemsToNodes.size( 1 );

    std::unique_ptr< FiniteElementBase >
    fe = feDiscretization->getFiniteElement( elementSubRegion.GetElementTypeString() );

    real64 const biotCoefficient = solid.getReference< real64 >( "BiotCoefficient" );

    arrayView2d< real64 const > const & density = fluid.density();

    int dim = 3;
    static constexpr int maxNumUDof = 24; // TODO: assuming linear HEX at most for the moment
    static constexpr int maxNumPDof = 1; // TODO: assuming piecewise constant (P0) only for the moment
    int nUDof = dim * numNodesPerElement;
    int nPDof = m_flowSolver->numDofPerCell();
    int numQuadraturePoints = fe->n_quadrature_points();

    forAll< serialPolicy >( elementSubRegion.size(), [=]( localIndex k )
    {
      stackArray1d< globalIndex, maxNumUDof > elementULocalDofIndex( nUDof );
      globalIndex elementPLocalDOfIndex;
      stackArray2d< real64, maxNumUDof * maxNumPDof > dRsdP( nUDof, nPDof );
      stackArray2d< real64, maxNumUDof * maxNumPDof > dRfdU( nPDof, nUDof );
      real64 Rf;

      dRsdP = 0.0;
      dRfdU = 0.0;
      Rf = 0.0;

      if( ghostRank[k] < 0 )
      {
        // Get dof local to global mapping
        for( localIndex a = 0; a < numNodesPerElement; ++a )
        {
          localIndex localNodeIndex = elemsToNodes[k][a];
          for( int i = 0; i < dim; ++i )
          {
            elementULocalDofIndex[static_cast< int >(a) * dim + i] = uDofNumber[localNodeIndex] + i;
          }
        }
        elementPLocalDOfIndex = pDofNumber[k];

        for( integer q = 0; q < numQuadraturePoints; ++q )
        {
          const realT detJq = detJ[k][q];

          for( integer a = 0; a < numNodesPerElement; ++a )
          {
            R1Tensor const dNdXa = dNdX[k][q][a];

            dRsdP( a * dim + 0, 0 ) += biotCoefficient * dNdXa[0] * detJq;
            dRsdP( a * dim + 1, 0 ) += biotCoefficient * dNdXa[1] * detJq;
            dRsdP( a * dim + 2, 0 ) += biotCoefficient * dNdXa[2] * detJq;
            dRfdU( 0, a * dim + 0 ) += density[k][0] * biotCoefficient * dNdXa[0] * detJq;
            dRfdU( 0, a * dim + 1 ) += density[k][0] * biotCoefficient * dNdXa[1] * detJq;
            dRfdU( 0, a * dim + 2 ) += density[k][0] * biotCoefficient * dNdXa[2] * detJq;

            localIndex localNodeIndex = elemsToNodes[k][a];

            real64 Rf_tmp = dNdXa[0] * incr_disp[localNodeIndex][0]
                            + dNdXa[1] * incr_disp[localNodeIndex][1]
                            + dNdXa[2] * incr_disp[localNodeIndex][2];
            Rf_tmp *= density[k][0] * biotCoefficient * detJq;
            Rf += Rf_tmp;
          }
        }

        matrix->add( elementULocalDofIndex.data(), &elementPLocalDOfIndex, dRsdP.data(), nUDof, nPDof );
        matrix->add( &elementPLocalDOfIndex, elementULocalDofIndex.data(), dRfdU.data(), nPDof, nUDof );
        rhs->add( &elementPLocalDOfIndex, &Rf, 1 );
      }
    } );
  } );

  matrix->close();
  rhs->close();
}

void PoroelasticSolver::ApplyBoundaryConditions( real64 const time_n,
                                                 real64 const dt,
                                                 DomainPartition * const domain,
                                                 DofManager const & dofManager,
                                                 ParallelMatrix & matrix,
                                                 ParallelVector & rhs )
{
  m_solidSolver->ApplyBoundaryConditions( time_n, dt, domain,
                                          dofManager,
                                          matrix,
                                          rhs );

  m_flowSolver->ApplyBoundaryConditions( time_n, dt, domain,
                                         dofManager,
                                         matrix,
                                         rhs );
}

real64 PoroelasticSolver::CalculateResidualNorm( DomainPartition const * const domain,
                                                 DofManager const & dofManager,
                                                 ParallelVector const & rhs )
{
  // compute norm of momentum balance residual equations
  real64 const momementumResidualNorm = m_solidSolver->CalculateResidualNorm( domain, dofManager, rhs );

  // compute norm of mass balance residual equations
  real64 const massResidualNorm = m_flowSolver->CalculateResidualNorm( domain, dofManager, rhs );

  return sqrt( momementumResidualNorm*momementumResidualNorm
               + massResidualNorm*massResidualNorm );
}

void PoroelasticSolver::CreatePreconditioner()
{
  if( m_linearSolverParameters.get().preconditionerType == "block" )
  {
    auto precond = std::make_unique< BlockPreconditioner< LAInterface > >( BlockShapeOption::UpperTriangular,
                                                                           SchurComplementOption::RowsumDiagonalProbing,
                                                                           BlockScalingOption::FrobeniusNorm );

    auto mechPrecond = LAInterface::createPreconditioner( m_solidSolver->getLinearSolverParameters() );
    precond->setupBlock( 0,
                         { { keys::TotalDisplacement, 0, 3 } },
                         std::make_unique< SeparateComponentPreconditioner< LAInterface > >( 3, std::move( mechPrecond ) ) );

    auto flowPrecond = LAInterface::createPreconditioner( m_flowSolver->getLinearSolverParameters() );
    precond->setupBlock( 1,
                         { { SinglePhaseBase::viewKeyStruct::pressureString, 0, 1 } },
                         std::move( flowPrecond ) );

    m_precond = std::move( precond );
  }
  else
  {
    m_precond = LAInterface::createPreconditioner( m_linearSolverParameters.get() );
  }
}

void PoroelasticSolver::SolveSystem( DofManager const & dofManager,
                                     ParallelMatrix & matrix,
                                     ParallelVector & rhs,
                                     ParallelVector & solution )
{
  solution.zero();

  SolverBase::SolveSystem( dofManager, matrix, rhs, solution );
}

void PoroelasticSolver::ApplySystemSolution( DofManager const & dofManager,
                                             ParallelVector const & solution,
                                             real64 const scalingFactor,
                                             DomainPartition * const domain )
{
  // update displacement field
  m_solidSolver->ApplySystemSolution( dofManager, solution, scalingFactor, domain );
  // update pressure field
  m_flowSolver->ApplySystemSolution( dofManager, solution, -scalingFactor, domain );
}

real64 PoroelasticSolver::SplitOperatorStep( real64 const & time_n,
                                             real64 const & dt,
                                             integer const cycleNumber,
                                             DomainPartition * const domain )
{
  real64 dtReturn = dt;
  real64 dtReturnTemporary;

  SolidMechanicsLagrangianFEM &
  solidSolver = *(this->getParent()->GetGroup( m_solidSolverName )->group_cast< SolidMechanicsLagrangianFEM * >());

  SinglePhaseBase &
  fluidSolver = *(this->getParent()->GetGroup( m_flowSolverName )->group_cast< SinglePhaseBase * >());

  fluidSolver.SetupSystem( domain,
                           fluidSolver.getDofManager(),
                           fluidSolver.getSystemMatrix(),
                           fluidSolver.getSystemRhs(),
                           fluidSolver.getSystemSolution() );

  solidSolver.SetupSystem( domain,
                           solidSolver.getDofManager(),
                           solidSolver.getSystemMatrix(),
                           solidSolver.getSystemRhs(),
                           solidSolver.getSystemSolution() );

  fluidSolver.ImplicitStepSetup( time_n, dt, domain,
                                 fluidSolver.getDofManager(),
                                 fluidSolver.getSystemMatrix(),
                                 fluidSolver.getSystemRhs(),
                                 fluidSolver.getSystemSolution() );

  solidSolver.ImplicitStepSetup( time_n, dt, domain,
                                 solidSolver.getDofManager(),
                                 solidSolver.getSystemMatrix(),
                                 solidSolver.getSystemRhs(),
                                 solidSolver.getSystemSolution() );

  this->ImplicitStepSetup( time_n, dt, domain, m_dofManager, m_matrix, m_rhs, m_solution );

  int iter = 0;
  while( iter < m_nonlinearSolverParameters.m_maxIterNewton )
  {
    if( iter == 0 )
    {
      // reset the states of all slave solvers if any of them has been reset
      fluidSolver.ResetStateToBeginningOfStep( domain );
      solidSolver.ResetStateToBeginningOfStep( domain );
      ResetStateToBeginningOfStep( domain );
    }

    GEOSX_LOG_LEVEL_RANK_0( 1, "\tIteration: " << iter+1  << ", FlowSolver: " );

    dtReturnTemporary = fluidSolver.NonlinearImplicitStep( time_n,
                                                           dtReturn,
                                                           cycleNumber,
                                                           domain,
                                                           fluidSolver.getDofManager(),
                                                           fluidSolver.getSystemMatrix(),
                                                           fluidSolver.getSystemRhs(),
                                                           fluidSolver.getSystemSolution() );

    if( dtReturnTemporary < dtReturn )
    {
      iter = 0;
      dtReturn = dtReturnTemporary;
      continue;
    }

    if( fluidSolver.getNonlinearSolverParameters().m_numNewtonIterations == 0 && iter > 0 )
    {
      GEOSX_LOG_LEVEL_RANK_0( 1, "***** The iterative coupling has converged in " << iter  << " iterations! *****\n" );
      break;
    }

    GEOSX_LOG_LEVEL_RANK_0( 1, "\tIteration: " << iter+1  << ", MechanicsSolver: " );

    solidSolver.ResetStressToBeginningOfStep( domain );
    dtReturnTemporary = solidSolver.NonlinearImplicitStep( time_n,
                                                           dtReturn,
                                                           cycleNumber,
                                                           domain,
                                                           solidSolver.getDofManager(),
                                                           solidSolver.getSystemMatrix(),
                                                           solidSolver.getSystemRhs(),
                                                           solidSolver.getSystemSolution() );

    solidSolver.updateStress( domain );

    if( dtReturnTemporary < dtReturn )
    {
      iter = 0;
      dtReturn = dtReturnTemporary;
      continue;
    }
    if( solidSolver.getNonlinearSolverParameters().m_numNewtonIterations > 0 )
    {
      this->UpdateDeformationForCoupling( domain );
    }
    ++iter;
  }

  fluidSolver.ImplicitStepComplete( time_n, dt, domain );
  solidSolver.ImplicitStepComplete( time_n, dt, domain );
  this->ImplicitStepComplete( time_n, dt, domain );

  return dtReturn;
}


REGISTER_CATALOG_ENTRY( SolverBase, PoroelasticSolver, std::string const &, Group * const )

} /* namespace geosx */
