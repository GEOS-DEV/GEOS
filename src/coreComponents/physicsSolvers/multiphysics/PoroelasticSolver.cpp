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
  m_couplingTypeOptionString( "FixedStress" ),
  m_couplingTypeOption()

{
  registerWrapper( viewKeyStruct::solidSolverNameString, &m_solidSolverName, 0 )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Name of the solid mechanics solver to use in the poroelastic solver" );

  registerWrapper( viewKeyStruct::fluidSolverNameString, &m_flowSolverName, 0 )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Name of the fluid mechanics solver to use in the poroelastic solver" );

  registerWrapper( viewKeyStruct::couplingTypeOptionStringString, &m_couplingTypeOptionString, 0 )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Coupling option: (FixedStress, TightlyCoupled)" );

}

void PoroelasticSolver::RegisterDataOnMesh( dataRepository::Group * const MeshBodies )
{
  for( auto & mesh : MeshBodies->GetSubGroups() )
  {
    ElementRegionManager * const elemManager = mesh.second->group_cast< MeshBody * >()->getMeshLevel( 0 )->getElemManager();
    elemManager->forElementSubRegions< CellElementSubRegion,
                                       FaceElementSubRegion >( [&]( auto * const elementSubRegion ) -> void
    {
      elementSubRegion->template registerWrapper< array1d< real64 > >( viewKeyStruct::totalMeanStressString )->
        setDescription( "Total Mean Stress" );
      elementSubRegion->template registerWrapper< array1d< real64 > >( viewKeyStruct::oldTotalMeanStressString )->
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

  localIndex const numLocalDof = dofManager.numLocalDofs();
  matrix.createWithLocalSize( numLocalDof, numLocalDof, 30, MPI_COMM_GEOSX );
  rhs.createWithLocalSize( numLocalDof, MPI_COMM_GEOSX );
  solution.createWithLocalSize( numLocalDof, MPI_COMM_GEOSX );

  dofManager.setSparsityPattern( m_matrix, true );

}

void PoroelasticSolver::ImplicitStepSetup( real64 const & time_n,
                                           real64 const & dt,
                                           DomainPartition * const domain,
                                           DofManager & dofManager,
                                           ParallelMatrix & matrix,
                                           ParallelVector & rhs,
                                           ParallelVector & solution )
{
  string ctOption = this->getReference< string >( viewKeyStruct::couplingTypeOptionStringString );

  if( ctOption == "FixedStress" )
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
  else if( ctOption == "TightlyCoupled" )
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
  if( m_couplingTypeOption == couplingTypeOption::TightlyCoupled )
  {
    m_solidSolver->updateStress( domain );
    m_solidSolver->ImplicitStepComplete( time_n, dt, domain );

    m_flowSolver->ImplicitStepComplete( time_n, dt, domain );
  }
}

void PoroelasticSolver::PostProcessInput()
{
  SolverBase::PostProcessInput();

  string ctOption = this->getReference< string >( viewKeyStruct::couplingTypeOptionStringString );

  if( ctOption == "FixedStress" )
  {
    this->m_couplingTypeOption = couplingTypeOption::FixedStress;

    // For this coupled solver the minimum number of Newton Iter should be 0 for both flow and solid solver otherwise it
    // will never converge.
    SolidMechanicsLagrangianFEM &
    solidSolver = *(this->getParent()->GetGroup( m_solidSolverName )->group_cast< SolidMechanicsLagrangianFEM * >());
    integer & minNewtonIterSolid = solidSolver.getNonlinearSolverParameters().m_minIterNewton;


    SinglePhaseBase &
    fluidSolver = *(this->getParent()->GetGroup( m_flowSolverName )->group_cast< SinglePhaseBase * >());
    integer & minNewtonIterFluid = fluidSolver.getNonlinearSolverParameters().m_minIterNewton;

    minNewtonIterSolid = 0;
    minNewtonIterFluid = 0;
  }
  else if( ctOption == "TightlyCoupled" )
  {
    this->m_couplingTypeOption = couplingTypeOption::TightlyCoupled;

    m_flowSolver  = this->getParent()->GetGroup< SinglePhaseBase >( m_flowSolverName );
    m_solidSolver = this->getParent()->GetGroup< SolidMechanicsLagrangianFEM >( m_solidSolverName );

    GEOSX_ERROR_IF( m_flowSolver == nullptr, "Flow solver not found or invalid type: " << m_flowSolverName );
    GEOSX_ERROR_IF( m_solidSolver == nullptr, "Solid solver not found or invalid type: " << m_solidSolverName );
  }
  else
  {
    GEOSX_ERROR( "invalid coupling type option" );
  }

}

void PoroelasticSolver::InitializePostInitialConditions_PreSubGroups( Group * const problemManager )
{
  if( m_couplingTypeOption == couplingTypeOption::FixedStress )
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
  if( m_couplingTypeOption == couplingTypeOption::FixedStress )
  {
    dt_return = SplitOperatorStep( time_n, dt, cycleNumber, domain->group_cast< DomainPartition * >() );
  }
  else if( m_couplingTypeOption == couplingTypeOption::TightlyCoupled )
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

  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup< MeshBody >( 0 )->getMeshLevel( 0 );

  ElementRegionManager * const elemManager = mesh->getElemManager();

  NodeManager * const nodeManager = domain->getMeshBody( 0 )->getMeshLevel( 0 )->getNodeManager();

  NumericalMethodsManager const * const numericalMethodManager =
    domain->getParent()->GetGroup< NumericalMethodsManager >( keys::numericalMethodsManager );

  FiniteElementDiscretizationManager const * const feDiscretizationManager =
    numericalMethodManager->GetGroup< FiniteElementDiscretizationManager >( keys::finiteElementDiscretizations );

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X = nodeManager->referencePosition();
  arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const & u = nodeManager->totalDisplacement();


  for( localIndex er=0; er<elemManager->numRegions(); ++er )
  {
    ElementRegionBase * const elemRegion = elemManager->GetRegion( er );

    FiniteElementDiscretization const * feDiscretization = feDiscretizationManager->GetGroup< FiniteElementDiscretization >( m_discretizationName );

    for( localIndex esr=0; esr<elemRegion->numSubRegions(); ++esr )
    {
      CellElementSubRegion * const cellElementSubRegion = elemRegion->GetSubRegion< CellElementSubRegion >( esr );

      arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes = cellElementSubRegion->nodeList();

      arrayView1d< real64 > const &
      totalMeanStress = cellElementSubRegion->getReference< array1d< real64 > >( viewKeyStruct::totalMeanStressString );

      arrayView1d< real64 > const &
      oldTotalMeanStress = cellElementSubRegion->getReference< array1d< real64 > >( viewKeyStruct::oldTotalMeanStressString );

      arrayView1d< real64 const > const &
      pres = cellElementSubRegion->getReference< array1d< real64 > >( FlowSolverBase::viewKeyStruct::pressureString );

      arrayView1d< real64 const > const &
      dPres = cellElementSubRegion->getReference< array1d< real64 > >( FlowSolverBase::viewKeyStruct::deltaPressureString );

      arrayView1d< real64 > const &
      poro = cellElementSubRegion->getReference< array1d< real64 > >( SinglePhaseBase::viewKeyStruct::porosityString );

      arrayView1d< real64 const > const &
      poroOld = cellElementSubRegion->getReference< array1d< real64 > >( SinglePhaseBase::viewKeyStruct::porosityOldString );

      arrayView1d< real64 const > const &
      volume = cellElementSubRegion->getReference< array1d< real64 > >( CellBlock::viewKeyStruct::elementVolumeString );

      arrayView1d< real64 > const &
      dVol = cellElementSubRegion->getReference< array1d< real64 > >( SinglePhaseBase::viewKeyStruct::deltaVolumeString );

      string const solidModelName = this->getSolidSolver()->getSolidMaterialName();

      arrayView1d< real64 const > const &
      bulkModulus = cellElementSubRegion->GetConstitutiveModels()->GetGroup( solidModelName )->getReference< array1d< real64 > >( "BulkModulus" );

      real64 const
      biotCoefficient = cellElementSubRegion->GetConstitutiveModels()->GetGroup( solidModelName )->getReference< real64 >( "BiotCoefficient" );

      arrayView3d< real64 const, solid::STRESS_USD > const &
      stress = cellElementSubRegion->GetConstitutiveModels()->GetGroup( solidModelName )->
                 getReference< array3d< real64, solid::STRESS_PERMUTATION > >( SolidBase::viewKeyStruct::stressString );


      localIndex const numNodesPerElement = elemsToNodes.size( 1 );
      localIndex const numQuadraturePoints = feDiscretization->m_finiteElement->n_quadrature_points();


      forall_in_range< parallelHostPolicy >( 0, cellElementSubRegion->size(),
                                             [=] GEOSX_HOST_DEVICE ( localIndex const ei )
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
    }
  }

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
  AssembleCouplingBlocks( domain,
                          dofManager,
                          matrix,
                          rhs );

}

void PoroelasticSolver::AssembleCouplingBlocks( DomainPartition * const domain,
                                                DofManager const & dofManager,
                                                ParallelMatrix & matrix,
                                                ParallelVector & rhs )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup< MeshBody >( 0 )->getMeshLevel( 0 );
  NodeManager * const nodeManager = mesh->getNodeManager();
//  ConstitutiveManager  * const constitutiveManager = domain->GetGroup<ConstitutiveManager
// >(keys::ConstitutiveManager);
  ElementRegionManager * const elemManager = mesh->getElemManager();
  NumericalMethodsManager const * numericalMethodManager = domain->getParent()->GetGroup< NumericalMethodsManager >( keys::numericalMethodsManager );
  FiniteElementDiscretizationManager const * feDiscretizationManager = numericalMethodManager->GetGroup< FiniteElementDiscretizationManager >(
    keys::finiteElementDiscretizations );

  string const uDofKey = dofManager.getKey( keys::TotalDisplacement );
  arrayView1d< globalIndex const > const & uDofNumber = nodeManager->getReference< globalIndex_array >( uDofKey );

  arrayView2d< real64 const, nodes::INCR_DISPLACEMENT_USD > const & incr_disp = nodeManager->incrementalDisplacement();

  string const pDofKey = dofManager.getKey( FlowSolverBase::viewKeyStruct::pressureString );

  matrix.open();
  rhs.open();

  // begin region loop
  for( localIndex er=0; er<elemManager->numRegions(); ++er )
  {
    ElementRegionBase * const elementRegion = elemManager->GetRegion( er );

    FiniteElementDiscretization const *
      feDiscretization = feDiscretizationManager->GetGroup< FiniteElementDiscretization >( m_discretizationName );

    elementRegion->forElementSubRegionsIndex< CellElementSubRegion >( [&]( localIndex const GEOSX_UNUSED_PARAM ( esr ),
                                                                           CellElementSubRegion const * const elementSubRegion )
    {
      arrayView3d< R1Tensor const > const &
      dNdX = elementSubRegion->getReference< array3d< R1Tensor > >( keys::dNdX );

      arrayView2d< real64 const > const & detJ = elementSubRegion->getReference< array2d< real64 > >( keys::detJ );

      arrayView1d< globalIndex const > const & pDofNumber = elementSubRegion->getReference< globalIndex_array >( pDofKey );

      arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes = elementSubRegion->nodeList();
      localIndex const numNodesPerElement = elemsToNodes.size( 1 );

      std::unique_ptr< FiniteElementBase >
      fe = feDiscretization->getFiniteElement( elementSubRegion->GetElementTypeString() );

      string const solidModelName = this->getSolidSolver()->getSolidMaterialName();
      real64 const
      biotCoefficient = elementSubRegion->GetConstitutiveModels()->GetGroup( solidModelName )->getReference< real64 >( "BiotCoefficient" );

      arrayView2d< real64 const > const &
      density = elementSubRegion->GetConstitutiveModels()
                  ->GetGroup( m_flowSolver->fluidIndex())
                  ->getReference< array2d< real64 > >( SingleFluidBase::viewKeyStruct::densityString );

      int dim = 3;
      int nUDof = dim * numNodesPerElement;
      int nPDof = m_flowSolver->numDofPerCell();

      array1d< globalIndex > elementULocalDofIndex( nUDof );
      globalIndex elementPLocalDOfIndex;
      array2d< real64 >      dRsdP( nUDof, nPDof );
      array2d< real64 >      dRfdU( nPDof, nUDof );
      real64 Rf;

      for( localIndex k=0; k<elementSubRegion->size(); ++k )
      {

        dRsdP = 0.0;
        dRfdU = 0.0;
        Rf = 0.0;

        if( elementSubRegion->m_ghostRank[k] < 0 )
        {
          // Get dof local to global mapping
          for( localIndex a=0; a<numNodesPerElement; ++a )
          {
            localIndex localNodeIndex = elemsToNodes[k][a];
            for( int i=0; i<dim; ++i )
            {
              elementULocalDofIndex[static_cast< int >(a)*dim+i] = uDofNumber[localNodeIndex]+i;
            }
          }
          elementPLocalDOfIndex = pDofNumber[k];

          R1Tensor dNdXa;

          for( integer q=0; q<fe->n_quadrature_points(); ++q )
          {
            const realT detJq = detJ[k][q];

            for( integer a=0; a<numNodesPerElement; ++a )
            {
              dNdXa = dNdX[k][q][a];

              dRsdP( a*dim+0, 0 ) += biotCoefficient * dNdXa[0] * detJq;
              dRsdP( a*dim+1, 0 ) += biotCoefficient * dNdXa[1] * detJq;
              dRsdP( a*dim+2, 0 ) += biotCoefficient * dNdXa[2] * detJq;
              dRfdU( 0, a*dim+0 ) += density[k][0] * biotCoefficient * dNdXa[0] * detJq;
              dRfdU( 0, a*dim+1 ) += density[k][0] * biotCoefficient * dNdXa[1] * detJq;
              dRfdU( 0, a*dim+2 ) += density[k][0] * biotCoefficient * dNdXa[2] * detJq;

              localIndex localNodeIndex = elemsToNodes[k][a];

              real64 Rf_tmp = dNdXa[0]*incr_disp[localNodeIndex][0]
                              + dNdXa[1]*incr_disp[localNodeIndex][1]
                              + dNdXa[2]*incr_disp[localNodeIndex][2];
              Rf_tmp *= density[k][0] * biotCoefficient * detJq;
              Rf += Rf_tmp;
            }
          }

          matrix.add( elementULocalDofIndex.data(), &elementPLocalDOfIndex, dRsdP.data(), nUDof, nPDof );
          matrix.add( &elementPLocalDOfIndex, elementULocalDofIndex.data(), dRfdU.data(), nPDof, nUDof );
          rhs.add( &elementPLocalDOfIndex, &Rf, 1 );
        }
      }
    } );

  }

  matrix.close();
  rhs.close();
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

void PoroelasticSolver::SolveSystem( DofManager const & dofManager,
                                     ParallelMatrix & matrix,
                                     ParallelVector & rhs,
                                     ParallelVector & solution )
{
  solution.zero();

  SolverBase::SolveSystem( dofManager, matrix, rhs, solution );

  // Debug for logLevel >= 2
  GEOSX_LOG_LEVEL_RANK_0( 2, "After ReservoirSolver::SolveSystem" );
  GEOSX_LOG_LEVEL_RANK_0( 2, "\nSolution:\n" << solution );
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
