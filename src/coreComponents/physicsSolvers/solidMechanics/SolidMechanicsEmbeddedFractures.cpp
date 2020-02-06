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

/*
 * SolidMechanicsEmbeddedFractures.cpp
 */

#include "SolidMechanicsEmbeddedFractures.hpp"

#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/contact/ContactRelationBase.hpp"
#include "constitutive/solid/LinearElasticIsotropic.hpp"
#include "managers/DomainPartition.hpp"
#include "managers/NumericalMethodsManager.hpp"
#include "mesh/EmbeddedSurfaceRegion.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"
#include "linearAlgebra/utilities/LAIHelperFunctions.hpp"

#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

SolidMechanicsEmbeddedFractures::SolidMechanicsEmbeddedFractures( const std::string& name,
                                                                  Group * const parent ):
      SolverBase(name,parent),
      m_solidSolverName(),
      m_solidSolver(nullptr)
{
  registerWrapper(viewKeyStruct::solidSolverNameString, &m_solidSolverName, 0)->
      setInputFlag(InputFlags::REQUIRED)->
      setDescription("Name of the solid mechanics solver in the rock matrix");

  registerWrapper(viewKeyStruct::contactRelationNameString, &m_contactRelationName, 0)->
      setInputFlag(InputFlags::REQUIRED)->
      setDescription("Name of contact relation to enforce constraints on fracture boundary.");

}

SolidMechanicsEmbeddedFractures::~SolidMechanicsEmbeddedFractures()
{
  // TODO Auto-generated destructor stub
}

void SolidMechanicsEmbeddedFractures::RegisterDataOnMesh( dataRepository::Group * const  MeshBodies )
{

  for( auto & mesh : MeshBodies->GetSubGroups() )
  {
    MeshLevel * meshLevel = Group::group_cast<MeshBody *>(mesh.second)->getMeshLevel(0);

    ElementRegionManager * const elemManager = meshLevel->getElemManager();
    {
      elemManager->forElementRegions<EmbeddedSurfaceRegion>( [&] ( EmbeddedSurfaceRegion * const region )
      {
        region->forElementSubRegions<EmbeddedSurfaceSubRegion>( [&]( EmbeddedSurfaceSubRegion * const subRegion )
         {
          subRegion->registerWrapper< array1d<real64> >( viewKeyStruct::dispJumpString )->setPlotLevel(PlotLevel::LEVEL_0);
          subRegion->registerWrapper< array1d<real64> >( viewKeyStruct::deltaDispJumpString );
          // subRegion->registerWrapper< array1d<R1Tensor> >( viewKeyStruct::dispJumpString )->setPlotLevel(PlotLevel::LEVEL_0); eventually it should be a tensor
          // subRegion->registerWrapper< array1d<R1Tensor> >( viewKeyStruct::deltaDispJumpString );
         });
       });
    }
  }
}



void SolidMechanicsEmbeddedFractures::ResetStateToBeginningOfStep( DomainPartition * const domain )
{
  m_solidSolver->ResetStateToBeginningOfStep(domain);
}

void SolidMechanicsEmbeddedFractures::ImplicitStepSetup( real64 const & time_n,
                                                         real64 const & dt,
                                                         DomainPartition * const domain,
                                                         DofManager & GEOSX_UNUSED_ARG( dofManager ),
                                                         ParallelMatrix & GEOSX_UNUSED_ARG( matrix ),
                                                         ParallelVector & GEOSX_UNUSED_ARG( rhs ),
                                                         ParallelVector & GEOSX_UNUSED_ARG( solution ) )
{
  m_solidSolver = this->getParent()->GetGroup<SolidMechanicsLagrangianFEM>(m_solidSolverName);

  m_solidSolver->ImplicitStepSetup( time_n, dt, domain,
                                    m_solidSolver->getDofManager(),
                                    m_solidSolver->getSystemMatrix(),
                                    m_solidSolver->getSystemRhs(),
                                    m_solidSolver->getSystemSolution() );
}

void SolidMechanicsEmbeddedFractures::ImplicitStepComplete( real64 const& time_n,
                                                            real64 const& dt,
                                                            DomainPartition * const domain)
{
  m_solidSolver->ImplicitStepComplete( time_n, dt, domain );
}

void SolidMechanicsEmbeddedFractures::PostProcessInput()
{

}

void SolidMechanicsEmbeddedFractures::InitializePostInitialConditions_PreSubGroups(Group * const GEOSX_UNUSED_ARG( problemManager ) )
{

}

real64 SolidMechanicsEmbeddedFractures::SolverStep( real64 const & time_n,
                                                    real64 const & dt,
                                                    int const cycleNumber,
                                                    DomainPartition * const domain )
{
  real64 dtReturn = dt;

  ImplicitStepSetup( time_n,
                     dt,
                     domain,
                     m_dofManager,
                     m_matrix,
                     m_rhs,
                     m_solution );

  SetupSystem( domain,
               m_dofManager,
               m_matrix,
               m_rhs,
               m_solution  );

  // currently the only method is implicit time integration
  dtReturn = this->NonlinearImplicitStep( time_n,
                                          dt,
                                          cycleNumber,
                                          domain,
                                          m_dofManager,
                                          m_matrix,
                                          m_rhs,
                                          m_solution );

  m_solidSolver->updateStress( domain );

  // final step for completion of timestep. typically secondary variable updates and cleanup.
  ImplicitStepComplete( time_n, dtReturn, domain );

  return dtReturn;
}

void SolidMechanicsEmbeddedFractures::SetupDofs( DomainPartition const * const domain,
                                                 DofManager & dofManager ) const
{
  GEOSX_MARK_FUNCTION;
  m_solidSolver->SetupDofs( domain, dofManager );

  MeshLevel const * const meshLevel = domain->getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager const * const elemManager = meshLevel->getElemManager();

  array1d<string> regions;
  elemManager->forElementRegions<EmbeddedSurfaceRegion>( [&]( EmbeddedSurfaceRegion const * const region ) {
    regions.push_back( region->getName() );
  } );

  dofManager.addField( viewKeyStruct::dispJumpString,
                       DofManager::Location::Elem,
                       1,
                       regions );

  dofManager.addCoupling( viewKeyStruct::dispJumpString,
                          viewKeyStruct::dispJumpString,
                          DofManager::Connectivity::Elem,
                          regions );
}

void SolidMechanicsEmbeddedFractures::SetupSystem( DomainPartition * const domain,
                                                   DofManager & GEOSX_UNUSED_ARG( dofManager ),
                                                   ParallelMatrix & GEOSX_UNUSED_ARG( matrix ),
                                                   ParallelVector & GEOSX_UNUSED_ARG( rhs ),
                                                   ParallelVector & GEOSX_UNUSED_ARG( solution ) )
{
  GEOSX_MARK_FUNCTION;

  m_solidSolver->SetupSystem( domain,
                              m_solidSolver->getDofManager(),
                              m_solidSolver->getSystemMatrix(),
                              m_solidSolver->getSystemRhs(),
                              m_solidSolver->getSystemSolution() );


  // setup coupled DofManager
  m_dofManager.setMesh( domain, 0, 0 );
  SetupDofs( domain, m_dofManager );

  // By not calling dofManager.reorderByRank(), we keep separate dof numbering for each field,
  // which allows constructing separate sparsity patterns for off-diagonal blocks of the matrix.
  // Once the solver moves to monolithic matrix, we can remove this method and just use SolverBase::SetupSystem.

  m_matrix11.createWithLocalSize( m_dofManager.numLocalDofs(viewKeyStruct::dispJumpString),
                                  m_dofManager.numLocalDofs(viewKeyStruct::dispJumpString),
                                  1,
                                  MPI_COMM_GEOSX);

  m_matrix01.createWithLocalSize( m_solidSolver->getSystemMatrix().localRows(),
                                  m_dofManager.numLocalDofs(viewKeyStruct::dispJumpString),
                                  9,
                                  MPI_COMM_GEOSX);

  m_matrix10.createWithLocalSize( m_dofManager.numLocalDofs(viewKeyStruct::dispJumpString),
                                  m_solidSolver->getSystemMatrix().localRows(),
                                  24,
                                  MPI_COMM_GEOSX);

  //dofManager.setSparsityPattern( m_matrix01, keys::TotalDisplacement, keys::DispJump ); I am guessing that this won't work coz coupling has not been created.
  //dofManager.setSparsityPattern( m_matrix10, keys::DispJump, keys::TotalDisplacement );

  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  NodeManager * const nodeManager = mesh->getNodeManager();
  ElementRegionManager * const elemManager = mesh->getElemManager();

  string const jumpDofKey = m_dofManager.getKey( viewKeyStruct::dispJumpString );
  string const dispDofKey = m_solidSolver->getDofManager().getKey( keys::TotalDisplacement );

  arrayView1d<globalIndex> const &
  dispDofNumber =  nodeManager->getReference<globalIndex_array>( dispDofKey );

  elemManager->forElementSubRegions<EmbeddedSurfaceSubRegion>([&]( EmbeddedSurfaceSubRegion const * const embeddedSurfaceSubRegion )
    {
    localIndex const numEmbeddedElems = embeddedSurfaceSubRegion->size();
    arrayView1d< localIndex const>  const & embeddedSurfaceToRegion    = embeddedSurfaceSubRegion->getSurfaceToRegionList();
    arrayView1d< localIndex const>  const & embeddedSurfaceToSubRegion = embeddedSurfaceSubRegion->getSurfaceToSubRegionList();
    arrayView1d< localIndex const>  const & embeddedSurfaceToCell      = embeddedSurfaceSubRegion->getSurfaceToCellList();

    arrayView1d<globalIndex> const &
    embeddedElementDofNumber = embeddedSurfaceSubRegion->getReference< array1d<globalIndex> >( jumpDofKey );

    for( localIndex k=0 ; k<numEmbeddedElems ; ++k )
    {
      CellBlock const * const subRegion = Group::group_cast<CellBlock const * const>(elemManager->GetRegion(embeddedSurfaceToRegion[k])->
          GetSubRegion(embeddedSurfaceToSubRegion[k]));
      array1d<globalIndex> activeDisplacementDOF(3 * subRegion->numNodesPerElement());
      array1d<globalIndex> activeJumpDOF(embeddedSurfaceSubRegion->numOfJumpEnrichments());
      array1d<real64> values( 3*subRegion->numNodesPerElement() );
      values = 1;

      for (localIndex i=0 ; i<embeddedSurfaceSubRegion->numOfJumpEnrichments(); ++i)
      {
        activeJumpDOF[i] = embeddedElementDofNumber[k]+i;
      }

      for( localIndex a=0 ; a<subRegion->numNodesPerElement() ; ++a )
      {
        const localIndex & node = subRegion->nodeList(embeddedSurfaceToCell[k], a);
        for( int d=0 ; d<3 ; ++d )
        {
          activeDisplacementDOF[a * 3 + d] = dispDofNumber[node] + d;
        }
      }

      m_matrix10.insert( activeJumpDOF.data(),
                         activeDisplacementDOF.data(),
                         values.data(),
                         activeJumpDOF.size(),
                         activeDisplacementDOF.size() );

      m_matrix01.insert( activeDisplacementDOF.data(),
                         activeJumpDOF.data(),
                         values.data(),
                         activeDisplacementDOF.size(),
                         activeJumpDOF.size() );
    }
    });

}

void SolidMechanicsEmbeddedFractures::AssembleSystem( real64 const time,
                                                      real64 const dt,
                                                      DomainPartition * const domain,
                                                      DofManager const & GEOSX_UNUSED_ARG( dofManager ),
                                                      ParallelMatrix & GEOSX_UNUSED_ARG( matrix ),
                                                      ParallelVector & GEOSX_UNUSED_ARG( rhs ) )
{
  GEOSX_MARK_FUNCTION;
  m_solidSolver->AssembleSystem( time,
                                 dt,
                                 domain,
                                 m_solidSolver->getDofManager(),
                                 m_solidSolver->getSystemMatrix(),
                                 m_solidSolver->getSystemRhs() );

   MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
   Group * const nodeManager = mesh->getNodeManager();
   ConstitutiveManager  * const constitutiveManager = domain->GetGroup<ConstitutiveManager >(keys::ConstitutiveManager);
   ElementRegionManager * const elemManager = mesh->getElemManager();
   NumericalMethodsManager const * numericalMethodManager = domain->getParent()->GetGroup<NumericalMethodsManager>(keys::numericalMethodsManager);
   FiniteElementDiscretizationManager const * feDiscretizationManager = numericalMethodManager->GetGroup<FiniteElementDiscretizationManager>(keys::finiteElementDiscretizations);

   ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> const fluidPres =
     elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>("pressure");

   ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> const dPres =
     elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>("deltaPressure");

   m_matrix11.zero();
   m_matrix10.zero();
   m_matrix01.zero();
   m_residual1.zero();

   m_matrix11.open();
   m_matrix10.open();
   m_matrix01.open();
   m_residual1.open();

   r1_array const& disp  = nodeManager->getReference<r1_array>(keys::TotalDisplacement);
   r1_array const& dDisp = nodeManager->getReference<r1_array>(keys::IncrementalDisplacement);

   r1_array const uhattilde;

   string const dofKey = m_dofManager.getKey( keys::TotalDisplacement );

   globalIndex_array const & globalDofNumber = nodeManager->getReference<globalIndex_array>( dofKey );

   ElementRegionManager::ConstitutiveRelationAccessor<ConstitutiveBase>
   constitutiveRelations = elemManager->ConstructFullConstitutiveAccessor<ConstitutiveBase>(constitutiveManager);

   ElementRegionManager::MaterialViewAccessor< real64 > const
   density = elemManager->ConstructFullMaterialViewAccessor< real64 >( "density0",
                                                                       constitutiveManager );

   constexpr int dim = 3;
   // begin region loop
   elemManager->forElementRegions<EmbeddedSurfaceRegion>( [&]( EmbeddedSurfaceRegion * const embeddedRegion )->void
   {
     FiniteElementDiscretization const *
     feDiscretization = feDiscretizationManager->GetGroup<FiniteElementDiscretization>(m_discretizationName);
     // loop of embeddeSubregions
     embeddedRegion->forElementSubRegions<EmbeddedSurfaceSubRegion>( [&]( EmbeddedSurfaceSubRegion * const embeddedSurfaceSubRegion )->void
     {
       localIndex const numEmbeddedElems = embeddedSurfaceSubRegion->size();
       arrayView1d< localIndex const>  const & embeddedSurfaceToRegion    = embeddedSurfaceSubRegion->getSurfaceToRegionList();
       arrayView1d< localIndex const>  const & embeddedSurfaceToSubRegion = embeddedSurfaceSubRegion->getSurfaceToSubRegionList();
       arrayView1d< localIndex const>  const & embeddedSurfaceToCell      = embeddedSurfaceSubRegion->getSurfaceToCellList();

       // arrayView1d<globalIndex> const &
       // embeddedElementDofNumber = embeddedSurfaceSubRegion->getReference< array1d<globalIndex> >( jumpDofKey );

       for( localIndex k=0 ; k<numEmbeddedElems ; ++k )
       {
         CellBlock const * const elementSubRegion = Group::group_cast<CellBlock const * const>(elemManager->GetRegion(embeddedSurfaceToRegion[k])->
             GetSubRegion(embeddedSurfaceToSubRegion[k]));

         // Basis functions derivatives
         array3d<R1Tensor> const &
         dNdX = elementSubRegion->getReference< array3d<R1Tensor> >(keys::dNdX);
         // transformation determinant
         arrayView2d<real64> const & detJ = elementSubRegion->getReference< array2d<real64> >(keys::detJ);

         arrayView2d< localIndex const, CellBlock::NODE_MAP_UNIT_STRIDE_DIM > const & elemsToNodes = elementSubRegion->nodeList();
         localIndex const numNodesPerElement = elemsToNodes.size(1);

         std::unique_ptr<FiniteElementBase>
         fe = feDiscretization->getFiniteElement( elementSubRegion->GetElementTypeString() );

         // Initialise local matrices and vectors
         Epetra_LongLongSerialDenseVector elementLocalDofIndex ( dim * numNodesPerElement );
         Epetra_SerialDenseVector         R1                   ( dim * numNodesPerElement );
         Epetra_SerialDenseMatrix         dR1dU                ( 3, dim * numNodesPerElement );
         Epetra_SerialDenseMatrix         dR0dw                ( dim * numNodesPerElement, 3 );
         Epetra_SerialDenseMatrix         dR1dw                ( 3, 3 );

         R1Tensor u_local[8];
         R1Tensor du_local[8];
         R1Tensor uhattilde_local[8];

         dR1dU.Scale(0);
         dR0dw.Scale(0);
         dR1dw.Scale(0);
         R1.Scale(0);

         // Get mechanical moduli tensor
         real64 c[6][6];
         LinearElasticIsotropic::KernelWrapper const & solidConstitutive =
             Group::group_cast<LinearElasticIsotropic const * const>(constitutiveRelations[embeddedSurfaceToRegion[k]][embeddedSurfaceToSubRegion[k]][0])->
             createKernelWrapper();
         solidConstitutive.GetStiffness( embeddedSurfaceToCell[k], c );

        // if(elemGhostRank[k] < 0)
         {
           for( localIndex a=0 ; a<numNodesPerElement ; ++a)
           {

             localIndex localNodeIndex = elemsToNodes[embeddedSurfaceToCell[k]][a];

             for( int i=0 ; i<dim ; ++i )
             {
               elementLocalDofIndex[static_cast<int>(a)*dim+i] = globalDofNumber[localNodeIndex]+i;
             }
           }

           CopyGlobalToLocal<8,R1Tensor>( elemsToNodes[embeddedSurfaceToCell[k]], disp, dDisp, u_local, du_local );

           R1Tensor dNdXa;
           for( integer q=0 ; q<fe->n_quadrature_points() ; ++q )
           {
             const realT detJq = detJ[embeddedSurfaceToCell[k]][q];

             dR1dw(0, 0) -= 1 * detJq;
             dR1dw(0, 1) -= 0 * detJq;
             dR1dw(0, 2) -= 0 * detJq;

             dR1dw(1, 0) -= 0 * detJq;
             dR1dw(1, 1) -= 1 * detJq;
             dR1dw(1, 2) -= 0 * detJq;

             dR1dw(2, 0) -= 0 * detJq;
             dR1dw(2, 1) -= 0 * detJq;
             dR1dw(2, 2) -= 1 * detJq;

             for( integer a=0 ; a<numNodesPerElement ; ++a ) // nodes loop
             {
               dNdXa = dNdX[embeddedSurfaceToCell[k]][q][a];
               dR1dU(0,a*dim+0)   -= 0 * detJq;
               dR1dU(1,a*dim+1)   -= 0 * detJq;
               dR1dU(2,a*dim+2)   -= 0 * detJq;

               dR0dw(a*dim+0, 0) -= 0;
               dR0dw(a*dim+1, 1) -= 0;
               dR0dw(a*dim+2, 2) -= 0;
             }
           }
         }
       }
     }); // subregion loop
   }); // region loop

    m_matrix11.close();
    m_matrix10.close();
    m_matrix01.close();
    m_residual1.close();
}

void SolidMechanicsEmbeddedFractures::ApplyBoundaryConditions( real64 const time,
                                                               real64 const dt,
                                                               DomainPartition * const domain,
                                                               DofManager const & GEOSX_UNUSED_ARG( dofManager ),
                                                               ParallelMatrix & GEOSX_UNUSED_ARG( matrix ),
                                                               ParallelVector & GEOSX_UNUSED_ARG( rhs ) )
{
  GEOSX_MARK_FUNCTION;
  m_solidSolver->ApplyBoundaryConditions( time,
                                          dt,
                                          domain,
                                          m_solidSolver->getDofManager(),
                                          m_solidSolver->getSystemMatrix(),
                                          m_solidSolver->getSystemRhs() );


  if( getLogLevel() == 2 )
    {
      // Before outputting anything generate permuation matrix and permute.
      MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
      NodeManager * const nodeManager = mesh->getNodeManager();

      LAIHelperFunctions::CreatePermutationMatrix(nodeManager,
                                                  m_solidSolver->getSystemMatrix().localRows(),
                                                  m_solidSolver->getSystemMatrix().localCols(),
                                                  3,
                                                  m_solidSolver->getDofManager().getKey( keys::TotalDisplacement ),
                                                  m_permutationMatrix0);

      GEOSX_LOG_RANK_0("***********************************************************");
      GEOSX_LOG_RANK_0("matrix00");
      GEOSX_LOG_RANK_0("***********************************************************");
      LAIHelperFunctions::PrintPermutedMatrix(m_solidSolver->getSystemMatrix(), m_permutationMatrix0, std::cout);
      MpiWrapper::Barrier();

      GEOSX_LOG_RANK_0("***********************************************************");
      GEOSX_LOG_RANK_0("residual0");
      GEOSX_LOG_RANK_0("***********************************************************");
      LAIHelperFunctions::PrintPermutedVector(m_solidSolver->getSystemRhs(), m_permutationMatrix0, std::cout);
      MpiWrapper::Barrier();
    }
}

real64
SolidMechanicsEmbeddedFractures::
CalculateResidualNorm( DomainPartition const * const domain,
                       DofManager const & GEOSX_UNUSED_ARG( dofManager ),
                       ParallelVector const & GEOSX_UNUSED_ARG( rhs ) )
{
  GEOSX_MARK_FUNCTION;

  real64 const solidResidualNorm = m_solidSolver->CalculateResidualNorm( domain,
                                                                     m_solidSolver->getDofManager(),
                                                                     m_solidSolver->getSystemRhs() );

  GEOSX_LOG_RANK_0("residual = "<< solidResidualNorm);

  return solidResidualNorm;
}



void
SolidMechanicsEmbeddedFractures::
ApplySystemSolution( DofManager const & GEOSX_UNUSED_ARG( dofManager ),
                     ParallelVector const & GEOSX_UNUSED_ARG( solution ),
                     real64 const scalingFactor,
                     DomainPartition * const domain )
{
  GEOSX_MARK_FUNCTION;
  m_solidSolver->ApplySystemSolution( m_solidSolver->getDofManager(),
                                      m_solidSolver->getSystemSolution(),
                                      scalingFactor,
                                      domain );

}

REGISTER_CATALOG_ENTRY( SolverBase, SolidMechanicsEmbeddedFractures, std::string const &, Group * const )
} /* namespace geosx */



