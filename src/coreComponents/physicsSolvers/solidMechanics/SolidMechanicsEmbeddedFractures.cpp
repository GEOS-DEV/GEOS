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
#include "mesh/NodeManager.hpp"
#include "mesh/EmbeddedSurfaceRegion.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"
#include "linearAlgebra/utilities/LAIHelperFunctions.hpp"
#include "linearAlgebra/interfaces/BlasLapackLA.hpp"


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
          //subRegion->registerWrapper< array1d<real64> >( viewKeyStruct::dispJumpString )->setPlotLevel(PlotLevel::LEVEL_0);
          // subRegion->registerWrapper< array1d<real64> >( viewKeyStruct::deltaDispJumpString );
          subRegion->registerWrapper< array1d<R1Tensor> >( viewKeyStruct::dispJumpString )->setPlotLevel(PlotLevel::LEVEL_0);
          subRegion->registerWrapper< array1d<R1Tensor> >( viewKeyStruct::deltaDispJumpString );
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
                       3,
                       regions );

  dofManager.addCoupling( viewKeyStruct::dispJumpString,
                          viewKeyStruct::dispJumpString,
                          DofManager::Connectivity::Elem,
                          regions );
}

void SolidMechanicsEmbeddedFractures::SetupSystem( DomainPartition * const domain,
                                                   DofManager & dofManager,
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
  dofManager.setMesh( domain, 0, 0 );
  SetupDofs( domain, dofManager );

  // By not calling dofManager.reorderByRank(), we keep separate dof numbering for each field,
  // which allows constructing separate sparsity patterns for off-diagonal blocks of the matrix.
  // Once the solver moves to monolithic matrix, we can remove this method and just use SolverBase::SetupSystem.

  m_matrix11.createWithLocalSize( dofManager.numLocalDofs(viewKeyStruct::dispJumpString),
                                  dofManager.numLocalDofs(viewKeyStruct::dispJumpString),
                                  1,
                                  MPI_COMM_GEOSX);

  m_matrix01.createWithLocalSize( m_solidSolver->getSystemMatrix().localRows(),
                                  m_dofManager.numLocalDofs(viewKeyStruct::dispJumpString),
                                  9,
                                  MPI_COMM_GEOSX);

  m_matrix10.createWithLocalSize( dofManager.numLocalDofs(viewKeyStruct::dispJumpString),
                                  m_solidSolver->getSystemMatrix().localRows(),
                                  24,
                                  MPI_COMM_GEOSX);

  m_residual1.createWithLocalSize(dofManager.numLocalDofs(viewKeyStruct::dispJumpString),
                                  MPI_COMM_GEOSX);



  //dofManager.setSparsityPattern( m_matrix01, keys::TotalDisplacement, keys::DispJump ); I am guessing that this won't work coz coupling has not been created.
  //dofManager.setSparsityPattern( m_matrix10, keys::DispJump, keys::TotalDisplacement );

  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  NodeManager * const nodeManager = mesh->getNodeManager();
  ElementRegionManager * const elemManager = mesh->getElemManager();

  string const jumpDofKey = dofManager.getKey( viewKeyStruct::dispJumpString );
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

  m_matrix10.close();
  m_matrix01.close();
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
   NodeManager * const nodeManager = mesh->getNodeManager();
   ConstitutiveManager  * const constitutiveManager = domain->GetGroup<ConstitutiveManager >(keys::ConstitutiveManager);
   ElementRegionManager * const elemManager = mesh->getElemManager();
   NumericalMethodsManager const * numericalMethodManager = domain->getParent()->GetGroup<NumericalMethodsManager>(keys::numericalMethodsManager);
   FiniteElementDiscretizationManager const * feDiscretizationManager = numericalMethodManager->GetGroup<FiniteElementDiscretizationManager>(keys::finiteElementDiscretizations);

   ParallelVector & residual0 = m_solidSolver->getSystemRhs();

   residual0.open();
   m_matrix11.open();
   m_matrix10.open();
   m_matrix01.open();
   m_residual1.open();

   m_matrix11.zero();
   m_matrix10.zero();
   m_matrix01.zero();
   m_residual1.zero();

   r1_array const& disp  = nodeManager->getReference<r1_array>(keys::TotalDisplacement);
   r1_array const& dDisp = nodeManager->getReference<r1_array>(keys::IncrementalDisplacement);
   array1d<R1Tensor> const & nodesCoord = nodeManager->referencePosition();

   r1_array const uhattilde;


   string const dofKey = m_dofManager.getKey( keys::TotalDisplacement );
   string const jumpDofKey = m_dofManager.getKey( viewKeyStruct::dispJumpString );

   globalIndex_array const & globalDofNumber = nodeManager->getReference<globalIndex_array>( dofKey );

   ElementRegionManager::ConstitutiveRelationAccessor<ConstitutiveBase>
   constitutiveRelations = elemManager->ConstructFullConstitutiveAccessor<ConstitutiveBase>(constitutiveManager);

   ElementRegionManager::MaterialViewAccessor< real64 > const
   density = elemManager->ConstructFullMaterialViewAccessor< real64 >( "density0",
                                                                       constitutiveManager );

   constexpr int dim = 3;
   static constexpr int nUdof = dim * 8; // this is hard-coded for now.
   // begin region loop
   elemManager->forElementRegions<EmbeddedSurfaceRegion>( [&]( EmbeddedSurfaceRegion * const embeddedRegion )->void
   {
     FiniteElementDiscretization const *
     feDiscretization = feDiscretizationManager->GetGroup<FiniteElementDiscretization>(m_solidSolver->getDiscretization());
     // loop of embeddeSubregions
     embeddedRegion->forElementSubRegions<EmbeddedSurfaceSubRegion>( [&]( EmbeddedSurfaceSubRegion * const embeddedSurfaceSubRegion )->void
     {
       localIndex const numEmbeddedElems = embeddedSurfaceSubRegion->size();
       arrayView1d< localIndex const>  const & embeddedSurfaceToRegion    = embeddedSurfaceSubRegion->getSurfaceToRegionList();
       arrayView1d< localIndex const>  const & embeddedSurfaceToSubRegion = embeddedSurfaceSubRegion->getSurfaceToSubRegionList();
       arrayView1d< localIndex const>  const & embeddedSurfaceToCell      = embeddedSurfaceSubRegion->getSurfaceToCellList();

       arrayView1d<globalIndex> const &
       embeddedElementDofNumber = embeddedSurfaceSubRegion->getReference< array1d<globalIndex> >( jumpDofKey );
       arrayView1d<R1Tensor const> const & w_global  = embeddedSurfaceSubRegion->getReference<array1d<R1Tensor> >(viewKeyStruct::dispJumpString);
       arrayView1d<R1Tensor const> const & dw_global = embeddedSurfaceSubRegion->getReference<array1d<R1Tensor> >(viewKeyStruct::deltaDispJumpString);

       arrayView1d< real64 > const & fractureSurfaceArea = embeddedSurfaceSubRegion->getElementArea();

       // loop over embedded surfaces
       for( localIndex k=0 ; k<numEmbeddedElems ; ++k )
       {
         // Get rock matrix element subregion
         CellBlock const * const elementSubRegion = Group::group_cast<CellBlock const * const>(elemManager->GetRegion(embeddedSurfaceToRegion[k])->
             GetSubRegion(embeddedSurfaceToSubRegion[k]));
         arrayView2d< localIndex const, CellBlock::NODE_MAP_UNIT_STRIDE_DIM > const & elemsToNodes = elementSubRegion->nodeList();
         // Get the number of nodes per element
         localIndex const numNodesPerElement = elemsToNodes.size(1);
         std::cout << "before" << std::endl;
         // Get finite element discretization info
         std::unique_ptr<FiniteElementBase>
         fe = feDiscretization->getFiniteElement( elementSubRegion->GetElementTypeString() );
         std::cout << "after" << std::endl;

         // Initialise local matrices and vectors
         array1d<globalIndex>             elementLocalDofIndex ( nUdof );
         array1d<globalIndex>             jumpLocalDofIndex    (   3   );

         array2d<real64>            Kwu_elem( 3, nUdof );
         array2d<real64>            Kuw_elem( nUdof, 3 );
         array2d<real64>            Kww_elem( 3, 3 );
         array1d<real64>            R1(3);
         array1d<real64>            R1wu(3);
         array1d<real64>            R0(nUdof);

         BlasLapackLA::matrixScale(0, Kwu_elem);
         BlasLapackLA::matrixScale(0, Kuw_elem);
         BlasLapackLA::matrixScale(0, Kww_elem);
         BlasLapackLA::vectorScale(0, R1);
         BlasLapackLA::vectorScale(0, R1wu);
         BlasLapackLA::vectorScale(0, R0);

         // Equilibrium and compatibility operators for the element
         // number of strain components x number of jump enrichments. The comp operator is different
         // at each Gauss point.
         array2d<real64>       eqMatrix(3, 6);
         array2d<real64>       compMatrix(6, 3);
         array2d<real64>       strainMatrix(6, numNodesPerElement * dim);

         R1Tensor u_local[8];
         R1Tensor du_local[8];

         array1d<real64>       u(numNodesPerElement * dim);
         array1d<real64>       w(3);

         // Get mechanical moduli tensor
         array2d<real64> dMatrix(6,6);
         LinearElasticIsotropic::KernelWrapper const & solidConstitutive =
             Group::group_cast<LinearElasticIsotropic const * const>(constitutiveRelations[embeddedSurfaceToRegion[k]][embeddedSurfaceToSubRegion[k]][0])->
             createKernelWrapper();
         solidConstitutive.GetStiffness( embeddedSurfaceToCell[k], dMatrix);

         // Basis functions derivatives
         array3d<R1Tensor> const &
         dNdX = elementSubRegion->getReference< array3d<R1Tensor> >(keys::dNdX);

         // transformation determinant
         arrayView2d<real64> const & detJ = elementSubRegion->getReference< array2d<real64> >(keys::detJ);

         // Fill in equilibrium operator
         array1d< real64 > const & cellVolume = elementSubRegion->getElementVolume();
         real64 hInv = fractureSurfaceArea[k] / cellVolume[embeddedSurfaceToCell[k]];  // AreaFrac / cellVolume
         AssembleEquilibriumOperator(eqMatrix, embeddedSurfaceSubRegion, k, hInv);

        // if(elemGhostRank[k] < 0)
         {
           CopyGlobalToLocal<8,R1Tensor>( elemsToNodes[embeddedSurfaceToCell[k]], disp, dDisp, u_local, du_local );

           // Dof number of jump enrichment
           for (int i= 0 ; i < dim ; i++ )
           {
             jumpLocalDofIndex[i] = embeddedElementDofNumber[k] + i;
             w(i) = w_global[k][i] + dw_global[k][i];
           }

           // Dof index of nodal displacements
           for( localIndex a=0 ; a<numNodesPerElement ; ++a)
           {

             localIndex localNodeIndex = elemsToNodes[embeddedSurfaceToCell[k]][a];

             for( int i=0 ; i<dim ; ++i )
             {
               elementLocalDofIndex[static_cast<int>(a)*dim+i] = globalDofNumber[localNodeIndex]+i;
               u(a*dim + i) = u_local[a][i] + du_local[a][i];
             }
           }

           // 1. Assembly of element matrices
           // local storage of contribution of each gauss point
           array2d<real64>            Kwu_gauss( 3, nUdof );
           array2d<real64>            Kuw_gauss( nUdof, 3 );
           array2d<real64>            Kww_gauss( 3, 3 );


           // intermediate objects to do BDC, EDB, EDC
           array2d<real64>            matBD( numNodesPerElement * dim, 6 );
           array2d<real64>            matED( 3, 6 );

           BlasLapackLA::matrixScale(0, matED);
           BlasLapackLA::matrixMatrixMultiply(eqMatrix, dMatrix, matED);

           for( integer q=0 ; q<fe->n_quadrature_points() ; ++q )
           {
             BlasLapackLA::matrixScale(0, Kwu_gauss);
             BlasLapackLA::matrixScale(0, Kuw_gauss);
             BlasLapackLA::matrixScale(0, Kww_gauss);


             const realT detJq = detJ[embeddedSurfaceToCell[k]][q];
             AssembleCompatibilityOperator(compMatrix,
                                           embeddedSurfaceSubRegion,
                                           k,
                                           q,
                                           elemsToNodes,
                                           nodesCoord,
                                           embeddedSurfaceToCell,
                                           numNodesPerElement,
                                           dNdX);

             AssembleStrainOperator(strainMatrix,
                                    embeddedSurfaceToCell[k],
                                    q,
                                    numNodesPerElement,
                                    dNdX);
            // transp(B)D
            BlasLapackLA::matrixTMatrixMultiply(strainMatrix, dMatrix, matBD);
            // EDC
            BlasLapackLA::matrixMatrixMultiply(matED,  compMatrix, Kww_gauss);
            // EDB
            BlasLapackLA::matrixMatrixMultiply(matED,  strainMatrix, Kwu_gauss);
            // transp(B)DB
            BlasLapackLA::matrixMatrixMultiply(matBD, compMatrix, Kuw_gauss);

            // multiply by determinant
            BlasLapackLA::matrixScale( detJq , Kwu_gauss);
            BlasLapackLA::matrixScale( detJq , Kuw_gauss);
            BlasLapackLA::matrixScale( detJq , Kww_gauss);

            // Add Gauss point contribution to element matrix
            BlasLapackLA::matrixMatrixAdd(Kww_gauss , Kww_elem);
            BlasLapackLA::matrixMatrixAdd(Kwu_gauss , Kwu_elem);
            BlasLapackLA::matrixMatrixAdd(Kuw_gauss , Kuw_elem);
           }

           BlasLapackLA::matrixVectorMultiply(Kww_elem, w, R1);
           BlasLapackLA::matrixVectorMultiply(Kwu_elem, u, R1wu);
           BlasLapackLA::vectorVectorAdd(R1wu, R1);

           BlasLapackLA::matrixVectorMultiply(Kuw_elem, w, R0);
         }

         // 2. Assembly into global system
         // fill in residuals
         residual0.add   (elementLocalDofIndex, R0);
         m_residual1.add ( jumpLocalDofIndex, R1 );

         // fill in matrices
         m_matrix11.add  ( jumpLocalDofIndex   , jumpLocalDofIndex,    Kww_elem);
         m_matrix10.add  ( jumpLocalDofIndex   , elementLocalDofIndex, Kwu_elem);
         m_matrix01.add  ( elementLocalDofIndex, jumpLocalDofIndex,    Kuw_elem);

       }   // loop over embedded surfaces
     }); // subregion loop
   }); // region loop

   // close all the objects
   residual0.close();
   m_residual1.close();
   m_matrix11.close();
   m_matrix10.close();
   m_matrix01.close();
}

void SolidMechanicsEmbeddedFractures::AssembleEquilibriumOperator(array2d<real64> & eqMatrix,
                                                                  EmbeddedSurfaceSubRegion * const embeddedSurfaceSubRegion,
                                                                  const localIndex k,
                                                                  const real64 hInv)
{
  GEOSX_MARK_FUNCTION;
  // Normal and tangent unit vectors
  R1Tensor const nVec  = embeddedSurfaceSubRegion->getNormalVector(k);
  R1Tensor const tVec1 = embeddedSurfaceSubRegion->getTangentVector1(k);
  R1Tensor const tVec2 = embeddedSurfaceSubRegion->getTangentVector2(k);

  BlasLapackLA::matrixScale(0, eqMatrix);

  int VoigtIndex;

  for (int i = 0; i < 3; ++i)
  {
    for (int j=0; j < 3; ++j)
    {
      if (i == j)
      {
        eqMatrix(0, i) += 0.5 * (nVec[i]  * nVec[j]);
        eqMatrix(1, i) += 0.5 * (tVec1[i] * nVec[j]);
        eqMatrix(2, i) += 0.5 * (tVec2[i] * nVec[j]);
      }else
      {
        VoigtIndex = 6 - i - j;
        eqMatrix(0, VoigtIndex) += nVec[i]  * nVec[j];
        eqMatrix(1, VoigtIndex) += tVec1[i] * nVec[j];
        eqMatrix(2, VoigtIndex) += tVec2[i] * nVec[j];
      }
    }
  }
  BlasLapackLA::matrixScale(hInv, eqMatrix);
}

void
SolidMechanicsEmbeddedFractures::
AssembleCompatibilityOperator(array2d<real64> & compMatrix,
                              EmbeddedSurfaceSubRegion * const embeddedSurfaceSubRegion,
                              localIndex const k,
                              localIndex const q,
                              arrayView2d< localIndex const, CellBlock::NODE_MAP_UNIT_STRIDE_DIM > const & elemsToNodes,
                              array1d<R1Tensor> const & nodesCoord,
                              arrayView1d< localIndex const> const & embeddedSurfaceToCell,
                              localIndex const numNodesPerElement,
                              array3d<R1Tensor> const & dNdX)
{
  GEOSX_MARK_FUNCTION;
  // Normal and tangent unit vectors
  R1Tensor const nVec  = embeddedSurfaceSubRegion->getNormalVector(k);
  R1Tensor const tVec1 = embeddedSurfaceSubRegion->getTangentVector1(k);
  R1Tensor const tVec2 = embeddedSurfaceSubRegion->getTangentVector2(k);

  // Fill in compatibility operator

  // 1. construct mvector sum(dNdX(a) * H(a)) value for each Gauss point
  R1Tensor mVec;
  R1Tensor dNdXa;
  real64 heavisideFun;
  mVec = 0.0;
  for (integer a=0 ; a<numNodesPerElement ; ++a)
  {
    // Shape function derivatives
    dNdXa = dNdX[embeddedSurfaceToCell[k]][q][a];
    // Heaviside
    heavisideFun = embeddedSurfaceSubRegion->
        ComputeHeavisideFunction(nodesCoord[ elemsToNodes[embeddedSurfaceToCell[k]][a] ], k);
    // sum contribution of each node
    mVec[0] += dNdXa[0] * heavisideFun;
    mVec[1] += dNdXa[1] * heavisideFun;
    mVec[2] += dNdXa[2] * heavisideFun;
  }

  BlasLapackLA::matrixScale(0, compMatrix);

  int VoigtIndex;

  // 2. fill in the operator itself
  for (int i = 0; i < 3; ++i)
  {
    for (int j=0; j < 3; ++j)
    {
      if (i == j)
      {
        compMatrix(i,0) -= 0.5 * (nVec[i]  * mVec[j]);
        compMatrix(i, 1) -= 0.5 * (tVec1[i] * mVec[j]);
        compMatrix(i,2) -= 0.5 * (tVec2[i] * mVec[j]);
      }else
      {
        VoigtIndex = 6 - i - j;
        compMatrix(VoigtIndex, 0) -= nVec[i]  * mVec[j];
        compMatrix(VoigtIndex, 1) -= tVec1[i] * mVec[j];
        compMatrix(VoigtIndex, 2) -= tVec2[i] * mVec[j];
      }
    }
  }
}

void SolidMechanicsEmbeddedFractures::AssembleStrainOperator(array2d<real64> & strainMatrix,
                                                             localIndex const elIndex,
                                                             localIndex const q,
                                                             localIndex const numNodesPerElement,
                                                             array3d<R1Tensor> const & dNdX)
{
  GEOSX_MARK_FUNCTION;
  BlasLapackLA::matrixScale(0, strainMatrix); // make 0

  R1Tensor dNdXa;

  for (integer a=0 ; a<numNodesPerElement ; ++a)
  {
    dNdXa = dNdX[elIndex][q][a];

    strainMatrix(0, a*3 + 0) = dNdXa[0];
    strainMatrix(1, a*3 + 1) = dNdXa[1];
    strainMatrix(2, a*3 + 2) = dNdXa[2];

    strainMatrix(3, a*3 + 1) = dNdXa[2];
    strainMatrix(3, a*3 + 2) = dNdXa[1];

    strainMatrix(4, a*3 + 0) = dNdXa[2];
    strainMatrix(4, a*3 + 2) = dNdXa[0];

    strainMatrix(5, a*3 + 0) = dNdXa[1];
    strainMatrix(5, a*3 + 1) = dNdXa[0];
  }


}

void SolidMechanicsEmbeddedFractures::ApplyBoundaryConditions( real64 const time,
                                                               real64 const dt,
                                                               DomainPartition * const domain,
                                                               DofManager const & GEOSX_UNUSED_ARG( dofManager ),
                                                               ParallelMatrix & GEOSX_UNUSED_ARG( matrix ),
                                                               ParallelVector & GEOSX_UNUSED_ARG( rhs ) )
{
  GEOSX_MARK_FUNCTION;

  std::cout << "I m here" << std::endl;
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
      GEOSX_LOG_RANK_0("matrixUU");
      GEOSX_LOG_RANK_0("***********************************************************");
      LAIHelperFunctions::PrintPermutedMatrix(m_solidSolver->getSystemMatrix(), m_permutationMatrix0, std::cout);
      MpiWrapper::Barrier();

      GEOSX_LOG_RANK_0("***********************************************************");
      GEOSX_LOG_RANK_0("matrixWW");
      GEOSX_LOG_RANK_0("***********************************************************");
      m_matrix11.print(std::cout);
      MpiWrapper::Barrier();

      GEOSX_LOG_RANK_0("***********************************************************");
      GEOSX_LOG_RANK_0("matrixUw");
      GEOSX_LOG_RANK_0("***********************************************************");
      m_matrix01.print(std::cout);
      MpiWrapper::Barrier();

      GEOSX_LOG_RANK_0("***********************************************************");
      GEOSX_LOG_RANK_0("matrixwU");
      GEOSX_LOG_RANK_0("***********************************************************");
      m_matrix10.print(std::cout);
      MpiWrapper::Barrier();

      GEOSX_LOG_RANK_0("***********************************************************");
      GEOSX_LOG_RANK_0("residual0");
      GEOSX_LOG_RANK_0("***********************************************************");
      LAIHelperFunctions::PrintPermutedVector(m_solidSolver->getSystemRhs(), m_permutationMatrix0, std::cout);
      MpiWrapper::Barrier();

      GEOSX_LOG_RANK_0("***********************************************************");
      GEOSX_LOG_RANK_0("residual1");
      GEOSX_LOG_RANK_0("***********************************************************");
      m_residual1.print(std::cout);
      MpiWrapper::Barrier();
    }
}

real64 SolidMechanicsEmbeddedFractures::CalculateResidualNorm( DomainPartition const * const domain,
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

void SolidMechanicsEmbeddedFractures::ApplySystemSolution( DofManager const & GEOSX_UNUSED_ARG( dofManager ),
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



