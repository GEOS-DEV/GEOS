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

SolidMechanicsEmbeddedFractures::SolidMechanicsEmbeddedFractures( const std::string & name,
                                                                  Group * const parent ):
  SolverBase( name, parent ),
  m_solidSolverName(),
  m_solidSolver( nullptr )
{
  registerWrapper( viewKeyStruct::solidSolverNameString, &m_solidSolverName )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Name of the solid mechanics solver in the rock matrix" );

  registerWrapper( viewKeyStruct::contactRelationNameString, &m_contactRelationName )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Name of contact relation to enforce constraints on fracture boundary." );

}

SolidMechanicsEmbeddedFractures::~SolidMechanicsEmbeddedFractures()
{
  // TODO Auto-generated destructor stub
}

void SolidMechanicsEmbeddedFractures::RegisterDataOnMesh( dataRepository::Group * const MeshBodies )
{

  for( auto & mesh : MeshBodies->GetSubGroups() )
  {
    MeshLevel * meshLevel = Group::group_cast< MeshBody * >( mesh.second )->getMeshLevel( 0 );

    ElementRegionManager * const elemManager = meshLevel->getElemManager();
    {
      elemManager->forElementRegions< EmbeddedSurfaceRegion >( [&] ( EmbeddedSurfaceRegion & region )
      {
        region.forElementSubRegions< EmbeddedSurfaceSubRegion >( [&]( EmbeddedSurfaceSubRegion & subRegion )
        {
          //subRegion->registerWrapper< array1d<real64> >( viewKeyStruct::dispJumpString
          // )->setPlotLevel(PlotLevel::LEVEL_0);
          // subRegion->registerWrapper< array1d<real64> >( viewKeyStruct::deltaDispJumpString );
          subRegion.registerWrapper< array1d< R1Tensor > >( viewKeyStruct::dispJumpString )->setPlotLevel( PlotLevel::LEVEL_0 );
          subRegion.registerWrapper< array1d< R1Tensor > >( viewKeyStruct::deltaDispJumpString );
        } );
      } );
    }
  }
}


void SolidMechanicsEmbeddedFractures::ResetStateToBeginningOfStep( DomainPartition * const domain )
{
  m_solidSolver->ResetStateToBeginningOfStep( domain );
}

void SolidMechanicsEmbeddedFractures::ImplicitStepSetup( real64 const & time_n,
                                                         real64 const & dt,
                                                         DomainPartition * const domain,
                                                         DofManager & GEOSX_UNUSED_PARAM( dofManager ),
                                                         ParallelMatrix & GEOSX_UNUSED_PARAM( matrix ),
                                                         ParallelVector & GEOSX_UNUSED_PARAM( rhs ),
                                                         ParallelVector & GEOSX_UNUSED_PARAM( solution ) )
{
  m_solidSolver = this->getParent()->GetGroup< SolidMechanicsLagrangianFEM >( m_solidSolverName );

  m_solidSolver->ImplicitStepSetup( time_n, dt, domain,
                                    m_solidSolver->getDofManager(),
                                    m_solidSolver->getSystemMatrix(),
                                    m_solidSolver->getSystemRhs(),
                                    m_solidSolver->getSystemSolution() );
}

void SolidMechanicsEmbeddedFractures::ImplicitStepComplete( real64 const & time_n,
                                                            real64 const & dt,
                                                            DomainPartition * const domain )
{
  m_solidSolver->ImplicitStepComplete( time_n, dt, domain );
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
               m_solution );

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

  array1d< string > regions;
  elemManager->forElementRegions< EmbeddedSurfaceRegion >( [&]( EmbeddedSurfaceRegion const & region ) {
    regions.push_back( region.getName() );
  } );

  dofManager.addField( viewKeyStruct::dispJumpString,
                       DofManager::Location::Elem,
                       3,
                       regions );

  dofManager.addCoupling( viewKeyStruct::dispJumpString,
                          viewKeyStruct::dispJumpString,
                          DofManager::Connector::Elem,
                          regions );
}

void SolidMechanicsEmbeddedFractures::SetupSystem( DomainPartition * const domain,
                                                   DofManager & dofManager,
                                                   ParallelMatrix & matrix,
                                                   ParallelVector & rhs,
                                                   ParallelVector & solution )
{
  GEOSX_MARK_FUNCTION;

  // By not calling dofManager.reorderByRank(), we keep separate dof numbering for each field,
  // which allows constructing separate sparsity patterns for off-diagonal blocks of the matrix.
  // Once the solver moves to monolithic matrix, we can remove this method and just use SolverBase::SetupSystem.
  // setup coupled DofManager
  dofManager.setMesh( domain, 0, 0 );
  SetupDofs( domain, dofManager );
  dofManager.reorderByRank();

  localIndex const numLocalDof = dofManager.numLocalDofs();

  // Monolithic system
  matrix.createWithLocalSize( numLocalDof, numLocalDof, 27, MPI_COMM_GEOSX );
  rhs.createWithLocalSize( numLocalDof, MPI_COMM_GEOSX );
  solution.createWithLocalSize( numLocalDof, MPI_COMM_GEOSX );

  dofManager.setSparsityPattern( matrix, false ); // don't close the matrix
  //dofManager.setSparsityPattern( matrix, keys::TotalDisplacement, keys::DispJump ); I am guessing that this won't work
  // coz coupling has not been created.
  //dofManager.setSparsityPattern( matrix, keys::DispJump, keys::TotalDisplacement );

  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup< MeshBody >( 0 )->getMeshLevel( 0 );
  NodeManager * const nodeManager = mesh->getNodeManager();
  ElementRegionManager * const elemManager = mesh->getElemManager();

  string const jumpDofKey = dofManager.getKey( viewKeyStruct::dispJumpString );
  string const dispDofKey = dofManager.getKey( keys::TotalDisplacement );

  arrayView1d< globalIndex const > const &
  dispDofNumber =  nodeManager->getReference< globalIndex_array >( dispDofKey );

  elemManager->forElementSubRegions< EmbeddedSurfaceSubRegion >( [&]( EmbeddedSurfaceSubRegion & embeddedSurfaceSubRegion )
  {
    localIndex const numEmbeddedElems = embeddedSurfaceSubRegion.size();
    arrayView1d< localIndex const >  const & embeddedSurfaceToRegion    = embeddedSurfaceSubRegion.getSurfaceToRegionList();
    arrayView1d< localIndex const >  const & embeddedSurfaceToSubRegion = embeddedSurfaceSubRegion.getSurfaceToSubRegionList();
    arrayView1d< localIndex const >  const & embeddedSurfaceToCell      = embeddedSurfaceSubRegion.getSurfaceToCellList();

    arrayView1d< globalIndex const > const &
    embeddedElementDofNumber = embeddedSurfaceSubRegion.getReference< array1d< globalIndex > >( jumpDofKey );


    for( localIndex k=0; k<numEmbeddedElems; ++k )
    {
      CellBlock const * const subRegion = Group::group_cast< CellBlock const * const >( elemManager->GetRegion( embeddedSurfaceToRegion[k] )->
                                                                                          GetSubRegion( embeddedSurfaceToSubRegion[k] ));
      array1d< globalIndex > activeDisplacementDOF( 3 * subRegion->numNodesPerElement());
      array1d< globalIndex > activeJumpDOF( embeddedSurfaceSubRegion.numOfJumpEnrichments());
      array1d< real64 > values( 3*subRegion->numNodesPerElement() );
      values = 1;

      for( localIndex i=0; i<embeddedSurfaceSubRegion.numOfJumpEnrichments(); ++i )
      {
        activeJumpDOF[i] = embeddedElementDofNumber[k]+i;
      }

      for( localIndex a=0; a<subRegion->numNodesPerElement(); ++a )
      {
        const localIndex & node = subRegion->nodeList( embeddedSurfaceToCell[k], a );
        for( int d=0; d<3; ++d )
        {
          activeDisplacementDOF[a * 3 + d] = dispDofNumber[node] + d;
        }
      }

      matrix.insert( activeJumpDOF.data(),
                     activeDisplacementDOF.data(),
                     values.data(),
                     activeJumpDOF.size(),
                     activeDisplacementDOF.size() );

      matrix.insert( activeDisplacementDOF.data(),
                     activeJumpDOF.data(),
                     values.data(),
                     activeDisplacementDOF.size(),
                     activeJumpDOF.size() );
    }
  } );

  matrix.close();
}

void SolidMechanicsEmbeddedFractures::AssembleSystem( real64 const time,
                                                      real64 const dt,
                                                      DomainPartition * const domain,
                                                      DofManager const & dofManager,
                                                      ParallelMatrix & matrix,
                                                      ParallelVector & rhs )
{
  GEOSX_MARK_FUNCTION;

  matrix.zero();
  rhs.zero();


  m_solidSolver->AssembleSystem( time,
                                 dt,
                                 domain,
                                 dofManager,
                                 matrix,
                                 rhs );

  matrix.open();
  rhs.open();

  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup< MeshBody >( 0 )->getMeshLevel( 0 );
  NodeManager * const nodeManager = mesh->getNodeManager();
  ConstitutiveManager * const constitutiveManager = domain->GetGroup< ConstitutiveManager >( keys::ConstitutiveManager );
  ElementRegionManager * const elemManager = mesh->getElemManager();
  NumericalMethodsManager const * numericalMethodManager = domain->getParent()->GetGroup< NumericalMethodsManager >( keys::numericalMethodsManager );
  FiniteElementDiscretizationManager const * feDiscretizationManager = numericalMethodManager->GetGroup< FiniteElementDiscretizationManager >(
    keys::finiteElementDiscretizations );

  arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const & disp  = nodeManager->totalDisplacement();
  arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const & dDisp = nodeManager->incrementalDisplacement();
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodesCoord = nodeManager->referencePosition();

  r1_array const uhattilde;

  string const dofKey     = dofManager.getKey( keys::TotalDisplacement );
  string const jumpDofKey = dofManager.getKey( viewKeyStruct::dispJumpString );

  globalIndex_array const & globalDofNumber = nodeManager->getReference< globalIndex_array >( dofKey );

  ElementRegionManager::ConstitutiveRelationAccessor< ConstitutiveBase >
  constitutiveRelations = elemManager->ConstructFullConstitutiveAccessor< ConstitutiveBase >( constitutiveManager );

  ElementRegionManager::MaterialViewAccessor< real64 > const
  density = elemManager->ConstructFullMaterialViewAccessor< real64 >( "density0",
                                                                      constitutiveManager );

  constexpr int dim = 3;
  static constexpr int maxNumUdof = dim * 8; // this is hard-coded for now.

  // Initialise local matrices and vectors
  array1d< globalIndex >             elementLocalDofIndex ( maxNumUdof );
  array1d< globalIndex >             jumpLocalDofIndex    ( 3 );

  array2d< real64 >            Kwu_elem( 3, maxNumUdof );
  array2d< real64 >            Kuw_elem( maxNumUdof, 3 );
  array2d< real64 >            Kww_elem( 3, 3 );
  array1d< real64 >            R1( 3 );
  array1d< real64 >            R0( maxNumUdof );
  array1d< real64 >            tractionVec( 3 );
  array2d< real64 >            dTdw( 3, 3 );

  // Equilibrium and compatibility operators for the element
  // number of strain components x number of jump enrichments. The comp operator is different
  // at each Gauss point.
  array2d< real64 >       eqMatrix( 3, 6 );
  array2d< real64 >       compMatrix( 6, 3 );
  array2d< real64 >       strainMatrix( 6, maxNumUdof );

  // local storage of contribution of each gauss point
  array2d< real64 >            Kwu_gauss( 3, maxNumUdof );
  array2d< real64 >            Kuw_gauss( maxNumUdof, 3 );
  array2d< real64 >            Kww_gauss( 3, 3 );

  // intermediate objects to do BDC, EDB, EDC
  array2d< real64 >            matBD( maxNumUdof * dim, 6 );
  array2d< real64 >            matED( 3, 6 );

  array1d< R1Tensor > u_local( 8 );
  array1d< R1Tensor > du_local( 8 );

  array1d< real64 >       u( maxNumUdof );
  array1d< real64 >       w( 3 );

  array2d< real64 > dMatrix( 6, 6 );

  // begin region loop
  elemManager->forElementRegions< EmbeddedSurfaceRegion >( [&]( EmbeddedSurfaceRegion & embeddedRegion )->void
  {
    FiniteElementDiscretization const *
    feDiscretization = feDiscretizationManager->GetGroup< FiniteElementDiscretization >( m_solidSolver->getDiscretization());
    // loop of embeddeSubregions
    embeddedRegion.forElementSubRegions< EmbeddedSurfaceSubRegion >( [&]( EmbeddedSurfaceSubRegion & embeddedSurfaceSubRegion )->void
    {
      localIndex const numEmbeddedElems = embeddedSurfaceSubRegion.size();
      arrayView1d< localIndex const >  const & embeddedSurfaceToRegion    = embeddedSurfaceSubRegion.getSurfaceToRegionList();
      arrayView1d< localIndex const >  const & embeddedSurfaceToSubRegion = embeddedSurfaceSubRegion.getSurfaceToSubRegionList();
      arrayView1d< localIndex const >  const & embeddedSurfaceToCell      = embeddedSurfaceSubRegion.getSurfaceToCellList();

      arrayView1d< globalIndex > const &
      embeddedElementDofNumber = embeddedSurfaceSubRegion.getReference< array1d< globalIndex > >( jumpDofKey );
      arrayView1d< R1Tensor const > const & w_global  = embeddedSurfaceSubRegion.getReference< array1d< R1Tensor > >( viewKeyStruct::dispJumpString );
      arrayView1d< R1Tensor const > const & dw_global = embeddedSurfaceSubRegion.getReference< array1d< R1Tensor > >( viewKeyStruct::deltaDispJumpString );

      arrayView1d< real64 > const & fractureSurfaceArea = embeddedSurfaceSubRegion.getElementArea();

      // loop over embedded surfaces
      for( localIndex k=0; k<numEmbeddedElems; ++k )
      {
        // Get rock matrix element subregion
        CellBlock const * const elementSubRegion = Group::group_cast< CellBlock const * const >( elemManager->GetRegion( embeddedSurfaceToRegion[k] )->
                                                                                                   GetSubRegion( embeddedSurfaceToSubRegion[k] ));
        CellBlock::NodeMapType const & elemsToNodes = elementSubRegion->nodeList();
        // Get the number of nodes per element
        localIndex const numNodesPerElement = elemsToNodes.size( 1 );
        // Get finite element discretization info
        std::unique_ptr< FiniteElementBase >
        fe = feDiscretization->getFiniteElement( elementSubRegion->GetElementTypeString() );

        // Resize based on numbe of dof of the subregion
        int nUdof = numNodesPerElement * 3;
        elementLocalDofIndex.resize( nUdof );
        Kwu_elem.resizeDimension< 1 >( nUdof );
        Kuw_elem.resizeDimension< 0 >( nUdof );
        R0.resize( nUdof );
        u.resize( nUdof );

        // Initialize
        Kwu_elem = 0.0;
        Kuw_elem = 0.0;
        Kww_elem = 0.0;
        R0 = 0.0;
        R1 = 0.0;
        dTdw = 0.0;
        eqMatrix = 0.0;
        compMatrix = 0.0;
        strainMatrix = 0.0;

        u_local.resize( numNodesPerElement );
        du_local.resize( numNodesPerElement );

        // Get mechanical moduli tensor
        LinearElasticIsotropic::KernelWrapper const & solidConstitutive =
          Group::group_cast< LinearElasticIsotropic * const >( constitutiveRelations[embeddedSurfaceToRegion[k]][embeddedSurfaceToSubRegion[k]][0] )->
            createKernelWrapper();
        solidConstitutive.GetStiffness( embeddedSurfaceToCell[k], dMatrix );

        // Basis functions derivatives
        arrayView3d< R1Tensor const > const &
        dNdX = elementSubRegion->getReference< array3d< R1Tensor > >( keys::dNdX );

        // transformation determinant
        arrayView2d< real64 const > const & detJ = elementSubRegion->getReference< array2d< real64 > >( keys::detJ );

        // Fill in equilibrium operator
        arrayView1d< real64 const > const & cellVolume = elementSubRegion->getElementVolume();
        real64 hInv = fractureSurfaceArea[k] / cellVolume[embeddedSurfaceToCell[k]];  // AreaFrac / cellVolume
        AssembleEquilibriumOperator( eqMatrix, embeddedSurfaceSubRegion, k, hInv );

        //if(elemGhostRank[k] < 0)

        // Dof index of nodal displacements
        for( localIndex a=0; a<numNodesPerElement; ++a )
        {
          localIndex localNodeIndex = elemsToNodes[embeddedSurfaceToCell[k]][a];

          for( int i=0; i<dim; ++i )
          {
            elementLocalDofIndex[static_cast< int >(a)*dim+i] = globalDofNumber[localNodeIndex]+i;
          }
        }

        for( localIndex i = 0; i < numNodesPerElement; ++i )
        {
          localIndex const nodeID = elemsToNodes( embeddedSurfaceToCell[k], i );
          u_local[ i ] = disp[ nodeID ];
          du_local[ i ] = dDisp[ nodeID ];
        }

        // Dof number of jump enrichment
        for( int i= 0; i < dim; i++ )
        {
          jumpLocalDofIndex[i] = embeddedElementDofNumber[k] + i;
          w( i ) = w_global[k][i] + 0*dw_global[k][i];
        }

        // copy values in the R1Tensor object to use in BlasLapack interface
        for( localIndex j=0; j<numNodesPerElement; ++j )
        {
          for( int i=0; i<dim; ++i )
          {
            u( j*dim + i ) = u_local[j][i] + 0*du_local[j][i];
          }
        }

        // 1. Assembly of element matrices

        // Resize local storages.
        Kwu_gauss.resizeDimension< 1 >( nUdof );
        Kuw_gauss.resizeDimension< 0 >( nUdof );
        matBD.resizeDimension< 0 >( nUdof );
        strainMatrix.resizeDimension< 1 >( nUdof );

        BlasLapackLA::matrixMatrixMultiply( eqMatrix, dMatrix, matED );

        // Compute traction
        ComputeTraction( constitutiveManager, w, tractionVec, dTdw );

        for( integer q=0; q<fe->n_quadrature_points(); ++q )
        {
          const realT detJq = detJ[embeddedSurfaceToCell[k]][q];
          AssembleCompatibilityOperator( compMatrix,
                                         embeddedSurfaceSubRegion,
                                         k,
                                         q,
                                         elemsToNodes,
                                         nodesCoord,
                                         embeddedSurfaceToCell,
                                         numNodesPerElement,
                                         dNdX );

          AssembleStrainOperator( strainMatrix,
                                  embeddedSurfaceToCell[k],
                                  q,
                                  numNodesPerElement,
                                  dNdX );
          // transp(B)D
          BlasLapackLA::matrixTMatrixMultiply( strainMatrix, dMatrix, matBD );
          // EDC
          BlasLapackLA::matrixMatrixMultiply( matED, compMatrix, Kww_gauss );
          // EDB
          BlasLapackLA::matrixMatrixMultiply( matED, strainMatrix, Kwu_gauss );
          // transp(B)DB
          BlasLapackLA::matrixMatrixMultiply( matBD, compMatrix, Kuw_gauss );

          // multiply by determinant
          BlasLapackLA::matrixScale( detJq, Kwu_gauss );
          BlasLapackLA::matrixScale( detJq, Kuw_gauss );
          BlasLapackLA::matrixScale( detJq, Kww_gauss );

          // Add Gauss point contribution to element matrix
          BlasLapackLA::matrixMatrixAdd( Kww_gauss, Kww_elem, -1 );
          BlasLapackLA::matrixMatrixAdd( Kwu_gauss, Kwu_elem, -1 );
          BlasLapackLA::matrixMatrixAdd( Kuw_gauss, Kuw_elem, -1 );
        }

        BlasLapackLA::matrixMatrixAdd( dTdw, Kww_elem, -1 );

        BlasLapackLA::matrixVectorMultiply( Kww_elem, w, R1 );
        BlasLapackLA::matrixVectorMultiply( Kwu_elem, u, R1, 1, 1 );

        BlasLapackLA::matrixVectorMultiply( Kuw_elem, w, R0, 1, 1 );

        BlasLapackLA::vectorVectorAdd( tractionVec, R1, 1 );

        // 2. Assembly into global system
        // fill in residual vector
        rhs.add ( elementLocalDofIndex, R0 );
        rhs.add ( jumpLocalDofIndex, R1 );

        // fill in matrix
        matrix.add  ( jumpLocalDofIndex, jumpLocalDofIndex, Kww_elem );
        matrix.add  ( jumpLocalDofIndex, elementLocalDofIndex, Kwu_elem );
        matrix.add  ( elementLocalDofIndex, jumpLocalDofIndex, Kuw_elem );

      }   // loop over embedded surfaces
    } ); // subregion loop
  } ); // region loop

  // close all the objects
  rhs.close();
  matrix.close();
}

void SolidMechanicsEmbeddedFractures::AssembleEquilibriumOperator( array2d< real64 > & eqMatrix,
                                                                   EmbeddedSurfaceSubRegion & embeddedSurfaceSubRegion,
                                                                   const localIndex k,
                                                                   const real64 hInv )
{
  GEOSX_MARK_FUNCTION;
  // Normal and tangent unit vectors
  R1Tensor const nVec  = embeddedSurfaceSubRegion.getNormalVector( k );
  R1Tensor const tVec1 = embeddedSurfaceSubRegion.getTangentVector1( k );
  R1Tensor const tVec2 = embeddedSurfaceSubRegion.getTangentVector2( k );

  BlasLapackLA::matrixScale( 0, eqMatrix );

  R2Tensor nDn, t1Dn, t2Dn, t1DnSym, t2DnSym;

  // n dyadic n
  nDn.dyadic_aa( nVec );

  // sym(n dyadic t1) and sym (n dyadic t2)
  t1DnSym.dyadic_ab( nVec, tVec1 );
  t2DnSym.dyadic_ab( nVec, tVec2 );
  t1Dn.dyadic_ab( tVec1, nVec );
  t2Dn.dyadic_ab( tVec2, nVec );

  t1DnSym += t1Dn;
  t2DnSym += t2Dn;
  t1DnSym *= 0.5;
  t2DnSym *= 0.5;

  int VoigtIndex;

  for( int i=0; i < 3; ++i )
  {
    for( int j=0; j < 3; ++j )
    {
      if( i == j )
      {
        VoigtIndex = 1;
      }
      else
      {
        VoigtIndex = 6 - i - j;
      }
      eqMatrix( 0, VoigtIndex ) += nDn ( i, j );
      eqMatrix( 1, VoigtIndex ) += t1DnSym( i, j );
      eqMatrix( 2, VoigtIndex ) += t2DnSym( i, j );
    }
  }
  BlasLapackLA::matrixScale( -hInv, eqMatrix );
}

void
SolidMechanicsEmbeddedFractures::
  AssembleCompatibilityOperator( array2d< real64 > & compMatrix,
                                 EmbeddedSurfaceSubRegion & embeddedSurfaceSubRegion,
                                 localIndex const k,
                                 localIndex const q,
                                 CellBlock::NodeMapType const & elemsToNodes,
                                 arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodesCoord,
                                 arrayView1d< localIndex const > const & embeddedSurfaceToCell,
                                 localIndex const numNodesPerElement,
                                 arrayView3d< R1Tensor const > const & dNdX )
{
  GEOSX_MARK_FUNCTION;
  // Normal and tangent unit vectors
  R1Tensor const nVec  = embeddedSurfaceSubRegion.getNormalVector( k );
  R1Tensor const tVec1 = embeddedSurfaceSubRegion.getTangentVector1( k );
  R1Tensor const tVec2 = embeddedSurfaceSubRegion.getTangentVector2( k );

  // Fill in compatibility operator

  // 1. construct mvector sum(dNdX(a) * H(a)) value for each Gauss point
  R1Tensor mVec;
  R1Tensor dNdXa;
  real64 heavisideFun;
  mVec = 0.0;
  for( integer a=0; a<numNodesPerElement; ++a )
  {
    // Shape function derivatives
    dNdXa = dNdX[embeddedSurfaceToCell[k]][q][a];
    // Heaviside
    heavisideFun = embeddedSurfaceSubRegion.
                     ComputeHeavisideFunction( nodesCoord[ elemsToNodes[embeddedSurfaceToCell[k]][a] ], k );
    // sum contribution of each node
    mVec[0] -= dNdXa[0] * heavisideFun;
    mVec[1] -= dNdXa[1] * heavisideFun;
    mVec[2] -= dNdXa[2] * heavisideFun;
  }

  BlasLapackLA::matrixScale( 0, compMatrix );

  // 2. fill in the operator itself

  R2Tensor nDm, t1Dm, t2Dm, nDmSym, t1DmSym, t2DmSym;

  // sym(n dyadic m)
  nDmSym.dyadic_ab( mVec, nVec );
  nDm.dyadic_ab( nVec, mVec );
  nDmSym += nDm;
  nDmSym *= 0.5;

  // sym(n dyadic t1) and sym (n dyadic t2)
  t1DmSym.dyadic_ab( mVec, tVec1 );
  t2DmSym.dyadic_ab( mVec, tVec2 );
  t1Dm.dyadic_ab( tVec1, mVec );
  t2Dm.dyadic_ab( tVec2, mVec );
  t1DmSym += t1Dm;
  t2DmSym += t2Dm;
  t1DmSym *= 0.5;
  t2DmSym *= 0.5;

  int VoigtIndex;

  for( int i=0; i < 3; ++i )
  {
    for( int j=0; j < 3; ++j )
    {
      if( i == j )
      {
        VoigtIndex = 1;
      }
      else
      {
        VoigtIndex = 6 - i - j;
      }
      compMatrix( VoigtIndex, 0 ) += nDmSym ( i, j );
      compMatrix( VoigtIndex, 1 ) += t1DmSym( i, j );
      compMatrix( VoigtIndex, 2 ) += t2DmSym( i, j );
    }
  }
}

void SolidMechanicsEmbeddedFractures::AssembleStrainOperator( array2d< real64 > & strainMatrix,
                                                              localIndex const elIndex,
                                                              localIndex const q,
                                                              localIndex const numNodesPerElement,
                                                              arrayView3d< R1Tensor const > const & dNdX )
{
  GEOSX_MARK_FUNCTION;
  strainMatrix = 0.0; // make 0

  R1Tensor dNdXa;

  for( integer a=0; a<numNodesPerElement; ++a )
  {
    dNdXa = dNdX[elIndex][q][a];

    strainMatrix( 0, a*3 + 0 ) = dNdXa[0];
    strainMatrix( 1, a*3 + 1 ) = dNdXa[1];
    strainMatrix( 2, a*3 + 2 ) = dNdXa[2];

    strainMatrix( 3, a*3 + 1 ) = dNdXa[2];
    strainMatrix( 3, a*3 + 2 ) = dNdXa[1];

    strainMatrix( 4, a*3 + 0 ) = dNdXa[2];
    strainMatrix( 4, a*3 + 2 ) = dNdXa[0];

    strainMatrix( 5, a*3 + 0 ) = dNdXa[1];
    strainMatrix( 5, a*3 + 1 ) = dNdXa[0];
  }


}

void SolidMechanicsEmbeddedFractures::ApplyBoundaryConditions( real64 const time,
                                                               real64 const dt,
                                                               DomainPartition * const domain,
                                                               DofManager const & dofManager,
                                                               ParallelMatrix & matrix,
                                                               ParallelVector & rhs )
{
  GEOSX_MARK_FUNCTION;

  m_solidSolver->ApplyBoundaryConditions( time,
                                          dt,
                                          domain,
                                          dofManager,
                                          matrix,
                                          rhs );
}

real64 SolidMechanicsEmbeddedFractures::CalculateResidualNorm( DomainPartition const * const domain,
                                                               DofManager const & dofManager,
                                                               ParallelVector const & rhs )
{
  GEOSX_MARK_FUNCTION;

  // Matrix residual
  real64 const solidResidualNorm = m_solidSolver->CalculateResidualNorm( domain,
                                                                         dofManager,
                                                                         rhs );

  // Fracture residual

  GEOSX_LOG_RANK_0( "residual = "<< solidResidualNorm );

  return solidResidualNorm;
}

void SolidMechanicsEmbeddedFractures::ApplySystemSolution( DofManager const & dofManager,
                                                           ParallelVector const & solution,
                                                           real64 const scalingFactor,
                                                           DomainPartition * const domain )
{
  GEOSX_MARK_FUNCTION;

  m_solidSolver->ApplySystemSolution( dofManager,
                                      solution,
                                      scalingFactor,
                                      domain );

  dofManager.addVectorToField( solution, viewKeyStruct::dispJumpString, viewKeyStruct::deltaDispJumpString, -scalingFactor );

  dofManager.addVectorToField( solution, viewKeyStruct::dispJumpString, viewKeyStruct::dispJumpString, -scalingFactor );


// TODO
//  std::map<string, string_array > fieldNames;
//  fieldNames["node"].push_back( viewKeyStruct::dispJumpString );
//  fieldNames["node"].push_back( viewKeyStruct::deltaDispJumpString );
//
//  CommunicationTools::SynchronizeFields( fieldNames,
//                                           domain->getMeshBody( 0 )->getMeshLevel( 0 ),
//                                           domain->getReference< array1d<NeighborCommunicator> >(
// domain->viewKeys.neighbors ) );

}

void SolidMechanicsEmbeddedFractures::ComputeTraction( ConstitutiveManager const * const constitutiveManager,
                                                       array1d< real64 >  const & dispJump,
                                                       array1d< real64 > & tractionVector,
                                                       array2d< real64 > & dTdw )
{
  // Compute traction vector on the fracture element
  ContactRelationBase const * const
  contactRelation = constitutiveManager->GetGroup< ContactRelationBase >( m_contactRelationName );

  // check if fracture is open
  bool open = dispJump[0] >= 0 ? true : false;

  if( open )
  {
    tractionVector[0] = 1e5;
    tractionVector[1] = 0.0;
    tractionVector[2] = 0.0;
    dTdw( 0, 0 ) = 0.0;
    dTdw( 0, 1 ) = 0.0;
    dTdw( 0, 2 ) = 0.0;
    dTdw( 1, 0 ) = 0.0;
    dTdw( 1, 1 ) = 0.0;
    dTdw( 1, 2 ) = 0.0;
    dTdw( 2, 0 ) = 0.0;
    dTdw( 2, 1 ) = 0.0;
    dTdw( 2, 2 ) = 0.0;
  }
  else
  {
    // Contact through penalty condition.
    tractionVector[0] = contactRelation->stiffness() * dispJump[0];
    tractionVector[1] = 0;
    tractionVector[2] = 0;
  }
}

REGISTER_CATALOG_ENTRY( SolverBase, SolidMechanicsEmbeddedFractures, std::string const &, Group * const )
} /* namespace geosx */
