/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file MultiResolutionHFSolver.cpp
 *
 */


#include "MultiResolutionHFSolver.hpp"

#include "common/GEOS_RAJA_Interface.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/contact/ContactSelector.hpp"
#include "constitutive/fluid/SingleFluidBase.hpp"
#include "discretizationMethods/NumericalMethodsManager.hpp"
#include "finiteElement/Kinematics.h"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/SurfaceElementRegion.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "mesh/utilities/ComputationalGeometry.hpp"
#include "mesh/mpiCommunications/NeighborCommunicator.hpp"
#include "physicsSolvers/contact/SolidMechanicsEFEMKernelsHelper.hpp"
#include "physicsSolvers/contact/SolidMechanicsEmbeddedFractures.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "physicsSolvers/surfaceGeneration/SurfaceGenerator.hpp"
#include "linearAlgebra/utilities/LAIHelperFunctions.hpp"



namespace geosx
{

namespace
{
//this function is not part of any class yet - just tries to find the number of a node given its coordinates
// localIndex getNodeNumberFromCoordinate( MeshLevel const & mesh, real64 const & x, real64 const & y, real64 const & z)
// {
//   real64 tolerance = 1e-10;
//   NodeManager const & meshNodeManager = mesh.getNodeManager();
//   arrayView2d<real64 const> const nodalPositions = meshNodeManager.referencePosition();
//   for (localIndex a=0; a < meshNodeManager.size(); a++)
//     {
//  real64 dist_sq = (nodalPositions( a,0 ) - x)*(nodalPositions( a,0 ) - x) + (nodalPositions( a,1 ) - y)*(nodalPositions( a,1 ) - y) +
// (nodalPositions( a,2 ) - z)*(nodalPositions( a,2 ) - z);
//  real64 dist = sqrt(dist_sq);
//  if (dist < tolerance)
//    {
//      return a;
//    }
//     }
//   return -1;
// }

//this function checks if an index is already in an array
bool isMember( localIndex const & ind, array1d< localIndex > const & arr )
{
  for( localIndex iter : arr )
  {
    if( iter == ind )
    {
      return true;
    }
  }
  return false;
}

}

using namespace dataRepository;
using namespace constitutive;

MultiResolutionHFSolver::MultiResolutionHFSolver( const string & name,
                                                  Group * const parent ):
  SolverBase( name, parent ),
  m_baseSolverName(),
  m_patchSolverName(),
  m_nodeFixDamage(),
  m_nodeFixDisp(),
  m_fixedDispList(),
  m_maxNumResolves( 10 )
{
  registerWrapper( viewKeyStruct::baseSolverNameString(), &m_baseSolverName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription(
    "Name of the EFEM (SolidMechanicsEmbeddedFractures) solver to be used as the base solver in the MultiResolution scheme" );

  registerWrapper( viewKeyStruct::patchSolverNameString(), &m_patchSolverName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription(
    "Name of the Phase-Field (PhaseFieldFracture) solver to be used as patch solver in the MultiResolution scheme" );

  registerWrapper( viewKeyStruct::maxNumResolvesString(), &m_maxNumResolves ).
    setApplyDefaultValue( 10 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Value to indicate how many resolves may be executed to perform surface generation after the execution of base and patch scale solvers. " );

}

void MultiResolutionHFSolver::RegisterDataOnMesh( dataRepository::Group & MeshBodies )
{

  // forMeshTarget(MeshBodies,[&](string const &, MeshLevel & mesh, ArrayView1d<string const> const & regionNames){

  //     //register fix displacement nodes
  //     NodeManager & nodes = mesh.getNodeManager();
  //     nodes.registerWrapper< array1d<localIndex> >(SolidMechanicsLagrangianFEM::viewKeyStruct::nameChosenString() )
  //  .setApplyDefaultValue(0.0)
  //  .setPlotLevel( PlotLevel::NOPLOT )
  //  .setDescription("").;

  //     //register fix displacement values
  //     NodeManager & nodes = mesh.getNodeManager();
  //     nodes.registerWrapper< array1d<localIndex> >(SolidMechanicsLagrangianFEM::viewKeyStruct::nameChosenString() )
  //  .setApplyDefaultValue(0.0)
  //  .setPlotLevel( PlotLevel::NOPLOT )
  //  .setDescription("").;

  //     //register fix damage nodes
  //     NodeManager & nodes = mesh.getNodeManager();
  //     nodes.registerWrapper< array1d<localIndex> >(PhaseFieldDamageFEM::viewKeyStruct::nameChosenString() )
  //  .setApplyDefaultValue(0.0)
  //  .setPlotLevel( PlotLevel::NOPLOT )
  //  .setDescription("").;


  // });
  // Andre - 03/14 - keep this function empty for later
  GEOSX_UNUSED_VAR( MeshBodies );
}

MultiResolutionHFSolver::~MultiResolutionHFSolver()
{
  // TODO Auto-generated destructor stub
}

void MultiResolutionHFSolver::implicitStepSetup( real64 const & time_n,
                                                 real64 const & dt,
                                                 DomainPartition & domain )
{
  //Andre - 14/03/22 - can we just call these functions without doing the boundary condition operations first?
  m_baseSolver->implicitStepSetup( time_n, dt, domain );
  m_patchSolver->implicitStepSetup( time_n, dt, domain );
}

void MultiResolutionHFSolver::postProcessInput()
{
  SolverBase::postProcessInput();
  m_baseSolver = &this->getParent().getGroup< SolidMechanicsEmbeddedFractures >( m_baseSolverName );
  m_patchSolver = &this->getParent().getGroup< PhaseFieldFractureSolver >( m_patchSolverName );
}

void MultiResolutionHFSolver::initializePostInitialConditionsPreSubGroups()
{}

void MultiResolutionHFSolver::resetStateToBeginningOfStep( DomainPartition & domain )
{
  m_baseSolver->resetStateToBeginningOfStep( domain );
  m_patchSolver->resetStateToBeginningOfStep( domain );
}

real64 MultiResolutionHFSolver::solverStep( real64 const & time_n,
                                            real64 const & dt,
                                            int const cycleNumber,
                                            DomainPartition & domain )
{
  real64 dtReturn = dt;
  dtReturn = splitOperatorStep( time_n, dt, cycleNumber, domain );
  return dtReturn;
}

//in the initial step and whenever the subproblem mesh moves, we must update the node map
// void MultiResolutionHFSolver::updateNodeMaps( MeshLevel const & base,
//            MeshLevel const & patch )
// {
//   // get list of nodes on the boundary of the patch
//   NodeManager const & patchNodeManager = patch.getNodeManager();
//   NodeManager const & baseNodeManager = base.getNodeManager();
//   SortedArrayView< localIndex const > const patchExternalSet = patchNodeManager.externalSet();
//   arrayView2d<real64 const> const patchPosition = patchNodeManager.referencePosition();
//   for (localIndex a : patchExternalSet )
//     {
//       localIndex A = getNodeNumberFromCoordinate(base, patchPosition(a,0), patchPosition(a,1), patchPosition(a,2));
//       m_nodeMapIndices( a ) = A;
//       //since the boundary points lie on element edges, these will have dimension 2. weights can be computed by taking the ratio between
// distances to 2 closes points from base mesh.
//       m_nodeMapWeights( a,1 ) = 1; //for now, we dont need this, so we just use 1.
//     }
// }

//Andre - 03/15 - this function will loop over all subdomain nodes and check their distance to the prescribed discrete crack. If the
// distance is smaller than 1 element size (subdomain), we set the damage in this node to be fixed at 1.
void MultiResolutionHFSolver::setInitialCrackDamageBCs( DofManager const & GEOSX_UNUSED_PARAM( dofManager ),
                                                        CRSMatrixView< real64, globalIndex const > const & GEOSX_UNUSED_PARAM( localMatrix ),
                                                        MeshLevel const & GEOSX_UNUSED_PARAM( patch ),
                                                        MeshLevel const & base )
{
  // get list of nodes on the boundary of the patch
  //NodeManager const & patchNodeManager = patch.getNodeManager();
  //arrayView1d< globalIndex const > const & dofIndex = patchNodeManager.getReference< array1d< globalIndex > >( dofManager.getKey( "Damage"
  // ) );
  //need dofManager from PhaseFieldDamage, not PhaseFieldFracture
  //arrayView1d< real64 const > const nodalDamage = patchNodeManager.getReference< array1d< real64 > >( "Damage" );
  ElementRegionManager const & baseElemManager = base.getElemManager();
  baseElemManager.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion const & cellElementSubRegion )
  {
    SortedArrayView< localIndex const > const fracturedElements = cellElementSubRegion.fracturedElementsList();
    //USE 1D ARRAY REGISTERED ON THE NODE MANAGER INSTEAD OF m_nodeFixDamage
    m_nodeFixDamage.resize( cellElementSubRegion.numNodesPerElement()*fracturedElements.size() );
    localIndex count = 0;
    for( localIndex a : fracturedElements )
    {
      //get all nodes of fracturedElements(a)
      for( localIndex b=0; b < cellElementSubRegion.numNodesPerElement(); b++ )
      {
        //if meshes arent identical, this should be nodeList(elem_mapped_to_patch(a),b) or check (2)
        localIndex c = cellElementSubRegion.nodeList( a, b );
        //append c = node b of element a
        if( !isMember( c, m_nodeFixDamage )) //avoid repetition
        {
          //USE REGISTERED ARRAY HERE
          //(2) alternatively, we could use node_mapped_to_patch(cellElementSubRegion.nodeList(a,b))
          m_nodeFixDamage( count ) = cellElementSubRegion.nodeList( a, b );
          ++count;
        }
      }
    }
    //USE REGISTERED ARRAY HERE
    m_nodeFixDamage.resize( count );
  } );
}

// this function will read the patch solution and locate the crack tip to uptade the crack geometry in the base solver
void MultiResolutionHFSolver::findPhaseFieldTip( R1Tensor & tip,
                                                 MeshLevel const & patch )
{
  //IMPORTANT PARTS MISSING
  //-best way to get quadrature point damage to take average - partially done
  //-verify constructor of array1ds - partially done
  //-verify how to get elemCenter - partially done
  //-verify dist function between R1Tensors - partially done
  //-verify MPI part - partially done
  //reference point must be prescribed in base coordinate system (usually injection source)
  R1Tensor m_referencePoint = {0.0, 0.0, 0.0}; //add this as input file parameter
  real64 threshold = 0.95;
  //get mpi communicator
  MPI_Comm const & comm = MPI_COMM_GEOSX;
  ElementRegionManager const & patchElemManager = patch.getElemManager();
   //this is risky, not sure FE_TYPE will come from patch
  //each rank has these arrays
  real64 rankMaxDist = 1e20;
  R1Tensor rankFarthestCenter = {0.0, 0.0, 0.0};
  //loop over all elements in patch and compute elemental average damage
  //observe that patch coordinates are relative, so, reference point should be added 
  //for elements in patch //nice loop to parallelize
  patchElemManager.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion const & cellElementSubRegion )
  {
    localIndex numQuadraturePointsPerElem = 0;
    finiteElement::FiniteElementBase const & fe = cellElementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );
    finiteElement::FiniteElementDispatchHandler< ALL_FE_TYPES >::dispatch3D( fe, [&] ( auto const finiteElement )
    {
      using FE_TYPE = TYPEOFREF( finiteElement );
      numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;
    });
    string const & damageModelName = cellElementSubRegion.getReference< string >( PhaseFieldDamageFEM::viewKeyStruct::solidModelNamesString());
    //this call may need to constitutive pass-thru loop to be generalized to multiple damage types
    const constitutive::Damage<ElasticIsotropic> & damageUpdates = cellElementSubRegion.getConstitutiveModel< Damage<ElasticIsotropic> >( damageModelName );
    arrayView2d< const real64 > allElemCenters = cellElementSubRegion.getElementCenter(); 
    const arrayView2d< real64 const > qp_damage = damageUpdates.getDamage();
    forAll< serialPolicy >( cellElementSubRegion.size(), [&] ( localIndex const k )
    {
      //compute elemental averaged damage
      real64 average_d = 0;
      R1Tensor elemCenter;
      elemCenter[0] = allElemCenters[k][0];
      elemCenter[1] = allElemCenters[k][1];
      elemCenter[2] = allElemCenters[k][2]; //this is trying to create a R1Tensor from a array2d< real64 >
      for(localIndex q=0; q<numQuadraturePointsPerElem; q++)
      {
        //get damage at quadrature point i
        average_d = average_d + qp_damage( k,q )/numQuadraturePointsPerElem;
      }
      //if elemental damage > 0.95 (or another threshold)
      if (average_d > threshold)
      {
        //check if this element is farther than current farthest
        R1Tensor elemVec = LVARRAY_TENSOROPS_INIT_LOCAL_3( m_referencePoint );
        LvArray::tensorOps::subtract< 3 >( elemVec, elemCenter );
        real64 dist = LvArray::tensorOps::l2Norm< 3 >( elemVec );
        if (dist > rankMaxDist)
        {
          rankMaxDist = dist;
          rankFarthestCenter = elemCenter;
        }
      }
    });
  });

  real64 globalMax = MpiWrapper::max<real64>( rankMaxDist,comm );
  if (std::abs(rankMaxDist-globalMax)<1e-12)
  {
    tip = rankFarthestCenter;
  }
  
}

//Andre - 03/15 - this function needs to access u_base (hopefully via node manager), the patch mesh (via node manager too),
//a map that translate patch node numbers to base node numbers (could be a MultiResolutionHFSolver member) and also the
//damage field in the patch domain (via node manage
//The function will read all this data and prepare a list of dofs and u values that will be used by the patch solver to set
//the boundary conditions
void MultiResolutionHFSolver::prepareSubProblemBCs( MeshLevel const & base,
                                                    MeshLevel & patch )
{

  // get list of nodes on the boundary of the patch
  FaceManager const & patchFaceManager = patch.getFaceManager();
  NodeManager & patchNodeManager = patch.getNodeManager();
  patchNodeManager.setIsExternal( patchFaceManager );
  SortedArray< localIndex > patchExternalSet = patchNodeManager.externalSet();
  //arrayView1d< integer const > const patchExternalSet = patchNodeManager.isExternal(); //option 2 remove if not needed
  //arrayView1d< globalIndex const > const patchLocalToGlobalMap = patchNodeManager.localToGlobalMap();
  arrayView1d< real64 const > const patchDamage = patchNodeManager.getReference< array1d< real64 > >( "Damage" );
  NodeManager const & baseNodeManager = base.getNodeManager();
  //arrayView2d< real64 const > const baseDisp = baseNodeManager.totalDisplacement();
  arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const baseDisp = baseNodeManager.getField< fields::solidMechanics::totalDisplacement >();
  //unordered_map< globalIndex, localIndex > const & patchGlobalToLocalMap = baseNodeManager.globalToLocalMap();

  //we need to zero the lists to erase data from previous step
  //this->eraseLists(); //maybe this can go inside resetToBeginningOfStep - already there()
  real64 damage_threshold = 0.3;

  //USE REGISTERED ARRAYS INSTEAD OF m_nodeFixDisp and m_fixedDispList
  m_nodeFixDisp.resize( patchExternalSet.size() );
  m_fixedDispList.resize( patchExternalSet.size(), 3 );
  localIndex count=0;
  for( localIndex a : patchExternalSet )
  {
    if( patchDamage[a] < damage_threshold )
    {
      // NOTE: there needs to be a translation between patch and base mesh for the indices and values/weights.

      //append node a
      //USE REGISTERED ARRAY HERE
      m_nodeFixDisp( count ) = a;
      localIndex const numBaseNodes = 1;  //this wont be 1 in other cases //m_nodeMapIndices.sizeOfArray( a );
      for( localIndex b=0; b<numBaseNodes; ++b )
      {
        //write displacements from the base domain to become a boundary condition in the patch domain
        // a should be replaced by mapped(a) if meshes arent identical
        // mapped(a) is the node in base that has the same physical coordinates of a in patch
        m_fixedDispList( count, 0 ) = baseDisp( a, 0 );
        m_fixedDispList( count, 1 ) = baseDisp( a, 1 );
        m_fixedDispList( count, 2 ) = baseDisp( a, 2 );
      }

      ++count;
    }
  }
  //USE REGISTERED ARRAY HERE
  m_nodeFixDisp.resize( count );
  m_fixedDispList.resize( count, 3 );

}

real64 MultiResolutionHFSolver::splitOperatorStep( real64 const & time_n,
                                                   real64 const & dt,
                                                   integer const cycleNumber,
                                                   DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;
  real64 dtReturn = dt;
  real64 dtReturnTemporary;

  SolidMechanicsEmbeddedFractures &
  baseSolver = this->getParent().getGroup< SolidMechanicsEmbeddedFractures >( m_baseSolverName );

  PhaseFieldFractureSolver &
  patchSolver = this->getParent().getGroup< PhaseFieldFractureSolver >( m_patchSolverName );

  PhaseFieldDamageFEM &
  patchDamageSolver = patchSolver.getParent().getGroup< PhaseFieldDamageFEM >( patchSolver.getDamageSolverName() );

  SolidMechanicsLagrangianFEM &
  patchSolidSolver = patchSolver.getParent().getGroup< SolidMechanicsLagrangianFEM >( patchSolver.getSolidSolverName() );

  baseSolver.setupSystem( domain,
                          baseSolver.getDofManager(),
                          baseSolver.getLocalMatrix(),
                          baseSolver.getSystemRhs(),
                          baseSolver.getSystemSolution(),
                          true );

  patchSolver.setupSystem( domain,
                           patchSolver.getDofManager(),
                           patchSolver.getLocalMatrix(),
                           patchSolver.getSystemRhs(),
                           patchSolver.getSystemSolution(),
                           true );
  //Do we need to modify anything here??

  baseSolver.implicitStepSetup( time_n, dt, domain );

  patchSolver.implicitStepSetup( time_n, dt, domain );

  this->implicitStepSetup( time_n, dt, domain );

  NonlinearSolverParameters & solverParams = getNonlinearSolverParameters();
  //although these iterations are not really Newton iterations, we will use this nomeclature to keep things consistent
  integer & iter = solverParams.m_numNewtonIterations;
  iter = 0;
  bool isConverged = false;
  while( iter < solverParams.m_maxIterNewton )
  {
    if( iter == 0 )
    {
      // reset the states of all slave solvers if any of them has been reset
      //this is potentially a code duplication since resetStateToBeginningOfStep(domain) already calls the slaves
      patchSolver.resetStateToBeginningOfStep( domain );
      baseSolver.resetStateToBeginningOfStep( domain );
      resetStateToBeginningOfStep( domain );
    }

    GEOSX_LOG_LEVEL_RANK_0( 1, "\tIteration: " << iter+1 << ", BaseSolver: " );

    //we probably want to run a phase-field solve in the patch problem at timestep 0 to get a smooth initial crack. Also, re-run this
    // anytime the base crack changes
    map< std::pair< string, string >, array1d< string > > const & baseTargets = baseSolver.getReference< map< std::pair< string, string >, array1d< string > > >(
      SolverBase::viewKeyStruct::meshTargetsString());
    auto const baseTarget = baseTargets.begin()->first;
    map< std::pair< string, string >, array1d< string > > const & patchTargets = patchSolver.getReference< map< std::pair< string, string >, array1d< string > > >(
      SolverBase::viewKeyStruct::meshTargetsString());
    auto const patchTarget = patchTargets.begin()->first;
    CRSMatrix< real64, globalIndex > & patchDamageLocalMatrix = patchDamageSolver.getLocalMatrix();
    this->setInitialCrackDamageBCs( patchDamageSolver.getDofManager(), patchDamageLocalMatrix.toViewConstSizes(), domain.getMeshBody( patchTarget.first ).getBaseDiscretization(),
                                    domain.getMeshBody( baseTarget.first ).getBaseDiscretization() );
    patchDamageSolver.setInitialCrackNodes( m_nodeFixDamage );
    //patchSolver.addCustomBCDamage(m_nodeFixDamage); //this function still doesnt exist
    //must prescribe the damage boundary conditions based on the location of the base crack relative to the subdomain mesh

    //now perform the subproblem run with no BCs on displacements, just to set the damage inital condition;
    if( iter == 0 )
    {
      //this->updateNodeMaps(); //this has to be done at every prop. step
      real64 dtUseless = patchSolver.solverStep( time_n,
                                                 dtReturn,
                                                 cycleNumber,
                                                 domain );
      GEOSX_UNUSED_VAR( dtUseless );
    }

    dtReturnTemporary = baseSolver.nonlinearImplicitStep( time_n,
                                                          dtReturn,
                                                          cycleNumber,
                                                          domain );

    if( dtReturnTemporary < dtReturn )
    {
      iter = 0;
      dtReturn = dtReturnTemporary;
      continue;
    }

    if( baseSolver.getNonlinearSolverParameters().m_numNewtonIterations >= 1 && iter > 0 )
    {
      GEOSX_LOG_LEVEL_RANK_0( 1, "***** The Global-Local iterative scheme has converged in " << iter << " iterations! *****\n" );
      isConverged = true;
      break;
    }

    //here, before calling the nonlinarImplicitStep of the patch solver, we must prescribe the displacement boundary conditions
    this->prepareSubProblemBCs( domain.getMeshBody( baseTarget.first ).getBaseDiscretization(), domain.getMeshBody( patchTarget.first ).getBaseDiscretization());

    //write disp BCs to local disp solver
    patchSolidSolver.setInternalBoundaryConditions( m_nodeFixDisp, m_fixedDispList );

    //finally, pass the shared boundary information to the patch solver
    //patchSolver.
    //patchSolver.addCustomBCDisp(m_nodeFixDisp, m_fixedDispList); //this function still doesnt exist

    GEOSX_LOG_LEVEL_RANK_0( 1, "\tIteration: " << iter+1 << ", PatchSolver: " );

    //probably, a nonlinearImplicitStep is not the right one for the Phase Field solver, maybe we should call its solverStep which calls
    // splitOperatorStep
    dtReturnTemporary = patchSolver.solverStep( time_n,
                                                dtReturn,
                                                cycleNumber,
                                                domain );

    if( dtReturnTemporary < dtReturn )
    {
      iter = 0;
      dtReturn = dtReturnTemporary;
      continue;
    }

    ++iter;
  }

  GEOSX_ERROR_IF( !isConverged, "MultiResolutionHFSolver::SplitOperatorStep() did not converge" );

  baseSolver.implicitStepComplete( time_n, dt, domain );
  patchSolver.implicitStepComplete( time_n, dt, domain );
  //this->implicitStepComplete( time_n, dt, domain );

  return dtReturn;
}

// Andre - 14/03/22 - This can probably be made trivial in the MR case

// void MultiResolutionHFSolver::setNextDt( real64 const & currentDt,
//                                          real64 & nextDt )
// {

//   if( m_numResolves[0] == 0 && m_numResolves[1] == 0 )
//   {
//     this->setNextDtBasedOnNewtonIter( currentDt, nextDt );
//   }

//   GEOSX_LOG_LEVEL_RANK_0( 3, this->getName() << ": nextDt request is "  << nextDt );
// }

REGISTER_CATALOG_ENTRY( SolverBase, MultiResolutionHFSolver, string const &, Group * const )
} /* namespace geosx */
