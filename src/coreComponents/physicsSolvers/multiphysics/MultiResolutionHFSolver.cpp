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
#include "physicsSolvers/solidMechanics/SolidMechanicsEFEMKernelsHelper.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsEmbeddedFractures.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "physicsSolvers/surfaceGeneration/SurfaceGenerator.hpp"
#include "linearAlgebra/utilities/LAIHelperFunctions.hpp"



namespace geosx
{

namespace
{
//this function is not part of any class yet - just tries to find the number of a node given its coordinates
  localIndex getNodeNumberFromCoordinate( MeshLevel const & mesh, real64 const & x, real64 const & y, real64 const & z)
  {
    real64 tolerance = 1e-10;
    NodeManager const & meshNodeManager = mesh.getNodeManager();
    arrayView1d<real64 const> const nodalPositions = meshNodeManager.referencePosition();
    for (localIndex a=0; a < meshNodeManager.size(); a++)
      {
	real64 dist_sq = (nodalPositions( a,0 ) - x)*(nodalPositions( a,0 ) - x) + (nodalPositions( a,1 ) - y)*(nodalPositions( a,1 ) - y) + (nodalPositions( a,2 ) - z)*(nodalPositions( a,2 ) - z);
	real64 dist = sqrt(dist_sq);
	if (dist < tolerance)
	  {
	    return a;
	  }
      }
    return -1;
  }

  //this function checks if an index is already in an array
  bool isMember( localIndex const & ind, array1d<localIndex> const & arr)
  {
    for ( localIndex iter : arr)
      {
	if (iter == ind)
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
  m_contactRelationName(),
  m_surfaceGeneratorName(),
  m_nodeFixDamage(),
  m_nodeFixDisp(),
  m_fixedDispList(),
  m_nodeMapIndices(),
  m_nodeMapWeights(),
  m_surfaceGenerator( nullptr ),
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

  registerWrapper( viewKeyStruct::surfaceGeneratorNameString(), &m_surfaceGeneratorName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the surface generator to use in the multiresolutionHF solver" );

  registerWrapper( viewKeyStruct::contactRelationNameString(), &m_contactRelationName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of contact relation to enforce constraints on fracture boundary." );

  registerWrapper( viewKeyStruct::maxNumResolvesString(), &m_maxNumResolves ).
    setApplyDefaultValue( 10 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Value to indicate how many resolves may be executed to perform surface generation after the execution of base and patch scale solvers. " );
 
}

void MultiResolutionHFSolver::RegisterDataOnMesh( dataRepository::Group & MeshBodies )
{
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
  m_surfaceGenerator = &this->getParent().getGroup< SurfaceGenerator >( m_surfaceGeneratorName );
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
void MultiResolutionHFSolver::updateNodeMaps( MeshLevel const & base,
						MeshLevel const & patch )
{
  // get list of nodes on the boundary of the patch
  NodeManager const & patchNodeManager = patch.getNodeManager();
  NodeManager const & baseNodeManager = base.getNodeManager();
  SortedArrayView< localIndex const > const patchExternalSet = patchNodeManager.externalSet();
  arrayView1d<real64 const> const patchPosition = patchNodeManager.referencePosition();
  for (localIndex a : patchExternalSet )  
    {
      localIndex A = getNodeNumberFromCoordinate(base, patchPosition(a,0), patchPosition(a,1), patchPosition(a,2));
      m_nodeMapIndices(a) = A;
      //since the boundary points lie on element edges, these will have dimension 2. weights can be computed by taking the ratio between distances to 2 closes points from base mesh.
      m_nodeMapWeights(a) = 1; //for now, we dont need this, so we just use 1.
    }  
}

//Andre - 03/15 - this function will loop over all subdomain nodes and check their distance to the prescribed discrete crack. If the distance is smaller than 1 element size (subdomain), we set the damage in this node to be fixed at 1. 
void MultiResolutionHFSolver::setInitialCrackDamageBCs( MeshLevel const & patch, MeshLevel const & base )
{

  // get list of nodes on the boundary of the patch
  NodeManager const & patchNodeManager = patch.getNodeManager();
  ElementRegionManager const & baseElemManager = base.getElemManager();
  baseElemManager.forElementSubRegions< CellElementSubRegion >( regionNames, [&](localIndex const, CellElementSubRegion const & cellElementSubRegion)
  {
    SortedArrayView< localIndex const > const fracturedElements = cellElementSubRegion.fracturedElementsList();
    m_nodeFixDamage.resize( numNodesPerElement()*fracturedElements.size() );
    localIndex count = 0;
    for (localIndex a : fracturedElements)
      {
	//get all nodes of fracturedElements(a)
	for (localIndex b=0; b < numNodesPerElement(); b++)
	  {
	    cellElementSubRegion.nodeList( a,b )      
	    //append dofs associated with node b of element a
	      if (!isMember(a,m_nodeFixDamage))
	      {
		m_nodeFixDamage(count) = a;
		++count;
	      }
	  } 
      }
  });
  m_nodeFixDamage.resize( count );
}  
  
//Andre - 03/15 - this function needs to access u_base (hopefully via node manager), the patch mesh (via node manager too),
//a map that translate patch node numbers to base node numbers (could be a MultiResolutionHFSolver member) and also the
//damage field in the patch domain (via node manager)
//The function will read all this data and prepare a list of dofs and u values that will be used by the patch solver to set
//the boundary conditions 
void MultiResolutionHFSolver::prepareSubProblemBCs( MeshLevel const & base,
                                                    MeshLevel const & patch )
{

  // get list of nodes on the boundary of the patch
  NodeManager const & patchNodeManager = patch.getNodeManager();
  SortedArrayView< localIndex const > const patchExternalSet = patchNodeManager.externalSet();
  //arrayView1d< globalIndex const > const patchLocalToGlobalMap = patchNodeManager.localToGlobalMap();
  arrayView1d< real64 const > const patchDamage = patchNodeManager.getReference<array1d<real64>( "damage" );
  NodeManager const & baseNodeManager = base.getNodeManager();
  arrayView1d<real64 const> const baseDisp = baseNodeManager.totalDisplacement();
  //unordered_map< globalIndex, localIndex > const & patchGlobalToLocalMap = baseNodeManager.globalToLocalMap();

  //we need to zero the lists to erase data from previous step
  //this->eraseLists(); //maybe this can go inside resetToBeginningOfStep - already there() 
  real64 damage_threshold = 0.3;

  m_nodeFixDisp.resize( patchExternalSet.size() );
  m_fixedDispList.resize( patchExternalSet.size(), 3 );
  localIndex count=0;
  for( localIndex a : patchExternalSet )
  {
      if( patchDamage[a] < damage_threshold )
      {
        // NOTE: there needs to be a translation between patch and base mesh for the indices and values/weights.

        //append dofs associated with node i
        m_nodeFixDisp(count) = a;
        localIndex const numBaseNodes = m_nodeMapIndices.sizeOfArray( a );
        for( localIndex b=0; b<numBaseNodes; ++b )
        {
          localIndex const B = m_nodeMapIndices[b]; // base node number associated to patch node i
          m_fixedDispList(count,0) += m_nodeMapWeights(a,b) * baseDisp(B,0);
          m_fixedDispList(count,1) += m_nodeMapWeights(a,b) * baseDisp(B,1);
          m_fixedDispList(count,2) += m_nodeMapWeights(a,b) * baseDisp(B,2);
        }

        ++count;
      }
    }

  m_nodeFixDisp.resize( count );
  m_fixedDispList.resize( count, 3 );
  
}  

real64 MultiResolutionHFSolver::splitOperatorStep( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                               real64 const & dt,
                                               integer const GEOSX_UNUSED_PARAM( cycleNumber ),
                                               DomainPartition & GEOSX_UNUSED_PARAM( domain ) )
{
  GEOSX_MARK_FUNCTION;
  real64 dtReturn = dt;
  real64 dtReturnTemporary;

  SolidMechanicsEmbeddedFractures &
  baseSolver = this->getParent().getGroup< SolidMechanicsEmbeddedFractures >( m_baseSolverName );

  PhaseFieldFractureSolver &
  patchSolver = this->getParent().getGroup< PhaseFieldFractureSolver >( m_patchSolverName );

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
                           patchSolver.getSystemSolution() );
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

    //we probably want to run a phase-field solve in the patch problem at timestep 0 to get a smooth initial crack. Also, re-run this anytime the base crack changes
    this->setInitialCrackDamageBCs();

    patchSolver.addCustomBCDamage(m_nodeFixDamage); //this function still doesnt exist
    //must prescribe the damage boundary conditions based on the location of the base crack relative to the subdomain mesh

    //now perform the subproblem run with no BCs on displacements, just to set the damage inital condition;
    if (iter == 1){
      this->updateNodeMap(); //this has to be done at every prop. step
      Real64 dtUseless = patchSolver.nonlinearImplicitStep(time_n,
							 dtReturn,
							 cycleNumber,
							 domain );
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

    if( baseSolver.getNonlinearSolverParameters().m_numNewtonIterations == 0 && iter > 0 )
    {
      GEOSX_LOG_LEVEL_RANK_0( 1, "***** The Global-Local iterative scheme has converged in " << iter << " iterations! *****\n" );
      isConverged = true;
      break;
    }

    //here, before calling the nonlinarImplicitStep of the patch solver, we must prescribe the displacement boundary conditions
    this->prepareSubProblemBCs();

    //finally, pass the shared boundary information to the patch solver
    //patchSolver.
    patchSolver.addCustomBCDisp(m_nodeFixDisp, m_fixedDispList); //this function still doesnt exist
    
    GEOSX_LOG_LEVEL_RANK_0( 1, "\tIteration: " << iter+1 << ", PatchSolver: " );

    //probably, a nonlinearImplicitStep is not the right one for the Phase Field solver, maybe we should call its solverStep which calls splitOperatorStep
    dtReturnTemporary = patchSolver.nonlinearImplicitStep(  time_n,
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
  this->implicitStepComplete( time_n, dt, domain );

  return dtReturn;
}

real64 MultiResolutionHFSolver::explicitStep( real64 const & time_n,
                                          real64 const & dt,
                                          const int cycleNumber,
                                          DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;
  m_baseSolver->solverStep( time_n, dt, cycleNumber, domain );
  m_patchSolver->solverStep( time_n, dt, cycleNumber, domain );

  return dt;
}


void MultiResolutionHFSolver::setupDofs( DomainPartition const & domain,
                                     DofManager & dofManager ) const
{
  GEOSX_MARK_FUNCTION;
  m_baseSolver->setupDofs( domain, dofManager );
  m_patchSolver->setupDofs( domain, dofManager );
}

// Andre - 14/03/22 - This can probably be made trivial in the MR case
void MultiResolutionHFSolver::setupSystem( DomainPartition & domain,
                                       DofManager & dofManager,
                                       CRSMatrix< real64, globalIndex > & localMatrix,
                                       ParallelVector & rhs,
                                       ParallelVector & solution,
                                       bool const setSparsity )
{
  GEOSX_MARK_FUNCTION;

  GEOSX_UNUSED_VAR( setSparsity );

  dofManager.setDomain( domain );

  setupDofs( domain, dofManager );
  dofManager.reorderByRank();

  localIndex const numLocalRows = dofManager.numLocalDofs();

  SparsityPattern< globalIndex > patternOriginal;
  dofManager.setSparsityPattern( patternOriginal );

  // Get the original row lengths (diagonal blocks only)
  array1d< localIndex > rowLengths( patternOriginal.numRows() );
  for( localIndex localRow = 0; localRow < patternOriginal.numRows(); ++localRow )
  {
    rowLengths[localRow] = patternOriginal.numNonZeros( localRow );
  }

  // Add the number of nonzeros induced by coupling
  //addFluxApertureCouplingNNZ( domain, dofManager, rowLengths.toView() );

  // Create a new pattern with enough capacity for coupled matrix
  SparsityPattern< globalIndex > pattern;
  pattern.resizeFromRowCapacities< parallelHostPolicy >( patternOriginal.numRows(),
                                                         patternOriginal.numColumns(),
                                                         rowLengths.data() );

  // Copy the original nonzeros
  for( localIndex localRow = 0; localRow < patternOriginal.numRows(); ++localRow )
  {
    globalIndex const * cols = patternOriginal.getColumns( localRow ).dataIfContiguous();
    pattern.insertNonZeros( localRow, cols, cols + patternOriginal.numNonZeros( localRow ) );
  }

  // Add the nonzeros from coupling
  //addFluxApertureCouplingSparsityPattern( domain, dofManager, pattern.toView() );

  localMatrix.assimilate< parallelDevicePolicy<> >( std::move( pattern ) );
  localMatrix.setName( this->getName() + "/matrix" );

  rhs.setName( this->getName() + "/rhs" );
  rhs.create( numLocalRows, MPI_COMM_GEOSX );

  solution.setName( this->getName() + "/solution" );
  solution.create( numLocalRows, MPI_COMM_GEOSX );

  //setUpDflux_dApertureMatrix( domain, dofManager, localMatrix );
}

//Andre - 03/16 - Does this function ever gets called?
void MultiResolutionHFSolver::assembleSystem( real64 const time,
                                          real64 const dt,
                                          DomainPartition & domain,
                                          DofManager const & dofManager,
                                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                          arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  m_baseSolver->assembleSystem( time,
                                 dt,
                                 domain,
                                 dofManager,
                                 localMatrix,
                                 localRhs );

  m_patchSolver->assembleSystem( time,
				 dt,
				 domain,
                                 dofManager,
                                 localMatrix,
                                 localRhs );
}

//Andre - 03/16 - Does this function ever gets called?
void MultiResolutionHFSolver::applyBoundaryConditions( real64 const time,
                                                   real64 const dt,
                                                   DomainPartition & domain,
                                                   DofManager const & dofManager,
                                                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                   arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  m_baseSolver->applyBoundaryConditions( time,
                                          dt,
                                          domain,
                                          dofManager,
                                          localMatrix,
                                          localRhs );

  m_patchSolver->applyBoundaryConditions( time,
                                         dt,
                                         domain,
                                         dofManager,
                                         localMatrix,
                                         localRhs );
}

//Andre - 03/16 - Does this function need to be implemented? It is just calling the base class function, which should be done automatically
void
MultiResolutionHFSolver::
  applySystemSolution( DofManager const & dofManager,
                       arrayView1d< real64 const > const & localSolution,
                       real64 const scalingFactor,
                       DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;
  SolverBase::applySystemSolution( dofManager, localSolution, scalingFactor, domain );
  
}

void MultiResolutionHFSolver::updateState( DomainPartition & domain )
{
  m_baseSolver->updateState( domain );
  m_patchSolver->updateState( domain );
}



// Andre Costa - 14/03/22 - we probably need one of these for the local solver as well, not sure what is the best design
real64
MultiResolutionHFSolver::scalingForSystemSolution( DomainPartition const & domain,
                                               DofManager const & dofManager,
                                               arrayView1d< real64 const > const & localSolution )
{
  return m_baseSolver->scalingForSystemSolution( domain,
                                                  dofManager,
                                                  localSolution );
}

void MultiResolutionHFSolver::setNextDt( real64 const & currentDt,
                                     real64 & nextDt )
{

  if( m_numResolves[0] == 0 && m_numResolves[1] == 0 )
  {
    this->setNextDtBasedOnNewtonIter( currentDt, nextDt );
  }
  else
  {
    SolverBase & surfaceGenerator = this->getParent().getGroup< SolverBase >( "SurfaceGen" );
    nextDt = surfaceGenerator.GetTimestepRequest() < 1e99 ? surfaceGenerator.GetTimestepRequest() : currentDt;
  }

  GEOSX_LOG_LEVEL_RANK_0( 3, this->getName() << ": nextDt request is "  << nextDt );
}

void MultiResolutionHFSolver::initializeNewFaceElements( DomainPartition const & )
{
//  m_flowSolver->
}


REGISTER_CATALOG_ENTRY( SolverBase, MultiResolutionHFSolver, string const &, Group * const )
} /* namespace geosx */
