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

using namespace dataRepository;
using namespace constitutive;

MultiResolutionHFSolver::MultiResolutionHFSolver( const string & name,
                                          Group * const parent ):
  SolverBase( name, parent ),
  m_contactRelationName(),
  m_surfaceGeneratorName(),
  m_surfaceGenerator( nullptr ),
  m_maxNumResolves( 10 )
{
  registerWrapper( viewKeyStruct::surfaceGeneratorNameString(), &m_surfaceGeneratorName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the surface generator to use in the multiresolutionHF solver" );

  registerWrapper( viewKeyStruct::contactRelationNameString(), &m_contactRelationName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of contact relation to enforce constraints on fracture boundary." );

  registerWrapper( viewKeyStruct::maxNumResolvesString(), &m_maxNumResolves ).
    setApplyDefaultValue( 10 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Value to indicate how many resolves may be executed to perform surface generation after the execution of global and local scale solvers. " );
 
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
  m_globalSolver->implicitStepSetup( time_n, dt, domain );
  m_localSolver->implicitStepSetup( time_n, dt, domain );

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
  m_globalSolver->resetStateToBeginningOfStep( domain );
  m_localSolver->resetStateToBeginningOfStep( domain );
  //reset lists of dofs to empty to look up new boundary conditions for subproblem
  m_dofListDisp = [];
  m_fixedDispList = []; 
  m_dofListDamage = [];
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
void MultiResolutionHFSolver::updateNodeMap()
{
  for i in innerMesh //accessible through nodeManager
    {
      xi,yi,zi = getXiYiZi(); //get cartesian coordinates of node i
      //sometimes, there is no outerMesh node at xi,yi,zi, in this case, we can just return -1 in the function below.
      int N = getNodeNumberFromCoordinate(outerMesh, xi, yi, zi);
      //WARNING: this function may need a tolerance, which should probably be proportional to element size
      m_nodeMap[i] = N;
    }  
}

//Andre - 03/15 - this function will loop over all subdomain nodes and check their distance to the prescribed discrete crack. If the distance is smaller than 1 element size (subdomain), we set the damage in this node to be fixed at 1. 
void MultiResolutionHFSolver::setInitialCrackDamageBCs()
{

  //we need to zero the lists to erase data from previous step
  //this->eraseLists(); //maybe this can go inside resetToBeginningOfStep() - already there! 
  for i=allNodesInPatch
    {
      Real64 dist_i = compute_dist_to_frac(i);
      if dist_i < h_patch //small element size - we must choose a good distance here
      {
      //append dofs associated with node i
      m_dofListDamage.append(i); 
      }
    }
}  
  
//Andre - 03/15 - this function needs to access u_global (hopefully via node manager), the local mesh (via node manager too),
//a map that translate local node numbers to global node numbers (could be a MultiResolutionHFSolver member) and also the
//damage field in the local domain (via node manager)
//The function will read all this data and prepare a list of dofs and u values that will be used by the local solver to set
//the boundary conditions 
void MultiResolutionHFSolver::prepareSubProblemBCs()
{

  //we need to zero the lists to erase data from previous step
  //this->eraseLists(); //maybe this can go inside resetToBeginningOfStep - already there() 
  Real64 damage_threshold = 0.3;
  Array1d<Real64> uGlobal = retrieve_global_disp_field(); //pseudo code
  Array1d<Real64> allNodesInLocalBoundary = getNodeNumbersOfAllNodesInPatchBoundary();
  for i=allNodesInLocalBoundary
    {
      Real64 damage_at_i = retrieve_patch_damage_field_at_i();  
      if damage_at_i > damage_threshold
      {
        //do nothing;
	continue;
      }
      //append dofs associated with node i
      m_dofListDisp.append(3*i-2); //ux //here, we can have a dof or node list (node list is shorter but require translating node to dof later to
      m_dofListDisp.append(3*i-1); //uy
      m_dofListDisp.append(3*i);   //uz
      //apply the BC
      int I = m_nodeMap[i]; // global node number associated to local node i
      m_fixedDispList.append(uGlobal(3*I-2)); 
      m_fixedDispList.append(uGlobal(3*I-1));
      m_fixedDispList.append(uGlobal(3*I));
    } 
  
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
  globalSolver = this->getParent().getGroup< SolidMechanicsEmbeddedFractures >( m_globalSolverName );

  PhaseFieldFractureSolver &
  localSolver = this->getParent().getGroup< PhaseFieldFractureSolver >( m_localSolverName );

  globalSolver.setupSystem( domain,
                            globalSolver.getDofManager(),
                            globalSolver.getLocalMatrix(),
                            globalSolver.getSystemRhs(),
                            globalSolver.getSystemSolution(),
                            true );

  localSolver.setupSystem( domain,
                           localSolver.getDofManager(),
                           localSolver.getLocalMatrix(),
                           localSolver.getSystemRhs(),
                           localSolver.getSystemSolution() );
  //Do we need to modify anything here??
  
  globalSolver.implicitStepSetup( time_n, dt, domain );

  localSolver.implicitStepSetup( time_n, dt, domain );

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
      localSolver.resetStateToBeginningOfStep( domain );
      globalSolver.resetStateToBeginningOfStep( domain );
      resetStateToBeginningOfStep( domain );
    }

    GEOSX_LOG_LEVEL_RANK_0( 1, "\tIteration: " << iter+1 << ", GlobalSolver: " );

    //we probably want to run a phase-field solve in the patch problem at timestep 0 to get a smooth initial crack. Also, re-run this anytime the global crack changes
    this->setInitialCrackDamageBCs();

    localSolver.addCustomBCDamage(m_dofListDamage); //this function still doesnt exist
    //must prescribe the damage boundary conditions based on the location of the global crack relative to the subdomain mesh

    //now perform the subproblem run with no BCs on displacements, just to set the damage inital condition;
    if (iter == 1){
      this->updateNodeMap(); //this has to be done at every prop. step
      Real64 dtUseless = localSolver.nonlinearImplicitStep(time_n,
							 dtReturn,
							 cycleNumber,
							 domain );
    }

    dtReturnTemporary = globalSolver.nonlinearImplicitStep( time_n,
                                                           dtReturn,
                                                           cycleNumber,
                                                           domain );

    if( dtReturnTemporary < dtReturn )
    {
      iter = 0;
      dtReturn = dtReturnTemporary;
      continue;
    }

    if( globalSolver.getNonlinearSolverParameters().m_numNewtonIterations == 0 && iter > 0 )
    {
      GEOSX_LOG_LEVEL_RANK_0( 1, "***** The Global-Local iterative scheme has converged in " << iter << " iterations! *****\n" );
      isConverged = true;
      break;
    }

    //here, before calling the nonlinarImplicitStep of the local solver, we must prescribe the displacement boundary conditions
    this->prepareSubProblemBCs();

    //finally, pass the shared boundary information to the local solver
    //localSolver.
    localSolver.addCustomBCDisp(m_dofListDisp, m_fixedDispList); //this function still doesnt exist
    
    GEOSX_LOG_LEVEL_RANK_0( 1, "\tIteration: " << iter+1 << ", LocalSolver: " );

    //probably, a nonlinearImplicitStep is not the right one for the Phase Field solver, maybe we should call its solverStep which calls splitOperatorStep
    dtReturnTemporary = localSolver.nonlinearImplicitStep(  time_n,
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

  globalSolver.implicitStepComplete( time_n, dt, domain );
  localSolver.implicitStepComplete( time_n, dt, domain );
  this->implicitStepComplete( time_n, dt, domain );

  return dtReturn;
}

real64 MultiResolutionHFSolver::explicitStep( real64 const & time_n,
                                          real64 const & dt,
                                          const int cycleNumber,
                                          DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;
  m_globalSolver->solverStep( time_n, dt, cycleNumber, domain );
  m_localSolver->solverStep( time_n, dt, cycleNumber, domain );

  return dt;
}


void MultiResolutionHFSolver::setupDofs( DomainPartition const & domain,
                                     DofManager & dofManager ) const
{
  GEOSX_MARK_FUNCTION;
  m_globalSolver->setupDofs( domain, dofManager );
  m_localSolver->setupDofs( domain, dofManager );
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


void MultiResolutionHFSolver::assembleSystem( real64 const time,
                                          real64 const dt,
                                          DomainPartition & domain,
                                          DofManager const & dofManager,
                                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                          arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  m_globalSolver->assembleSystem( time,
                                 dt,
                                 domain,
                                 dofManager,
                                 localMatrix,
                                 localRhs );

  m_localSolver->assembleSystem( time,
				 dt,
				 domain,
                                 dofManager,
                                 localMatrix,
                                 localRhs );
}

void MultiResolutionHFSolver::applyBoundaryConditions( real64 const time,
                                                   real64 const dt,
                                                   DomainPartition & domain,
                                                   DofManager const & dofManager,
                                                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                   arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  m_globalSolver->applyBoundaryConditions( time,
                                          dt,
                                          domain,
                                          dofManager,
                                          localMatrix,
                                          localRhs );

  m_localSolver->applyBoundaryConditions( time,
                                         dt,
                                         domain,
                                         dofManager,
                                         localMatrix,
                                         localRhs );
}

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
  m_globalSolver->updateState( domain );
  m_localSolver->updateState( domain );
}



// Andre Costa - 14/03/22 - we probably need one of these for the local solver as well, not sure what is the best design
real64
MultiResolutionHFSolver::scalingForSystemSolution( DomainPartition const & domain,
                                               DofManager const & dofManager,
                                               arrayView1d< real64 const > const & localSolution )
{
  return m_globalSolver->scalingForSystemSolution( domain,
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
