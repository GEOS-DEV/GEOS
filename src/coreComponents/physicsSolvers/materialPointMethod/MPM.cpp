/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file MPM.cpp
 */

// Source includes
#include <physicsSolvers/materialPointMethod/MPM.hpp>
#include "MPMKernels.hpp"

#include "mpiCommunications/CommunicationTools.hpp"
#include "common/TimingMacros.hpp"
#include "common/DataTypes.hpp"
#include "finiteElement/FiniteElementDiscretizationManager.hpp"
#include "managers/DomainPartition.hpp"
#include "managers/NumericalMethodsManager.hpp"
#include "managers/GeosxState.hpp"

namespace geosx
{

namespace dataRepository
{
namespace keys
{}
}

using namespace dataRepository;

/*----------------------------------------------------------------------------------
 * MPM: Solid mechanics solver using the material point method.
 * ---------------------------------------------------------------------------------
 *
 *
 *
 * ---------------------------------------------------------------------------------
 */


/* CONSTRUCTOR
   First, let us inspect the constructor of a "MPM" object.
   This constructor does three important things:
   1 - It constructs an instance of the MPM class (here: using the SolverBase constructor and passing through the arguments).
   2 - It sets some default values for the MPM-specific private variables (here: m_fieldName and m_timeIntegrationOption).
   3 - It creates and activates a "registerWrapper" for each private variable.
   This is where the private variables are declared either as REQUIRED or OPTIONAL.
   An error is thrown if a REQUIRED variable is not specified in the XML file,
   along with the description of this variable and possible enum values if relevant.
   The description that is set is used in auto-generated documentation and console error messages.
 */

//START_SPHINX_INCLUDE_01
MPM::MPM( const string & name,
                        Group * const parent ):
  SolverBase( name, parent ),
  m_fieldName( "primaryField" ),
  m_timeIntegrationOption( TimeIntegrationOption::ImplicitTransient )
{
  registerWrapper( mpmViewKeys.timeIntegrationOption.key(), &m_timeIntegrationOption ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Time integration method. Options are:\n* " + EnumStrings< TimeIntegrationOption >::concat( "\n* " ) );

  registerWrapper( mpmViewKeys.fieldVarName.key(), &m_fieldName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of field variable" );
}
//END_SPHINX_INCLUDE_01


// Destructor
MPM::~MPM()
{
  // TODO Auto-generated destructor stub
}

/* REGISTER THE PDE SOLUTION DATA ON THE MESH
   In the MPM solver, we compute the solution of the partial differential equation "numerically".
   This means that we are not solving this PDE everywhere,
   we are computing the solution at specific locations in space.
   To do that, we have to register the Laplace solver so that it works on nodes of a mesh.
   This registration process is done here, in three steps:
   1 - for each mesh body (if the mesh is split), we collect he "nodes" (nodes carry their location information),
   2 - On nodes, we register a new property called m_fieldName and give it a type (here, the type is an array of real64)
   3 - We set some information for this property on the nodes: what is their "PlotLevel"? how can they be described?
     The PlotLevel is a flag that instructs GEOSX to export values of this property for instance so that they can be plotted.
     All properties mounted on nodes carry a certain PlotLevel value. This way, every time GEOSX triggers an
     output event (a request to "print out" data), all properties at or above a certain PlotLevel are automatically exported.
     The description here is simply an additional metadata for the newly mounted property.
 */
//START_SPHINX_INCLUDE_02
void MPM::registerDataOnMesh( Group & meshBodies )
{
  meshBodies.forSubGroups< MeshBody >( [&] ( MeshBody & meshBody )
  {
    NodeManager & nodes = meshBody.getMeshLevel( 0 ).getNodeManager();

    nodes.registerWrapper< real64_array >( m_fieldName ).
      setApplyDefaultValue( 0.0 ).
      setPlotLevel( PlotLevel::LEVEL_0 ).
      setDescription( "Primary field variable" );
  } );
}
//END_SPHINX_INCLUDE_02


/* STEPPING IN TIME
   Here, we decide how we march in time in the resolutions based on the possible
   three options set in the XML file (Steady state, Implicit transient, or Explicit transient).
   Based on these options, we can either perform an Explicit Step (forward Euler),
   or an Implicit Step (backward Euler).
   The implementation of the Explicit or Implicit Steps are found in the SolverBase.
   From now on, we oscillate between specific Laplace solver operations if implemented and more generic SolverBase operations.
   The initial values of the solver step are all at time_n, and the solver attempts to advance by a time step of dt.
   This dt time step size is specified initially by the user; and unfortunately, it can sometimes be too large for convergence.
   The SolverStep method thus returns the time step value that is was actually capable of solving for with good convergence.
 */

real64 MPM::solverStep( real64 const & time_n,
                               real64 const & dt,
                               const int cycleNumber,
                               DomainPartition & domain )
{
  real64 dtReturn = dt;
  if( m_timeIntegrationOption == TimeIntegrationOption::ExplicitTransient )
  {
    dtReturn = explicitStep( time_n, dt, cycleNumber, domain );
  }
  else if( m_timeIntegrationOption == TimeIntegrationOption::ImplicitTransient ||
           m_timeIntegrationOption == TimeIntegrationOption::SteadyState )
  {
	std::cout<<"MPM solver currently supports only ExplicitTransient Integration Option! Sorry."<<std::endl;
    dtReturn = this->linearImplicitStep( time_n, dt, cycleNumber, domain );
  }
  return dtReturn;
}

real64 MPM::explicitStep(
  real64 const & GEOSX_UNUSED_PARAM( time_n ),
  real64 const & dt,
  const int GEOSX_UNUSED_PARAM( cycleNumber ),
  DomainPartition & GEOSX_UNUSED_PARAM( domain ) )
{
  return dt;
}


/*
   IMPLICIT STEP SETUP
   This method uses the system setup from MPM (see below).
   It "deactivates" the time variables (with the GEOSX_UNUSED_PARAM macro) and does a steady state system set-up.
 */
void MPM::implicitStepSetup( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                    real64 const & GEOSX_UNUSED_PARAM( dt ),
                                    DomainPartition & domain )
{
  // Computation of the sparsity pattern
  setupSystem( domain, m_dofManager, m_localMatrix, m_localRhs, m_localSolution );
}

void MPM::implicitStepComplete( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                       real64 const & GEOSX_UNUSED_PARAM( dt ),
                                       DomainPartition & GEOSX_UNUSED_PARAM( domain ) )
{}

void MPM::setupDofs( DomainPartition const & GEOSX_UNUSED_PARAM( domain ),
                            DofManager & dofManager ) const
{
  dofManager.addField( m_fieldName,
                       DofManager::Location::Node );

  dofManager.addCoupling( m_fieldName,
                          m_fieldName,
                          DofManager::Connector::Elem );
}

/* SETUP SYSTEM
   Setting up the system using the base class method
 */

void MPM::setupSystem( DomainPartition & domain,
                              DofManager & dofManager,
                              CRSMatrix< real64, globalIndex > & localMatrix,
                              array1d< real64 > & localRhs,
                              array1d< real64 > & localSolution,
                              bool const setSparsity )
{
  GEOSX_MARK_FUNCTION;
  SolverBase::setupSystem( domain, dofManager, localMatrix, localRhs, localSolution, setSparsity );

  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  NodeManager const & nodeManager = mesh.getNodeManager();
  arrayView1d< globalIndex const > const &
  dofIndex = nodeManager.getReference< globalIndex_array >( dofManager.getKey( m_fieldName ) );

  SparsityPattern< globalIndex > sparsityPattern( dofManager.numLocalDofs(),
                                                  dofManager.numGlobalDofs(),
                                                  8*8*3 );

  finiteElement::fillSparsity< CellElementSubRegion,
                               MPMKernel >( mesh,
                                                   targetRegionNames(),
                                                   this->getDiscretizationName(),
                                                   dofIndex,
                                                   dofManager.rankOffset(),
                                                   sparsityPattern );

  sparsityPattern.compress();
  localMatrix.assimilate< parallelDevicePolicy<> >( std::move( sparsityPattern ) );

}


/*
   ASSEMBLE SYSTEM
   This is the most important method to assemble the matrices needed before sending them to our solver.
   For a system A.x = B (with x the unknown), here, we use:
   - A : "localMatrix" this represents a Compressed Row Storage (optimized for sparse) matrix of real64 values associated with their index,
   - B : "localRhs" this represents a vector (1d array) of real64 numbers specified at the equation's right-hand side.
   The "local" prefix indicates that we are working on a local problem here, and the parallelization is performed at a higher level.
   This assembly step collects all the information needed to create the matrices localMatrix and localRhs, and the computation of values
   is done in a specific Laplace kernel optimized for parallel performance. Here we:
   1 - identify and point to the mesh of this domain,
   2 - find the node manager of this mesh,
   3 - extract the indices of the nodes that will be solved for (ie. the degrees of freedom or "dof")
   4 - pass all this information to a Laplace-specific finite element computation kernel.
   The call to the kernel is a templated call designed for performance (we will not explain the kernel here).
   See the implementation in MPMKernel.cpp.
 */
//START_SPHINX_INCLUDE_04
void MPM::assembleSystem( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                 real64 const GEOSX_UNUSED_PARAM( dt ),
                                 DomainPartition & domain,
                                 DofManager const & dofManager,
                                 CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                 arrayView1d< real64 > const & localRhs )
{
  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  NodeManager & nodeManager = mesh.getNodeManager();

  arrayView1d< globalIndex const > const &
  dofIndex =  nodeManager.getReference< array1d< globalIndex > >( dofManager.getKey( m_fieldName ) );


  finiteElement::
    regionBasedKernelApplication< parallelDevicePolicy< 32 >,
                                  constitutive::NullModel,
                                  CellElementSubRegion,
                                  MPMKernel >( mesh,
                                                      targetRegionNames(),
                                                      this->getDiscretizationName(),
                                                      arrayView1d< string const >(),
                                                      dofIndex,
                                                      dofManager.rankOffset(),
                                                      localMatrix,
                                                      localRhs,
                                                      m_fieldName );



  //END_SPHINX_INCLUDE_04
}

void MPM::applySystemSolution( DofManager const & dofManager,
                                      arrayView1d< real64 const > const & localSolution,
                                      real64 const scalingFactor,
                                      DomainPartition & domain )
{
  dofManager.addVectorToField( localSolution,
                               m_fieldName,
                               m_fieldName,
                               scalingFactor );

  // Synchronize ghost nodes
  std::map< string, string_array > fieldNames;
  fieldNames["node"].emplace_back( m_fieldName );

  getGlobalState().getCommunicationTools().synchronizeFields( fieldNames,
                                                              domain.getMeshBody( 0 ).getMeshLevel( 0 ),
                                                              domain.getNeighbors(),
                                                              true );
}

/*
   APPLY BOUNDARY CONDITIONS
   Here, this call is the generic call from SolverBase.
   All it does is to call a specific Dirichlet boundary condition implemented for this solver
 */
void MPM::applyBoundaryConditions( real64 const time_n,
                                          real64 const dt,
                                          DomainPartition & domain,
                                          DofManager const & dofManager,
                                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                          arrayView1d< real64 > const & localRhs )
{
  applyDirichletBCImplicit( time_n + dt, dofManager, domain, localMatrix, localRhs );
}

/*
   SOLVE SYSTEM
   This method is simply initiating the solution and right-hand side
   and pass is to the base class solver.
 */
void MPM::solveSystem( DofManager const & dofManager,
                              ParallelMatrix & matrix,
                              ParallelVector & rhs,
                              ParallelVector & solution )
{
  rhs.scale( -1.0 ); // TODO decide if we want this here
  solution.zero();
  SolverBase::solveSystem( dofManager, matrix, rhs, solution );
}

/*
   DIRICHLET BOUNDARY CONDITIONS
   This is the boundary condition method applied for this particular solver.
   It is called by the more generic "applyBoundaryConditions" method.
 */
void MPM::applyDirichletBCImplicit( real64 const time,
                                           DofManager const & dofManager,
                                           DomainPartition & domain,
                                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                           arrayView1d< real64 > const & localRhs )
{
  FieldSpecificationManager const & fsManager = getGlobalState().getFieldSpecificationManager();

  fsManager.apply( time,
                   domain,
                   "nodeManager",
                   m_fieldName,
                   [&]( FieldSpecificationBase const & bc,
                        string const &,
                        SortedArrayView< localIndex const > const & targetSet,
                        Group & targetGroup,
                        string const & GEOSX_UNUSED_PARAM( fieldName ) )
  {
    bc.applyBoundaryConditionToSystem< FieldSpecificationEqual, parallelDevicePolicy< 32 > >( targetSet,
                                                                                              time,
                                                                                              targetGroup,
                                                                                              m_fieldName,
                                                                                              dofManager.getKey( m_fieldName ),
                                                                                              dofManager.rankOffset(),
                                                                                              localMatrix,
                                                                                              localRhs );
  } );
}

void MPM::resetStateToBeginningOfStep( DomainPartition & GEOSX_UNUSED_PARAM( domain ) )
{}

//START_SPHINX_INCLUDE_00
REGISTER_CATALOG_ENTRY( SolverBase, MPM, string const &, Group * const )
//END_SPHINX_INCLUDE_00
} /* namespace ANST */
