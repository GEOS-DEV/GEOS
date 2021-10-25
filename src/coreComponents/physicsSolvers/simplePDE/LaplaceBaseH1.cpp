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
 * @file LaplaceBaseH1.cpp
 */

#include "LaplaceBaseH1.hpp"

#include "dataRepository/InputFlags.hpp"
#include "mainInterface/GeosxState.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"

namespace geosx
{

using namespace dataRepository;

//START_SPHINX_INCLUDE_CONSTRUCTOR
LaplaceBaseH1::LaplaceBaseH1( const string & name,
                              Group * const parent ):
  SolverBase( name, parent ),
  m_fieldName( "primaryField" ),
  m_timeIntegrationOption( TimeIntegrationOption::ImplicitTransient )
{
  this->registerWrapper( viewKeyStruct::timeIntegrationOption(), &m_timeIntegrationOption ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Time integration method. Options are:\n* " + EnumStrings< TimeIntegrationOption >::concat( "\n* " ) );

  this->registerWrapper( viewKeyStruct::fieldVarName(), &m_fieldName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of field variable" );

}
//END_SPHINX_INCLUDE_CONSTRUCTOR

LaplaceBaseH1::~LaplaceBaseH1()
{
  // TODO Auto-generated destructor stub
}

/* REGISTER THE PDE SOLUTION DATA ON THE MESH
   In the Laplace solver, we compute the solution of the partial differential equation
   "numerically".  This means that we are not solving this PDE everywhere, we are computing the
   solution at specific locations in space.  To do that, we have to register the Laplace solver so
   that it works on nodes of a mesh.  This registration process is done here, in three steps:
   1 - for each mesh body (if the mesh is split), we collect the "nodes" (nodes carry their location
   information);
   2 - On nodes, we register a new property called m_fieldName and give it a type (here, the type is
   an array of real64);
   3 - We set some information for this property on the nodes: what is their "PlotLevel"? how can
   they be described?  The PlotLevel is a flag that instructs GEOSX to export values of this
   property for instance so that they can be plotted.  All properties mounted on nodes carry a
   certain PlotLevel value. This way, every time GEOSX triggers an output event (a request to "print
   out" data), all properties at or above a certain PlotLevel are automatically exported.  The
   description here is simply an additional metadata for the newly mounted property.
 */
//START_SPHINX_INCLUDE_REGISTERDATAONMESH
void LaplaceBaseH1::registerDataOnMesh( Group & meshBodies )
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
//END_SPHINX_INCLUDE_REGISTERDATAONMESH


/* STEPPING IN TIME
   Here, we decide how we march in time in the resolutions based on the possible two options set in
   the XML file (Steady state or Implicit transient). In the case of Implicit transient, we perform
   an implicit step (backward Euler). The implementation of the Implicit Step is found in the
   SolverBase.  From now on, we oscillate between specific Laplace solver operations if implemented
   and more generic SolverBase operations.  The initial values of the solver step are all at time_n,
   and the solver attempts to advance by a time step of dt.  This dt time step size is specified
   initially by the user and the solverStep method also returns its value.
 */

real64 LaplaceBaseH1::solverStep( real64 const & time_n,
                                  real64 const & dt,
                                  const int cycleNumber,
                                  DomainPartition & domain )
{
  return this->linearImplicitStep( time_n, dt, cycleNumber, domain );
}

/*
   IMPLICIT STEP SETUP
   This method uses the system setup from SolverBase (see below).
   The current system of this class does not use the time variable. The macro GEOSX_UNUSED_PARAM is
   therefore used here to avoid a compilation error.
 */
void LaplaceBaseH1::implicitStepSetup( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                       real64 const & GEOSX_UNUSED_PARAM( dt ),
                                       DomainPartition & domain )
{
  // Computation of the sparsity pattern
  setupSystem( domain, m_dofManager, m_localMatrix, m_localRhs, m_localSolution );
}

void LaplaceBaseH1::implicitStepComplete( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                          real64 const & GEOSX_UNUSED_PARAM( dt ),
                                          DomainPartition & GEOSX_UNUSED_PARAM( domain ) )
{}

void LaplaceBaseH1::setupDofs( DomainPartition const & GEOSX_UNUSED_PARAM( domain ),
                               DofManager & dofManager ) const
{
  dofManager.addField( m_fieldName,
                       DofManager::Location::Node,
                       1,
                       targetRegionNames() );

  dofManager.addCoupling( m_fieldName,
                          m_fieldName,
                          DofManager::Connector::Elem );
}

void LaplaceBaseH1::applySystemSolution( DofManager const & dofManager,
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

void LaplaceBaseH1::updateState( DomainPartition & domain )
{
  GEOSX_UNUSED_VAR( domain );
}

/*
   APPLY BOUNDARY CONDITIONS
   Here, this call is the generic call from SolverBase.
   All it does is to call a specific Dirichlet boundary condition implemented for this solver
 */
void LaplaceBaseH1::applyBoundaryConditions( real64 const time_n,
                                             real64 const dt,
                                             DomainPartition & domain,
                                             DofManager const & dofManager,
                                             CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                             arrayView1d< real64 > const & localRhs )
{
  applyDirichletBCImplicit( time_n + dt, dofManager, domain, localMatrix, localRhs );
}

/*
   DIRICHLET BOUNDARY CONDITIONS
   This is the boundary condition method applied for this particular solver.
   It is called by the more generic "applyBoundaryConditions" method.
 */
void LaplaceBaseH1::
  applyDirichletBCImplicit( real64 const time,
                            DofManager const & dofManager,
                            DomainPartition & domain,
                            CRSMatrixView< real64, globalIndex const > const & localMatrix,
                            arrayView1d< real64 > const & localRhs )
{
  FieldSpecificationManager const & fsManager = FieldSpecificationManager::getInstance();

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
    bc.applyBoundaryConditionToSystem< FieldSpecificationEqual,
                                       parallelDevicePolicy< 32 > >( targetSet,
                                                                     time,
                                                                     targetGroup,
                                                                     m_fieldName,
                                                                     dofManager.getKey( m_fieldName ),
                                                                     dofManager.rankOffset(),
                                                                     localMatrix,
                                                                     localRhs );
  } );
}

/*
   SOLVE SYSTEM
   This method is simply initiating the solution and right-hand side
   and pass it to the base class solver.
 */
void LaplaceBaseH1::solveSystem( DofManager const & dofManager,
                                 ParallelMatrix & matrix,
                                 ParallelVector & rhs,
                                 ParallelVector & solution )
{
  rhs.scale( -1.0 ); // TODO decide if we want this here
  solution.zero();
  SolverBase::solveSystem( dofManager, matrix, rhs, solution );
}

void LaplaceBaseH1::resetStateToBeginningOfStep( DomainPartition & GEOSX_UNUSED_PARAM( domain ) )
{}

} // namespace geosx
