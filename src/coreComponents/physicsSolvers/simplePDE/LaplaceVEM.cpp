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
 * @file LaplaceVEM.cpp
 */

// Source includes
#include "LaplaceVEM.hpp"

#include "common/TimingMacros.hpp"
#include "common/DataTypes.hpp"
#include "discretizationMethods/NumericalMethodsManager.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "mainInterface/GeosxState.hpp"
#include "virtualElement/ConformingVirtualElementOrder1.hpp"

namespace geosx
{

namespace dataRepository
{
namespace keys
{}
}

using namespace dataRepository;

/*----------------------------------------------------------------------------------
 * LaplaceVEM: Solving Laplace's partial differential equation with virtual elements
 * ---------------------------------------------------------------------------------
 *
 * What does this solver do?
 * --------------------------
 *
 * This solver finds a solution f(x,y,z) to the Laplace equation: div ( grad ( f )) = 0
 * This common elliptic PDE represents the solution of a steady-state heat transfer, for instance.
 *
 * Where can I find an example of what it does?
 * --------------------------------------------
 *
 * Integrated tests associated to this solver are found in the ./integratedTests/ folder
 * These tests consist of computing the steady-state temperature profile in a simple cube-shaped domain
 * with fixed temperatures applied on two opposite cube faces ("Dirichlet" boundary conditions: imposing a value).
 * Feel free to run these tests cases, check out the XML input files, and inspect the output.
 *
 * Implementation: before we start:
 * ---------------------------------
 * In this implementation, the solution function (called above f) is called m_fieldName.
 * The variable m_fieldName is a string that points to a data container (an array) that
 * holds the numerical values of the PDE solution for each location at which f is evaluated.
 *
 * Let's take a look at the implementation step by step.
 *
 * ---------------------------------------------------------------------------------
 */


/* CONSTRUCTOR
   First, let us inspect the constructor of a "LaplaceVEM" object.
   This constructor does three important things:
   1 - It constructs an instance of the LaplaceVEM class (here: using the SolverBase constructor and passing through the arguments).
   2 - It sets some default values for the LaplaceVEM-specific private variables (here: m_fieldName and m_timeIntegrationOption).
   3 - It creates and activates a "registerWrapper" for each private variable.
   This is where the private variables are declared either as REQUIRED or OPTIONAL.
   An error is thrown if a REQUIRED variable is not specified in the XML file,
   along with the description of this variable and possible enum values if relevant.
   The description that is set is used in auto-generated documentation and console error messages.
 */

//START_SPHINX_INCLUDE_01
LaplaceVEM::LaplaceVEM( const string & name,
                        Group * const parent ):
  SolverBase( name, parent ),
  m_fieldName( "primaryField" ),
  m_timeIntegrationOption( TimeIntegrationOption::ImplicitTransient )
{
  registerWrapper( laplaceFEMViewKeys.timeIntegrationOption.key(), &m_timeIntegrationOption ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Time integration method. Options are:\n* " + EnumStrings< TimeIntegrationOption >::concat( "\n* " ) );

  registerWrapper( laplaceFEMViewKeys.fieldVarName.key(), &m_fieldName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of field variable" );
}
//END_SPHINX_INCLUDE_01


// Destructor
LaplaceVEM::~LaplaceVEM()
{
  // TODO Auto-generated destructor stub
}

/* REGISTER THE PDE SOLUTION DATA ON THE MESH
   In the LaplaceVEM solver, we compute the solution of the partial differential equation "numerically".
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
void LaplaceVEM::registerDataOnMesh( Group & meshBodies )
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

real64 LaplaceVEM::solverStep( real64 const & time_n,
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
    dtReturn = this->linearImplicitStep( time_n, dt, cycleNumber, domain );
  }
  return dtReturn;
}

/*
   IMPLICIT STEP SETUP
   This method uses the system setup from LaplaceVEM (see below).
   It "deactivates" the time variables (with the GEOSX_UNUSED_PARAM macro) and does a steady state system set-up.
 */
void LaplaceVEM::implicitStepSetup( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                    real64 const & GEOSX_UNUSED_PARAM( dt ),
                                    DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  // Note: here we cannot use SolverBase::setupSystem, because it does:
  //       m_localMatrix.assimilate< parallelDevicePolicy<> >( std::move( pattern ) );
  //       and that creates problems (integratedTests failures) on Lassen for CPU-only implicit simulations

  m_dofManager.setMesh( domain, 0, 0 );

  setupDofs( domain, m_dofManager );
  m_dofManager.reorderByRank();

  localIndex const numLocalRows = m_dofManager.numLocalDofs();

  SparsityPattern< globalIndex > pattern;
  m_dofManager.setSparsityPattern( pattern );
  m_localMatrix.assimilate< serialPolicy >( std::move( pattern ) );

  m_localRhs.resize( numLocalRows );
  m_localSolution.resize( numLocalRows );

  m_localMatrix.setName( this->getName() + "/localMatrix" );
  m_localRhs.setName( this->getName() + "/localRhs" );
  m_localSolution.setName( this->getName() + "/localSolution" );
}

void LaplaceVEM::implicitStepComplete( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                       real64 const & GEOSX_UNUSED_PARAM( dt ),
                                       DomainPartition & GEOSX_UNUSED_PARAM( domain ) )
{}

void LaplaceVEM::setupDofs( DomainPartition const & GEOSX_UNUSED_PARAM( domain ),
                            DofManager & dofManager ) const
{
  dofManager.addField( m_fieldName,
                       DofManager::Location::Node );

  dofManager.addCoupling( m_fieldName,
                          m_fieldName,
                          DofManager::Connector::Elem );
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
   See the implementation in LaplaceFEMKernel.cpp.
 */
//START_SPHINX_INCLUDE_04
void LaplaceVEM::assembleSystem( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                 real64 const GEOSX_UNUSED_PARAM( dt ),
                                 DomainPartition & domain,
                                 DofManager const & dofManager,
                                 CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                 arrayView1d< real64 > const & localRhs )
{
  using VEM = virtualElement::ConformingVirtualElementOrder1< m_maxCellNodes, m_maxFaceNodes >;

  MeshLevel & mesh = domain.getMeshBodies().getGroup< MeshBody >( 0 ).getMeshLevel( 0 );
  NodeManager & nodeManager = mesh.getNodeManager();
  FaceManager const & faceManager = mesh.getFaceManager();
  EdgeManager const & edgeManager = mesh.getEdgeManager();

  // Get geometric properties used to compute projectors.
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > nodesCoords =
    nodeManager.referencePosition();
  FaceManager::NodeMapType const & faceToNodeMap = faceManager.nodeList();
  FaceManager::EdgeMapType const & faceToEdgeMap = faceManager.edgeList();
  EdgeManager::NodeMapType const & edgeToNodeMap = edgeManager.nodeList();
  arrayView2d< real64 const > const faceCenters = faceManager.faceCenter();
  arrayView2d< real64 const > const faceNormals = faceManager.faceNormal();
  arrayView1d< real64 const > const faceAreas = faceManager.faceArea();
  arrayView1d< globalIndex const > const & dofIndex =
    nodeManager.getReference< array1d< globalIndex > >( dofManager.getKey( m_fieldName ) );

  real64 const diffusion = 1.0;
  globalIndex const rankOffset = dofManager.rankOffset();

  forTargetRegionsComplete( mesh, [&]( localIndex const,
                                       localIndex const,
                                       ElementRegionBase & elementRegion )
  {
    elementRegion.forElementSubRegions< CellElementSubRegion >
      ( [&]( CellElementSubRegion const & elemSubRegion )
    {
      CellBlock::NodeMapType const & elemToNodeMap = elemSubRegion.nodeList();
      CellBlock::FaceMapType const & elementToFaceMap = elemSubRegion.faceList();
      arrayView2d< real64 const > elemCenters = elemSubRegion.getElementCenter();
      arrayView1d< real64 const > elemVolumes = elemSubRegion.getElementVolume();
      arrayView1d< integer const > const & elemGhostRank = elemSubRegion.ghostRank();
      localIndex const numCells = elemSubRegion.size();
      forAll< serialPolicy >( numCells, [=] ( localIndex const cellIndex )
      {
        real64 basisDerivativesIntegralMean[VEM::maxSupportPoints][3];
        globalIndex elemDofIndex[VEM::maxSupportPoints];
        real64 element_matrix[VEM::maxSupportPoints][VEM::maxSupportPoints];

        if( elemGhostRank[cellIndex] < 0 )
        {
          VEM::computeProjectors( cellIndex, nodesCoords, elemToNodeMap, elementToFaceMap,
                                  faceToNodeMap, faceToEdgeMap, edgeToNodeMap,
                                  faceCenters, faceNormals, faceAreas,
                                  elemCenters[cellIndex], elemVolumes[cellIndex] );
          localIndex const numSupportPoints = VEM::getNumSupportPoints();
          for( localIndex a = 0; a < numSupportPoints; ++a )
          {
            elemDofIndex[a] = dofIndex[ elemToNodeMap( cellIndex, a ) ];
            for( localIndex b = 0; b < numSupportPoints; ++b )
            {
              element_matrix[a][b] = VEM::calcStabilizationValue( a, b );
            }
            localRhs[a] = 0.0;
          }
          for( localIndex q = 0; q < VEM::getNumQuadraturePoints(); ++q )
          {
            VEM::calcGradN( q, basisDerivativesIntegralMean );
            for( localIndex a = 0; a < numSupportPoints; ++a )
            {
              for( localIndex b = 0; b < numSupportPoints; ++b )
              {
                element_matrix[a][b] += diffusion * VEM::transformedQuadratureWeight( q ) *
                                        (basisDerivativesIntegralMean[a][0] * basisDerivativesIntegralMean[b][0] +
                                         basisDerivativesIntegralMean[a][1] * basisDerivativesIntegralMean[b][1] +
                                         basisDerivativesIntegralMean[a][2] * basisDerivativesIntegralMean[b][2] );
              }
            }
          }
          for( localIndex a = 0; a < numSupportPoints; ++a )
          {
            localIndex const dof = LvArray::integerConversion< localIndex >
                                     ( elemDofIndex[a] - rankOffset );
            if( dof < 0 || dof >= localMatrix.numRows() )
              continue;
            localMatrix.template addToRowBinarySearchUnsorted< serialAtomic >
              ( dof, elemDofIndex, element_matrix[ a ], numSupportPoints );
          }
        }
      } );
    } );
  } );

  //END_SPHINX_INCLUDE_04
}

void LaplaceVEM::applySystemSolution( DofManager const & dofManager,
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
void LaplaceVEM::applyBoundaryConditions( real64 const time_n,
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
void LaplaceVEM::solveSystem( DofManager const & dofManager,
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
void LaplaceVEM::applyDirichletBCImplicit( real64 const time,
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
    bc.applyBoundaryConditionToSystem< FieldSpecificationEqual, serialPolicy >( targetSet,
                                                                                time,
                                                                                targetGroup,
                                                                                m_fieldName,
                                                                                dofManager.getKey( m_fieldName ),
                                                                                dofManager.rankOffset(),
                                                                                localMatrix,
                                                                                localRhs );
  } );
}

void LaplaceVEM::resetStateToBeginningOfStep( DomainPartition & GEOSX_UNUSED_PARAM( domain ) )
{}

//START_SPHINX_INCLUDE_00
REGISTER_CATALOG_ENTRY( SolverBase, LaplaceVEM, string const &, Group * const )
//END_SPHINX_INCLUDE_00
} /* namespace ANST */
