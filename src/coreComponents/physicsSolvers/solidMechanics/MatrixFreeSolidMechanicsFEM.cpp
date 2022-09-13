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
 * @file MatrixFreeSolidMechanicsFEM.cpp
 */

// Source includes
#include "MatrixFreeSolidMechanicsFEM.hpp"
#include "TeamSolidMechanicsFEMKernels.hpp"
#include "finiteElement/kernelInterface/KernelBase.hpp"
#include "mesh/MeshBody.hpp"
#include "linearAlgebra/solvers/CgSolver.hpp"
#include "linearAlgebra/solvers/PreconditionerIdentity.hpp"
#include "linearAlgebra/common/LinearOperatorWithBC.hpp"
#include "constitutive/solid/ElasticIsotropic.hpp"

namespace geosx
{

namespace dataRepository
{
namespace keys
{}
}

using namespace dataRepository;
using namespace constitutive;

MatrixFreeSolidMechanicsFEMOperator::
  MatrixFreeSolidMechanicsFEMOperator( DomainPartition & domain,
                                       map< std::pair< string, string >, array1d< string > > const & meshTargets,
                                       DofManager & dofManager,
                                       string const & finiteElementName ):
    m_meshBodies( domain.getMeshBodies() ),
    m_meshTargets( meshTargets ),
    m_dofManager( dofManager ),
    m_finiteElementName( finiteElementName )
{ }

MatrixFreeSolidMechanicsFEMOperator::
  MatrixFreeSolidMechanicsFEMOperator( dataRepository::Group & meshBodies,
                                       map< std::pair< string, string >, array1d< string > > const & meshTargets,
                                       DofManager & dofManager,
                                       string const & finiteElementName ):
    m_meshBodies( meshBodies ),
    m_meshTargets( meshTargets ),
    m_dofManager( dofManager ),
    m_finiteElementName( finiteElementName )
{ }

void MatrixFreeSolidMechanicsFEMOperator::apply( ParallelVector const & src, ParallelVector & dst ) const
{
  dst.zero();
  arrayView1d< real64 const > const localSrc = src.values();
  arrayView1d< real64 > const localDst = dst.open();

  for( auto const & target: m_meshTargets )
  {
    string const meshBodyName = target.first.first;
    string const meshLevelName = target.first.second;
    arrayView1d< string const > const & regionNames = target.second.toViewConst();
    MeshBody & meshBody = m_meshBodies.getGroup< MeshBody >( meshBodyName );

    MeshLevel * meshLevelPtr = meshBody.getMeshLevels().getGroupPointer< MeshLevel >( meshLevelName );
    if( meshLevelPtr==nullptr )
    {
      meshLevelPtr = meshBody.getMeshLevels().getGroupPointer< MeshLevel >( MeshBody::groupStructKeys::baseDiscretizationString() );
    }
    MeshLevel & mesh = *meshLevelPtr;
      
    auto const & totalDisplacement = mesh.getNodeManager().totalDisplacement();
    arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > localSrc2d( totalDisplacement.dimsArray(), totalDisplacement.stridesArray(), 0, localSrc.dataBuffer() );
    arrayView2d< real64, nodes::TOTAL_DISPLACEMENT_USD > localDst2d( totalDisplacement.dimsArray(), totalDisplacement.stridesArray(), 0, localDst.dataBuffer() );
    TeamSolidMechanicsFEMKernelFactory kernelFactory( localSrc2d, localDst2d );
  
    finiteElement::
      regionBasedKernelApplication< team_launch_policy,
                                    constitutive::SolidBase,
                                    CellElementSubRegion >( mesh,
                                                            regionNames,
                                                            m_finiteElementName,
                                                            "solidMaterialNames",
                                                            kernelFactory );
  }
  dst.close();
}

void MatrixFreeSolidMechanicsFEMOperator::computeDiagonal( ParallelVector & diagonal ) const
{
  GEOSX_ERROR( "computeDiagonal: operation not yet implemented" );
}

globalIndex MatrixFreeSolidMechanicsFEMOperator::numGlobalRows() const
{
  return m_dofManager.numGlobalDofs();
}

globalIndex MatrixFreeSolidMechanicsFEMOperator::numGlobalCols() const
{
  return m_dofManager.numGlobalDofs();
}

localIndex MatrixFreeSolidMechanicsFEMOperator::numLocalRows() const
{
  return m_dofManager.numLocalDofs();
}

localIndex MatrixFreeSolidMechanicsFEMOperator::numLocalCols() const
{
  return m_dofManager.numLocalDofs();
}

MPI_Comm MatrixFreeSolidMechanicsFEMOperator::comm() const
{
  return MPI_COMM_GEOSX;
}

/*----------------------------------------------------------------------------------
 * LaplaceFEM: Solving Laplace's partial differential equation with finite elements
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
   First, let us inspect the constructor of a "MatrixFreeSolidMechanicsFEM" object.
   This constructor does three important things:
   1 - It constructs an instance of the MatrixFreeSolidMechanicsFEM class (here: using the SolverBase constructor and passing through the arguments).
   2 - It sets some default values for the MatrixFreeSolidMechanicsFEM-specific private variables (here: m_fieldName and m_timeIntegrationOption).
   3 - It creates and activates a "registerWrapper" for each private variable.
   This is where the private variables are declared either as REQUIRED or OPTIONAL.
   An error is thrown if a REQUIRED variable is not specified in the XML file,
   along with the description of this variable and possible enum values if relevant.
   The description that is set is used in auto-generated documentation and console error messages.
 */

//START_SPHINX_INCLUDE_CONSTRUCTOR
MatrixFreeSolidMechanicsFEM::MatrixFreeSolidMechanicsFEM( const string & name,
                        Group * const parent ):
  SolverBase( name, parent ),
  m_fieldName( "TotalDisplacement" )
{}
//END_SPHINX_INCLUDE_CONSTRUCTOR

MatrixFreeSolidMechanicsFEM::~MatrixFreeSolidMechanicsFEM()
{
  // TODO Auto-generated destructor stub
}

real64 MatrixFreeSolidMechanicsFEM::solverStep( real64 const & time_n,
                                                real64 const & dt,
                                                const int cycleNumber,
                                                DomainPartition & domain )
{
  setupSystem( domain,
               m_dofManager,
               m_localMatrix,
               m_rhs,
               m_solution,
               false );

  MatrixFreeSolidMechanicsFEMOperator unconstrained_solid_mechanics(
    domain,
    getMeshTargets(),
    m_dofManager,
    this->getDiscretizationName() );

  LinearOperatorWithBC< ParallelVector, FieldType > constrained_solid_mechanics(
    *this,
    unconstrained_solid_mechanics,
    domain,
    m_dofManager,
    m_fieldName,
    time_n+dt,
    LinearOperatorWithBC< ParallelVector, FieldType >::DiagPolicy::DiagonalOne );

  constrained_solid_mechanics.computeConstrainedRHS( m_rhs, m_solution );

  std::cout<< "rhs: " << m_rhs << std::endl;

  MatrixFreePreconditionerIdentity< HypreInterface > identity( m_dofManager );

  auto & params = m_linearSolverParameters.get();
  params.isSymmetric = true;

  CgSolver< ParallelVector > solver( params, constrained_solid_mechanics, identity );
  
  solver.solve( m_rhs, m_solution );

  std::cout << "m_solution: " << m_solution << std::endl;

  applySystemSolution( m_dofManager, m_solution.values(), 1.0, domain );

  return dt;
}

void MatrixFreeSolidMechanicsFEM::setupDofs( DomainPartition const & GEOSX_UNUSED_PARAM( domain ),
                                             DofManager & dofManager ) const
{
  GEOSX_MARK_FUNCTION;
  dofManager.addField( keys::TotalDisplacement,
                       FieldLocation::Node,
                       3,
                       getMeshTargets() );

  dofManager.addCoupling( keys::TotalDisplacement,
                          keys::TotalDisplacement,
                          DofManager::Connector::Elem );
}

void MatrixFreeSolidMechanicsFEM::registerDataOnMesh( Group & meshBodies )
{
  SolverBase::registerDataOnMesh( meshBodies );

  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & meshLevel,
                                                    arrayView1d< string const > const & regionNames )
  {
    NodeManager & nodes = meshLevel.getNodeManager();

    nodes.registerWrapper< array2d< real64, nodes::TOTAL_DISPLACEMENT_PERM > >( keys::TotalDisplacement ).
      setPlotLevel( PlotLevel::LEVEL_0 ).
      setRegisteringObjects( this->getName()).
      setDescription( "An array that holds the total displacements on the nodes." ).
      reference().resizeDimension< 1 >( 3 );
  });

}

void MatrixFreeSolidMechanicsFEM::setConstitutiveNamesCallSuper( ElementSubRegionBase & subRegion ) const
{
  SolverBase::setConstitutiveNamesCallSuper( subRegion );

  subRegion.registerWrapper< string >( "solidMaterialNames" ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setSizedFromParent( 0 );

  string & solidMaterialName = subRegion.getReference< string >( "solidMaterialNames" );
  solidMaterialName = SolverBase::getConstitutiveName< SolidBase >( subRegion );
  GEOSX_ERROR_IF( solidMaterialName.empty(), GEOSX_FMT( "SolidBase model not found on subregion {}", subRegion.getName() ) );
}


void
MatrixFreeSolidMechanicsFEM::applySystemSolution( DofManager const & dofManager,
                                                  arrayView1d< real64 const > const & localSolution,
                                                  real64 const scalingFactor,
                                                  DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;
  dofManager.addVectorToField( localSolution,
                               keys::TotalDisplacement,
                               keys::TotalDisplacement,
                               scalingFactor );
                               

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {
    auto const & disp = mesh.getNodeManager().totalDisplacement();
    std::cout<<disp<<std::endl;

  } );

  // forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
  //                                                               MeshLevel & mesh,
  //                                                               arrayView1d< string const > const & )

  // {
  //   FieldIdentifiers fieldsToBeSync;

  //   fieldsToBeSync.addFields( FieldLocation::Node, { keys::TotalDisplacement } );

  //   CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync,
  //                                                        mesh,
  //                                                        domain.getNeighbors(),
  //                                                        true );
  // } );
}



//START_SPHINX_INCLUDE_REGISTER
REGISTER_CATALOG_ENTRY( SolverBase, MatrixFreeSolidMechanicsFEM, string const &, Group * const )
//END_SPHINX_INCLUDE_REGISTER
} /* namespace geosx */
