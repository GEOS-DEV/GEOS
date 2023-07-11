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
#define SELECTED_FE_TYPES H1_Hexahedron_Lagrange1_GaussLegendre2

// Source includes
#include "MatrixFreeSolidMechanicsFEM.hpp"
#include "kernels/SmallStrainResidual.hpp"
#include "finiteElement/kernelInterface/KernelBase.hpp"
#include "mesh/MeshBody.hpp"
#include "linearAlgebra/solvers/CgSolver.hpp"
#include "linearAlgebra/solvers/PreconditionerIdentity.hpp"
#include "linearAlgebra/common/LinearOperatorWithBC.hpp"
#include "constitutive/solid/ElasticIsotropic.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsFields.hpp"

#include "SolidMechanicsFields.hpp"
#include "LvArray/src/output.hpp"

namespace geos
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
  GEOS_MARK_FUNCTION;

  arrayView1d< real64 const > const localSrc = src.values();
  arrayView1d< real64 > const localDst = dst.open();
  // We do it by hand to avoid hypre call
  using POLICY = parallelDeviceAsyncPolicy<>;
  forAll< POLICY >( localDst.size(), [localDst] GEOS_HOST_DEVICE ( localIndex const i )
  {
    localDst[ i ] = 0.0;
  } );

    // {
    // std::cout<<"MatrixFreeSolidMechanicsFEMOperator::apply - bp1"<<std::endl;
    // LvArray::print< parallelDevicePolicy< 32 > >( localSrc );
    // std::cout<<"MatrixFreeSolidMechanicsFEMOperator::apply - bp2"<<std::endl;
    // }



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
      
    auto const & totalDisplacement = mesh.getNodeManager().getField<fields::solidMechanics::totalDisplacement>();
    arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > localSrc2d( totalDisplacement.dimsArray(), totalDisplacement.stridesArray(), 0, localSrc.dataBuffer() );
    arrayView2d< real64, nodes::TOTAL_DISPLACEMENT_USD > localDst2d( totalDisplacement.dimsArray(), totalDisplacement.stridesArray(), 0, localDst.dataBuffer() );

#if 0
    TeamSolidMechanicsFEMKernelFactory kernelFactory( localSrc2d, localDst2d );  
    finiteElement::
      regionBasedKernelApplication< team_launch_policy,
                                    constitutive::SolidBase,
                                    CellElementSubRegion >( mesh,
                                                            regionNames,
                                                            m_finiteElementName,
                                                            "solidMaterialNames",
                                                            kernelFactory );
#else
    auto kernelFactory = solidMechanicsLagrangianFEMKernels::SmallStrainResidualFactory( localSrc2d, 
                                                                                         localDst2d,
                                                                                         0, 
                                                                                         "" );

    finiteElement::
    regionBasedKernelApplication< parallelDevicePolicy< 32 >,
                                  constitutive::ElasticIsotropic,
                                  CellElementSubRegion >( mesh,
                                                          regionNames,
                                                          m_finiteElementName,
                                                          "solidMaterialNames",
                                                          kernelFactory );
#endif
  }

  // {
  //   std::cout<<"MatrixFreeSolidMechanicsFEMOperator::apply - bp3"<<std::endl;
  //   LvArray::print< parallelDevicePolicy< 32 > >( localDst );
  //   std::cout<<"MatrixFreeSolidMechanicsFEMOperator::apply - bp4"<<std::endl;
  // }

  dst.close();
}

void MatrixFreeSolidMechanicsFEMOperator::computeDiagonal( ParallelVector & GEOS_UNUSED_PARAM(diagonal) ) const
{
  GEOS_ERROR( "computeDiagonal: operation not yet implemented" );
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
  m_fieldName( "totalDisplacement" )
{}
//END_SPHINX_INCLUDE_CONSTRUCTOR

MatrixFreeSolidMechanicsFEM::~MatrixFreeSolidMechanicsFEM()
{
  // TODO Auto-generated destructor stub
}

real64 MatrixFreeSolidMechanicsFEM::solverStep( real64 const & time_n,
                                                real64 const & dt,
                                                const int GEOS_UNUSED_PARAM(cycleNumber),
                                                DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

//  std::cout<<"MatrixFreeSolidMechanicsFEM::solverStep - begin"<<std::endl;
  m_dofManager.setDomain( domain );
  setupDofs( domain, m_dofManager );
  m_dofManager.reorderByRank();
  m_rhs.setName( this->getName() + "/rhs" );
  m_rhs.create( m_dofManager.numLocalDofs(), MPI_COMM_GEOSX );
  m_solution.setName( this->getName() + "/solution" );
  m_solution.create( m_dofManager.numLocalDofs(), MPI_COMM_GEOSX );

//  std::cout<<"     MatrixFreeSolidMechanicsFEM::solverStep - bp1"<<std::endl;

  MatrixFreeSolidMechanicsFEMOperator unconstrained_solid_mechanics(
    domain,
    getMeshTargets(),
    m_dofManager,
    this->getDiscretizationName() );

//  std::cout<<"     MatrixFreeSolidMechanicsFEM::solverStep - bp2"<<std::endl;

  LinearOperatorWithBC< ParallelVector, FieldType > constrained_solid_mechanics(
    *this,
    unconstrained_solid_mechanics,
    domain,
    m_dofManager,
    m_fieldName,
    time_n+dt,
    LinearOperatorWithBC< ParallelVector, FieldType >::DiagPolicy::DiagonalOne );

//  std::cout<<"     MatrixFreeSolidMechanicsFEM::solverStep - bp3"<<std::endl;

// std::cout<< "m_rhs0: "<<std::endl<< m_rhs << std::endl;

  constrained_solid_mechanics.computeConstrainedRHS( m_rhs, m_solution );

//  std::cout<< "m_rhs: "<<std::endl<< m_rhs << std::endl;
//  std::cout<< "solution: " <<std::endl<< m_solution << std::endl;

//   std::cout<<"     MatrixFreeSolidMechanicsFEM::solverStep - bp4"<<std::endl;


  // MatrixFreePreconditionerIdentity< HypreInterface > identity( m_dofManager );

// std::cout<<"     MatrixFreeSolidMechanicsFEM::solverStep - bp5"<<std::endl;

  auto & params = m_linearSolverParameters.get();
  params.isSymmetric = true;

// std::cout<<"     MatrixFreeSolidMechanicsFEM::solverStep - bp6"<<std::endl;

  // CgSolver< ParallelVector > solver( params, constrained_solid_mechanics, identity );
  UnprecCgSolver< ParallelVector > solver( params, constrained_solid_mechanics );
  
//   std::cout<<"     MatrixFreeSolidMechanicsFEM::solverStep - bp7"<<std::endl;

  solver.solve( m_rhs, m_solution );
  const real64 elapsed_seconds = solver.result().solveTime;
  std::cout << "solve time: " << elapsed_seconds << "s\n";
  const integer num_iter = solver.result().numIterations;
  const double time_per_iter = elapsed_seconds / num_iter;
  std::cout << "Time per CG iteration: " << time_per_iter << "s\n";
  const size_t num_dofs = m_dofManager.numLocalDofs();
  std::cout << "Number of local dofs: " << num_dofs << "dofs\n";
  std::cout << "Throughput: " << num_dofs/time_per_iter*1e-6 << "MDofs/s\n";

// std::cout<<"     MatrixFreeSolidMechanicsFEM::solverStep - bp8"<<std::endl;

//  std::cout << "m_solution: " << m_solution << std::endl;

  applySystemSolution( m_dofManager, m_solution.values(), 1.0, domain );

//  std::cout<<"MatrixFreeSolidMechanicsFEM::solverStep - end"<<std::endl;

  return dt;
}

void MatrixFreeSolidMechanicsFEM::setupDofs( DomainPartition const & GEOS_UNUSED_PARAM( domain ),
                                             DofManager & dofManager ) const
{
  GEOS_MARK_FUNCTION;
//  std::cout<<"MatrixFreeSolidMechanicsFEM::setupDofs - begin"<<std::endl;

  dofManager.addField( fields::solidMechanics::totalDisplacement::key(),
                       FieldLocation::Node,
                       3,
                       getMeshTargets() );

  dofManager.addCoupling( fields::solidMechanics::totalDisplacement::key(),
                          fields::solidMechanics::totalDisplacement::key(),
                          DofManager::Connector::Elem );

//  std::cout<<"MatrixFreeSolidMechanicsFEM::setupDofs - end"<<std::endl;

}

void MatrixFreeSolidMechanicsFEM::registerDataOnMesh( Group & meshBodies )
{
  SolverBase::registerDataOnMesh( meshBodies );

  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & meshLevel,
                                                    arrayView1d< string const > const & GEOS_UNUSED_PARAM(regionNames) )
  {
    NodeManager & nodes = meshLevel.getNodeManager();

    nodes.registerField< fields::solidMechanics::totalDisplacement >( getName() ).
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
  GEOS_ERROR_IF( solidMaterialName.empty(), GEOS_FMT( "SolidBase model not found on subregion {}", subRegion.getName() ) );
}


void
MatrixFreeSolidMechanicsFEM::applySystemSolution( DofManager const & dofManager,
                                                  arrayView1d< real64 const > const & localSolution,
                                                  real64 const scalingFactor,
                                                  DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;
//  std::cout<<"MatrixFreeSolidMechanicsFEM::applySystemSolution - begin"<<std::endl;
  dofManager.addVectorToField( localSolution,
                               fields::solidMechanics::totalDisplacement::key(),
                               fields::solidMechanics::totalDisplacement::key(),
                               scalingFactor );
                               

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & GEOS_UNUSED_PARAM(mesh),
                                                                arrayView1d< string const > const & )
  {
    // auto const & disp = mesh.getNodeManager().getField<fields::solidMechanics::totalDisplacement>();
//    std::cout<<disp<<std::endl;

  } );

  // forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
  //                                                               MeshLevel & mesh,
  //                                                               arrayView1d< string const > const & )

  // {
  //   FieldIdentifiers fieldsToBeSync;

  //   fieldsToBeSync.addFields( FieldLocation::Node, { fields::solidMechanics::totalDisplacement } )::key;

  //   CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync,
  //                                                        mesh,
  //                                                        domain.getNeighbors(),
  //                                                        true );
  // } );
//    std::cout<<"MatrixFreeSolidMechanicsFEM::applySystemSolution - end"<<std::endl;

}



//START_SPHINX_INCLUDE_REGISTER
REGISTER_CATALOG_ENTRY( SolverBase, MatrixFreeSolidMechanicsFEM, string const &, Group * const )
//END_SPHINX_INCLUDE_REGISTER
} /* namespace geos */
#undef SELECTED_FE_TYPES
