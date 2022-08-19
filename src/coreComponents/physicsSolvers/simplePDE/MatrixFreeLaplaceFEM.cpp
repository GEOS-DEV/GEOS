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
 * @file MatrixFreeLaplaceFEM.cpp
 */

// Source includes
#include "MatrixFreeLaplaceFEM.hpp"
#include "LaplaceFEMKernels.hpp"
#include "TeamLaplaceFEMKernels.hpp"
#include "finiteElement/kernelInterface/KernelBase.hpp"
#include "mesh/MeshBody.hpp"
#include "linearAlgebra/solvers/CgSolver.hpp"
#include "linearAlgebra/solvers/PreconditionerIdentity.hpp"

namespace geosx
{

namespace dataRepository
{
namespace keys
{}
}

using namespace dataRepository;

MatrixFreeLaplaceFEMOperator::
  MatrixFreeLaplaceFEMOperator( DomainPartition & domain,
                                map< string, array1d< string > > & meshTargets,
                                DofManager & dofManager,
                                string const & finiteElementName ):
    m_meshBodies( domain.getMeshBodies() ),
    m_meshTargets( meshTargets ),
    m_dofManager( dofManager ),
    m_finiteElementName( finiteElementName )
{ }

MatrixFreeLaplaceFEMOperator::
  MatrixFreeLaplaceFEMOperator( dataRepository::Group & meshBodies,
                                map< string, array1d< string > > & meshTargets,
                                DofManager & dofManager,
                                string const & finiteElementName ):
    m_meshBodies( meshBodies ),
    m_meshTargets( meshTargets ),
    m_dofManager( dofManager ),
    m_finiteElementName( finiteElementName )
{ }

void MatrixFreeLaplaceFEMOperator::apply( ParallelVector const & src, ParallelVector & dst ) const
{
  dst.zero();
  arrayView1d< real64 const > const localSrc = src.values();
  arrayView1d< real64 > const localDst = dst.open();
  for( auto const & target: m_meshTargets )
  {
    string const meshBodyName = target.first;
    arrayView1d< string const > const & regionNames = target.second.toViewConst();
    MeshBody & meshBody = m_meshBodies.getGroup< MeshBody >( meshBodyName );
    meshBody.forMeshLevels( [&]( MeshLevel & mesh )
    {

      TeamLaplaceFEMKernelFactory2 kernelFactory( localSrc, localDst );

      string const dummyString = "dummy";
      finiteElement::
        regionBasedKernelApplication< team_launch_policy,
                                      constitutive::NullModel,
                                      CellElementSubRegion >( mesh,
                                                              regionNames,
                                                              m_finiteElementName,
                                                              dummyString,
                                                              kernelFactory );

    } );
  }
  dst.close();
}

globalIndex MatrixFreeLaplaceFEMOperator::numGlobalRows() const
{
  return m_dofManager.numGlobalDofs();
}

globalIndex MatrixFreeLaplaceFEMOperator::numGlobalCols() const
{
  return m_dofManager.numGlobalDofs();
}

localIndex MatrixFreeLaplaceFEMOperator::numLocalRows() const
{
  return m_dofManager.numLocalDofs();
}

localIndex MatrixFreeLaplaceFEMOperator::numLocalCols() const
{
  return m_dofManager.numLocalDofs();
}

MPI_Comm MatrixFreeLaplaceFEMOperator::comm() const
{
  return MPI_COMM_GEOSX;
}

MatrixFreePreconditionerIdentity::MatrixFreePreconditionerIdentity( DofManager & dofManager )
: m_dofManager( dofManager )
{ }

globalIndex MatrixFreePreconditionerIdentity::numGlobalRows() const
{
  return m_dofManager.numGlobalDofs();
}

globalIndex MatrixFreePreconditionerIdentity::numGlobalCols() const
{
  return m_dofManager.numGlobalDofs();
}

localIndex MatrixFreePreconditionerIdentity::numLocalRows() const
{
  return m_dofManager.numLocalDofs();
}

localIndex MatrixFreePreconditionerIdentity::numLocalCols() const
{
  return m_dofManager.numLocalDofs();
}

MPI_Comm MatrixFreePreconditionerIdentity::comm() const
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
   First, let us inspect the constructor of a "MatrixFreeLaplaceFEM" object.
   This constructor does three important things:
   1 - It constructs an instance of the MatrixFreeLaplaceFEM class (here: using the SolverBase constructor and passing through the arguments).
   2 - It sets some default values for the MatrixFreeLaplaceFEM-specific private variables (here: m_fieldName and m_timeIntegrationOption).
   3 - It creates and activates a "registerWrapper" for each private variable.
   This is where the private variables are declared either as REQUIRED or OPTIONAL.
   An error is thrown if a REQUIRED variable is not specified in the XML file,
   along with the description of this variable and possible enum values if relevant.
   The description that is set is used in auto-generated documentation and console error messages.
 */

//START_SPHINX_INCLUDE_CONSTRUCTOR
MatrixFreeLaplaceFEM::MatrixFreeLaplaceFEM( const string & name,
                        Group * const parent ):
  LaplaceBaseH1( name, parent )
{}
//END_SPHINX_INCLUDE_CONSTRUCTOR

MatrixFreeLaplaceFEM::~MatrixFreeLaplaceFEM()
{
  // TODO Auto-generated destructor stub
}

real64 MatrixFreeLaplaceFEM::solverStep( real64 const & time_n,
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

  MatrixFreeLaplaceFEMOperator unconstrained_laplace( domain, m_meshTargets, m_dofManager, this->getDiscretizationName() );
  LinearOperatorWithBC< ParallelVector > constrained_laplace( *this,
                                                              unconstrained_laplace,
                                                              domain,
                                                              m_dofManager,
                                                              m_fieldName,
                                                              time_n+dt,
                                                              LinearOperatorWithBC< ParallelVector >::
                                                                DiagPolicy::
                                                                  DiagonalOne );
  constrained_laplace.computeConstrainedRHS( m_rhs );
  MatrixFreePreconditionerIdentity identity( m_dofManager );
  auto & params = m_linearSolverParameters.get();
  params.isSymmetric = true;
  CgSolver< ParallelVector > solver( params, constrained_laplace, identity );
  solver.solve( m_solution, m_rhs );
  return dt;
}

real64 MatrixFreeLaplaceFEM::explicitStep( real64 const & time_n,
                                           real64 const & dt,
                                           const int GEOSX_UNUSED_PARAM( cycleNumber ),
                                           DomainPartition & domain )
{
  m_rhs.zero();
  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel & mesh,
                                                arrayView1d< string const > const & regionNames )
  {

    arrayView1d< real64 > const localRhs = m_rhs.open();
  
    TeamLaplaceFEMKernelFactory kernelFactory( localRhs, m_fieldName );

    string const dummyString = "dummy";
    finiteElement::
      regionBasedKernelApplication< team_launch_policy,
                                    constitutive::NullModel,
                                    CellElementSubRegion >( mesh,
                                                            regionNames,
                                                            this->getDiscretizationName(),
                                                            dummyString,
                                                            kernelFactory );

    m_rhs.close();
  } );
  return dt;
}

//START_SPHINX_INCLUDE_REGISTER
REGISTER_CATALOG_ENTRY( SolverBase, MatrixFreeLaplaceFEM, string const &, Group * const )
//END_SPHINX_INCLUDE_REGISTER
} /* namespace geosx */
