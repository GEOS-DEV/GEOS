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

namespace geosx
{

namespace dataRepository
{
namespace keys
{}
}

using namespace dataRepository;

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
               m_solution );

  return this->explicitStep( time_n, dt, cycleNumber, domain );
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
  return 1; // TODO: is this what we want?
}

//START_SPHINX_INCLUDE_REGISTER
REGISTER_CATALOG_ENTRY( SolverBase, MatrixFreeLaplaceFEM, string const &, Group * const )
//END_SPHINX_INCLUDE_REGISTER
} /* namespace geosx */
