/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistrubute it and/or modify it under
 * the terms of the GNU Lesser General Public Liscense (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

#include "SolverBase.hpp"
#include "managers/DomainPartition.hpp"
#include "mesh/MeshBody.hpp"
#include "systemSolverInterface/LinearSolverWrapper.hpp"
#include "systemSolverInterface/EpetraBlockSystem.hpp"
#include "systemSolverInterface/SystemSolverParameters.hpp"

namespace geosx
{

using namespace dataRepository;

SolverBase::SolverBase( std::string const & name,
                        ManagedGroup * const parent ):
  ManagedGroup( name, parent ),
  m_linearSolverWrapper(nullptr),
  m_linearSystem(nullptr),
  m_verboseLevel(0),
  m_gravityVector( R1Tensor(0.0) )
{
  // register group with repository. Have Repository own object.
  this->RegisterGroup<SystemSolverParameters>( groupKeys.systemSolverParameters.Key() );

  this->RegisterViewWrapper( viewKeyStruct::verboseLevelString, &m_verboseLevel, 0 );
  this->RegisterViewWrapper( viewKeyStruct::gravityVectorString, &m_gravityVector, 0 );

  if( this->globalGravityVector() != nullptr )
  {
    m_gravityVector=*globalGravityVector();
  }

  m_linearSolverWrapper = new systemSolverInterface::LinearSolverWrapper();
  m_linearSystem = new systemSolverInterface::EpetraBlockSystem();
}

SolverBase::~SolverBase()
{
  delete m_linearSolverWrapper;
  delete m_linearSystem;
}

SolverBase::CatalogInterface::CatalogType& SolverBase::GetCatalog()
{
  static SolverBase::CatalogInterface::CatalogType catalog;
  return catalog;
}

void SolverBase::FillDocumentationNode()
{


  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();
  docNode->setName(this->CatalogName());    // If this method lived in Managed
                                            // groups, this could be done
                                            // automatically
  docNode->setSchemaType("Node");

  docNode->AllocateChildNode( keys::courant,
                              keys::courant,
                              -1,
                              "real64",
                              "real64",
                              "courant Number",
                              "courant Number",
                              "0.7",
                              "",
                              1,
                              1,
                              0 );

  docNode->AllocateChildNode( keys::maxDt,
                              keys::maxDt,
                              -1,
                              "real64",
                              "real64",
                              "Maximum Stable Timestep",
                              "Maximum Stable Timestep",
                              "0.0",
                              "",
                              0,
                              1,
                              0 );

  docNode->AllocateChildNode( viewKeyStruct::verboseLevelString,
                              viewKeyStruct::verboseLevelString,
                              -1,
                              "integer",
                              "integer",
                              "verbosity level",
                              "verbosity level",
                              "0",
                              "",
                              0,
                              1,
                              0 );


}



void SolverBase::SolverStep( real64 const& time_n,
                           real64 const& dt,
                           const int cycleNumber,
                           ManagedGroup * domain )
{
}


real64 SolverBase::LinearImplicitStep( real64 const & time_n,
                                       real64 const & dt,
                                       integer const cycleNumber,
                                       DomainPartition * const domain )
{
  real64 dt_return = dt;

  // call setup for physics solver. Pre step allocations etc.
  ImplicitStepSetup( time_n, dt, domain );

  // call assemble to fill the matrix and the rhs
  AssembleSystem( domain, m_linearSystem, time_n+dt, dt );

  // call the default linear solver on the system
  SolveSystem( m_linearSystem,
               getSystemSolverParameters() );

  // apply the system solution to the fields/variables
  ApplySystemSolution( m_linearSystem, 1.0, 0, domain );

  // final step for completion of timestep. typically secondary variable updates and cleanup.
  ImplicitStepComplete( time_n, dt,  domain );

  // return the achieved timestep
  return dt_return;
}

real64 SolverBase::NonlinearImplicitStep( real64 const & time_n,
                                           real64 const & dt,
                                           integer const cycleNumber,
                                           DomainPartition * const domain )
{
  real64 stepDt = dt;

  // call setup for physics solver. Pre step allocations etc.
  ImplicitStepSetup( time_n, stepDt, domain );

  SystemSolverParameters const * const solverParams = getSystemSolverParameters();
  real64 const newtonTol = solverParams->newtonTol();
  integer const maxNewtonIter = solverParams->maxIterNewton();
  integer const maxNumberDtCuts = 2;
  integer const maxNumberLineSearchCuts = 5;

  integer isConverged = 0;
  for( int dtAttempt = 0 ; dtAttempt<maxNumberDtCuts ; ++dtAttempt )
  {
    real64 lastResidual = 1e99;
    for( int k=0 ; k<maxNewtonIter ; ++k )
    {

      // call assemble to fill the matrix and the rhs
      real64 residualNorm = AssembleSystem( domain, m_linearSystem, time_n, stepDt );


      if( residualNorm < newtonTol )
      {
        isConverged = 1;
        break;
      }

      // do line search in case residual has increased
      if( residualNorm > lastResidual )
      {
        real64 scaleFactor = -1.0;
        for( int lineSearchIteration=0 ; lineSearchIteration<maxNumberLineSearchCuts ; ++lineSearchIteration )
        {
          scaleFactor *= 0.5;
          ApplySystemSolution( m_linearSystem, scaleFactor, 0, domain );
          residualNorm = AssembleSystem( domain, m_linearSystem, time_n, stepDt );
          if( residualNorm < lastResidual )
          {
            break;
          }
        }
      }

      // call the default linear solver on the system
      SolveSystem( m_linearSystem,
                   getSystemSolverParameters() );

      // apply the system solution to the fields/variables
      ApplySystemSolution( m_linearSystem, 1.0, 0, domain );

      lastResidual = residualNorm;
    }
    if( isConverged )
    {
      // final step for completion of timestep. typically secondary variable updates and cleanup.
      ImplicitStepComplete( time_n, stepDt,  domain );
      break;
    }
    else
    {
      // cut timestep and reset state
      stepDt *= 0.5;
      ResetStateToBeginningOfStep( domain );
    }
  }


  // return the achieved timestep
  return stepDt;
}


real64 SolverBase::ExplicitStep( real64 const & time_n,
                               real64 const & dt,
                               integer const cycleNumber,
                               DomainPartition * const domain )
{
  GEOS_ERROR( "SolverBase::ExplicitStep called!. Should be overridden.");
  return 0;
}

void SolverBase::ImplicitStepSetup( real64 const& time_n,
                        real64 const& dt,
                        DomainPartition * const domain )
{
  GEOS_ERROR( "SolverBase::ImplicitStepSetup called!. Should be overridden.");
}

real64 SolverBase::AssembleSystem( DomainPartition * const domain,
                             systemSolverInterface::EpetraBlockSystem * const blockSystem,
                             real64 const time,
                             real64 const dt )
{
  GEOS_ERROR( "SolverBase::Assemble called!. Should be overridden.");
  return 0;
}

void SolverBase::SolveSystem( systemSolverInterface::EpetraBlockSystem * const blockSystem,
                              SystemSolverParameters const * const params )
{
  GEOS_ERROR( "SolverBase::SolveSystem called!. Should be overridden.");
}

void SolverBase::ApplySystemSolution( systemSolverInterface::EpetraBlockSystem const * const blockSystem,
                          real64 const scalingFactor,
                          localIndex const dofOffset,
                          DomainPartition * const  )
{
  GEOS_ERROR( "SolverBase::ApplySystemSolution called!. Should be overridden.");
}

void SolverBase::ResetStateToBeginningOfStep( DomainPartition * const  )
{
  GEOS_ERROR( "SolverBase::ResetStateToBeginningOfStep called!. Should be overridden.");
}

void SolverBase::ImplicitStepComplete( real64 const & time,
                           real64 const & dt,
                           DomainPartition * const domain )
{
  GEOS_ERROR( "SolverBase::ImplicitStepComplete called!. Should be overridden.");
}




void SolverBase::CreateChild( string const & childKey, string const & childName )
{
  if( CatalogInterface::hasKeyName(childKey) )
  {
    std::cout << "Adding Solver of type " << childKey << ", named " << childName << std::endl;
    this->RegisterGroup( childName, CatalogInterface::Factory( childKey, childName, this ) );
  }
}


R1Tensor const * SolverBase::globalGravityVector() const
{
  R1Tensor const * rval = nullptr;
  if( getParent()->getName() == "Solvers" )
  {
    rval = &(getParent()->
             group_cast<SolverBase const *>()->
             getGravityVector());
  }

  return rval;
}

SystemSolverParameters * SolverBase::getSystemSolverParameters()
{
  return this->GetGroup<SystemSolverParameters>(groupKeys.systemSolverParameters);
}


} /* namespace ANST */
