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
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

#include "SolverBase.hpp"
#include "PhysicsSolverManager.hpp"
#include "managers/DomainPartition.hpp"
#include "mesh/MeshBody.hpp"
#include "systemSolverInterface/EpetraBlockSystem.hpp"

namespace geosx
{

using namespace dataRepository;

SolverBase::SolverBase( std::string const & name,
                        ManagedGroup * const parent ):
  ExecutableGroup( name, parent ),
  m_linearSolverWrapper(),
  m_verboseLevel(0),
  m_gravityVector( R1Tensor(0.0) ),
  m_systemSolverParameters( groupKeyStruct::systemSolverParametersString, this )//,
//  m_blockLocalDofNumber()
{
  // register group with repository. Have Repository own object.
  this->RegisterGroup( groupKeyStruct::systemSolverParametersString, &m_systemSolverParameters, 0 );

  this->RegisterViewWrapper( viewKeyStruct::verboseLevelString, &m_verboseLevel, 0 );
  this->RegisterViewWrapper( viewKeyStruct::gravityVectorString, &m_gravityVector, 0 );
//  this->RegisterViewWrapper( viewKeyStruct::blockLocalDofNumberString, &m_blockLocalDofNumber, 0 );

  if( this->globalGravityVector() != nullptr )
  {
    m_gravityVector=*globalGravityVector();
  }

//  m_linearSolverWrapper = new systemSolverInterface::LinearSolverWrapper();

}

SolverBase::~SolverBase()
{
//  delete m_linearSolverWrapper;
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

void SolverBase::FillOtherDocumentationNodes( dataRepository::ManagedGroup * const rootGroup )
{

}


real64 SolverBase::SolverStep( real64 const& time_n,
                           real64 const& dt,
                           const int cycleNumber,
                           DomainPartition * domain )
{
  return 0;
}


void SolverBase::Execute( real64 const& time_n,
                          real64 const& dt,
                          const int cycleNumber,
                          ManagedGroup * domain )
{
  if ( dt > 0 )
  {
    SolverStep(time_n, dt, cycleNumber, domain->group_cast<DomainPartition*>());
  }
}


real64 SolverBase::LinearImplicitStep( real64 const & time_n,
                                       real64 const & dt,
                                       integer const cycleNumber,
                                       DomainPartition * const domain,
                                       systemSolverInterface::EpetraBlockSystem * const blockSystem )
{
  // call setup for physics solver. Pre step allocations etc.
  ImplicitStepSetup( time_n, dt, domain, blockSystem );

  // call assemble to fill the matrix and the rhs
  AssembleSystem( domain, blockSystem, time_n+dt, dt );

  // apply boundary conditions to system
  ApplyBoundaryConditions( domain, blockSystem, time_n, dt );

  // call the default linear solver on the system
  SolveSystem( blockSystem,
               getSystemSolverParameters() );

  // apply the system solution to the fields/variables
  ApplySystemSolution( blockSystem, 1.0, domain );

  // final step for completion of timestep. typically secondary variable updates and cleanup.
  ImplicitStepComplete( time_n, dt,  domain );

  // return the achieved timestep
  return dt;
}

real64 SolverBase::NonlinearImplicitStep( real64 const & time_n,
                                           real64 const & dt,
                                           integer const cycleNumber,
                                           DomainPartition * const domain,
                                           systemSolverInterface::EpetraBlockSystem * const blockSystem )
{
  // dt may be cut during the course of this step, so we are keeping a local
  // value to track the achieved dt for this step.
  real64 stepDt = dt;

  // call setup for physics solver. Pre step allocations etc.
  ImplicitStepSetup( time_n, stepDt, domain, blockSystem );

  SystemSolverParameters const * const solverParams = getSystemSolverParameters();
  real64 const newtonTol = solverParams->newtonTol();
  integer const maxNewtonIter = solverParams->maxIterNewton();
  integer const maxNumberDtCuts = 2;
  integer const maxNumberLineSearchCuts = 5;

  // a flag to denote whether we have converged
  integer isConverged = 0;

  // outer loop attempts to apply full timestep, and managed the cutting of the timestep if
  // required.
  for( int dtAttempt = 0 ; dtAttempt<maxNumberDtCuts ; ++dtAttempt )
  {
    // main Newton loop
    // keep residual from previous iteration in case we need to do a line search
    real64 lastResidual = 1e99;
    for( int k=0 ; k<maxNewtonIter ; ++k )
    {

      // call assemble to fill the matrix and the rhs
      AssembleSystem( domain, blockSystem, time_n, stepDt );

      // apply boundary conditions to system
      ApplyBoundaryConditions( domain, blockSystem, time_n, stepDt );

      // get residual norm
      real64 residualNorm = CalculateResidualNorm(blockSystem, domain);

      if (m_verboseLevel >= 1)
      {
        std::cout << "Attempt: " << dtAttempt  << ", Newton: " << k
                  << ", R = " << residualNorm << std::endl;
      }

      // if the residual norm is less than the Newton tolerance we denote that we have
      // converged and break from the Newton loop immediately.
      if( residualNorm < newtonTol )
      {
        isConverged = 1;
        break;
      }


      // do line search in case residual has increased
      if( residualNorm > lastResidual )
      {

        // flag to determine if we should solve the system and apply the solution. If the line
        // search fails we just bail.
        int lineSearchSuccess = 0;

        // scale factor is value applied to the previous solution. In this case we want to
        // subtract a portion of the previous solution.
        real64 scaleFactor = -1.0;

        // main loop for the line search.
        for( int lineSearchIteration=0 ; lineSearchIteration<maxNumberLineSearchCuts ; ++lineSearchIteration )
        {
          // cut the scale factor by half. This means that the scale factors will
          // have values of -0.5, -0.25, -0.125, ...
          scaleFactor *= 0.5;
          ApplySystemSolution( blockSystem, scaleFactor, domain );

          // re-assemble system
          AssembleSystem( domain, blockSystem, time_n, stepDt );

          // apply boundary conditions to system
          ApplyBoundaryConditions( domain, blockSystem, time_n, stepDt );

          // get residual norm
          residualNorm = CalculateResidualNorm(blockSystem, domain);

          if (m_verboseLevel >= 1)
          {
            std::cout << "Attempt: " << dtAttempt + 1 << ", Newton: " << k + 1
                      << ", Line search: " << lineSearchIteration + 1
                      << ", R = " << residualNorm << std::endl;
          }

          // if the residual norm is less than the last residual, we can proceed to the
          // solution step
          if( residualNorm < lastResidual )
          {
            lineSearchSuccess = 1;
            break;
          }
        }

        // if line search failed, then break out of the main Newton loop. Timestep will be cut.
        if( !lineSearchSuccess )
        {
          break;
        }
      }

      // call the default linear solver on the system
      SolveSystem( blockSystem,
                   getSystemSolverParameters() );

      // apply the system solution to the fields/variables
      ApplySystemSolution( blockSystem, 1.0, domain );

      lastResidual = residualNorm;
    }
    if( isConverged )
    {
      // break out of outer loop
      break;
    }
    else
    {
      // cut timestep and reset state to beginning of step, and restart the Newton loop.
      stepDt *= 0.5;
      ResetStateToBeginningOfStep( domain );
    }
  }

  if( !isConverged )
  {
    std::cout<<"Convergence not achieved.";
  }

  // final step for completion of timestep. typically secondary variable updates and cleanup.
  ImplicitStepComplete( time_n, stepDt,  domain );

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
                        DomainPartition * const domain,
                        systemSolverInterface::EpetraBlockSystem * const blockSystem )
{
  GEOS_ERROR( "SolverBase::ImplicitStepSetup called!. Should be overridden.");
}

void SolverBase::AssembleSystem( DomainPartition * const domain,
                                 systemSolverInterface::EpetraBlockSystem * const blockSystem,
                                 real64 const time,
                                 real64 const dt )
{
  GEOS_ERROR( "SolverBase::Assemble called!. Should be overridden.");
}


void SolverBase::ApplyBoundaryConditions( DomainPartition * const domain,
                                          systemSolverInterface::EpetraBlockSystem * const blockSystem,
                                          real64 const time,
                                          real64 const dt )
{
  GEOS_ERROR( "SolverBase::SolveSystem called!. Should be overridden.");
}

real64
SolverBase::
CalculateResidualNorm( systemSolverInterface::EpetraBlockSystem const *const blockSystem,
                       DomainPartition * const domain )
{
  GEOS_ERROR( "SolverBase::CalculateResidualNorm called!. Should be overridden.");
  return 0;
}

void SolverBase::SolveSystem( systemSolverInterface::EpetraBlockSystem * const blockSystem,
                              SystemSolverParameters const * const params )
{
  GEOS_ERROR( "SolverBase::SolveSystem called!. Should be overridden.");
}

void SolverBase::ApplySystemSolution( systemSolverInterface::EpetraBlockSystem const * const blockSystem,
                                      real64 const scalingFactor,
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


void SolverBase::SolveSystem( systemSolverInterface::EpetraBlockSystem * const blockSystem,
                              SystemSolverParameters const * const params,
                              systemSolverInterface::BlockIDs const blockID )
{
  Epetra_FEVector * const
  solution = blockSystem->GetSolutionVector( blockID );

  Epetra_FEVector * const
  residual = blockSystem->GetResidualVector( blockID );
  residual->Scale(-1.0);

  solution->Scale(0.0);

  m_linearSolverWrapper.SolveSingleBlockSystem( blockSystem,
                                                 params,
                                                 blockID );

  if( verboseLevel() >= 2 )
  {
    solution->Print(std::cout);
  }

}

//void SolverBase::CreateChild( string const & childKey, string const & childName )
//{
//  if( CatalogInterface::hasKeyName(childKey) )
//  {
//    std::cout << "Adding Solver of type " << childKey << ", named " << childName << std::endl;
//    this->RegisterGroup( childName, CatalogInterface::Factory( childKey, childName, this ) );
//  }
//}


R1Tensor const * SolverBase::globalGravityVector() const
{
  R1Tensor const * rval = nullptr;
  if( getParent()->getName() == "Solvers" )
  {
    rval = &(getParent()->getReference<R1Tensor>( viewKeyStruct::gravityVectorString ));
  }

  return rval;
}


systemSolverInterface::EpetraBlockSystem const * SolverBase::getLinearSystemRepository() const
{
  return &( getParent()->
            getReference<systemSolverInterface::EpetraBlockSystem>( PhysicsSolverManager::
                                                                    viewKeyStruct::
                                                                    blockSystemRepositoryString ) );
}

systemSolverInterface::EpetraBlockSystem * SolverBase::getLinearSystemRepository()
{
  return &( getParent()->
            getReference<systemSolverInterface::EpetraBlockSystem>( PhysicsSolverManager::
                                                                    viewKeyStruct::
                                                                    blockSystemRepositoryString ) );
}

void SolverBase::FinalInitialization(ManagedGroup * const group)
{

}

} /* namespace ANST */
