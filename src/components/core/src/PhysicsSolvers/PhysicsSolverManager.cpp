/*
 * PhysicsSolverManager.cpp
 *
 *  Created on: Sep 7, 2016
 *      Author: rrsettgast
 */

#include "PhysicsSolverManager.hpp"
#include "SolverBase.hpp"

namespace geosx
{

using namespace dataRepository;

PhysicsSolverManager::PhysicsSolverManager( std::string const & name,
                                            ManagedGroup * const parent ):
  ManagedGroup( name, parent)
{
}

PhysicsSolverManager::~PhysicsSolverManager()
{
}



void PhysicsSolverManager::CreateSolver( string const & solverCatalogKey, string const & solverName )
{
  std::unique_ptr<SolverBase> solver = SolverBase::CatalogInterface::Factory( solverCatalogKey, solverName, this );
  this->RegisterGroup( solverName, std::move(solver) );
}


} /* namespace geosx */
