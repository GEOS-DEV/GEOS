/*
 * PhysicsSolverManager.cpp
 *
 *  Created on: Sep 7, 2016
 *      Author: rrsettgast
 */

#include "PhysicsSolverManager.hpp"

#include "../../../cxx-utilities/src/src/DocumentationNode.hpp"
#include "SolverBase.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace cxx_utilities;

PhysicsSolverManager::PhysicsSolverManager( std::string const & name,
                                            ManagedGroup * const parent ):
  ManagedGroup( name, parent)
{}

PhysicsSolverManager::~PhysicsSolverManager()
{}


void PhysicsSolverManager::FillDocumentationNode( dataRepository::ManagedGroup * const /*group*/ )
{
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();
  docNode->setName("Solvers");
  docNode->setSchemaType("UniqueNode");
  docNode->setShortDescription("Solver manager");

}


void PhysicsSolverManager::CreateChild( string const & childKey, string const & childName )
{
  std::cout << "Adding Solver: " << childKey << ", " << childName << std::endl;
  std::unique_ptr<SolverBase> solver = SolverBase::CatalogInterface::Factory( childKey, childName, this );
  this->RegisterGroup<SolverBase>( childName, std::move(solver) );
}



} /* namespace geosx */
