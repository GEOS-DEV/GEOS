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
{
}

PhysicsSolverManager::~PhysicsSolverManager()
{
}



// SolverBase * PhysicsSolverManager::CreateSolver( string const & solverCatalogKey, string const & solverName )
// {
//   std::unique_ptr<SolverBase> solver = SolverBase::CatalogInterface::Factory( solverCatalogKey, solverName, this );
//   SolverBase * rval = this->RegisterGroup<SolverBase>( solverName, std::move(solver) );

//   return rval;
// }

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


// void PhysicsSolverManager::ReadXML( dataRepository::ManagedGroup * domain,
//                                     xmlWrapper::xmlNode const & problemNode )
// {

//   // Solvers
//   xmlWrapper::xmlNode topLevelNode = problemNode.child("Solvers");
//   std::cout << "Solvers:" << std::endl;
//   if (topLevelNode == NULL)
//   {
//     throw std::invalid_argument("Solver block not present in input xml file!");
//   }
//   else
//   {
//     for (xmlWrapper::xmlNode solverNode=topLevelNode.first_child(); solverNode; solverNode=solverNode.next_sibling())
//     {
//       std::cout << "   " << solverNode.name() << std::endl;

//       // Register the new solver
//       std::string solverID = solverNode.attribute("name").value();
//       SolverBase * newSolver = CreateSolver( solverNode.name(), solverID );

//       // Set the documentation node
//       newSolver->SetDocumentationNodes( domain );

//       // Register fields in the solver and parse options
//       newSolver->BuildDataStructure( domain );

//       newSolver->ReadXML(solverNode );
//     }
//   }

// }


} /* namespace geosx */
