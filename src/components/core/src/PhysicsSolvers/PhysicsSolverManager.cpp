/*
 * PhysicsSolverManager.cpp
 *
 *  Created on: Sep 7, 2016
 *      Author: rrsettgast
 */

#include "PhysicsSolverManager.hpp"
#include "SolverBase.hpp"
#include <pugixml.hpp>

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



SolverBase & PhysicsSolverManager::CreateSolver( string const & solverCatalogKey, string const & solverName )
{
  std::unique_ptr<SolverBase> solver = SolverBase::CatalogInterface::Factory( solverCatalogKey, solverName, this );
  SolverBase & rval = this->RegisterGroup( solverName, std::move(solver) );

  return rval;
}

void PhysicsSolverManager::ReadXML( pugi::xml_node const & problemNode )
{
  // Solvers
  pugi::xml_node topLevelNode = problemNode.child("Solvers");
  std::cout << "Solvers:" << std::endl;
  if (topLevelNode == NULL)
  {
    throw std::invalid_argument("Solver block not present in input xml file!");
  }
  else
  {
    // Determine the number of active solvers, resize the solver collection
    std_size_t nSolvers = std::distance(topLevelNode.children().begin(), topLevelNode.children().end());

    int ii = 0;

    for (pugi::xml_node solverNode=topLevelNode.first_child(); solverNode; solverNode=solverNode.next_sibling())
    {
      std::cout << "   " << solverNode.name() << std::endl;

      // Register the new solver
      std::string solverID = solverNode.attribute("name").value();
      SolverBase & newSolver = CreateSolver( solverNode.name(), solverID );


      // Register fields in the solver and parse options
//      newSolver.Registration( &domain );
      newSolver.ReadXML(solverNode);
      ii++;
    }
  }

}


} /* namespace geosx */
