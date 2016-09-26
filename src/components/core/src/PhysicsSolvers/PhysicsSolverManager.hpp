/*
 * PhysicsSolverManager.hpp
 *
 *  Created on: Sep 7, 2016
 *      Author: rrsettgast
 */

#ifndef SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_PHYSICSSOLVERMANAGER_HPP_
#define SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_PHYSICSSOLVERMANAGER_HPP_

#include "dataRepository/ManagedGroup.hpp"

namespace pugi
{
class xml_node;
}

namespace cxx_utilities
{
struct DocumentationNode;
}

namespace geosx
{
class SolverBase;

class PhysicsSolverManager : public dataRepository::ManagedGroup
{
public:
  PhysicsSolverManager( std::string const & name,
                        ManagedGroup * const parent );

  virtual ~PhysicsSolverManager();

  SolverBase & CreateSolver( string const & solverCatalogKey, string const & solverName );


  void ReadXML( dataRepository::ManagedGroup& domain, pugi::xml_node const & problemNode, cxx_utilities::DocumentationNode & docNode  );

};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_PHYSICSSOLVERMANAGER_HPP_ */
