/*
 * PhysicsSolverManager.hpp
 *
 *  Created on: Sep 7, 2016
 *      Author: rrsettgast
 */

#ifndef SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_PHYSICSSOLVERMANAGER_HPP_
#define SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_PHYSICSSOLVERMANAGER_HPP_

#include "dataRepository/ManagedGroup.hpp"

namespace geosx
{

class PhysicsSolverManager : public dataRepository::ManagedGroup
{
public:
  PhysicsSolverManager( std::string const & name,
                        ManagedGroup * const parent );

  virtual ~PhysicsSolverManager();

  void CreateSolver( string const & solverCatalogKey, string const & solverName );

};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_PHYSICSSOLVERMANAGER_HPP_ */
