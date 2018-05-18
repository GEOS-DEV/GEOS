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

namespace geosx
{
class SolverBase;

class PhysicsSolverManager : public dataRepository::ManagedGroup
{
public:
  PhysicsSolverManager( std::string const & name,
                        ManagedGroup * const parent );

  virtual ~PhysicsSolverManager() override;

  virtual void FillDocumentationNode() override;

  virtual void CreateChild( string const & childKey, string const & childName ) override;

private:
  PhysicsSolverManager() = delete;

};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_PHYSICSSOLVERMANAGER_HPP_ */
