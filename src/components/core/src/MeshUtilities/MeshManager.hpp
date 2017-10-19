/*
 * MeshManager.hpp
 *
 *  Created on: Sep 7, 2016
 *      Author: rrsettgast
 */

#ifndef SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_MESHMANAGER_HPP_
#define SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_MESHMANAGER_HPP_

#include "dataRepository/ManagedGroup.hpp"


namespace geosx
{
class SolverBase;

class MeshManager : public dataRepository::ManagedGroup
{
public:
  MeshManager( std::string const & name,
                        ManagedGroup * const parent );

  virtual ~MeshManager();

  virtual void FillDocumentationNode( dataRepository::ManagedGroup * const group ) override;

  virtual void CreateChild( string const & childKey, string const & childName ) override;

  void GenerateMeshes( dataRepository::ManagedGroup * const domain );

private:
  MeshManager() = delete;

};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_MESHMANAGER_HPP_ */
