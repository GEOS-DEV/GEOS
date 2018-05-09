/*
 * GeometricObjectManager.hpp
 *
 *  Created on: Sep 7, 2016
 *      Author: rrsettgast
 */

#ifndef SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_GEOMETRICOBJECTMANAGER_HPP_
#define SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_GEOMETRICOBJECTMANAGER_HPP_

#include "dataRepository/ManagedGroup.hpp"


namespace geosx
{

class GeometricObjectManager : public dataRepository::ManagedGroup
{
public:
  GeometricObjectManager( std::string const & name,
                          ManagedGroup * const parent );

  virtual ~GeometricObjectManager() override;

  virtual void FillDocumentationNode() override;

  virtual void CreateChild( string const & childKey, string const & childName ) override;

private:
  GeometricObjectManager() = delete;

};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_GEOMETRICOBJECTMANAGER_HPP_ */
