/*
 * FiniteElementSpaceManager.hpp
 *
 *  Created on: Dec 5, 2017
 *      Author: sherman
 */

#ifndef SRC_COMPONENTS_CORE_SRC_FINITEELEMENT_FINITEELEMENTSPACEMANAGER_HPP_
#define SRC_COMPONENTS_CORE_SRC_FINITEELEMENT_FINITEELEMENTSPACEMANAGER_HPP_

#include "dataRepository/ManagedGroup.hpp"

namespace geosx
{
namespace dataRepository
{
namespace keys
{
}
}


class FiniteElementSpaceManager : public dataRepository::ManagedGroup
{
public:
  FiniteElementSpaceManager() = delete;
  FiniteElementSpaceManager(string const & name, ManagedGroup * const parent);
  virtual ~FiniteElementSpaceManager();

  virtual void FillDocumentationNode() override;
  virtual void CreateChild( string const & childKey, string const & childName ) override;


};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_FINITEELEMENT_FINITEELEMENTSPACEMANAGER_HPP_ */
