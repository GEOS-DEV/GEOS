/*
 * BasisFunctionManager.hpp
 *
 *  Created on: Dec 5, 2017
 *      Author: sherman
 */

#ifndef SRC_COMPONENTS_CORE_SRC_FINITEELEMENT_BASISFUNCTIONMANAGER_HPP_
#define SRC_COMPONENTS_CORE_SRC_FINITEELEMENT_BASISFUNCTIONMANAGER_HPP_

#include "dataRepository/ManagedGroup.hpp"


namespace geosx
{
namespace dataRepository
{
namespace keys
{
}
}


class BasisFunctionManager : public dataRepository::ManagedGroup
{
public:
  BasisFunctionManager() = delete;
  BasisFunctionManager(string const & name, ManagedGroup * const parent);
  virtual ~BasisFunctionManager();

  virtual void FillDocumentationNode() override;
  virtual void CreateChild( string const & childKey, string const & childName ) override;


};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_FINITEELEMENT_BASISFUNCTIONMANAGER_HPP_ */
