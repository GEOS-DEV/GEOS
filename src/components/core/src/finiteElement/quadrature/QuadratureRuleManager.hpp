/*
 * QuadratureRuleManager.hpp
 *
 *  Created on: Apr 18, 2017
 *      Author: rrsettgast
 */

#ifndef SRC_COMPONENTS_CORE_SRC_FINITEELEMENT_QUADRATURERULEMANAGER_HPP_
#define SRC_COMPONENTS_CORE_SRC_FINITEELEMENT_QUADRATURERULEMANAGER_HPP_

#include "dataRepository/ManagedGroup.hpp"
#include "fileIO/xmlWrapper.hpp"

namespace geosx
{
namespace dataRepository
{
namespace keys
{
}
}


class QuadratureRuleManager : public dataRepository::ManagedGroup
{
public:
  QuadratureRuleManager() = delete;
  QuadratureRuleManager(string const & name, ManagedGroup * const parent);
  virtual ~QuadratureRuleManager();

  virtual void FillDocumentationNode() override;
  virtual void CreateChild( string const & childKey, string const & childName ) override;
  virtual void ReadXMLsub( xmlWrapper::xmlNode const & targetNode ) override;

};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_FINITEELEMENT_QUADRATURERULEMANAGER_HPP_ */
