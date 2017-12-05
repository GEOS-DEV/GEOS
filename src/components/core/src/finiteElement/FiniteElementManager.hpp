/*
 * FiniteElementManager.hpp
 *
 *  Created on: Apr 18, 2017
 *      Author: rrsettgast
 */

#ifndef SRC_COMPONENTS_CORE_SRC_FINITEELEMENT_FINITEELEMENTMANAGER_HPP_
#define SRC_COMPONENTS_CORE_SRC_FINITEELEMENT_FINITEELEMENTMANAGER_HPP_

#include "FiniteElementSpace.hpp"

namespace geosx
{
namespace dataRepository
{
namespace keys
{
string const finiteElementManager = "finiteElementManager";
string const basisFunctions = "BasisFunctions";
string const quadratureRules = "QuadratureRules";
string const finiteElements = "FiniteElements";
}
}


class FiniteElementManager : public dataRepository::ManagedGroup
{
public:
  FiniteElementManager() = delete;
  FiniteElementManager(string const & name, ManagedGroup * const parent);
  virtual ~FiniteElementManager();

  void FillDocumentationNode() override final;
  virtual void CreateChild( string const & childKey, string const & childName ) override;


};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_FINITEELEMENT_FINITEELEMENTMANAGER_HPP_ */
