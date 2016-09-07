/*
 * FiniteElementSpace.hpp
 *
 *  Created on: Aug 4, 2016
 *      Author: rrsettgast
 */

#ifndef SRC_COMPONENTS_CORE_SRC_FINITEELEMENT_FINITEELEMENTSPACE_HPP_
#define SRC_COMPONENTS_CORE_SRC_FINITEELEMENT_FINITEELEMENTSPACE_HPP_
#include "../dataRepository/ManagedGroup.hpp"

namespace geosx
{


class FiniteElementSpace : public dataRepository::ManagedGroup
{
public:
  FiniteElementSpace() = delete;

  explicit FiniteElementSpace( std::string const & name, ManagedGroup * const parent );

  ~FiniteElementSpace();

  virtual void Registration( dataRepository::ManagedGroup * const parent );

  virtual dataRepository::ManagedGroup & getNodeManager();
  virtual dataRepository::ManagedGroup & getEdgeManager();
  virtual dataRepository::ManagedGroup & getFaceManager();
  virtual dataRepository::ManagedGroup & getElementManager();




};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_FINITEELEMENT_FINITEELEMENTSPACE_HPP_ */
