/*
 * FiniteElementSpace.hpp
 *
 *  Created on: Aug 4, 2016
 *      Author: rrsettgast
 */

#ifndef SRC_COMPONENTS_CORE_SRC_FINITEELEMENT_FINITEELEMENTSPACE_HPP_
#define SRC_COMPONENTS_CORE_SRC_FINITEELEMENT_FINITEELEMENTSPACE_HPP_
#include "../dataRepository/SynchronizedGroup.hpp"

namespace geosx
{


class FiniteElementSpace : public dataRepository::SynchronizedGroup
{
public:
  FiniteElementSpace() = delete;

  explicit FiniteElementSpace( std::string const & name, SynchronizedGroup * const parent );

  ~FiniteElementSpace();

  virtual void Registration( dataRepository::SynchronizedGroup * const parent );

  virtual dataRepository::SynchronizedGroup & getNodeManager();
  virtual dataRepository::SynchronizedGroup & getEdgeManager();
  virtual dataRepository::SynchronizedGroup & getFaceManager();
  virtual dataRepository::SynchronizedGroup & getElementManager();




};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_FINITEELEMENT_FINITEELEMENTSPACE_HPP_ */
