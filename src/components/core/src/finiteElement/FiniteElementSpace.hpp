/*
 * FiniteElementSpace.hpp
 *
 *  Created on: Aug 4, 2016
 *      Author: rrsettgast
 */

#ifndef SRC_COMPONENTS_CORE_SRC_FINITEELEMENT_FINITEELEMENTSPACE_HPP_
#define SRC_COMPONENTS_CORE_SRC_FINITEELEMENT_FINITEELEMENTSPACE_HPP_
#include "dataRepository/WrapperCollection.hpp"

namespace geosx
{


class FiniteElementSpace : public dataRepository::WrapperCollection
{
public:
  FiniteElementSpace() = delete;

  explicit FiniteElementSpace( std::string const & name, WrapperCollection * const parent );

  ~FiniteElementSpace();

  virtual void Registration( dataRepository::WrapperCollection * const parent );

  virtual dataRepository::WrapperCollection & getNodeManager();
  virtual dataRepository::WrapperCollection & getEdgeManager();
  virtual dataRepository::WrapperCollection & getFaceManager();
  virtual dataRepository::WrapperCollection & getElementManager();




};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_FINITEELEMENT_FINITEELEMENTSPACE_HPP_ */
