/*
 * DataObject.h
 *
 *  Created on: Jun 8, 2016
 *      Author: settgast
 */

#ifndef CORE_SRC_DATASTRUCTURES_DATAOBJECT_HPP_
#define CORE_SRC_DATASTRUCTURES_DATAOBJECT_HPP_

#include "sidre/sidre.hpp"

namespace geosx
{

  class DataObject
  {
  public:
    DataObject();
    virtual ~DataObject();
  };

} /* namespace geosx */

#endif /* CORE_SRC_DATASTRUCTURES_DATAOBJECT_HPP_ */
