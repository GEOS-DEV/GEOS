/*
 * DataObjectBase.cpp
 *
 *  Created on: Jun 17, 2016
 *      Author: rrsettgast
 */

#include "DataObjectBase.hpp"

namespace geosx {

DataObjectBase::DataObjectBase( const std::string& name ):
  m_fieldName(name)
{}

DataObjectBase::DataObjectBase() {
	// TODO Auto-generated constructor stub

}

DataObjectBase::~DataObjectBase() {
	// TODO Auto-generated destructor stub
}

} /* namespace geosx */
