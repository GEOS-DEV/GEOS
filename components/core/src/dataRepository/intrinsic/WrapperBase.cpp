/*
 * DataObjectBase.cpp
 *
 *  Created on: Jun 17, 2016
 *      Author: rrsettgast
 */

#include "WrapperBase.hpp"

namespace geosx {
namespace dataRepository
{

WrapperBase::WrapperBase( std::string const & name ):
  m_name(name)
{}


WrapperBase::~WrapperBase()
{}


WrapperBase::WrapperBase( WrapperBase&& source ):
    m_name( std::move(source.m_name) )
{}

}
} /* namespace geosx */
