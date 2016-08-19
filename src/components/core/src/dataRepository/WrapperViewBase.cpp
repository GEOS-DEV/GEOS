/*
 * DataObjectBase.cpp
 *
 *  Created on: Jun 17, 2016
 *      Author: rrsettgast
 */

#include "WrapperViewBase.hpp"

#include "slic/slic.hpp"
#include "SynchronizedGroup.hpp"

namespace geosx
{
namespace dataRepository
{

WrapperViewBase::WrapperViewBase( std::string const & name,
                          SynchronizedGroup * const parent ) :
  m_name(name),
  m_parent(parent),
  m_sizedFromParent(0),
  m_sidreView(nullptr)
{
  SLIC_ERROR_IF(parent==nullptr,"parameter WrapperCollection * const parent must not be nullptr");

  if( parent->getSidreGroup()->hasView(name) )
  {
    m_sidreView = parent->getSidreGroup()->getView(name);
  }
  else
  {
    m_sidreView = parent->getSidreGroup()->createView(name);
  }
}


WrapperViewBase::~WrapperViewBase()
{}


WrapperViewBase::WrapperViewBase( WrapperViewBase&& source ) :
  m_name( std::move(source.m_name) ),
  m_parent( source.m_parent),
  m_sizedFromParent( source.m_sizedFromParent),
  m_sidreView( source.m_sidreView )
{}

}
} /* namespace geosx */
