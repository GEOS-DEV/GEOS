/*
 * DataObjectBase.cpp
 *
 *  Created on: Jun 17, 2016
 *      Author: rrsettgast
 */

#include "ViewWrapperBase.hpp"

#include "ManagedGroup.hpp"
#include "RestartFlags.hpp"

#ifdef USE_ATK
#include "slic/slic.hpp"
#endif

namespace geosx
{
namespace dataRepository
{


ViewWrapperBase::ViewWrapperBase( std::string const & name,
                                  ManagedGroup * const parent):
  m_name(name),
  m_parent(parent),
  m_sizedFromParent(1),
  m_restart_flags(RestartFlags::WRITE_AND_READ)
#ifdef USE_ATK
  ,m_sidreView(nullptr)
#endif
{
#ifdef USE_ATK
  SLIC_ERROR_IF(parent==nullptr,"parameter WrapperCollection * const parent must not be nullptr");

  if( parent->getSidreGroup()->hasView(name) )
  {
    m_sidreView = parent->getSidreGroup()->getView(name);
  }
  else
  {
    m_sidreView = parent->getSidreGroup()->createView(name);
  }
#endif
}


ViewWrapperBase::~ViewWrapperBase()
{}


ViewWrapperBase::ViewWrapperBase( ViewWrapperBase&& source ):
  m_name( std::move(source.m_name) ),
  m_parent( source.m_parent),
  m_sizedFromParent( source.m_sizedFromParent),
  m_restart_flags(RestartFlags::WRITE_AND_READ)
#ifdef USE_ATK
  ,m_sidreView( source.m_sidreView )
#endif

{}

void ViewWrapperBase::resize()
{
  resize(m_parent->size());
}


}
} /* namespace geosx */
