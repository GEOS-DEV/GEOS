/*
 * DataObjectBase.cpp
 *
 *  Created on: Jun 17, 2016
 *      Author: rrsettgast
 */

#include "ViewWrapperBase.hpp"

#include "ManagedGroup.hpp"

#if ATK_FOUND
#include "slic/slic.hpp"
#endif

namespace geosx
{
namespace dataRepository
{

ViewWrapperBase::ViewWrapperBase( std::string const & name,
                                  ManagedGroup * const parent ) :
#if ATK_FOUND
  m_sidreView(nullptr),
#endif
  m_name(name),
  m_parent(parent),
  m_sizedFromParent(1)

{
#if ATK_FOUND
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


ViewWrapperBase::ViewWrapperBase( ViewWrapperBase&& source ) :
  #if ATK_FOUND
  m_sidreView( source.m_sidreView ),
#endif
  m_name( std::move(source.m_name) ),
  m_parent( source.m_parent),
  m_sizedFromParent( source.m_sizedFromParent)
{}

void ViewWrapperBase::resize()
{
  resize(m_parent->size());
}


}
} /* namespace geosx */
