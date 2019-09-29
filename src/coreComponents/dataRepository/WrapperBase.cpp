/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/** @file */

#include "WrapperBase.hpp"

#include "Group.hpp"
#include "RestartFlags.hpp"


namespace geosx
{
namespace dataRepository
{


WrapperBase::WrapperBase( std::string const & name,
                          Group * const parent ):
  m_name( name ),
  m_parent( parent ),
  m_sizedFromParent( 1 ),
  m_restart_flags( RestartFlags::WRITE_AND_READ ),
  m_plotLevel( PlotLevel::LEVEL_3 ),
  m_inputFlag( InputFlags::INVALID ),
  m_description(),
  m_registeringObjects()
#ifdef GEOSX_USE_ATK
  , m_sidreView( nullptr )
#endif
{
#ifdef GEOSX_USE_ATK
  GEOS_ERROR_IF( parent==nullptr, "parameter WrapperCollection * const parent must not be nullptr" );

  if( parent->getSidreGroup()->hasView( name ) )
  {
    m_sidreView = parent->getSidreGroup()->getView( name );
  }
  else
  {
    m_sidreView = parent->getSidreGroup()->createView( name );
  }
#endif
}


WrapperBase::~WrapperBase()
{}


WrapperBase::WrapperBase( WrapperBase && source ):
  m_name( std::move( source.m_name ) ),
  m_parent( source.m_parent ),
  m_sizedFromParent( source.m_sizedFromParent ),
  m_restart_flags( source.m_restart_flags )
#ifdef GEOSX_USE_ATK
  , m_sidreView( source.m_sidreView )
#endif

{}

void WrapperBase::resize()
{
  resize( m_parent->size());
}

void WrapperBase::CopyWrapperAttributes( WrapperBase const & source )
{
  m_name = source.m_name;
  m_sizedFromParent = source.m_sizedFromParent;
  m_restart_flags = source.m_restart_flags;
}

#if defined(USE_TOTALVIEW_OUTPUT)
int WrapperBase::setTotalviewDisplay() const
{
  //std::cout<<"exectuing WrapperBase::setTotalviewDisplay()"<<std::endl;
//  TV_ttf_add_row("TYPE", TV_ttf_type_ascii_string, type.c_str() );
  TV_ttf_add_row( "m_name", totalview::typeName< string >().c_str(), &m_name );
  TV_ttf_add_row( "m_parent", totalview::typeName< Group >().c_str(), m_parent );
  TV_ttf_add_row( "m_sizedFromParent", "int", &m_sizedFromParent );
  TV_ttf_add_row( "m_restart_flags", totalview::typeName< RestartFlags >().c_str(), &m_restart_flags );
  TV_ttf_add_row( "m_plotLevel", totalview::typeName< PlotLevel >().c_str(), &m_plotLevel );
  TV_ttf_add_row( "m_inputFlag", totalview::typeName< InputFlags >().c_str(), &m_inputFlag );
  TV_ttf_add_row( "m_description", totalview::typeName< string >().c_str(), &m_description );
  size_t junk = m_registeringObjects.size();
  TV_ttf_add_row( "m_registeringObjects",
                  totalview::format< string, size_t >( 1, &junk ).c_str(),
                  m_registeringObjects.data() );

  return 0;
}
#endif


}
} /* namespace geosx */

#if defined(USE_TOTALVIEW_OUTPUT)
/**
 * @brief Global function correlated with WrapperBase to be called by Totalview when displaying
 *        a WrapperBase as a VieWrapper<T>
 * @param wrapper A pointer to the wrapper that will be displayed.
 * @return 0
 */
int TV_ttf_display_type( const geosx::dataRepository::WrapperBase * wrapper )
{
  if( wrapper!=nullptr )
  {
    //std::cout<<"displaying WrapperBase "<<wrapper->getName()<<" as "<<wrapper->totalviewTypeName()<<std::endl;
// keep this and try to make it work later on.
//    rval = TV_ttf_add_row( "casted_this", wrapper->totalviewTypeName().c_str(), wrapper );
    wrapper->setTotalviewDisplay();
  }
  return 0;
}
#endif
