/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/** @file */

#include "ViewWrapperBase.hpp"

#include "ManagedGroup.hpp"
#include "RestartFlags.hpp"


namespace geosx
{
namespace dataRepository
{


ViewWrapperBase::ViewWrapperBase( std::string const & name,
                                  ManagedGroup * const parent ):
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


ViewWrapperBase::~ViewWrapperBase()
{}


ViewWrapperBase::ViewWrapperBase( ViewWrapperBase && source ):
  m_name( std::move( source.m_name ) ),
  m_parent( source.m_parent ),
  m_sizedFromParent( source.m_sizedFromParent ),
  m_restart_flags( source.m_restart_flags )
#ifdef GEOSX_USE_ATK
  , m_sidreView( source.m_sidreView )
#endif

{}

void ViewWrapperBase::resize()
{
  resize( m_parent->size());
}

void ViewWrapperBase::CopyWrapperAttributes( ViewWrapperBase const & source )
{
  m_name = source.m_name;
  m_sizedFromParent = source.m_sizedFromParent;
  m_restart_flags = source.m_restart_flags;
}

#if defined(USE_TOTALVIEW_OUTPUT)
int ViewWrapperBase::setTotalviewDisplay() const
{
  //std::cout<<"exectuing ViewWrapperBase::setTotalviewDisplay()"<<std::endl;
//  TV_ttf_add_row("TYPE", TV_ttf_type_ascii_string, type.c_str() );
  TV_ttf_add_row( "m_name", totalview::typeName< string >().c_str(), &m_name );
  TV_ttf_add_row( "m_parent", totalview::typeName< ManagedGroup >().c_str(), m_parent );
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
 * @brief Global function correlated with ViewWrapperBase to be called by Totalview when displaying
 *        a ViewWrapperBase as a VieWrapper<T>
 * @param wrapper A pointer to the wrapper that will be displayed.
 * @return 0
 */
int TV_ttf_display_type( const geosx::dataRepository::ViewWrapperBase * wrapper )
{
  if( wrapper!=nullptr )
  {
    //std::cout<<"displaying ViewWrapperBase "<<wrapper->getName()<<" as "<<wrapper->totalviewTypeName()<<std::endl;
// keep this and try to make it work later on.
//    rval = TV_ttf_add_row( "casted_this", wrapper->totalviewTypeName().c_str(), wrapper );
    wrapper->setTotalviewDisplay();
  }
  return 0;
}
#endif
