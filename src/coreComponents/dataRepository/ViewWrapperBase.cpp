/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
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

/*
 * DataObjectBase.cpp
 *
 *  Created on: Jun 17, 2016
 *      Author: rrsettgast
 */

#include "ViewWrapperBase.hpp"

#include "ManagedGroup.hpp"
#include "RestartFlags.hpp"


namespace geosx
{
namespace dataRepository
{


ViewWrapperBase::ViewWrapperBase( std::string const & name,
                                  ManagedGroup * const parent):
  m_name(name),
  m_parent(parent),
  m_sizedFromParent(1),
  m_restart_flags(RestartFlags::WRITE_AND_READ),
  m_plotLevel(PlotLevel::LEVEL_3),
  m_inputFlag(InputFlags::INVALID),
  m_description()
#ifdef GEOSX_USE_ATK
  ,m_sidreView(nullptr)
#endif
{
#ifdef GEOSX_USE_ATK
  GEOS_ERROR_IF(parent==nullptr,"parameter WrapperCollection * const parent must not be nullptr");

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
#ifdef GEOSX_USE_ATK
  ,m_sidreView( source.m_sidreView )
#endif

{}

void ViewWrapperBase::resize()
{
  resize(m_parent->size());
}


}
} /* namespace geosx */
