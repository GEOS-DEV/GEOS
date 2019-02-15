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

#include "WellElementSubRegion.hpp"

#include "../wells/WellElementManager.hpp"
#include "../wells/WellElement.hpp"
#include "../wells/Well.hpp"

namespace geosx
{

WellElementSubRegion::WellElementSubRegion( string const & name, ManagedGroup * const parent ):
  ElementSubRegionBase( name, parent )
{
  RegisterViewWrapper( viewKeyStruct::wellElementIndexString, &m_wellElementIndex, false );
}


WellElementSubRegion::~WellElementSubRegion()
{}

WellElement const * WellElementSubRegion::getWellElement( localIndex iwelem ) const
{
  Well const * parent = getParent()->group_cast<Well const *>();
  WellElementManager const * wellElementManager
    = parent->GetGroup<WellElementManager>( Well::groupKeyStruct::wellElementsString );
  WellElement const * wellElement
    = wellElementManager->getWellElement( m_wellElementIndex[iwelem] );
  return wellElement;
}

WellElement * WellElementSubRegion::getWellElement( localIndex iwelem )
{
  Well * parent = getParent()->group_cast<Well *>();
  WellElementManager * wellElementManager
    = parent->GetGroup<WellElementManager>( Well::groupKeyStruct::wellElementsString );
  WellElement * wellElement
    = wellElementManager->getWellElement( m_wellElementIndex[iwelem] );
  return wellElement;
}

void WellElementSubRegion::InitializePreSubGroups( ManagedGroup * const problemManager )
{
  // todo later: MPI partitioning
  
  // for now, don't bother with that
  WellElementManager const * wellElementManager
    = getParent()->GetGroup<WellElementManager>( Well::groupKeyStruct::wellElementsString );

  // set the size to the number of global elements
  resize( wellElementManager->numWellElementsGlobal() );
  // dummy map from local to global
  for (localIndex iwelem = 0; iwelem < size(); ++iwelem)
  {
    m_wellElementIndex[iwelem] = iwelem;
  }
}

  
}
