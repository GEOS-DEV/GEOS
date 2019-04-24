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
  RegisterViewWrapper( viewKeyStruct::nextWellElementIndexString, &m_nextWellElementIndex, false );

  RegisterViewWrapper( viewKeyStruct::gravityDepthString, &m_gravityDepth, false );
  RegisterViewWrapper( ElementSubRegionBase::viewKeyStruct::elementCenterString, &m_elementCenter, false );
  RegisterViewWrapper( ElementSubRegionBase::viewKeyStruct::elementVolumeString, &m_elementVolume, false );
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
  // TODO: make this work in parallel
  // TODO: this function is temporary and will be rewritten entirely
    
  // dummy map from local to global
  for (localIndex iwelem = 0; iwelem < size(); ++iwelem)
  {
    m_wellElementIndex[iwelem] = iwelem;

    if (iwelem == 0)
    {
      m_nextWellElementIndex[iwelem] = -1;
    }
    else
    {
      string const nextWellElementName = getWellElement( iwelem )->getNextWellElementName();
      // this is a temporary hack
      for (localIndex iwelemNext = 0; iwelemNext < size(); ++iwelemNext)
      {
        if (getWellElement( iwelemNext )->getName() == nextWellElementName)
        {
          m_nextWellElementIndex[iwelem] = iwelemNext;
          break;
        }
      }
    }


    m_elementVolume[iwelem] = 1.;
    m_elementCenter[iwelem] = getWellElement( iwelem )->getLocation();
  }
}

void WellElementSubRegion::InitializePostInitialConditions_PreSubGroups( ManagedGroup * const problemManager )
{
}
  
}
