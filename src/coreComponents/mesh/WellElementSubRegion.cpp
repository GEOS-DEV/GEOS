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
  RegisterViewWrapper( viewKeyStruct::nextWellElementIndexString, &m_nextWellElementIndex, false );

  RegisterViewWrapper( viewKeyStruct::gravityDepthString, &m_gravityDepth, false );
  RegisterViewWrapper( ElementSubRegionBase::viewKeyStruct::elementCenterString, &m_elementCenter, false );
  RegisterViewWrapper( ElementSubRegionBase::viewKeyStruct::elementVolumeString, &m_elementVolume, false );
}


WellElementSubRegion::~WellElementSubRegion()
{}

void WellElementSubRegion::InitializePreSubGroups( ManagedGroup * const problemManager )
{
  // TODO: this function is temporary and needs to be improved
  //       depending on the parallelization strategy 
   
  Well const * const parent = getParent()->group_cast<Well *>();

  WellElementManager const * const wellElementManager = 
    parent->GetGroup<WellElementManager>( Well::groupKeyStruct::wellElementsString );

  // initialize the attribute of the well elements
  for (localIndex iwelem = 0; iwelem < size(); ++iwelem)
  {
 
    // get the current well element
    WellElement const * const wellElement = 
      wellElementManager->getWellElement( iwelem );

    // initialize the next well elem index
    m_nextWellElementIndex[iwelem] = -2;

    if (iwelem == 0) // well head, set next element to -1
    {
      m_nextWellElementIndex[iwelem] = -1;
    }
    else // else find the next well element by name
    {
      // get the next well element name
      string const nextWellElementName = wellElement->getNextWellElementName();

      // find the index of the next well element and save it
      for (localIndex iwelemNext = 0; iwelemNext < size(); ++iwelemNext)
      {
        if (iwelemNext == iwelem)
        {
          continue;
        }

        // get a possible next well elem
        WellElement const * const nextWellElement = 
          wellElementManager->getWellElement( iwelemNext );

        // if the names match, save the index
        if (nextWellElement->getName() == nextWellElementName)
        {
          m_nextWellElementIndex[iwelem] = iwelemNext;
          break;
        }
      }
    }

    // error message if the well element is not found
    if (m_nextWellElementIndex[iwelem] == -2)
    {
      GEOS_ERROR("Invalid next well element name: " << wellElement->getNextWellElementName()
                 << " for well element " << wellElement->getName() );
    }    
    
    // save volume and location
    m_elementVolume[iwelem] = 1.; 
    m_elementCenter[iwelem] = wellElement->getLocation();
  }
}

void WellElementSubRegion::InitializePostInitialConditions_PreSubGroups( ManagedGroup * const problemManager )
{
}
  
}
