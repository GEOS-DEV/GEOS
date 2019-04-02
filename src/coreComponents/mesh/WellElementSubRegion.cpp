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
  RegisterViewWrapper( viewKeyStruct::wellElementVolumeString, &m_wellElementVolume, false );
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
  // @Francois: do not resize here!!!! otherwise you resize phaseNames, etc
  
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
      m_nextWellElementIndex[iwelem] = iwelem - 1;
    }
    m_wellElementVolume[iwelem] = 1.;
  }
}

void WellElementSubRegion::InitializePostInitialConditions_PreSubGroups( ManagedGroup * const problemManager )
{
  R1Tensor const & gravity = getParent()->group_cast<Well *>()->getGravityVector();
  arrayView1d<real64> & gravDepth = getReference<array1d<real64>>( viewKeyStruct::gravityDepthString );

  for (localIndex iwelem = 0; iwelem < size(); ++iwelem)
  {
    WellElement const * wellElement = getWellElement( iwelem );
    gravDepth[iwelem] = Dot( wellElement->getLocation(), gravity );
  }
}
  
}
