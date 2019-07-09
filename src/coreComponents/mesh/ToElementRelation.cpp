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

/*
 * ToElementRelation.cpp
 *
 *  Created on: Jun 6, 2018
 *      Author: settgast
 */

#include "ToElementRelation.hpp"

#include "ElementRegionManager.hpp"

namespace geosx
{


void erase( OrderedVariableToManyElementRelation & relation,
            localIndex const firstIndex,
            localIndex const er,
            localIndex const esr,
            localIndex const ei )
{
  for( localIndex a=relation.m_toElementRegion.sizeOfArray(firstIndex)-1 ; a>=0 ; --a )
  {
    if( er==relation.m_toElementRegion[firstIndex][a] &&
        esr==relation.m_toElementSubRegion[firstIndex][a] &&
        ei==relation.m_toElementIndex[firstIndex][a] )
    {
      relation.m_toElementRegion.eraseFromArray(firstIndex,a);
      relation.m_toElementSubRegion.eraseFromArray( firstIndex,a);
      relation.m_toElementIndex.eraseFromArray(firstIndex,a);
    }
  }
}

void insert( OrderedVariableToManyElementRelation & relation,
             localIndex const firstIndex,
             localIndex const er,
             localIndex const esr,
             localIndex const ei )
{
  bool alreadyPresent = false;
  for( localIndex a=0 ; a<relation.m_toElementRegion.sizeOfArray(firstIndex) ; ++a )
  {
    if( er==relation.m_toElementRegion[firstIndex][a] &&
        esr==relation.m_toElementSubRegion[firstIndex][a] &&
        ei==relation.m_toElementIndex[firstIndex][a] )
    {
      alreadyPresent = true;
    }
  }
  if( !alreadyPresent )
  {
    relation.m_toElementRegion.appendToArray( firstIndex, er );
    relation.m_toElementSubRegion.appendToArray( firstIndex, esr );
    relation.m_toElementIndex.appendToArray( firstIndex, ei );
  }
}

} /* namespace geosx */
