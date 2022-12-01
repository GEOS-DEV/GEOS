/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file ToParticleRelation.cpp
 */

#include "ToParticleRelation.hpp"

#include "ParticleManager.hpp"

namespace geosx
{


void erase( OrderedVariableToManyParticleRelation & relation,
            localIndex const firstIndex,
            localIndex const er,
            localIndex const esr,
            localIndex const ei )
{
  for( localIndex a=relation.m_toParticleRegion.sizeOfArray( firstIndex )-1; a>=0; --a )
  {
    if( er==relation.m_toParticleRegion[firstIndex][a] &&
        esr==relation.m_toParticleSubRegion[firstIndex][a] &&
        ei==relation.m_toParticleIndex[firstIndex][a] )
    {
      relation.m_numParticles[firstIndex]--;
      relation.m_toParticleRegion.eraseFromArray( firstIndex, a );
      relation.m_toParticleSubRegion.eraseFromArray( firstIndex, a );
      relation.m_toParticleIndex.eraseFromArray( firstIndex, a );
    }
  }
}

void insert( OrderedVariableToManyParticleRelation & relation,
             localIndex const firstIndex,
             localIndex const er,
             localIndex const esr,
             localIndex const ei )
{
  bool alreadyPresent = false;
  for( localIndex a=0; a<relation.m_toParticleRegion.sizeOfArray( firstIndex ); ++a )
  {
    if( er==relation.m_toParticleRegion[firstIndex][a] &&
        esr==relation.m_toParticleSubRegion[firstIndex][a] &&
        ei==relation.m_toParticleIndex[firstIndex][a] )
    {
      alreadyPresent = true;
    }
  }
  if( !alreadyPresent )
  {
    relation.m_numParticles[firstIndex]++;
    relation.m_toParticleRegion.emplaceBack( firstIndex, er );
    relation.m_toParticleSubRegion.emplaceBack( firstIndex, esr );
    relation.m_toParticleIndex.emplaceBack( firstIndex, ei );
  }
}

} /* namespace geosx */
