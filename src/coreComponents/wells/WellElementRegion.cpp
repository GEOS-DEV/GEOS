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

/**
 * @file WellElementRegion.cpp
 */

#include "MPI_Communications/CommunicationTools.hpp"

#include "WellElementRegion.hpp"
#include "WellElementSubRegion.hpp"

namespace geosx
{
using namespace dataRepository;

WellElementRegion::WellElementRegion( string const & name, Group * const parent ):
  ElementRegionBase( name, parent ),
  m_subRegionName(name+"uniqueSubRegion"),
  m_wellControlsName(""),
  m_wellGeneratorName("")
{
  registerWrapper( viewKeyStruct::wellControlsString,  &m_wellControlsName,  false );
  registerWrapper( viewKeyStruct::wellGeneratorString, &m_wellGeneratorName, false );

  this->GetGroup( viewKeyStruct::elementSubRegions )
      ->RegisterGroup<WellElementSubRegion>( m_subRegionName );

}

WellElementRegion::~WellElementRegion()
{}


void WellElementRegion::GenerateWell( MeshLevel & mesh, 
                                      InternalWellGenerator const & wellGeometry,
                                      globalIndex nodeOffsetGlobal,
                                      globalIndex elemOffsetGlobal )
{
  // get the (unique) subregion
  WellElementSubRegion * const 
  subRegion = this->GetGroup( ElementRegionBase::viewKeyStruct::elementSubRegions )
                  ->GetGroup<WellElementSubRegion>( m_subRegionName );

  GEOS_ERROR_IF( subRegion == nullptr, 
                 "Well subRegion " << this->m_subRegionName << " not found in well region " << getName() );
  subRegion->SetWellControlsName( m_wellControlsName );

  PerforationData * const perforationData = subRegion->GetPerforationData();
  perforationData->SetNumPerforationsGlobal( wellGeometry.GetNumPerforations() );

  globalIndex const numElemsGlobal        = wellGeometry.GetNumElements();
  globalIndex const numPerforationsGlobal = wellGeometry.GetNumPerforations();


  // 1) select the local perforations based on connectivity to the local reservoir elements
  perforationData->ConnectToMeshElements( mesh, wellGeometry );

  globalIndex const matchedPerforations = CommunicationTools::Sum( perforationData->size() );
  GEOS_ERROR_IF( matchedPerforations != numPerforationsGlobal, 
                 "Invalid mapping perforation-to-element in well " << this->getName() );


  // 2) classify well elements based on connectivity to local mesh partition
  array1d<integer> elemStatusGlobal( numElemsGlobal );
  elemStatusGlobal = WellElementSubRegion::WellElemStatus::UNOWNED;

  array1d<globalIndex const> const & perfElemIdGlobal = wellGeometry.GetPerfElemIndex();

  for (localIndex iperfGlobal = 0; iperfGlobal < numPerforationsGlobal; ++iperfGlobal)
  {
    globalIndex const iwelemGlobal = perfElemIdGlobal[iperfGlobal];

    if (perforationData->m_globalToLocalMap.count( iperfGlobal ) > 0)
    {
      elemStatusGlobal[iwelemGlobal] |= WellElementSubRegion::WellElemStatus::LOCAL;
    }
    else
    {
      elemStatusGlobal[iwelemGlobal] |= WellElementSubRegion::WellElemStatus::REMOTE;
    }
  }


  // 3) select the local well elements and mark boundary nodes (for ghosting)
  subRegion->Generate( mesh, 
                       wellGeometry, 
                       elemStatusGlobal, 
                       nodeOffsetGlobal, 
                       elemOffsetGlobal );


  // 4) find out which rank is the owner of the top segment
  localIndex const refElemIdLocal = subRegion->GetTopWellElementIndex();

  array1d<localIndex> allRankTopElem;
  CommunicationTools::allGather( refElemIdLocal, allRankTopElem );
  int topRank = -1;
  for (int irank = 0; irank < allRankTopElem.size(); ++irank)
  {
    if (allRankTopElem[irank] >= 0)
    {
      GEOS_ASSERT( topRank < 0 );
      topRank = irank;
    }
  }
  GEOS_ASSERT( topRank >= 0 );
  subRegion->SetTopRank( topRank );


  // 5) construct the local perforation to well element map
  perforationData->ConnectToWellElements( wellGeometry, 
                                          subRegion->m_globalToLocalMap, 
                                          elemOffsetGlobal );  

}

REGISTER_CATALOG_ENTRY( ObjectManagerBase, WellElementRegion, std::string const &, Group * const )

} /* namespace geosx */
