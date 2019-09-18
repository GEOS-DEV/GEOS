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

#include "MPI_Communications/CommunicationTools.hpp"

#include  "WellGeneratorBase.hpp"

#include "PerforationData.hpp"
#include "Perforation.hpp"

namespace geosx
{

using namespace dataRepository;

WellGeneratorBase::WellGeneratorBase( string const & name, Group * const parent ):
  MeshGeneratorBase( name, parent ),
  m_numElemsPerSegment(0),
  m_crossSectionArea(0),
  m_wellRegionName(""),
  m_wellControlsName(""),
  m_meshBodyName(""),
  m_numElems(0),
  m_numNodesPerElem(2),
  m_numNodes(0),
  m_numPerforations(0),
  m_nDims(3),
  m_polylineHeadNodeId(-1)
{
  registerWrapper(keys::nElems, &m_numElemsPerSegment, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setSizedFromParent(0)->
    setDescription("number of well elements per polyline segment");

  registerWrapper(keys::wellRegionName, &m_wellRegionName, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setSizedFromParent(0)->
    setDescription("name of the well element region");

  registerWrapper(keys::wellControlsName, &m_wellControlsName, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setSizedFromParent(0)->
    setDescription("name of the set of constraints associated with this well");

  registerWrapper(keys::meshBodyName, &m_meshBodyName, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setSizedFromParent(0)->
    setDescription("name of the reservoir mesh associated with this well");
}

WellGeneratorBase::~WellGeneratorBase()
{
}

void WellGeneratorBase::PostProcessInput()
{
}

Group * WellGeneratorBase::CreateChild( string const & childKey, string const & childName )
{
  if ( childKey == keys::perforation )
  {
    ++m_numPerforations;

    // keep track of the perforations that have been added
    m_perforationList.push_back( childName );

    return RegisterGroup<Perforation>( childName );
  }
  else
  {
    GEOS_ERROR( "Unrecognized node: " << childKey );
  }
  return nullptr;
}

void WellGeneratorBase::DebugWellGeometry() const 
{
  if (CommunicationTools::MPI_Rank( MPI_COMM_GEOSX ) != 0)
  {
    return;
  } 

  std::cout << std::endl;
  std::cout << "++++++++++++++++++++++++++" << std::endl;
  std::cout << "InternalWellGenerator = " << getName() << std::endl;
  std::cout << "MPI rank = " << CommunicationTools::MPI_Rank( MPI_COMM_GEOSX ) << std::endl;
  std::cout << "Number of well elements = " << m_numElems << std::endl;
  
  for (globalIndex iwelem = 0; iwelem < m_numElems; ++iwelem)
  {
    std::cout << "m_elemCenterCoords[" << iwelem << "] = " << m_elemCenterCoords[iwelem] 
              << std::endl;
    std::cout << "m_nextElemId[" << iwelem << "] = " << m_nextElemId[iwelem] 
              << std::endl;
    std::cout << "m_prevElemId[" << iwelem << "] = " << m_prevElemId[iwelem][0] 
              << std::endl;
    for (globalIndex inode = 0; inode < m_numNodesPerElem; ++inode)
    {
      std::cout << "m_elemToNodesMap[" << iwelem << "][" << inode << "] = " << m_elemToNodesMap[iwelem][inode]
                << std::endl;
    }
  }

  std::cout << "Number of well nodes = " << m_numNodes << std::endl;
  
  for (globalIndex inode = 0; inode < m_numNodes; ++inode)
  {
    std::cout << "m_nodeCoords[" << inode << "] = " << m_nodeCoords[inode] 
              << std::endl;
    std::cout << "m_nodeDistFromHead[" << inode << "] = " << m_nodeDistFromHead[inode]
              << std::endl;
  }

  std::cout << "Number of perforations = " << m_numPerforations << std::endl;

  for (globalIndex iperf = 0; iperf < m_numPerforations; ++iperf) 
  {
    std::cout << "m_perfCoords[" << iperf << "] = " << m_perfCoords[iperf] 
              << std::endl;
    std::cout << "m_perfTrans[" << iperf << "] = " << m_perfTrans[iperf] 
              << std::endl;
    std::cout << "m_perfElemId[" << iperf << "] = " << m_perfElemId[iperf] 
              << std::endl;
  }

}

} // namespace
