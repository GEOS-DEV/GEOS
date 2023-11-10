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

/*
 * @file InternalWellGenerator.cpp
 *
 */

#include "InternalWellGenerator.hpp"

#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "mesh/Perforation.hpp"
#include "mesh/generators/LineBlockABC.hpp"
#include "LvArray/src/genericTensorOps.hpp"
#include "physicsSolvers/fluidFlow/wells/WellControls.hpp"

namespace geos
{
using namespace dataRepository;

InternalWellGenerator::InternalWellGenerator( string const & name, Group * const parent ):
  WellGeneratorBase( name, parent ),
  m_numElemsPerSegment( 0 ),
  m_minSegmentLength( 1e-2 ),
  m_minElemLength( 1e-3 ),
  m_radius( 0 ),
  m_wellRegionName( "" ),
  m_wellControlsName( "" ),
  m_numElems( 0 ),
  m_numNodesPerElem( 2 ),
  m_numNodes( 0 ),
  m_nDims( 3 ),
  m_polylineHeadNodeId( -1 )
{
  enableLogLevelInput();

  registerWrapper( viewKeyStruct::polylineNodeCoordsString(), &m_polyNodeCoords ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Physical coordinates of the well polyline nodes" );

  registerWrapper( viewKeyStruct::polylineSegmentConnString(), &m_segmentToPolyNodeMap ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Connectivity of the polyline segments" );

  registerWrapper( viewKeyStruct::radiusString(), &m_radius ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Radius of the well [m]" );

  registerWrapper( viewKeyStruct::numElementsPerSegmentString(), &m_numElemsPerSegment ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Number of well elements per polyline segment" );

  registerWrapper( viewKeyStruct::minSegmentLengthString(), &m_minSegmentLength ).
    setApplyDefaultValue( 1e-2 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Minimum length of a well segment [m]" );

  registerWrapper( viewKeyStruct::minElementLengthString(), &m_minElemLength ).
    setApplyDefaultValue( 1e-3 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Minimum length of a well element, computed as (segment length / number of elements per segment ) [m]" );

  registerWrapper( viewKeyStruct::wellRegionNameString(), &m_wellRegionName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Name of the well element region" );

  registerWrapper( viewKeyStruct::wellControlsNameString(), &m_wellControlsName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Name of the set of constraints associated with this well" );
}

void InternalWellGenerator::postProcessInput()
{
  GEOS_THROW_IF( m_polyNodeCoords.size( 1 ) != m_nDims,
                 "InternalWell " << getWrapperDataContext( viewKeyStruct::polylineNodeCoordsString() ) <<
                 ": Invalid number of physical coordinates.",
                 InputError );

  GEOS_THROW_IF( m_segmentToPolyNodeMap.size( 1 ) != 2,
                 "InternalWell " << getWrapperDataContext( viewKeyStruct::polylineSegmentConnString() ) <<
                 ": Invalid size.",
                 InputError );

  GEOS_THROW_IF( m_polyNodeCoords.size( 0 )-1 != m_segmentToPolyNodeMap.size( 0 ),
                 "Incompatible sizes of " << getWrapperDataContext( viewKeyStruct::polylineNodeCoordsString() ) <<
                 " and " << getWrapperDataContext( viewKeyStruct::polylineSegmentConnString() ),
                 InputError );

  GEOS_THROW_IF( m_radius <= 0,
                 "InternalWell " << getWrapperDataContext( viewKeyStruct::radiusString() ) <<
                 ": Radius value must be greater that 0.",
                 InputError );

  GEOS_THROW_IF( m_wellRegionName.empty(),
                 "InternalWell " << getWrapperDataContext( viewKeyStruct::wellRegionNameString() ) <<
                 ": Empty well region name.",
                 InputError );

  GEOS_THROW_IF( m_wellControlsName.empty(),
                 "InternalWell " << getWrapperDataContext( viewKeyStruct::wellControlsNameString() ) <<
                 ": Empty well constraint name.",
                 InputError );

  // TODO: add more checks here
  // TODO: check that the connectivity of the well is valid
  // TODO: check that with no branching we can go from top to bottom and touch all the elements
}

void InternalWellGenerator::generateWellGeometry( )
{
  // count the number of well elements to create
  m_numElems = m_numElemsPerSegment * m_segmentToPolyNodeMap.size( 0 );
  m_numNodes = m_numElems + 1;

  // resize the well element, node, and perforation arrays
  m_elemCenterCoords.resize( m_numElems, 3 );
  m_nextElemId.resize( m_numElems );
  m_prevElemId.resize( m_numElems );
  m_elemToNodesMap.resizeDimension< 0 >( m_numElems );
  m_elemToNodesMap.resizeDimension< 1 >( m_numNodesPerElem );
  m_elemVolume.resize( m_numElems );

  m_nodeDistFromHead.resize( m_numNodes );
  m_nodeCoords.resize( m_numNodes, 3 );

  m_perfCoords.resize( m_numPerforations, 3 );
  m_perfDistFromHead.resize( m_numPerforations );
  m_perfTransmissibility.resize( m_numPerforations );
  m_perfElemId.resize( m_numPerforations );

  // construct a reverse map from the polyline nodes to the segments
  constructPolylineNodeToSegmentMap();

  // detect the head polyline node based on depth
  findPolylineHeadNodeIndex();

  // compute the location and distance from top of the well elements
  discretizePolyline();

  // map the perforations to the well elements
  connectPerforationsToWellElements();

  // make sure that the perforation locations are valid
  checkPerforationLocationsValidity();

  if( getLogLevel() >= 1 )
  {
    debugWellGeometry();
  }

}

void InternalWellGenerator::constructPolylineNodeToSegmentMap()
{
  m_polyNodeToSegmentMap.resize( m_polyNodeCoords.size( 0 ) );

  bool foundSmallElem = false;

  // loop over the segments
  for( globalIndex iseg = 0; iseg < m_segmentToPolyNodeMap.size( 0 ); ++iseg )
  {
    globalIndex const ipolyNode_a = m_segmentToPolyNodeMap[iseg][0];
    globalIndex const ipolyNode_b = m_segmentToPolyNodeMap[iseg][1];

    real64 vSeg[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( m_polyNodeCoords[ipolyNode_a] );
    LvArray::tensorOps::subtract< 3 >( vSeg, m_polyNodeCoords[ipolyNode_b] );

    real64 const segmentLength = LvArray::tensorOps::l2Norm< 3 >( vSeg );

    // various checks and warnings on the segment and element length
    GEOS_THROW_IF( segmentLength < m_minSegmentLength,
                   getDataContext() << ": Error in the topology of the well. We detected a" <<
                   " polyline segment measuring less than " << m_minSegmentLength << "m. \n" <<
                   "You can change the minimum segment length using the field " << viewKeyStruct::minSegmentLengthString(),
                   InputError );

    GEOS_THROW_IF( m_polyNodeCoords[ipolyNode_a][2] < m_polyNodeCoords[ipolyNode_b][2],
                   getDataContext() << ": Error in the topology of the well. In the polyline"<<
                   ", each segment must be going down. \n" <<
                   "This is not the case between polyline nodes " << m_polyNodeCoords[ipolyNode_a] << " and " << m_polyNodeCoords[ipolyNode_b],
                   InputError );

    if( segmentLength / m_numElemsPerSegment < m_minElemLength )
    {
      foundSmallElem = true;
    }

    // map the polyline node ids to the polyline segment ids
    m_polyNodeToSegmentMap[ipolyNode_a].insert( iseg );
    m_polyNodeToSegmentMap[ipolyNode_b].insert( iseg );
  }

  if( foundSmallElem )
  {
    GEOS_LOG_RANK_0( "\nWarning for well " << getDataContext() << ": the chosen number of well" <<
                     " elements per polyline segment (" << m_numElemsPerSegment << ") leads to well" <<
                     " elements measuring less than " << m_minElemLength << "m in the topology.\n" <<
                     "The simulation can proceed like that, but small well elements may cause" <<
                     " numerical issues, so it is something to keep an eye on.\n" <<
                     "You can get rid of this message by changing the field " << viewKeyStruct::minElementLengthString() );
  }
}

void InternalWellGenerator::findPolylineHeadNodeIndex()
{
  // we assume that the first *polyline segment* in the XML input is the well head segment
  globalIndex const wellHeadSegId = 0;

  // get the corresponding node indices
  globalIndex const ipolyNode_a = m_segmentToPolyNodeMap[wellHeadSegId][0];
  globalIndex const ipolyNode_b = m_segmentToPolyNodeMap[wellHeadSegId][1];

  // we determine which node is the head node based on depth
  // therefore here we throw an error if the well head segment is horizontal
  GEOS_THROW_IF( !(m_polyNodeCoords[ipolyNode_a][2] < m_polyNodeCoords[ipolyNode_b][2])
                 && !(m_polyNodeCoords[ipolyNode_a][2] > m_polyNodeCoords[ipolyNode_b][2]),
                 getDataContext() << ": The head polyline segment cannot be horizontal in the well" <<
                 " since we use depth to determine which of its nodes is to head node of the well.\n" <<
                 "If you are trying to set up a horizontal well, please simply add a non-horizontal" <<
                 " segment at the top of the well, and this error will go away.",
                 InputError );

  // detect the top node, assuming z oriented upwards
  m_polylineHeadNodeId =
    ( m_polyNodeCoords[ipolyNode_a][2] > m_polyNodeCoords[ipolyNode_b][2] )
    ? ipolyNode_a
    : ipolyNode_b;

  real64 const headZcoord = m_polyNodeCoords[m_polylineHeadNodeId][2];
  for( globalIndex inode = 0; inode < m_polyNodeCoords.size( 0 ); ++inode )
  {
    if( inode == m_polylineHeadNodeId )
    {
      continue;
    }
    real64 const currentZcoord = m_polyNodeCoords[inode][2];

    GEOS_THROW_IF( !(currentZcoord < headZcoord),
                   getDataContext() << ": Error in the topology since we found a well node that" <<
                   " is above the head node",
                   InputError );
  }
}

void InternalWellGenerator::discretizePolyline()
{
  // initialize well elements and node ids
  globalIndex ipolyNodeTop  = m_polylineHeadNodeId;
  globalIndex iwelemCurrent = 0;
  globalIndex isegCurrent   = 0;

  // set the location of the first well node and distance from well head
  LvArray::tensorOps::copy< 3 >( m_nodeCoords[0], m_polyNodeCoords[ipolyNodeTop] );
  m_nodeDistFromHead[0] = 0.0;

  // note: this part of the code does not support well branching
  // TODO: check that there is only one branch
  // TODO: read wells with branching (already supported elsewhere in the code)

  // go through the well from top to bottom
  for( globalIndex is = 0; is < m_segmentToPolyNodeMap.size( 0 ); ++is )
  {

    GEOS_THROW_IF( isegCurrent == -1,
                   getWrapperDataContext( viewKeyStruct::polylineSegmentConnString() ) <<
                   ": Invalid map.",
                   InputError );

    globalIndex const ipolyNodeBottom = ( ipolyNodeTop == m_segmentToPolyNodeMap[isegCurrent][0] )
                                      ? m_segmentToPolyNodeMap[isegCurrent][1]
                                      : m_segmentToPolyNodeMap[isegCurrent][0];

    real64 vPoly[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( m_polyNodeCoords[ipolyNodeBottom] );
    LvArray::tensorOps::subtract< 3 >( vPoly, m_polyNodeCoords[ipolyNodeTop] );

    // add the well elements and well nodes corresponding to this polyline segment
    for( localIndex iw = 0; iw < m_numElemsPerSegment; ++iw )
    {

      // 1) set the element location, connectivity
      real64 const scaleCenter = (iw + 0.5) / static_cast< real64 >(m_numElemsPerSegment);
      LvArray::tensorOps::copy< 3 >( m_elemCenterCoords[iwelemCurrent], vPoly );
      LvArray::tensorOps::scale< 3 >( m_elemCenterCoords[iwelemCurrent], scaleCenter );
      LvArray::tensorOps::add< 3 >( m_elemCenterCoords[iwelemCurrent], m_polyNodeCoords[ipolyNodeTop] );

      GEOS_THROW_IF( iwelemCurrent >= m_numElems,
                     getDataContext() << ": Invalid well topology",
                     InputError );

      globalIndex const iwelemTop    = iwelemCurrent - 1;
      globalIndex const iwelemBottom = iwelemCurrent + 1;
      m_nextElemId[iwelemCurrent] = iwelemTop;
      m_prevElemId[iwelemCurrent].resize( 1 );
      m_prevElemId[iwelemCurrent][0] = (iwelemBottom < m_numElems)
                                     ? iwelemBottom
                                     : -1;

      // 2) set the node properties
      globalIndex const iwellNodeTop    = iwelemCurrent;
      globalIndex const iwellNodeBottom = iwelemCurrent+1;
      m_elemToNodesMap[iwelemCurrent][LineBlockABC::NodeLocation::TOP]    = iwellNodeTop;
      m_elemToNodesMap[iwelemCurrent][LineBlockABC::NodeLocation::BOTTOM] = iwellNodeBottom;

      real64 const scaleBottom = (iw + 1.0) / m_numElemsPerSegment;
      LvArray::tensorOps::copy< 3 >( m_nodeCoords[iwellNodeBottom], vPoly );
      LvArray::tensorOps::scale< 3 >( m_nodeCoords[iwellNodeBottom], scaleBottom );
      LvArray::tensorOps::add< 3 >( m_nodeCoords[iwellNodeBottom], m_polyNodeCoords[ipolyNodeTop] );

      // 3) increment the distance from the well head to the bottom of current element
      real64 vWellElem[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( m_nodeCoords[iwellNodeBottom] );
      LvArray::tensorOps::subtract< 3 >( vWellElem, m_nodeCoords[iwellNodeTop] );
      m_nodeDistFromHead[iwellNodeBottom]  = LvArray::tensorOps::l2Norm< 3 >( vWellElem );
      m_nodeDistFromHead[iwellNodeBottom] += m_nodeDistFromHead[iwellNodeTop];

      // 4) set element volume
      m_elemVolume[iwelemCurrent] = LvArray::tensorOps::l2Norm< 3 >( vWellElem ) * M_PI * m_radius * m_radius;

      // 4) increment the element counter
      ++iwelemCurrent;
    }

    // then consider the next polyline segment
    ipolyNodeTop = ipolyNodeBottom;
    isegCurrent  = getNextSegmentIndex( isegCurrent, ipolyNodeTop );
  }

  // set the previous index for the bottom segment
  m_prevElemId[iwelemCurrent-1].resize( 1 );
  m_prevElemId[iwelemCurrent-1][0] = -1;

}


void InternalWellGenerator::connectPerforationsToWellElements()
{

  // assign a well element to each perforation
  for( globalIndex iperf = 0; iperf < m_numPerforations; ++iperf )
  {

    // get the perforation and its properties
    Perforation const & perf = this->getGroup< Perforation >( m_perforationList[iperf] );
    m_perfDistFromHead[iperf]  = perf.getDistanceFromWellHead();
    m_perfTransmissibility[iperf] = perf.getWellTransmissibility();

    // search in all the elements of this well between head and bottom
    globalIndex iwelemTop    = 0;
    globalIndex iwelemBottom = m_numElems - 1;

    // check the validity of the perforation before starting
    real64 const wellLength = m_nodeDistFromHead[m_elemToNodesMap[iwelemBottom][LineBlockABC::NodeLocation::BOTTOM]];

    GEOS_THROW_IF( m_perfDistFromHead[iperf] > wellLength,
                   perf.getWrapperDataContext( Perforation::viewKeyStruct::distanceFromHeadString() ) <<
                   ": Distance from well perforation to head (" <<
                   Perforation::viewKeyStruct::distanceFromHeadString() << " = " << m_perfDistFromHead[iperf] <<
                   ") is larger than well polyline length (" << wellLength <<
                   ")\n \n You should check the following values:" <<
                   "\n   1 - " << perf.getWrapperDataContext( Perforation::viewKeyStruct::distanceFromHeadString() ) <<
                   "\n   2 - " << getWrapperDataContext( viewKeyStruct::polylineNodeCoordsString() ) << ", Z values",
                   InputError );

    // start binary search
    const globalIndex maxNumSteps = m_numElems + 1;
    globalIndex currentNumSteps = 0;
    while( iwelemTop < iwelemBottom )
    {
      globalIndex iwelemMid =
        static_cast< globalIndex >(floor( static_cast< real64 >(iwelemTop + iwelemBottom) / 2.0 ));
      real64 const headToBottomDist =
        m_nodeDistFromHead[m_elemToNodesMap[iwelemMid][LineBlockABC::NodeLocation::BOTTOM]];

      if( headToBottomDist < m_perfDistFromHead[iperf] )
      {
        iwelemTop = iwelemMid + 1;
      }
      else
      {
        iwelemBottom = iwelemMid;
      }

      GEOS_THROW_IF( currentNumSteps > maxNumSteps,
                     perf.getDataContext() << ": Perforation cannot be mapped to a well element.",
                     InputError );

      currentNumSteps++;
    }

    // set the index of the matched element
    globalIndex iwelemMatched = iwelemTop;
    GEOS_THROW_IF( iwelemMatched >= m_numElems,
                   getDataContext() << ": Invalid well topology.",
                   InputError );

    m_perfElemId[iperf] = iwelemMatched;

    // compute the physical location of the perforation
    globalIndex const inodeTop    = m_elemToNodesMap[iwelemMatched][LineBlockABC::NodeLocation::TOP];
    globalIndex const inodeBottom = m_elemToNodesMap[iwelemMatched][LineBlockABC::NodeLocation::BOTTOM];
    real64 const elemLength       = m_nodeDistFromHead[inodeBottom] - m_nodeDistFromHead[inodeTop];
    real64 const topToPerfDist    = m_perfDistFromHead[iperf] - m_nodeDistFromHead[inodeTop];

    GEOS_THROW_IF( (elemLength <= 0) || (topToPerfDist < 0),
                   getDataContext() << ": Invalid well topology.",
                   InputError );

    LvArray::tensorOps::copy< 3 >( m_perfCoords[iperf], m_nodeCoords[inodeBottom] );
    LvArray::tensorOps::subtract< 3 >( m_perfCoords[iperf], m_nodeCoords[inodeTop] );
    LvArray::tensorOps::scale< 3 >( m_perfCoords[iperf], topToPerfDist / elemLength );
    LvArray::tensorOps::add< 3 >( m_perfCoords[iperf], m_nodeCoords[inodeTop] );

  }
}


globalIndex InternalWellGenerator::getNextSegmentIndex( globalIndex topSegId,
                                                        globalIndex currentPolyNodeId ) const
{
  globalIndex nextSegId = -1;

  // get the index of the two segments sharing this node
  for( auto is : m_polyNodeToSegmentMap[currentPolyNodeId] )
  {
    // if this is not equal to the current segId, we found the next segId
    if( is != topSegId )
    {
      nextSegId = is;
      break;
    }
  }

  return nextSegId;
}

void InternalWellGenerator::checkPerforationLocationsValidity()
{
  array1d< array1d< localIndex > > elemToPerfMap;
  elemToPerfMap.resize( m_numElems );

  for( globalIndex iperf = 0; iperf < m_numPerforations; ++iperf )
  {
    elemToPerfMap[m_perfElemId[iperf]].emplace_back( iperf );
  }

  // merge perforations to make sure that no well element is shared between two MPI domains
  // TODO: instead of merging perforations, split the well elements and do not change the physical location of the
  // perforation
  int const mpiSize = MpiWrapper::commSize( MPI_COMM_GEOSX );
  if( mpiSize > 1 )
  {
    mergePerforations( elemToPerfMap );
  }

  for( globalIndex iwelem = 0; iwelem < m_numElems; ++iwelem )
  {
    // check that there is always a perforation in the last well element (otherwise, the problem is not well posed)
    for( localIndex iwelemPrev = 0; iwelemPrev < m_prevElemId[iwelem].size(); ++iwelemPrev )
    {
      GEOS_THROW_IF( m_prevElemId[iwelem][iwelemPrev] == -1 && elemToPerfMap[iwelem].size() == 0,
                     getDataContext() << "The bottom element of the well does not have a perforation. " <<
                     "This is needed to have a well-posed problem. \n\n" <<
                     "Here are the two possible ways to solve this problem: \n\n" <<
                     "1) Adding a perforation located close to the bottom of the well. " <<
                     "To do that, compute the total length of the well polyline (by summing the length of the well segments defined by the keywords \"polylineNodeCoords\" and \"polylineSegmentConn\") " <<
                     "and place a perforation whose \"distanceFromHead\" is slightly smaller than this total length. \n \n" <<
                     "2) Shorten the well polyline. " <<
                     "To do that, reduce the length of the well polyline by shortening the segments defined by the keywords \"polylineNodeCoords\" and \"polylineSegmentConn\", or by removing a segment.",
                     InputError );
    }
  }
}

void InternalWellGenerator::mergePerforations( array1d< array1d< localIndex > > const & elemToPerfMap )
{

  for( globalIndex iwelem = 0; iwelem < m_numElems; ++iwelem )
  {
    // collect the indices of the elems with more that one perforation
    if( elemToPerfMap[iwelem].size() > 1 )
    {
      // find the perforation with the largest Peaceman index and keep its location
      globalIndex iperfMaxTransmissibility = elemToPerfMap[iwelem][0];
      real64 maxTransmissibility = m_perfTransmissibility[iperfMaxTransmissibility];
      for( localIndex ip = 1; ip < elemToPerfMap[iwelem].size(); ++ip )
      {
        if( m_perfTransmissibility[elemToPerfMap[iwelem][ip]] > maxTransmissibility )
        {
          iperfMaxTransmissibility = elemToPerfMap[iwelem][ip];
          maxTransmissibility = m_perfTransmissibility[iperfMaxTransmissibility];
        }
      }

      // assign the coordinates of the perf with the largest trans to the other perfs on this elem
      for( localIndex ip = 0; ip < elemToPerfMap[iwelem].size(); ++ip )
      {
        if( elemToPerfMap[iwelem][ip] == iperfMaxTransmissibility )
        {
          continue;
        }

        GEOS_LOG_RANK_0( "\n \nThe GEOSX wells currently have the following limitation in parallel: \n" <<
                         "We cannot allow an element of the well mesh to have two or more perforations associated with it. \n" <<
                         "So, in the present simulation, perforation #" << elemToPerfMap[iwelem][ip] <<
                         " of well " << getDataContext() <<
                         " is moved from " << m_perfCoords[elemToPerfMap[iwelem][ip]] <<
                         " to " << m_perfCoords[iperfMaxTransmissibility] <<
                         " to make sure that no element of the well mesh has two perforations associated with it. \n" <<
                         "To circumvent this issue, please increase the value of \"numElementsPerSegment\" that controls the number of (uniformly distributed) well mesh elements per segment of the well polyline. \n" <<
                         "Our recommendation is to choose \"numElementsPerSegment\" such that each well mesh element has at most one perforation. \n\n" );
        LvArray::tensorOps::copy< 3 >( m_perfCoords[elemToPerfMap[iwelem][ip]], m_perfCoords[iperfMaxTransmissibility] );
      }
    }
  }
}

void InternalWellGenerator::debugWellGeometry() const
{
  if( MpiWrapper::commRank( MPI_COMM_GEOSX ) != 0 )
  {
    return;
  }

  std::cout << std::endl;
  std::cout << "++++++++++++++++++++++++++" << std::endl;
  std::cout << "InternalWellGenerator = " << getName() << std::endl;
  std::cout << "MPI rank = " << MpiWrapper::commRank( MPI_COMM_GEOSX ) << std::endl << std::endl;
  std::cout << "Number of well elements = " << m_numElems << std::endl;

  for( globalIndex iwelem = 0; iwelem < m_numElems; ++iwelem )
  {
    std::cout << "Well element #" << iwelem << std::endl;
    std::cout << "Coordinates of the element center: " << m_elemCenterCoords[iwelem] << std::endl;
    if( m_nextElemId[iwelem] < 0 )
    {
      std::cout << "No next well element" << std::endl;
    }
    else
    {
      std::cout << "Next well element # = " << m_nextElemId[iwelem] << std::endl;
    }
    if( m_prevElemId[iwelem][0] < 0 )
    {
      std::cout << "No previous well element" << std::endl;
    }
    else
    {
      std::cout << "Previous well element #" << m_prevElemId[iwelem][0] << std::endl;
    }
    for( globalIndex inode = 0; inode < m_numNodesPerElem; ++inode )
    {
      if( inode == 0 )
      {
        std::cout << "First well node: #" << m_elemToNodesMap[iwelem][inode] << std::endl;
      }
      else
      {
        std::cout << "Second well node: #" << m_elemToNodesMap[iwelem][inode] << std::endl;
      }
    }
  }

  std::cout << std::endl << "Number of perforations = " << m_numPerforations << std::endl;

  for( globalIndex iperf = 0; iperf < m_numPerforations; ++iperf )
  {
    std::cout << "Perforation #" << iperf << std::endl;
    std::cout << "Coordinates of the perforation: " << m_perfCoords[iperf] << std::endl;
    std::cout << "Is connected to well element #" << m_perfElemId[iperf] << std::endl;
  }
  std::cout << std::endl;

}

REGISTER_CATALOG_ENTRY( WellGeneratorBase, InternalWellGenerator, string const &, Group * const )
}
