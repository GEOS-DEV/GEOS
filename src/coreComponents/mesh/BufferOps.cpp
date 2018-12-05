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


#include "BufferOps.hpp"
#include "dataRepository/BufferOps_inline.hpp"
#include "ToElementRelation.hpp"
#include "ElementRegionManager.hpp"
#include "codingUtilities/Utilities.hpp"

namespace geosx
{
namespace bufferOps
{


template< bool DO_PACKING >
localIndex Pack( char*& buffer,
                 OrderedVariableToManyElementRelation const & var,
                 arrayView1d<localIndex const> const & packList,
                 ElementRegionManager const * const elementRegionManager )
{
  localIndex sizeOfPackedChars = 0;

  sizeOfPackedChars += bufferOps::Pack<DO_PACKING>( buffer, packList.size() );
  for( localIndex a=0 ; a<packList.size() ; ++a )
  {
    localIndex index = packList[a];
    sizeOfPackedChars += bufferOps::Pack<DO_PACKING>( buffer, var.m_toElementRegion[index].size() );
    for( localIndex b=0 ; b<var.m_toElementRegion[index].size() ; ++b )
    {
      localIndex elemRegionIndex             = var.m_toElementRegion[index][b];
      ElementRegion const * const elemRegion = elementRegionManager->GetRegion(elemRegionIndex);

      localIndex elemSubRegionIndex                  = var.m_toElementSubRegion[index][b];
      CellBlockSubRegion const * const elemSubRegion = elemRegion->GetSubRegion(elemSubRegionIndex);

      localIndex elemIndex = var.m_toElementIndex[index][b];

      sizeOfPackedChars += bufferOps::Pack<DO_PACKING>( buffer, elemRegionIndex );
      sizeOfPackedChars += bufferOps::Pack<DO_PACKING>( buffer, elemSubRegionIndex );
      sizeOfPackedChars += bufferOps::Pack<DO_PACKING>( buffer, elemSubRegion->m_localToGlobalMap[elemIndex] );

    }
  }

  return sizeOfPackedChars;
}
template localIndex Pack<true>( char*&,
                                OrderedVariableToManyElementRelation const &,
                                arrayView1d<localIndex const> const &,
                                ElementRegionManager const * const );
template localIndex Pack<false>( char*&,
                                 OrderedVariableToManyElementRelation const &,
                                 arrayView1d<localIndex const> const &,
                                 ElementRegionManager const * const );


localIndex Unpack( char const * & buffer,
                   OrderedVariableToManyElementRelation & var,
                   arrayView1d<localIndex const> const & packList,
                   ElementRegionManager const * const elementRegionManager,
                   bool const clearFlag )
{
  localIndex sizeOfUnpackedChars = 0;

  localIndex numIndicesUnpacked;
  sizeOfUnpackedChars += bufferOps::Unpack( buffer, numIndicesUnpacked );
  GEOS_ERROR_IF( numIndicesUnpacked != packList.size(), "");

  for( localIndex a=0 ; a<packList.size() ; ++a )
  {
    localIndex const index = packList[a];
    sizeOfUnpackedChars += bufferOps::Unpack( buffer, numIndicesUnpacked );

    set< std::tuple<localIndex, localIndex, localIndex> > values;
    if( !clearFlag )
    {
      for( localIndex b=0 ; b<var.m_toElementRegion[index].size() ; ++b )
      {
        values.insert( std::make_tuple( var.m_toElementRegion[index][b],
                                        var.m_toElementSubRegion[index][b],
                                        var.m_toElementIndex[index][b] ) );
      }
    }

    for( localIndex b=0 ; b<numIndicesUnpacked ; ++b )
    {

      localIndex elemRegionIndex;
      sizeOfUnpackedChars += bufferOps::Unpack( buffer, elemRegionIndex );
      ElementRegion const * const elemRegion = elementRegionManager->GetRegion(elemRegionIndex);

      localIndex elemSubRegionIndex;
      sizeOfUnpackedChars += bufferOps::Unpack( buffer, elemSubRegionIndex );

      globalIndex globalElementIndex;
      sizeOfUnpackedChars += bufferOps::Unpack( buffer, globalElementIndex );
      CellBlockSubRegion const * const elemSubRegion = elemRegion->GetSubRegion(elemSubRegionIndex);

      localIndex localElementIndex = softMapLookup( elemSubRegion->m_globalToLocalMap,
                                                    globalElementIndex,
                                                    localIndex(-1) );

      values.insert( std::make_tuple( elemRegionIndex,
                                      elemSubRegionIndex,
                                      localElementIndex ) );
    }

    localIndex const newSize = values.size();

    var.m_toElementRegion[index].resize( newSize );
    var.m_toElementSubRegion[index].resize( newSize );
    var.m_toElementIndex[index].resize( newSize );

    {
      localIndex b=0;
      for( auto & value : values )
      {
        var.m_toElementRegion[index][b]    = std::get<0>(value);
        var.m_toElementSubRegion[index][b] = std::get<1>(value);
        var.m_toElementIndex[index][b]     = std::get<2>(value);
        ++b;
      }
    }
  }

  return sizeOfUnpackedChars;
}


template< bool DO_PACKING >
localIndex Pack( char*& buffer,
                 FixedToManyElementRelation const & var,
                 arrayView1d<localIndex const> const & packList,
                 ElementRegionManager const * const elementRegionManager )
{
  localIndex sizeOfPackedChars = 0;

  sizeOfPackedChars += bufferOps::Pack<DO_PACKING>( buffer, packList.size() );
  for( localIndex a=0 ; a<packList.size() ; ++a )
  {
    localIndex index = packList[a];
    sizeOfPackedChars += bufferOps::Pack<DO_PACKING>( buffer, var.m_toElementRegion.size(1) );
    for( localIndex b=0 ; b<var.m_toElementRegion.size(1) ; ++b )
    {
      localIndex elemRegionIndex             = var.m_toElementRegion[index][b];

      if( elemRegionIndex == -1 )
      {
        sizeOfPackedChars += bufferOps::Pack<DO_PACKING>( buffer, localIndex(-1) );
        sizeOfPackedChars += bufferOps::Pack<DO_PACKING>( buffer, localIndex(-1) );
        sizeOfPackedChars += bufferOps::Pack<DO_PACKING>( buffer, localIndex(-1) );
      }
      else
      {
        ElementRegion const * const elemRegion = elementRegionManager->GetRegion(elemRegionIndex);

        localIndex elemSubRegionIndex                  = var.m_toElementSubRegion[index][b];
        CellBlockSubRegion const * const elemSubRegion = elemRegion->GetSubRegion(elemSubRegionIndex);

        localIndex elemIndex = var.m_toElementIndex[index][b];

        sizeOfPackedChars += bufferOps::Pack<DO_PACKING>( buffer, elemRegionIndex );
        sizeOfPackedChars += bufferOps::Pack<DO_PACKING>( buffer, elemSubRegionIndex );
        sizeOfPackedChars += bufferOps::Pack<DO_PACKING>( buffer, elemSubRegion->m_localToGlobalMap[elemIndex] );
      }
    }
  }

  return sizeOfPackedChars;
}
template localIndex Pack<true>( char*&,
                                FixedToManyElementRelation const &,
                                arrayView1d<localIndex const> const &,
                                ElementRegionManager const * const );
template localIndex Pack<false>( char*&,
                                 FixedToManyElementRelation const &,
                                 arrayView1d<localIndex const> const &,
                                 ElementRegionManager const * const );


localIndex Unpack( char const * & buffer,
                   FixedToManyElementRelation & var,
                   arrayView1d<localIndex const> const & packList,
                   ElementRegionManager const * const elementRegionManager,
                   bool const clearFlag )
{
  localIndex sizeOfUnpackedChars = 0;

  localIndex numIndicesUnpacked;
  sizeOfUnpackedChars += bufferOps::Unpack( buffer, numIndicesUnpacked );
  GEOS_ERROR_IF( numIndicesUnpacked != packList.size(), "");

  for( localIndex a=0 ; a<packList.size() ; ++a )
  {
    localIndex index = packList[a];
    sizeOfUnpackedChars += bufferOps::Unpack( buffer, numIndicesUnpacked );
    GEOS_ERROR_IF( numIndicesUnpacked != var.m_toElementRegion.size(1), "");

    for( localIndex b=0 ; b<numIndicesUnpacked ; ++b )
    {
      localIndex recvElemRegionIndex;
      localIndex recvElemSubRegionIndex;
      globalIndex globalElementIndex;

      sizeOfUnpackedChars += bufferOps::Unpack( buffer, recvElemRegionIndex );
      sizeOfUnpackedChars += bufferOps::Unpack( buffer, recvElemSubRegionIndex );
      sizeOfUnpackedChars += bufferOps::Unpack( buffer, globalElementIndex );

      if( recvElemRegionIndex!=-1 && recvElemSubRegionIndex!=-1 && globalElementIndex!=-11 )
      {
        ElementRegion const * const
        elemRegion = elementRegionManager->GetRegion(recvElemRegionIndex);

        CellBlockSubRegion const * const
        elemSubRegion = elemRegion->GetSubRegion(recvElemSubRegionIndex);

        localIndex const recvElemIndex = softMapLookup( elemSubRegion->m_globalToLocalMap,
                                                        globalElementIndex,
                                                        localIndex(-1) );

        for( localIndex c=0; c<var.m_toElementRegion.size(1) ; ++c )
        {
          localIndex & elemRegionIndex = var.m_toElementRegion[index][c];
          localIndex & elemSubRegionIndex = var.m_toElementSubRegion[index][c];
          localIndex & elemIndex = var.m_toElementIndex[index][c];
          if( ( elemRegionIndex==-1 || elemSubRegionIndex==-1 || elemIndex==-1 ) )
          {
            elemRegionIndex = recvElemRegionIndex;
            elemSubRegionIndex = recvElemSubRegionIndex;
            elemIndex = recvElemIndex;
          }
        }
      }
    }
  }

  return sizeOfUnpackedChars;
}

}
}
