/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */


#include "BufferOps.hpp"
#include "dataRepository/BufferOps.hpp"
#include "ToElementRelation.hpp"
#include "codingUtilities/Utilities.hpp"
#include "ElementRegionManager.hpp"

namespace geosx
{
namespace bufferOps
{


template< bool DO_PACKING >
localIndex Pack( buffer_unit_type *& buffer,
                 OrderedVariableToManyElementRelation const & var,
                 arrayView1d<localIndex const> const & packList,
                 ElementRegionManager const * const elementRegionManager )
{
  localIndex sizeOfPackedChars = 0;

  sizeOfPackedChars += bufferOps::Pack<DO_PACKING>( buffer, packList.size() );
  for( localIndex a=0 ; a<packList.size() ; ++a )
  {
    localIndex index = packList[a];
    sizeOfPackedChars += bufferOps::Pack<DO_PACKING>( buffer, var.m_toElementRegion.sizeOfArray(index) );
    for( localIndex b=0 ; b<var.m_toElementRegion.sizeOfArray(index) ; ++b )
    {
      localIndex elemRegionIndex    = var.m_toElementRegion[index][b];
      localIndex elemSubRegionIndex = var.m_toElementSubRegion[index][b];
      localIndex elemIndex          = var.m_toElementIndex[index][b];

      sizeOfPackedChars += bufferOps::Pack<DO_PACKING>( buffer, elemRegionIndex );
      sizeOfPackedChars += bufferOps::Pack<DO_PACKING>( buffer, elemSubRegionIndex );

      if( elemRegionIndex!=-1 && elemSubRegionIndex!=-1 && elemIndex!=-1 )
      {
        ElementRegionBase const * const elemRegion = elementRegionManager->GetRegion(elemRegionIndex);
        ElementSubRegionBase const * const elemSubRegion = elemRegion->GetSubRegion(elemSubRegionIndex);
        sizeOfPackedChars += bufferOps::Pack<DO_PACKING>( buffer, elemSubRegion->m_localToGlobalMap[elemIndex] );
      }
      else
      {
        sizeOfPackedChars += bufferOps::Pack<DO_PACKING>( buffer, elemIndex );
      }
    }
  }

  return sizeOfPackedChars;
}
template localIndex Pack<true>( buffer_unit_type *&,
                                OrderedVariableToManyElementRelation const &,
                                arrayView1d<localIndex const> const &,
                                ElementRegionManager const * const );
template localIndex Pack<false>( buffer_unit_type *&,
                                 OrderedVariableToManyElementRelation const &,
                                 arrayView1d<localIndex const> const &,
                                 ElementRegionManager const * const );


localIndex Unpack( buffer_unit_type const * & buffer,
                   OrderedVariableToManyElementRelation & var,
                   arrayView1d<localIndex const> const & packList,
                   ElementRegionManager const * const elementRegionManager,
                   bool const clearFlag )
{
  localIndex sizeOfUnpackedChars = 0;

  localIndex numIndicesUnpacked;
  sizeOfUnpackedChars += bufferOps::Unpack( buffer, numIndicesUnpacked );
  GEOSX_ERROR_IF( numIndicesUnpacked != packList.size(), "");

  for( localIndex a=0 ; a<packList.size() ; ++a )
  {
    localIndex const index = packList[a];
    sizeOfUnpackedChars += bufferOps::Unpack( buffer, numIndicesUnpacked );

    std::set< std::tuple<localIndex, localIndex, localIndex> > values;
    if( !clearFlag )
    {
      for( localIndex b=0 ; b<var.m_toElementRegion.sizeOfArray(index) ; ++b )
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

      localIndex elemSubRegionIndex;
      sizeOfUnpackedChars += bufferOps::Unpack( buffer, elemSubRegionIndex );

      globalIndex globalElementIndex;
      sizeOfUnpackedChars += bufferOps::Unpack( buffer, globalElementIndex );

      if( elemRegionIndex!=-1 && elemSubRegionIndex!=-1 )
      {
        ElementRegionBase const * const elemRegion = elementRegionManager->GetRegion(elemRegionIndex);
        ElementSubRegionBase const * const elemSubRegion = elemRegion->GetSubRegion(elemSubRegionIndex);

        localIndex localElementIndex = softMapLookup( elemSubRegion->m_globalToLocalMap,
                                                      globalElementIndex,
                                                      localIndex(-1) );
        if( localElementIndex!=-1 )
        {
          values.insert( std::make_tuple( elemRegionIndex,
                                          elemSubRegionIndex,
                                          localElementIndex ) );
        }
      }
      else
      {
        values.insert( std::make_tuple( -1, -1, -1 ) );
      }
    }

    localIndex const newSize = values.size();

    var.m_toElementRegion.resizeArray( index, newSize );
    var.m_toElementSubRegion.resizeArray( index, newSize );
    var.m_toElementIndex.resizeArray( index, newSize );

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
localIndex Pack( buffer_unit_type *& buffer,
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
      localIndex elemRegionIndex    = var.m_toElementRegion[index][b];
      localIndex elemSubRegionIndex = var.m_toElementSubRegion[index][b];
      localIndex elemIndex          = var.m_toElementIndex[index][b];

      sizeOfPackedChars += bufferOps::Pack<DO_PACKING>( buffer, elemRegionIndex );
      sizeOfPackedChars += bufferOps::Pack<DO_PACKING>( buffer, elemSubRegionIndex );

      if( elemRegionIndex!=-1 && elemSubRegionIndex!=-1 && elemIndex!=-1 )
      {
        ElementRegionBase const * const elemRegion = elementRegionManager->GetRegion(elemRegionIndex);
        ElementSubRegionBase const * const elemSubRegion = elemRegion->GetSubRegion(elemSubRegionIndex);
        sizeOfPackedChars += bufferOps::Pack<DO_PACKING>( buffer, elemSubRegion->m_localToGlobalMap[elemIndex] );
      }
      else
      {
        sizeOfPackedChars += bufferOps::Pack<DO_PACKING>( buffer, elemIndex );
      }
    }
  }

  return sizeOfPackedChars;
}
template localIndex Pack<true>( buffer_unit_type *&,
                                FixedToManyElementRelation const &,
                                arrayView1d<localIndex const> const &,
                                ElementRegionManager const * const );
template localIndex Pack<false>( buffer_unit_type *&,
                                 FixedToManyElementRelation const &,
                                 arrayView1d<localIndex const> const &,
                                 ElementRegionManager const * const );


localIndex Unpack( buffer_unit_type const * & buffer,
                   FixedToManyElementRelation & var,
                   arrayView1d<localIndex const> const & packList,
                   ElementRegionManager const * const elementRegionManager,
                   bool const clearFlag )
{
  localIndex sizeOfUnpackedChars = 0;

  localIndex numIndicesUnpacked;
  sizeOfUnpackedChars += bufferOps::Unpack( buffer, numIndicesUnpacked );
  GEOSX_ERROR_IF( numIndicesUnpacked != packList.size(), "");

  for( localIndex a=0 ; a<packList.size() ; ++a )
  {
    localIndex index = packList[a];
    localIndex numSubIndicesUnpacked;
    sizeOfUnpackedChars += bufferOps::Unpack( buffer, numSubIndicesUnpacked );
    GEOSX_ERROR_IF( numSubIndicesUnpacked != var.m_toElementRegion.size(1), "");

    for( localIndex b=0 ; b<numSubIndicesUnpacked ; ++b )
    {
      localIndex recvElemRegionIndex;
      localIndex recvElemSubRegionIndex;
      globalIndex globalElementIndex;

      sizeOfUnpackedChars += bufferOps::Unpack( buffer, recvElemRegionIndex );
      sizeOfUnpackedChars += bufferOps::Unpack( buffer, recvElemSubRegionIndex );
      sizeOfUnpackedChars += bufferOps::Unpack( buffer, globalElementIndex );

      if( recvElemRegionIndex!=-1 && recvElemSubRegionIndex!=-1 && globalElementIndex!=-1 )
      {
        ElementRegionBase const * const
        elemRegion = elementRegionManager->GetRegion(recvElemRegionIndex);

        ElementSubRegionBase const * const
        elemSubRegion = elemRegion->GetSubRegion(recvElemSubRegionIndex);

        localIndex const recvElemIndex = softMapLookup( elemSubRegion->m_globalToLocalMap,
                                                        globalElementIndex,
                                                        localIndex(-1) );

        for( localIndex c=0; c<var.m_toElementRegion.size(1) ; ++c )
        {
          localIndex & elemRegionIndex = var.m_toElementRegion[index][c];
          localIndex & elemSubRegionIndex = var.m_toElementSubRegion[index][c];
          localIndex & elemIndex = var.m_toElementIndex[index][c];
          if( ( elemRegionIndex==recvElemRegionIndex &&
                elemSubRegionIndex==recvElemSubRegionIndex &&
                elemIndex==recvElemIndex ) )
          {
            break;
          }
          else if( ( elemRegionIndex==-1 || elemSubRegionIndex==-1 || elemIndex==-1 ) )
          {
            elemRegionIndex = recvElemRegionIndex;
            elemSubRegionIndex = recvElemSubRegionIndex;
            elemIndex = recvElemIndex;
            break;
          }
          else
          {
            //TODO need a better criteria and an error check here
          }
        }
      }
      else if( clearFlag )
      {
        var.m_toElementRegion[index][b] = -1;
        var.m_toElementSubRegion[index][b] = -1;
        var.m_toElementIndex[index][b] = -1;
      }
    }
  }

  return sizeOfUnpackedChars;
}

}
}
