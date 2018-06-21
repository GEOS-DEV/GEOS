
#include "BufferOps.hpp"
#include "dataRepository/BufferOps.hpp"
#include "ToElementRelation.hpp"
#include "ElementRegionManager.hpp"
#include "codingUtilities/Utilities.hpp"

namespace geosx
{
namespace bufferOps
{


template< bool DO_PACKING >
localIndex Pack( char*& buffer,
                 UnorderedVariableToManyElementRelation const & var,
                 array<localIndex> const & packList,
                 ElementRegionManager const * const elementRegionManager )
{
  localIndex sizeOfPackedChars = 0;

  sizeOfPackedChars += Pack<DO_PACKING>( buffer, packList.size() );
  for( localIndex a=0 ; a<packList.size() ; ++a )
  {
    localIndex index = packList[a];
    sizeOfPackedChars += Pack<DO_PACKING>( buffer, var.m_toElementRegion[index].size() );
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
                                UnorderedVariableToManyElementRelation const &,
                                array<localIndex> const &,
                                ElementRegionManager const * const );
template localIndex Pack<false>( char*&,
                                 UnorderedVariableToManyElementRelation const &,
                                 array<localIndex> const &,
                                 ElementRegionManager const * const );


localIndex Unpack( char const * & buffer,
                   UnorderedVariableToManyElementRelation & var,
                   array<localIndex> const & packList,
                   ElementRegionManager const * const elementRegionManager )
{
  localIndex sizeOfUnpackedChars = 0;

  localIndex numIndicesUnpacked;
  sizeOfUnpackedChars += Unpack( buffer, numIndicesUnpacked );
  GEOS_ASSERT( numIndicesUnpacked==packList.size(), "")

  for( localIndex a=0 ; a<packList.size() ; ++a )
  {
    localIndex index = packList[a];
    sizeOfUnpackedChars += Unpack( buffer, numIndicesUnpacked );
    var.m_toElementRegion[index].resize( numIndicesUnpacked );
    var.m_toElementSubRegion[index].resize( numIndicesUnpacked );
    var.m_toElementIndex[index].resize( numIndicesUnpacked );
//    GEOS_ASSERT( numIndicesUnpacked==var.m_toElementRegion[index].size(), "")

    for( localIndex b=0 ; b<var.m_toElementRegion[index].size() ; ++b )
    {
      localIndex elemRegionIndex             = var.m_toElementRegion[index][b];
      ElementRegion const * const elemRegion = elementRegionManager->GetRegion(elemRegionIndex);

      localIndex elemSubRegionIndex                  = var.m_toElementSubRegion[index][b];
      CellBlockSubRegion const * const elemSubRegion = elemRegion->GetSubRegion(elemSubRegionIndex);


      sizeOfUnpackedChars += bufferOps::Unpack( buffer, var.m_toElementRegion[index][b] );
      sizeOfUnpackedChars += bufferOps::Unpack( buffer, var.m_toElementSubRegion[index][b] );

      globalIndex globalElementIndex;
      sizeOfUnpackedChars += bufferOps::Unpack( buffer, globalElementIndex );

      var.m_toElementIndex[index][b] = softMapLookup( elemSubRegion->m_globalToLocalMap,
                                                      globalElementIndex,
                                                      localIndex(-1) );
    }
  }

  return sizeOfUnpackedChars;
}













template< bool DO_PACKING >
localIndex Pack( char*& buffer,
                 FixedToManyElementRelation const & var,
                 array<localIndex> const & packList,
                 ElementRegionManager const * const elementRegionManager )
{
  localIndex sizeOfPackedChars = 0;

  sizeOfPackedChars += Pack<DO_PACKING>( buffer, packList.size() );
  for( localIndex a=0 ; a<packList.size() ; ++a )
  {
    localIndex index = packList[a];
    sizeOfPackedChars += Pack<DO_PACKING>( buffer, var.m_toElementRegion[index].size() );
    for( localIndex b=0 ; b<var.m_toElementRegion[index].size() ; ++b )
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
                                array<localIndex> const &,
                                ElementRegionManager const * const );
template localIndex Pack<false>( char*&,
                                 FixedToManyElementRelation const &,
                                 array<localIndex> const &,
                                 ElementRegionManager const * const );


localIndex Unpack( char const * & buffer,
                   FixedToManyElementRelation & var,
                   array<localIndex> const & packList,
                   ElementRegionManager const * const elementRegionManager )
{
  localIndex sizeOfUnpackedChars = 0;

  localIndex numIndicesUnpacked;
  sizeOfUnpackedChars += Unpack( buffer, numIndicesUnpacked );
  GEOS_ASSERT( numIndicesUnpacked==packList.size(), "")

  for( localIndex a=0 ; a<packList.size() ; ++a )
  {
    localIndex index = packList[a];
    sizeOfUnpackedChars += Unpack( buffer, numIndicesUnpacked );
    GEOS_ASSERT( numIndicesUnpacked==var.m_toElementRegion[index].size(), "")

    for( localIndex b=0 ; b<var.m_toElementRegion[index].size() ; ++b )
    {
      localIndex & elemRegionIndex = var.m_toElementRegion[index][b];
      sizeOfUnpackedChars += bufferOps::Unpack( buffer, elemRegionIndex );

      if( elemRegionIndex==-1 )
      {
        sizeOfUnpackedChars += bufferOps::Unpack( buffer, var.m_toElementSubRegion[index][b] );
        sizeOfUnpackedChars += bufferOps::Unpack( buffer, var.m_toElementIndex[index][b] );
      }
      else
      {
        ElementRegion const * const elemRegion = elementRegionManager->GetRegion(elemRegionIndex);

        localIndex & elemSubRegionIndex = var.m_toElementSubRegion[index][b];
        sizeOfUnpackedChars += bufferOps::Unpack( buffer, elemSubRegionIndex );

        CellBlockSubRegion const * const elemSubRegion = elemRegion->GetSubRegion(elemSubRegionIndex);

        globalIndex globalElementIndex;
        sizeOfUnpackedChars += bufferOps::Unpack( buffer, globalElementIndex );
        var.m_toElementIndex[index][b] = softMapLookup( elemSubRegion->m_globalToLocalMap,
                                                        globalElementIndex,
                                                        localIndex(-1) );

      }


    }
  }

  return sizeOfUnpackedChars;
}







}
}
