
#include "BufferOps.hpp"
#include "dataRepository/BufferOps.hpp"
#include "ToElementRelation.hpp"
#include "ElementRegionManager.hpp"

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

  for( localIndex a=0 ; a<packList.size() ; ++a )
  {
    localIndex index = packList[a];
    for( localIndex b=0 ; b<var.m_toElementRegion[index].size() ; ++b )
    {
      localIndex elemRegionIndex             = var.m_toElementRegion[index][b];
      ElementRegion const * const elemRegion = elementRegionManager->GetRegion(elemRegionIndex);

      localIndex elemSubRegionIndex                  = var.m_toElementSubRegion[index][b];
      CellBlockSubRegion const * const elemSubRegion = elemRegion->GetSubRegion(elemSubRegionIndex);

      localIndex elemIndex = var.m_toElementRegion[index][b];

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
}
}
