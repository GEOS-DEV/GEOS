
#ifndef MESH_BUFFEROPS_H_
#define MESH_BUFFEROPS_H_

#include "common/DataTypes.hpp"

class UnorderedVariableToManyElementRelation;
class ElementRegionManager;

namespace geosx
{
namespace bufferOps
{

template< bool DO_PACKING >
localIndex Pack( char*& buffer,
                 UnorderedVariableToManyElementRelation const & var,
                 array<localIndex> const & packList,
                 ElementRegionManager const * const elementRegionManager );




}
}
#endif
