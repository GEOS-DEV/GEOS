#include "mesh/ExtrinsicMeshData.hpp"

// FIXME C++17 Remove this file when GEOSX is using c++17 (which introduces inline variables)
//             Otherwise you get some `undefined reference to` when linking.
namespace geosx
{
namespace extrinsicMeshData
{

decltype( ChildIndex::key ) constexpr ChildIndex::key;
decltype( ParentIndex::key ) constexpr ParentIndex::key;

} // end of namespace extrinsicMeshData
} // end of namespace geosx
