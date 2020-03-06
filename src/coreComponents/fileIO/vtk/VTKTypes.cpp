#include "VTKTypes.hpp"

#include "common/DataTypes.hpp"

#include <vtkType.h>

#include <unordered_map>

namespace geosx
{
namespace vtk
{

static std::unordered_map< std::type_index, int > geosxToVTKTypeMap =
{
  {std::type_index( typeid( integer ) ), VTK_INT},
  {std::type_index( typeid( localIndex ) ), VTK_INT},
  {std::type_index( typeid( globalIndex ) ), VTK_LONG},
  {std::type_index( typeid( real32 ) ), VTK_FLOAT},
  {std::type_index( typeid( real64 ) ), VTK_DOUBLE}
};

int GetVTKType( std::type_index const & type )
{
  GEOSX_ERROR_IF( geosxToVTKTypeMap.count( type ) == 0, 
                  "Type " << type.name() << "not recognized for a VTK output" );
  return geosxToVTKTypeMap.at( type );
}
}
}
