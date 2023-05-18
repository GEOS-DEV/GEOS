#ifndef GEOS_PYTHON_PYVTKOUTPUTTYPE_HPP_
#define GEOS_PYTHON_PYVTKOUTPUTTYPE_HPP_

#include "LvArray/src/python/pythonForwardDeclarations.hpp"
#include "mesh/DomainPartition.hpp"

namespace geos
{
namespace python
{

PyTypeObject * getPyVTKOutputType();

} // namespace python
} // namespace geos

#endif
