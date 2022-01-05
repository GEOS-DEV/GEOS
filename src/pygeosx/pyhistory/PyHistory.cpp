#define PY_SSIZE_T_CLEAN
#include <Python.h>

// Source includes
#include "../pygeosx.hpp"
#include "PyHistory.hpp"

#include "PyHistoryCollectionType.hpp"
#include "PyHistoryOutputType.hpp"

namespace geosx
{
extern std::unique_ptr< GeosxState > g_state;

extern bool g_alreadyInitialized;

namespace python
{


static PyModuleDef pyhistorymodule = {
  PyModuleDef_HEAD_INIT,
  "pyhistory",
  "pyhistory module for HistoryCollection and HistoryOutput",
  -1,
  NULL, NULL, NULL, NULL, NULL
};


PyMODINIT_FUNC
PyInit_pyhistory( void )
{
  LvArray::python::PyObjectRef<> module{ PyModule_Create( &pyhistorymodule ) };
  if( module == nullptr )
  {
    return nullptr;
  }

  if( !LvArray::python::addTypeToModule( module, geosx::python::getPyHistoryCollectionType(), "HistoryCollection" ) )
  {
    return nullptr;
  }

  if( !LvArray::python::addTypeToModule( module, geosx::python::getPyHistoryOutputType(), "HistoryOutput" ) )
  {
    return nullptr;
  }

  return module.release();
}

}
}
