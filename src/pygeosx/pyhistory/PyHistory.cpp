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

static PyMethodDef HistoryFuncs[] = {
  { nullptr, nullptr, 0, nullptr }        /* Sentinel */
};

BEGIN_ALLOW_DESIGNATED_INITIALIZERS

static PyModuleDef pyhistorymodule = {
  PyModuleDef_HEAD_INIT,
  .m_name = "pyhistory",
  .m_doc = "pyhistory module for HistoryCollection and HistoryOutput",
  .m_size = -1,
  .m_methods = HistoryFuncs
};

END_ALLOW_DESIGNATED_INITIALIZERS


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
