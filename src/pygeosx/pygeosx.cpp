// Python.h must be included first.
#define PY_SSIZE_T_CLEAN
#include <Python.h>

// Source includes
#include "managers/initialization.hpp"
#include "managers/GeosxState.hpp"

// System includes
#include <chrono>

namespace geosx
{

// static std::unique_ptr< GeosxState > state;

// bool initialize( int argc, char ** argv )
// {
//   basicSetup( argc, argv );
//   state = std::make_unique< GeosxState >();
//   state.initializeDataRepository();
// }

// void run()
// {
//   GEOSX_ERROR_IF( state == nullptr, "" );
//   state.run()
// }

// PyObject * get( std::string const & path, std::string const & wrapper, bool const modify )
// {
//   GEOSX_ERROR_IF( state == nullptr, "" );

//   Group * group = m_state.m_problemManager.GetGroupByPath( path );
//   GEOSX_ERROR_IF( group == nullptr, "" );
//   Wrapper * wrapper = group->getWrapperBase( wrapper );
//   GEOSX_ERROR_IF( wrapper == nullptr, "" );
//   return wrapper->createPythonObject( modify );
// }

// void finalize()
// {
//   state = nullptr;
//   basicCleanup();
// }

static PyObject *
run( PyObject *self, PyObject *args )
{
  GEOSX_UNUSED_VAR( self );
  GEOSX_UNUSED_VAR( args );

  GEOSX_LOG( "In c++!" );

  constexpr int argc = 3;
  std::string argStrings[ argc ] = { "", "-i", "integratedTests/solidMechanicsSSLE/SSLE-sedov.xml" };
  char * argv[ argc + 1 ] = { const_cast< char * >( argStrings[ 0 ].data() ),
                              const_cast< char * >( argStrings[ 1 ].data() ),
                              const_cast< char * >( argStrings[ 2 ].data() ),
                              nullptr };

  basicSetup( argc, argv, true );

  {
    GeosxState state;

    bool const problemToRun = state.initializeDataRepository();
    if ( problemToRun )
    { state.run(); }
  }

  basicCleanup();

  Py_RETURN_NONE;
}

} // namespace geosx

static PyMethodDef pygeosxFuncs[] = {
  { "run", geosx::run, METH_VARARGS, 
    "run the simulation" },
  { NULL, NULL, 0, NULL }        /* Sentinel */
};

/**
 * Initialize the module object for Python with the exported functions
 */
static struct PyModuleDef pygeosxModuleFunctions = {
  PyModuleDef_HEAD_INIT,
  "pygeosx",   /* name of module */
  "Module for testing numpy views of LvArray::Array objects", /* module documentation, may be NULL */
  -1,       /* size of per-interpreter state of the module or -1 if the module keeps state in global variables. */
  pygeosxFuncs,
  NULL,
  NULL,
  NULL,
  NULL,
};


PyMODINIT_FUNC
PyInit_pygeosx(void)
{
  return PyModule_Create(&pygeosxModuleFunctions);
}
