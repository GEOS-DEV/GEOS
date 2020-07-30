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


// Python.h must be included first.
#define PY_SSIZE_T_CLEAN
#include <Python.h>

// Source includes
#include "pygeosx.hpp"
#include "PyGroup.hpp"
#include "PyWrapper.hpp"
#include "managers/initialization.hpp"
#include "managers/GeosxState.hpp"

// System includes
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcast-qual"
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
#pragma GCC diagnostic pop
#include <chrono>

namespace geosx
{

std::unique_ptr< GeosxState > g_state;

static constexpr char const * initializeDocString =
"initialize(rank, args)\n"
"--\n\n"
"Initialize GEOSX, this must be the first module call.\n"
"\n"
"Parameters\n"
"__________\n"
"rank : int\n"
"    The MPI rank of the process.\n"
"args : list of strings\n"
"    The list of command line arguments to pass to GEOSX.\n"
"\n"
"Returns\n"
"_______\n"
"Group\n"
"    The ProblemManager.";
static PyObject * initialize( PyObject * self, PyObject * args )
{
  GEOSX_UNUSED_VAR( self );

  if ( g_state != nullptr )
  {
    PyErr_SetString( PyExc_RuntimeError, "state already initialized" );
    return nullptr;
  }

  long pythonMPIRank;
  PyObject * list;
  if ( !PyArg_ParseTuple( args, "lO", &pythonMPIRank, &list ) )
  { return nullptr; }

  PyObjectRef iterator{ PyObject_GetIter( list ) };
  if ( iterator == nullptr )
  { return nullptr; }

  std::vector< std::string > stringArgs;
  for( PyObjectRef item{ PyIter_Next( iterator ) }; item != nullptr; item = PyIter_Next( iterator ) )
  {
    PyObjectRef ascii { PyUnicode_AsASCIIString( item ) };
    if ( ascii == nullptr )
    { return nullptr; }

    char const * const stringValue = PyBytes_AsString( ascii );
    if ( stringValue == nullptr )
    { return nullptr; }

    stringArgs.push_back( stringValue );
  }

  if ( PyErr_Occurred() )
  { return nullptr; }

  std::vector< char * > argv( stringArgs.size() + 1 );
  for ( std::size_t i = 0; i < stringArgs.size(); ++i )
  { argv[ i ] = const_cast< char * >( stringArgs[ i ].data() ); }

  g_state = std::make_unique< GeosxState >( basicSetup( argv.size() - 1, argv.data(), true ) );

  if ( pythonMPIRank != MpiWrapper::Comm_rank() )
  {
    PyErr_SetString( PyExc_ValueError, "Python MPI rank does not align with GEOSX MPI rank" );
    return nullptr;
  }

  g_state->initializeDataRepository();

  return python::createNewPyGroup( g_state->getProblemManagerAsGroup() );
}

static constexpr char const * reinitDocString = "";
static PyObject * reinit( PyObject * self, PyObject * args )
{
  GEOSX_UNUSED_VAR( self );

  if ( g_state == nullptr || g_state->getState() != State::COMPLETED )
  {
    PyErr_SetString( PyExc_RuntimeError, "State must be COMPLETED" );
    return nullptr;
  }

  PyObject * list;
  if ( !PyArg_ParseTuple( args, "O", &list ) )
  { return nullptr; }

  PyObjectRef iterator{ PyObject_GetIter( list ) };
  if ( iterator == nullptr )
  { return nullptr; }

  std::vector< std::string > stringArgs;
  for( PyObjectRef item{ PyIter_Next( iterator ) }; item != nullptr; item = PyIter_Next( iterator ) )
  {
    PyObjectRef ascii { PyUnicode_AsASCIIString( item ) };
    if ( ascii == nullptr )
    { return nullptr; }

    char const * const stringValue = PyBytes_AsString( ascii );
    if ( stringValue == nullptr )
    { return nullptr; }

    stringArgs.push_back( stringValue );
  }

  if ( PyErr_Occurred() )
  { return nullptr; }

  std::vector< char * > argv( stringArgs.size() + 1 );
  for ( std::size_t i = 0; i < stringArgs.size(); ++i )
  { argv[ i ] = const_cast< char * >( stringArgs[ i ].data() ); }

  // Must first delete the existing state.
  g_state = nullptr;
  g_state = std::make_unique< GeosxState >( parseCommandLineOptions( argv.size() - 1, argv.data() ) );

  g_state->initializeDataRepository();

  return python::createNewPyGroup( g_state->getProblemManagerAsGroup() );
}

static constexpr char const * applyInitialConditionsDocString =
"applyInitialConditions()\n"
"--\n\n"
"Apply the initial conditions.\n"
"\n"
"Returns\n"
"_______\n"
"None\n";
static PyObject * applyInitialConditions( PyObject * self, PyObject * args )
{
  GEOSX_UNUSED_VAR( self, args );

  if ( g_state == nullptr ){
    PyErr_SetString( PyExc_RuntimeError, "state must be initialized" );
    return nullptr;
  }
  g_state->applyInitialConditions();
  Py_RETURN_NONE;
}

static constexpr char const * runDocString =
"run()\n"
"--\n\n"
"Enter the event loop.\n"
"\n"
"Returns\n"
"_______\n"
"int\n"
"    The state of the simulation. If the simulation has ended the value is `COMPLETED`. If the "
"simulation still has steps left to run the value is `READY_TO_RUN`.";
static PyObject * run( PyObject * self, PyObject * args )
{
  GEOSX_UNUSED_VAR( self, args );

  if ( g_state == nullptr ){
    PyErr_SetString( PyExc_RuntimeError, "state must be initialized" );
    return nullptr;
  }
  g_state->run();
  return PyLong_FromLong( static_cast< int >( g_state->getState() ) );
}

static constexpr char const * finalizeDocString =
"finalize()\n"
"--\n\n"
"Finalize GEOSX. After this no calls into pygeosx or to MPI are allowed.\n"
"\n"
"Returns\n"
"_______\n"
"None\n";
static PyObject * finalize( PyObject * self, PyObject * args )
{
  GEOSX_UNUSED_VAR( self, args );

  if ( g_state == nullptr )
  {
    PyErr_SetString( PyExc_RuntimeError, "State either not initialized or already finalized." );
    return nullptr;
  }

  g_state = nullptr;
  basicCleanup();

  Py_RETURN_NONE;
}

} // namespace geosx


/**
 * Add geosx::State enums to the given module. Return the module, or nullptr on failure
 */
static bool addConstants( PyObject * module )
{
  std::array< std::pair< long, char const * >, 4 > const constants = { {
    { static_cast< long >( geosx::State::COMPLETED ), "COMPLETED" },
    { static_cast< long >( geosx::State::INITIALIZED ), "INITIALIZED" },
    { static_cast< long >( geosx::State::UNINITIALIZED ), "UNINITIALIZED" },
    { static_cast< long >( geosx::State::READY_TO_RUN ), "READY_TO_RUN" }
  } };

  for ( std::pair< long, char const * > const & pair : constants )
  {
    if ( PyModule_AddIntConstant( module, pair.second, pair.first ) )
    {
      PyErr_SetString( PyExc_RuntimeError, "couldn't add constant" );
      return false;
    }
  }

  return true;
}

/**
 * Add exit handler to a module. Register the module's 'finalize' function,
 * which should take no arguments, with the `atexit` standard library module.
 */
static bool addExitHandler( PyObject * module ){
  geosx::PyObjectRef atexit_module { PyImport_ImportModule( "atexit" ) };
  
  if ( atexit_module == nullptr )
  { return 0; }

  geosx::PyObjectRef atexit_register_pyfunc { PyObject_GetAttrString( atexit_module, "register" ) };
  if ( atexit_register_pyfunc == nullptr )
  { return 0; }

  geosx::PyObjectRef finalize_pyfunc { PyObject_GetAttrString( module, "finalize" ) };
  if ( finalize_pyfunc == nullptr )
  { return 0; }

  if ( !PyCallable_Check( atexit_register_pyfunc ) || !PyCallable_Check( finalize_pyfunc ) )
  { return 0; }

  geosx::PyObjectRef returnval { PyObject_CallFunctionObjArgs( atexit_register_pyfunc, finalize_pyfunc.get(), nullptr ) };
  
  return returnval != nullptr;
}

// Allow mixing designated and non-designated initializers in the same initializer list.
// I don't like the pragmas but the designated initializers is the only sane way to do this stuff.
// The other option is to put this in a `.c` file and compile with the C compiler, but that seems like more work.
#pragma GCC diagnostic push
#if defined( __clang_version__ )
  #pragma GCC diagnostic ignored "-Wc99-designator"
#else
  #pragma GCC diagnostic ignored "-Wpedantic"
  #pragma GCC diagnostic ignored "-Wmissing-field-initializers"
#endif

/**
 *
 */
static PyMethodDef pygeosxFuncs[] = {
  { "initialize", geosx::initialize, METH_VARARGS, geosx::initializeDocString },
  { "reinit", geosx::reinit, METH_VARARGS, geosx::reinitDocString },
  { "applyInitialConditions", geosx::applyInitialConditions, METH_NOARGS, geosx::applyInitialConditionsDocString },
  { "run", geosx::run, METH_NOARGS, geosx::runDocString },
  { "finalize", geosx::finalize, METH_NOARGS, geosx::finalizeDocString },
  { nullptr, nullptr, 0, nullptr }        /* Sentinel */
};

static constexpr char const * pygeosxDocString =
"Python driver for GEOSX.";

/**
 * Initialize the module object for Python with the exported functions
 */
static struct PyModuleDef pygeosxModuleFunctions = {
  PyModuleDef_HEAD_INIT,
  .m_name = "pygeosx",
  .m_doc = pygeosxDocString,
  .m_size = -1,
  .m_methods = pygeosxFuncs
};

#pragma GCC diagnostic pop

/**
 * Initialize the module with functions, constants, and exit handler
 */
PyMODINIT_FUNC
PyInit_pygeosx(void)
{
  import_array();

  if( PyType_Ready( geosx::python::getPyGroupType() ) < 0 )
  { return nullptr; }

  if( PyType_Ready( geosx::python::getPyWrapperType() ) < 0 )
  { return nullptr; }

  geosx::PyObjectRef module{ PyModule_Create( &pygeosxModuleFunctions ) };
  if ( module == nullptr )
  { return nullptr; }

  if ( !addExitHandler( module ) )
  {
    if ( PyErr_Occurred() == nullptr )
    { PyErr_SetString( PyExc_RuntimeError, "couldn't add exit handler" ); }

    return nullptr;
  }

  if( !addConstants( module ) )
  { return nullptr; }

  Py_INCREF( geosx::python::getPyGroupType() );
  if ( PyModule_AddObject( module, "Group", reinterpret_cast< PyObject * >( geosx::python::getPyGroupType() ) ) < 0 )
  {
    Py_DECREF( geosx::python::getPyGroupType() );
    return nullptr;
  }

  Py_INCREF( geosx::python::getPyWrapperType() );
  if ( PyModule_AddObject( module, "Wrapper", reinterpret_cast< PyObject * >( geosx::python::getPyWrapperType() ) ) < 0 )
  {
    Py_DECREF( geosx::python::getPyWrapperType() );
    return nullptr;
  }

  // Since we return module we don't want to decrease the reference count.
  return module.release();
}
