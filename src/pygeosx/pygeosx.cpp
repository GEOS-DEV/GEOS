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
#define NPY_NO_DEPRECATED_API NPY_1_15_API_VERSION
#include <numpy/arrayobject.h>
#pragma GCC diagnostic pop
#include <chrono>

namespace geosx
{

std::unique_ptr< GeosxState > & getState()
{
  static std::unique_ptr< GeosxState > s_state;
  return s_state;
};

/**
 *
 */
static PyObject * initialize( PyObject * self, PyObject * args )
{
  GEOSX_UNUSED_VAR( self );

  if ( getState() != nullptr ){
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

  basicSetup( argv.size() - 1, argv.data(), true );

  if ( pythonMPIRank != MpiWrapper::Comm_rank() ){
    PyErr_SetString( PyExc_ValueError, "Python MPI rank does not align with GEOSX MPI rank" );
    return nullptr;
  }

  getState() = std::make_unique< GeosxState >();
  getState()->initializeDataRepository();

  return python::createNewPyGroup( getState()->getProblemManager() );
}

/**
 *
 */
static PyObject * applyInitialConditions( PyObject * self, PyObject * args )
{
  GEOSX_UNUSED_VAR( self, args );

  if ( getState() == nullptr ){
    PyErr_SetString( PyExc_RuntimeError, "state must be initialized" );
    return nullptr;
  }
  getState()->applyInitialConditions();
  Py_RETURN_NONE;
}

/**
 *
 */
static PyObject * run( PyObject * self, PyObject * args )
{
  GEOSX_UNUSED_VAR( self, args );

  if ( getState() == nullptr ){
    PyErr_SetString( PyExc_RuntimeError, "state must be initialized" );
    return nullptr;
  }
  getState()->run();
  return PyLong_FromLong( static_cast< int >( getState()->getState() ) );
}

/**
 *
 */
static PyObject * finalize( PyObject * self, PyObject * args )
{
  GEOSX_UNUSED_VAR( self, args );

  if ( getState() == nullptr ){
    PyErr_SetString( PyExc_RuntimeError, "state already finalized" );
    return nullptr;
  }

  getState() = nullptr;
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
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wc99-designator"

/**
 *
 */
static PyMethodDef pygeosxFuncs[] = {
  { "initialize", geosx::initialize, METH_VARARGS,
    "" },
  { "applyInitialConditions", geosx::applyInitialConditions, METH_NOARGS,
    "" },
  { "run", geosx::run, METH_NOARGS,
    "" },
  { "finalize", geosx::finalize, METH_NOARGS,
    "" },
  { nullptr, nullptr, 0, nullptr }        /* Sentinel */
};

/**
 * Initialize the module object for Python with the exported functions
 */
static struct PyModuleDef pygeosxModuleFunctions = {
  PyModuleDef_HEAD_INIT,
  .m_name = "pygeosx",
  .m_doc = "Module for testing numpy views of LvArray::Array objects",
  .m_size = -1,
  .m_methods = pygeosxFuncs
};

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

#pragma GCC diagnostic pop
