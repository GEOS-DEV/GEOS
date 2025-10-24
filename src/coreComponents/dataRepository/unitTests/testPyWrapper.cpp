/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

// Python.h must be included first.
#define PY_SSIZE_T_CLEAN
#include <Python.h>

// Source includes
#include "dataRepository/Group.hpp"
#include "dataRepository/Wrapper.hpp"
#include "dataRepository/python/PyWrapper.hpp"
#include "dataRepository/python/PyGroup.hpp"

namespace geos
{
namespace python
{
namespace testing
{

// Static test data - initialized once
static conduit::Node s_node;
static dataRepository::Group s_rootGroup( "root", s_node );
static dataRepository::Group * s_testGroup = nullptr;
static bool s_initialized = false;

/**
 * @brief Initialize the test group with wrappers (called once).
 */
static void initializeTestGroup()
{
  if( s_initialized )
  {
    return;
  }

  s_testGroup = &s_rootGroup.registerGroup( "testGroup" );

  // Register various wrapper types for testing
  s_testGroup->registerWrapper< std::string >( "stringWrapper" ).
    setDefaultValue( "initial_string" );

  s_testGroup->registerWrapper< int >( "intWrapper" ).
    setDefaultValue( 0 );

  s_testGroup->registerWrapper< long >( "longWrapper" ).
    setDefaultValue( 0L );

  s_testGroup->registerWrapper< double >( "doubleWrapper" ).
    setDefaultValue( 0.0 );

  s_initialized = true;
}


static PyObject * getStringWrapper( PyObject * self, PyObject * args )
{
  GEOS_UNUSED_VAR( self, args );
  initializeTestGroup();
  return createNewPyWrapper( s_testGroup->getWrapper< std::string >( "stringWrapper" ) );
}

static PyObject * getIntWrapper( PyObject * self, PyObject * args )
{
  GEOS_UNUSED_VAR( self, args );
  initializeTestGroup();
  return createNewPyWrapper( s_testGroup->getWrapper< int >( "intWrapper" ) );
}

static PyObject * getLongWrapper( PyObject * self, PyObject * args )
{
  GEOS_UNUSED_VAR( self, args );
  initializeTestGroup();
  return createNewPyWrapper( s_testGroup->getWrapper< long >( "longWrapper" ) );
}

static PyObject * getDoubleWrapper( PyObject * self, PyObject * args )
{
  GEOS_UNUSED_VAR( self, args );
  initializeTestGroup();
  return createNewPyWrapper( s_testGroup->getWrapper< double >( "doubleWrapper" ) );
}

BEGIN_ALLOW_DESIGNATED_INITIALIZERS

static PyMethodDef TestPyWrapperMethods[] = {
  { "get_string_wrapper", getStringWrapper, METH_NOARGS, "Get a string wrapper for testing." },
  { "get_int_wrapper", getIntWrapper, METH_NOARGS, "Get an int wrapper for testing." },
  { "get_long_wrapper", getLongWrapper, METH_NOARGS, "Get a long wrapper for testing." },
  { "get_double_wrapper", getDoubleWrapper, METH_NOARGS, "Get a double wrapper for testing." },
  { nullptr, nullptr, 0, nullptr }
};

static struct PyModuleDef testPyWrapperModule = {
  PyModuleDef_HEAD_INIT,
  .m_name = "testPyWrapper",
  .m_doc = "",
  .m_size = -1,
  .m_methods = TestPyWrapperMethods
};

END_ALLOW_DESIGNATED_INITIALIZERS

PyMODINIT_FUNC PyInit_testPyWrapper( void )
{
  PyObject * module = PyModule_Create( &testPyWrapperModule );
  if( module == nullptr )
  {
    return nullptr;
  }

  // Ensure PyWrapper type is ready
  PyTypeObject * wrapperType = getPyWrapperType();
  if( PyType_Ready( wrapperType ) < 0 )
  {
    Py_DECREF( module );
    return nullptr;
  }

  return module;
}

} // namespace testing
} // namespace python
} // namespace geos
