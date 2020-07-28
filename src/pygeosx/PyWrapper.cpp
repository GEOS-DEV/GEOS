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
#include "PyWrapper.hpp"
#include "pygeosx.hpp"

#define VERIFY_NON_NULL_AND_RETURN( obj ) \
  if ( obj == nullptr ) \
  { \
    PyErr_SetString( PyExc_RuntimeError, "Passed a nullptr as an argument" ); \
    return nullptr; \
  } \

#define VERIFY_NON_NULL_WRAPPER_AND_RETURN( pywrapper ) \
  if ( pywrapper->wrapper == nullptr ) \
  { \
    PyErr_SetString( PyExc_AttributeError, "Wrapper has not been initialized" ); \
    return nullptr; \
  } \

namespace geosx
{
namespace python
{

/**
 *
 */
struct PyWrapper
{
  PyObject_HEAD
  dataRepository::WrapperBase * wrapper;
};

/**
 *
 */
static dataRepository::WrapperBase * getWrapperBase( PyObject * const obj )
{
  VERIFY_NON_NULL_AND_RETURN( obj );

  int isInstanceOfPyWrapper = PyObject_IsInstance( obj, reinterpret_cast< PyObject * >( getPyWrapperType() ) );
  if ( isInstanceOfPyWrapper < 0 )
  { return nullptr; }

  if ( isInstanceOfPyWrapper == 0 )
  {
    PyErr_SetString( PyExc_AttributeError, "Expect an argument of type Wrapper." );
    return nullptr;
  }

  PyWrapper * wrapper = reinterpret_cast< PyWrapper * >( obj );
  VERIFY_NON_NULL_WRAPPER_AND_RETURN( wrapper );

  return wrapper->wrapper;
}

/**
 *
 */
static PyObject * PyWrapper_repr( PyObject * const obj )
{
  dataRepository::WrapperBase const * const wrapper = getWrapperBase( obj );
  if ( wrapper == nullptr )
  { return nullptr; }

  std::string const path = wrapper->getPath();
  std::string const type = LvArray::demangle( typeid( *wrapper ).name() );
  std::string const repr = path + " ( " + type + " )";
  return PyUnicode_FromString( repr.c_str() );
}

/**
 *
 */
static PyObject * PyWrapper_value( PyWrapper * const self, PyObject * const args )
{
  VERIFY_NON_NULL_WRAPPER_AND_RETURN( self );

  int modify;
  if ( !PyArg_ParseTuple( args, "p", &modify ) )
  { return nullptr; }

  PyObject * const ret = self->wrapper->createPythonObject( modify );
  if ( ret == nullptr && PyErr_Occurred() == nullptr )
  { Py_RETURN_NONE; }

  return ret;
}


// Allow mixing designated and non-designated initializers in the same initializer list.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wc99-designator"

static PyMethodDef PyWrapperMethods[] = {
  { "value", (PyCFunction) PyWrapper_value, METH_VARARGS,
    "" },
  { nullptr, nullptr, 0, nullptr }        /* Sentinel */
};

static PyTypeObject PyWrapperType = {
  PyVarObject_HEAD_INIT( nullptr, 0 )
  .tp_name = "Wrapper",
  .tp_basicsize = sizeof( PyWrapper ),
  .tp_itemsize = 0,
  .tp_repr = PyWrapper_repr,
  .tp_flags = Py_TPFLAGS_DEFAULT,
  .tp_doc = "",
  .tp_methods = PyWrapperMethods,
  .tp_new = PyType_GenericNew,
};

#pragma GCC diagnostic pop

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
PyObject * createNewPyWrapper( dataRepository::WrapperBase & wrapper )
{
  // Create a new Group and set the dataRepository::Group it points to.
  PyObject * const ret = PyObject_CallFunction( reinterpret_cast< PyObject * >( getPyWrapperType() ), "" );
  PyWrapper * const retWrapper = reinterpret_cast< PyWrapper * >( ret );
  if ( retWrapper == nullptr )
  { return nullptr; }

  retWrapper->wrapper = &wrapper;

  return ret;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
PyTypeObject * getPyWrapperType()
{ return &PyWrapperType; }

} // namespace python
} // namespace geosx
