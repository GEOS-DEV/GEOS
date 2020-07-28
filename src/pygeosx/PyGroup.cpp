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
#include "PyGroup.hpp"
#include "pygeosx.hpp"
#include "PyWrapper.hpp"

#define VERIFY_NON_NULL_AND_RETURN( obj, retvalue ) \
  if ( obj == nullptr ) \
  { \
    PyErr_SetString( PyExc_RuntimeError, "Passed a nullptr as an argument" ); \
    return retvalue; \
  } \

#define VERIFY_NON_NULL_GROUP_AND_RETURN( pygroup ) \
  if ( pygroup->group == nullptr ) \
  { \
    PyErr_SetString( PyExc_AttributeError, "group has not been initialized" ); \
    return nullptr; \
  } \

namespace geosx
{
namespace python
{

/**
 *
 */
struct PyGroup
{
  PyObject_HEAD
  dataRepository::Group * group;
};

/**
 *
 */
static dataRepository::Group * getGroup( PyObject * const obj )
{
  VERIFY_NON_NULL_AND_RETURN( obj, nullptr );

  int isInstanceOfPyGroup = PyObject_IsInstance( obj, reinterpret_cast< PyObject * >( getPyGroupType() ) );
  if ( isInstanceOfPyGroup < 0 )
  { return nullptr; }

  if ( isInstanceOfPyGroup == 0 )
  {
    PyErr_SetString( PyExc_AttributeError, "Expect an argument of type Group." );
    return nullptr;
  }

  PyGroup * group = reinterpret_cast< PyGroup * >( obj );
  VERIFY_NON_NULL_GROUP_AND_RETURN( group );

  return group->group;
}

/**
 *
 */
static PyObject * PyGroup_repr( PyObject * const obj )
{
  dataRepository::Group const * const group = getGroup( obj );
  if ( group == nullptr )
  { return nullptr; }

  std::string const path = group->getPath();
  std::string const type = LvArray::demangle( typeid( *group ).name() );
  std::string const repr = path + " ( " + type + " )";
  return PyUnicode_FromString( repr.c_str() );
}

/**
 *
 */
static PyObject * PyGroup_groups( PyGroup * const self, PyObject * const args )
{
  GEOSX_UNUSED_VAR( args );

  VERIFY_NON_NULL_AND_RETURN( self, nullptr );
  VERIFY_NON_NULL_GROUP_AND_RETURN( self );

  localIndex const numSubGroups = self->group->numSubGroups();

  PyObject * pyList = PyList_New( numSubGroups );
  Py_ssize_t i = 0;
  bool error = false;
  self->group->forSubGroups( [pyList, &i, &error] ( dataRepository::Group & subGroup )
  {
    PyObject * const retGroup = createNewPyGroup( subGroup );
    error = error || retGroup == nullptr;

    PyList_SET_ITEM( pyList, i, retGroup );
    ++i;
  } );

  if ( error )
  {
    Py_DECREF( pyList );
    return nullptr;
  }

  return pyList;
}

/**
 *
 */
static PyObject * PyGroup_wrappers( PyGroup * const self, PyObject * const args )
{
  GEOSX_UNUSED_VAR( args );

  VERIFY_NON_NULL_AND_RETURN( self, nullptr );
  VERIFY_NON_NULL_GROUP_AND_RETURN( self );

  localIndex const numWrappers = self->group->numWrappers();

  PyObject * pyList = PyList_New( numWrappers );
  Py_ssize_t i = 0;
  bool error = false;
  self->group->forWrappers( [pyList, &i, &error] ( dataRepository::WrapperBase & wrapper )
  {
    PyObject * const retWrapper = createNewPyWrapper( wrapper );
    error = error || retWrapper == nullptr;

    PyList_SET_ITEM( pyList, i, retWrapper );
    ++i;
  } );

  if ( error )
  {
    Py_DECREF( pyList );
    return nullptr;
  }

  return pyList;
}

/**
 *
 */
static PyObject * PyGroup_getGroup( PyGroup * const self, PyObject * const args )
{
  VERIFY_NON_NULL_AND_RETURN( self, nullptr );
  VERIFY_NON_NULL_GROUP_AND_RETURN( self );

  PyObject * unicodePath;
  if ( !PyArg_ParseTuple( args, "U", &unicodePath ) )
  { return nullptr; }

  PyObjectRef asciiPath { PyUnicode_AsASCIIString( unicodePath ) };
  if ( asciiPath == nullptr )
  { return nullptr; }

  char const * const path = PyBytes_AsString( asciiPath );
  if ( path == nullptr )
  { return nullptr; }

  dataRepository::Group * const result = self->group->GetGroupByPath( path );
  if ( result == nullptr )
  {
    GEOSX_LOG_RANK( "Group " << self->group->getPath() << "/" << path << " does not exist." );
    Py_RETURN_NONE;
  }

  // Create a new Group and set the dataRepository::Group it points to.
  return createNewPyGroup( *result );
}


// Allow mixing designated and non-designated initializers in the same initializer list.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wc99-designator"

static PyMethodDef PyGroup_methods[] = {
  { "groups", (PyCFunction) PyGroup_groups, METH_NOARGS,
    "" },
  { "wrappers", (PyCFunction) PyGroup_wrappers, METH_NOARGS,
    "" },
  { "getGroup", (PyCFunction) PyGroup_getGroup, METH_VARARGS,
    "" },
  { nullptr, nullptr, 0, nullptr }        /* Sentinel */
};

static PyTypeObject PyGroupType = {
  PyVarObject_HEAD_INIT( nullptr, 0 )
  .tp_name = "Group",
  .tp_basicsize = sizeof( PyGroup ),
  .tp_itemsize = 0,
  .tp_repr = PyGroup_repr,
  .tp_flags = Py_TPFLAGS_DEFAULT,
  .tp_doc = "",
  .tp_methods = PyGroup_methods,
  .tp_new = PyType_GenericNew,
};

#pragma GCC diagnostic pop

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
PyObject * createNewPyGroup( dataRepository::Group & group )
{
  // Create a new Group and set the dataRepository::Group it points to.
  PyObject * const ret = PyObject_CallFunction( reinterpret_cast< PyObject * >( getPyGroupType() ), "" );
  PyGroup * const retGroup = reinterpret_cast< PyGroup * >( ret );
  if ( retGroup == nullptr )
  { return nullptr; }

  retGroup->group = &group;

  return ret;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
PyTypeObject * getPyGroupType()
{ return &PyGroupType; }

} // namespace python
} // namespace geosx
