/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#pragma once

// Source includes
#include "dataRepository/Group.hpp"
#include "LvArray/src/python/pythonForwardDeclarations.hpp"
#include "PyWrapper.hpp"

namespace geosx
{
namespace python
{

/**
 *
 */


static constexpr char const * PyGroup_getWrapperDocString =
  "get_wrapper(self, path, default, /)\n"
  "--\n\n"
  "Return the `Wrapper` at the relative path ``path``; ``default`` is optional.\n\n"
  "If no wrapper exists and ``default`` is not given, raise a ``ValueError``;\n"
  "otherwise return ``default``.\n"
  "\n"
  "Parameters\n"
  "__________\n"
  "path : str\n"
  "    The relative path of the wrapper to return.\n"
  "\n"
  "Returns\n"
  "_______\n"
  "Wrapper\n"
  "    The wrapper at the relative path.";


template< typename T >
static PyObject * PyGroup_getWrapper( T * const self, PyObject * const args ) noexcept
{
  PYTHON_ERROR_IF( self == nullptr, PyExc_RuntimeError, "Passed a nullptr as self.", nullptr );

  PyObject * unicodePath;
  PyObject * defaultReturnValue = nullptr;
  if( !PyArg_ParseTuple( args, "U|O", &unicodePath, &defaultReturnValue ) )
  {
    return nullptr;
  }

  LvArray::python::PyObjectRef<> asciiPath { PyUnicode_AsASCIIString( unicodePath ) };
  if( asciiPath == nullptr )
  {
    return nullptr;
  }

  char const * const path = PyBytes_AsString( asciiPath );
  if( path == nullptr )
  {
    return nullptr;
  }

  string groupPath, wrapperName;
  std::tie( groupPath, wrapperName ) = splitPath( path );

  try
  {
    dataRepository::Group & group = self->group->getGroupByPath( groupPath );
    dataRepository::WrapperBase & wrapper = group.getWrapperBase( wrapperName );
    return createNewPyWrapper( wrapper );
  }
  // If the path isn't valid then either getGroupByPath or getWrapperBase will throw a std::domain_error
  catch( std::domain_error const & e )
  {
    // If no default return value was specified then this results in a Python exception.
    if( defaultReturnValue == nullptr )
    {
      PyErr_SetString( PyExc_KeyError, e.what() );
      return nullptr;
    }

    // Otherwise we return the default value.
    Py_INCREF( defaultReturnValue );
    return defaultReturnValue;
  }
}




PyObject * createNewPyGroup( dataRepository::Group & group );

/**
 *
 */
PyTypeObject * getPyGroupType();

} // namespace python
} // namespace geosx
