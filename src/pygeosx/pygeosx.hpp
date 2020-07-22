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

// Source includes
#include "PyGroup.hpp"
#include "LvArray/src/python/pythonForwardDeclarations.hpp"
#include "managers/GeosxState.hpp"

namespace geosx
{

// TODO: corbett This should go in LvArray
class PyObjectRef
{
public:

  /**
   * @brief Create an uninitialized (nullptr) reference.
   */
  PyObjectRef() = default;

  /**
   * @brief Take ownership of a reference to @p src.
   * @p src The object to be referenced.
   */
  explicit PyObjectRef( PyObject * const src ):
    m_object( src )
  {}

  // Theese could be implemented but I don't have a use case yet.
  PyObjectRef( PyObjectRef const & ) = delete;
  PyObjectRef( PyObjectRef && ) = delete;

  /**
   * @brief Decrease the reference count to the current object.
   */
  ~PyObjectRef()
  {
    if ( m_object != nullptr )
    { Py_DECREF( m_object ); }
  }

  PyObjectRef & operator=( PyObjectRef const & ) = delete;
  PyObjectRef & operator=( PyObjectRef && ) = delete;

  /**
   * @brief Decrease the reference count to the current object and take ownership
   *   of a new reference.
   * @p src The new object to be referenced.
   * @return *this.
   */
  PyObjectRef & operator=( PyObject * src )
  {
    if ( m_object != nullptr )
    { Py_DECREF( m_object ); }

    m_object = src;
    return *this;
  }

  operator PyObject*()
  { return m_object; }

  PyObject * get() const
  { return m_object; }

  PyObject * release()
  { 
    PyObject * const ret = m_object;
    m_object = nullptr;
    return ret;
  }

private:
  PyObject * m_object = nullptr;
};

std::unique_ptr< GeosxState > & getState();

} // namespace geosx
