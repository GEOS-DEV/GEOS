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

#ifndef GEOS_PYTHON_PYWRAPPER_HPP_
#define GEOS_PYTHON_PYWRAPPER_HPP_

// Source includes
#include "dataRepository/WrapperBase.hpp"

namespace geos
{
namespace python
{


/**
 * @brief Python wrapper object for geos::dataRepository::WrapperBase.
 *
 * This struct defines the Python-facing interface for accessing and modifying
 * wrapped C++ values from Python.
 */
struct PyWrapper
{
  PyObject_HEAD

  /// Documentation string for the PyWrapper Python type.
  static constexpr char const * docString =
    "A Python interface to geos::dataRepository::WrapperBase.";

  /// Pointer to the underlying C++ WrapperBase instance.
  dataRepository::WrapperBase * wrapper;
};


/**
 *
 */
PyObject * createNewPyWrapper( dataRepository::WrapperBase & wrapper );

/**
 *
 */
PyTypeObject * getPyWrapperType();

/**
 * @brief Set the value of a scalar wrapper from Python.
 *
 * This function allows assigning a scalar value from Python into the underlying
 * C++ Wrapper<T> instance. Supported types include:
 * - std::string
 * - int
 * - long
 * - double
 *
 * If the type of the Python object does not match the wrapped C++ type,
 * a Python TypeError is raised. Same for unsupported types.
 *
 * @param self Pointer to the PyWrapper instance.
 * @param args A Python tuple containing a single value to assign.
 * @return Py_None on success, or nullptr on failure with a Python exception set.
 */
PyObject * PyWrapper_setValue( PyWrapper * self, PyObject * args );


} // namespace python
} // namespace geos

#endif
