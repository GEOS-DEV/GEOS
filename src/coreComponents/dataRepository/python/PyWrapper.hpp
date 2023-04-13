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

#ifndef GEOSX_PYTHON_PYWRAPPER_HPP_
#define GEOSX_PYTHON_PYWRAPPER_HPP_

// Source includes
#include "dataRepository/WrapperBase.hpp"

namespace geos
{
namespace python
{

/**
 *
 */
PyObject * createNewPyWrapper( dataRepository::WrapperBase & wrapper );

/**
 *
 */
PyTypeObject * getPyWrapperType();

} // namespace python
} // namespace geos

#endif
