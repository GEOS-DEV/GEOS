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

#ifndef GEOS_PYTHON_PYGROUP_HPP_
#define GEOS_PYTHON_PYGROUP_HPP_

// Source includes
#include "dataRepository/Group.hpp"
#include "LvArray/src/python/pythonForwardDeclarations.hpp"

#include "PyWrapper.hpp"

namespace geos
{
namespace python
{



PyObject * createNewPyGroup( dataRepository::Group & group );



} // namespace python
} // namespace geos

#endif
