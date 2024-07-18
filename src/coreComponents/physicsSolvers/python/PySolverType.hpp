/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 Chevron 
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#ifndef GEOS_PYTHON_PYSOLVERTYPE_HPP_
#define GEOS_PYTHON_PYSOLVERTYPE_HPP_

#include "LvArray/src/python/pythonForwardDeclarations.hpp"

namespace geos
{
namespace python
{


PyTypeObject * getPySolverType();

} // namespace python
} // namespace geos

#endif
