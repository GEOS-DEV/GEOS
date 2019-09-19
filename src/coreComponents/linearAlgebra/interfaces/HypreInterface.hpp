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

/**
 * @file HypreInterface.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_HYPREINTERFACE_HPP_
#define GEOSX_LINEARALGEBRA_HYPREINTERFACE_HPP_

//#include "HypreSparseMatrix.hpp"
//#include "HypreVector.hpp"
//#include "HypreSolver.hpp"

namespace geosx
{

/**
 * \class HypreInterface
 * \brief This class holds aliases based on the Hypre library.
 */

struct HypreInterface
{
  static void initialize( int & GEOSX_UNUSED_ARG( argc ), char ** & GEOSX_UNUSED_ARG( argv ) ) {}

  static void finalize() {}
};

} /* namespace geosx */

#endif /*GEOSX_LINEARALGEBRA_HYPREINTERFACE_HPP_*/
