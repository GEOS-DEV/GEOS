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
 * @file Coarsening.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_MULTISCALE_COARSENING_HPP_
#define GEOSX_LINEARALGEBRA_MULTISCALE_COARSENING_HPP_

#include "common/DataTypes.hpp"
#include "utilities/LinearSolverParameters.hpp"

namespace geosx
{
namespace multiscale
{

class MeshLevel;

namespace coarsening
{

void buildCoarseMesh( MeshLevel & fineMesh,
                      MeshLevel & coarseMesh,
                      LinearSolverParameters::Multiscale::Coarsening const & coarse_params,
                      array1d< string > const & boundaryNodeSets );

} // namespace coarsening
} // namespace multiscale
} // namespace geosx

#endif //GEOSX_LINEARALGEBRA_MULTISCALE_COARSENING_HPP_
