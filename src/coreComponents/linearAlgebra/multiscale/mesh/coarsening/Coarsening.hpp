/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file Coarsening.hpp
 */

#ifndef GEOS_LINEARALGEBRA_MULTISCALE_COARSENING_HPP_
#define GEOS_LINEARALGEBRA_MULTISCALE_COARSENING_HPP_

#include "common/DataTypes.hpp"
#include "utilities/LinearSolverParameters.hpp"

namespace geos
{
namespace multiscale
{

class MeshLevel;

namespace coarsening
{

/**
 * @brief Perform the mesh coarsening and populate the coarse mesh data structure.
 * @param fineMesh fine-level mesh
 * @param coarseMesh coarse-level mesh to be filled
 * @param params coarsening parameters
 * @param boundaryNodeSets a list of global domain boundary node set names, used to improve coarse node selection
 */
void buildCoarseMesh( MeshLevel & fineMesh,
                      MeshLevel & coarseMesh,
                      LinearSolverParameters::Multiscale::Coarsening const & params,
                      string_array const & boundaryNodeSets );

} // namespace coarsening
} // namespace multiscale
} // namespace geos

#endif //GEOS_LINEARALGEBRA_MULTISCALE_COARSENING_HPP_
