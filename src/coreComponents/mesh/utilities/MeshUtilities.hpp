/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file MeshUtilities.hpp
 */

#ifndef GEOSX_MESH_UTILITIES_MESHUTILITIES_HPP
#define GEOSX_MESH_UTILITIES_MESHUTILITIES_HPP

#include "common/DataTypes.hpp"

namespace geosx
{

namespace dataRepository
{
class Group;
}
class NodeManager;


/**
 * @class MeshUtilities
 * @brief This class is used to generate the utilities for the mesh.
 */
namespace MeshUtilities
{
/**
 * @brief Build all the node sets from a geometric object in the DomainPartition.
 * @param[in] geometry a pointer to the group in the data repository
 * @param[in] nodeManager pointer to the NodeManager object in the DomainPartition
 */
void generateNodesets( dataRepository::Group const & geometry,
                       NodeManager & nodeManager );

} // namespace MeshUtilities

} // namespace geosx

#endif /* GEOSX_MESH_UTILITIES_MESHUTILITIES_HPP */
