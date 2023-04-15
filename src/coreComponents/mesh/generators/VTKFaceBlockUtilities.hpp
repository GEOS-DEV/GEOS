/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2020-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#ifndef GEOS_VTKFACEBLOCKUTILITIES_HPP
#define GEOS_VTKFACEBLOCKUTILITIES_HPP

#include "common/DataTypes.hpp"
#include "CellBlockManager.hpp"

#include <vtkDataSet.h>
#include <vtkSmartPointer.h>

namespace geos
{

/**
 * @brief Import face block @p faceBlockName from @p vtkMesh into the @p cellBlockManager.
 * @param[in] filePath Path to the multi-block vtk file.
 * @param[in] faceBlockName The face block name to include. (It's the name of the block in the multi-block file.)
 * @param[in] mesh The 3d vtk mesh.
 * @param[inout] cellBlockManager The face block instance (with name @p faceBlockName) will be attached to the @p cellBlockManager
 * @return The vtk mesh instance for the considered face block.
 */
vtkSmartPointer< vtkDataSet > importFractureNetwork( Path const & filePath,
                                                     string const & faceBlockName,
                                                     vtkSmartPointer< vtkDataSet > mesh,
                                                     CellBlockManager & cellBlockManager );

}

#endif // include guard
