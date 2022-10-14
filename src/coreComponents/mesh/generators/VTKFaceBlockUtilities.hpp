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

#ifndef GEOSX_VTKFACEBLOCKUTILITIES_HPP
#define GEOSX_VTKFACEBLOCKUTILITIES_HPP

#include "common/DataTypes.hpp"
#include "CellBlockManager.hpp"

#include <vtkDataSet.h>
#include <vtkSmartPointer.h>

namespace geosx {

/**
 * @brief Import face block @p faceBlockName from @p vtkMesh into the @p cellBlockManager.
 * @param[in] faceBlockName The face block name to include. It's both the name in the vtk file, and it will be the face block name.
 * @param[in] vtkMesh The vtk mesh.
 * @param[inout] cellBlockManager The face block instance will be attached to the @p cellBlockManager
 */
void importFracture( string const & faceBlockName,
                     vtkSmartPointer< vtkDataSet > vtkMesh,
                     CellBlockManager & cellBlockManager );

}

#endif // include guard
