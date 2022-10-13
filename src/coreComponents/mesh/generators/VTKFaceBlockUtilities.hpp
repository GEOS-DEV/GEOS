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

namespace geosx {

/**
 * @brief Import face block @p faceBlockName from @p vtkMesh into the @p cellBlockManager.
 * @param faceBlockName
 * @param vtkMesh
 * @param cellBlockManager
 * @deprecated This is most likely not entirely the place to do this.
 */
void ImportFracture( string const & faceBlockName,
                     vtkSmartPointer< vtkDataSet > vtkMesh,
                     CellBlockManager & cellBlockManager );

}

#endif // include guard
