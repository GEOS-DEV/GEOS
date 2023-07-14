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

namespace geos::vtk
{

/**
 * @brief Attach the face block information to the cell block manager.
 * @param faceBlockName The name of the face block.
 * @param faceMesh The vtk mesh for the face block.
 * @param mesh The vtk volumic mesh.
 * @param cellBlockManager The cell block manager that will receive the face block information.
 */
void importFractureNetwork( string const & faceBlockName,
                            vtkSmartPointer< vtkDataSet > faceMesh,
                            vtkSmartPointer< vtkDataSet > mesh,
                            CellBlockManager & cellBlockManager );
}

#endif // include guard
