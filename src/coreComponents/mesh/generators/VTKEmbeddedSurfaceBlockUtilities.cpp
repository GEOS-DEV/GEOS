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


#include "VTKEmbeddedSurfaceBlockUtilities.hpp"
#include "CellBlockManager.hpp"

namespace geos::vtk
{

void importFractureNetwork( string const & embeddedSurfaceBlockName,
                            vtkSmartPointer< vtkDataSet > embeddedSurfaceMesh,
                            vtkSmartPointer< vtkDataSet > mesh,
                            CellBlockManager & cellBlockManager ){
    
// prepare the data
vtkIdType const numEdfmSurf = embeddedSurfaceMesh->GetNumberOfCells();


EmbeddedSurfaceBlock & embeddedSurfBlock = cellBlockManager.registerEmbeddedSurfaceBlock(embeddedSurfaceBlockName);                         


embeddedSurfBlock.setNumEmbeddedSurfElem(2828);
embeddedSurfBlock.setEmbeddedSurfElemNodes(2828);
embeddedSurfBlock.setEmbeddedSurfElemNormalVectors(2828);
embeddedSurfBlock.setEmbeddedSurfElemTangentialLengthVectors(2828);
embeddedSurfBlock.setEmbeddedSurfElemTangentialWidthVectors(2828);
embeddedSurfBlock.setEmbeddedSurfElemTo3dElem(2828);
                           
                           
                           
                           }
}

