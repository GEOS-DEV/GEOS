/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file RESQMLUtilities.hpp
 */

#ifndef GEOSX_MESH_GENERATORS_RESQMLUTILITIES_HPP
#define GEOSX_MESH_GENERATORS_RESQMLUTILITIES_HPP

#include "common/DataTypes.hpp"

#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDataSet.h>

#include "fesapi/resqml2/UnstructuredGridRepresentation.h"
namespace geosx
{

/**
 * @brief Load a RESQML UnstructuredGriRepresentation in a vtkUnstructuredGrid
 *
 * @param[in] grid The RESQML UnstructuredGriRepresentation
 * @return The loaded dataset
 */
vtkSmartPointer< vtkUnstructuredGrid >
loadUnstructuredGridRepresentation( RESQML2_NS::UnstructuredGridRepresentation *grid );

/**
 * @brief Load a Property in an existing dataset
 *
 * @param[in] dataset The existing dataset
 * @param[in] valuesProperty The RESQML Property
 * @param[in] fieldNameInGEOSX The name of property in GEOSX
 * @return The dataset with the loaded property
 */
vtkSmartPointer< vtkDataSet >
loadProperty( vtkSmartPointer< vtkDataSet > dataset, RESQML2_NS::AbstractValuesProperty *valuesProperty, string fieldNameInGEOSX );

/**
 * @brief Create a cell array of regions with an array of RESQML SubRepresentations
 *
 * @param dataset The existing dataset
 * @param regions The array of RESQML SubRepresentations
 * @param attributeName The name of the vtk cell array
 * @return The dataset with the loaded regions
 */
vtkSmartPointer< vtkDataSet >
createRegions( vtkSmartPointer< vtkDataSet > dataset, std::vector< RESQML2_NS::SubRepresentation * > regions, string attributeName );

/**
 * @brief Create as many surfaces in a dataset as RESQML SubRepresentations
 *
 * @param dataset The existing dataset
 * @param surfaces The array of RESQML SubRepresentations
 * @param regionAttributeName The name of the region array name
 * @return The dataset with the loaded surfaces
 */
vtkUnstructuredGrid *
createSurfaces( vtkSmartPointer< vtkDataSet > dataset, std::vector< RESQML2_NS::SubRepresentation * > surfaces, string regionAttributeName );


} // namespace geosx

#endif /* GEOSX_MESH_GENERATORS_RESQMLUTILITIES_HPP */
