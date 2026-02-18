/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file ComputationalGeometry.cpp
 */

#include "ComputationalGeometry.hpp"

namespace geos
{

namespace computationalGeometry
{

array2d<real64> computeVelocity(arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & normals, arrayView1d< real64 const > const & fluxes, ElementType& elem /*subRegion.getElementType()*/) 
{

switch (elem) {

  case geos::ElementType::Triangle:
    //TODO pre-post rotate
    GEOS_ERROR("Not Implemented Yet !");
    break;
    // return computeVelocities_<3>(normals,fluxes, rotation);
  case geos::ElementType::Quadrilateral:
  //TODO pre-post rotate
    GEOS_ERROR("Not Implemented Yet !");
    break;
    // return computeVelocities_<4>(normals,fluxes, rotation);
  case geos::ElementType::Tetrahedron:
    return computeVelocities_<4>(normals,fluxes);
  case geos::ElementType::Pyramid:
  case geos::ElementType::Prism5:
  case geos::ElementType::Wedge:
    return computeVelocities_<5>(normals,fluxes);
  case geos::ElementType::Hexahedron:
  case geos::ElementType::Prism6:
    return computeVelocities_<6>(normals,fluxes);
  case geos::ElementType::Prism7:
    return computeVelocities_<7>(normals,fluxes);
  case geos::ElementType::Prism8:
    return computeVelocities_<8>(normals,fluxes);
  case geos::ElementType::Prism9:
    return computeVelocities_<9>(normals,fluxes);
  case geos::ElementType::Prism10:
    return computeVelocities_<10>(normals,fluxes);
  case geos::ElementType::Prism11:
    return computeVelocities_<11>(normals,fluxes);
    

  case geos::ElementType::Polygon:
  case geos::ElementType::Polyhedron:
  default:
    GEOS_ERROR("Velocities are not computed on 2D Polygons cell type");
    break;

}

  return {};
}


} /* namespace computationalGeometry */

} /* namespace geos */
