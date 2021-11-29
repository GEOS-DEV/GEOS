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
 * @file ElementType.hpp
 */

#ifndef GEOSX_MESH_ELEMENTTYPE_HPP
#define GEOSX_MESH_ELEMENTTYPE_HPP

#include "codingUtilities/EnumStrings.hpp"

namespace geosx
{

/**
 * @brief Denotes type of cell/element shape
 */
enum class ElementType : integer
{
  Line,          ///< Two-node line segment
  Triangle,      ///< Three-node triangle
  Quadrilateral, ///< Four-node quadrilateral
  Polygon,       ///< General polygonal element
  Tetrahedron,   ///< Four-node tetrahedral element
  Pyramid,       ///< Five-node pyramid element
  Prism,         ///< Six-node wedge element
  Hexahedron,    ///< Eight-node hexahedral element
  Polyhedron     ///< General polyhedral element
};

/// Strings for ElementType
ENUM_STRINGS( ElementType,
              "BEAM",
              "C2D3",
              "C2D4",
              "Polygon",
              "C3D4",
              "C3D5",
              "C3D6",
              "C3D8",
              "Polyhedron" );

} // namespace geosx

#endif //GEOSX_MESH_ELEMENTTYPE_HPP
