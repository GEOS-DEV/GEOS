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
  Vertex,        ///< Single-node vertex element
  Line,          ///< Two-node line segment
  Triangle,      ///< Three-node triangle
  Quadrilateral, ///< Four-node quadrilateral
  Polygon,       ///< General polygonal element
  Tetrahedron,   ///< Four-node tetrahedral element
  Pyramid,       ///< Five-node pyramid element
  Prism,         ///< Six-node wedge element
  Hexahedron,    ///< Eight-node hexahedral element
  Polyhedron,    ///< General polyhedral element
  // NOTE: If you add anything below Polyhedron,
  // don't forget to update numElementTypes() below.
};

/**
 * @brief @return number of supported element types
 * @note this MUST be updated if a new element type is inserted after Polyhedron
 */
inline constexpr integer numElementTypes()
{
  return static_cast< integer >( ElementType::Polyhedron ) + 1;
}

/**
 * @brief Get number of spatial dimensions of element type
 * @param elementType type of element
 * @return number of spatial dimensions (1-3)
 */
inline int getElementDim( ElementType const elementType )
{
  switch( elementType )
  {
    case ElementType::Vertex:        return 0;
    case ElementType::Line:          return 1;
    case ElementType::Triangle:
    case ElementType::Quadrilateral:
    case ElementType::Polygon:       return 2;
    case ElementType::Tetrahedron:
    case ElementType::Pyramid:
    case ElementType::Prism:
    case ElementType::Hexahedron:
    case ElementType::Polyhedron:    return 3;
  }
  return 0;
}

/// Strings for ElementType
ENUM_STRINGS( ElementType,
              "Vertex",
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
