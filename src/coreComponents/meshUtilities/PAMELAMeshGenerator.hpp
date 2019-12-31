/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file PAMELAMeshGenerator.cpp
 */

#pragma once

#include "dataRepository/Group.hpp"
#include "codingUtilities/Utilities.hpp"

//This is an include of PAMELA
#include "Mesh/Mesh.hpp"
#include "MeshDataWriters/Writer.hpp"

#include "MeshGeneratorBase.hpp"

namespace geosx
{

class PAMELAMeshGenerator : public MeshGeneratorBase
{
public:
  PAMELAMeshGenerator( const std::string& name,
                       Group * const parent );

  virtual ~PAMELAMeshGenerator() override;

  static string CatalogName() { return "PAMELAMeshGenerator"; }

  struct viewKeyStruct
  {
    constexpr static auto filePathString = "file";
    constexpr static auto scaleString = "scale";
    constexpr static auto fieldsToImportString = "fieldsToImport";
    constexpr static auto fieldNamesInGEOSXString = "fieldNamesInGEOSX";
    constexpr static auto reverseZString = "reverseZ";
  };

  virtual void GenerateElementRegions( DomainPartition& domain ) override;

  virtual Group * CreateChild( string const & childKey, string const & childName ) override;

  virtual void GenerateMesh( DomainPartition * const domain ) override;

  virtual void GetElemToNodesRelationInBox ( const std::string& elementType,
                                             const int index[],
                                             const int& iEle,
                                             int nodeIDInBox[],
                                             const int size ) override;

  virtual void RemapMesh ( dataRepository::Group * const domain ) override;

protected:
  void PostProcessInput() override final;

private:

  /// Mesh in the data structure of PAMELA.
  std::unique_ptr< PAMELA::Mesh >  m_pamelaMesh;

  /// Names of the fields to be copied from PAMELA to GEOSX data structure
  string_array m_fieldsToImport;

  /// Path to the mesh file
  Path m_filePath;

  /// Scale factor that will be applied to the point coordinates
  real64 m_scale;

  string_array m_fieldNamesInGEOSX;

  int m_isZReverse;

  const std::unordered_map<PAMELA::ELEMENTS::TYPE, string, PAMELA::ELEMENTS::EnumClassHash> ElementToLabel
    =
    {
    { PAMELA::ELEMENTS::TYPE::VTK_VERTEX, "VERTEX"},
    { PAMELA::ELEMENTS::TYPE::VTK_LINE, "LINE"  },
    { PAMELA::ELEMENTS::TYPE::VTK_TRIANGLE, "TRIANGLE" },
    { PAMELA::ELEMENTS::TYPE::VTK_QUAD, "QUAD" },
    { PAMELA::ELEMENTS::TYPE::VTK_TETRA, "TETRA" },
    { PAMELA::ELEMENTS::TYPE::VTK_HEXAHEDRON, "HEX" },
    { PAMELA::ELEMENTS::TYPE::VTK_WEDGE, "WEDGE" },
    { PAMELA::ELEMENTS::TYPE::VTK_PYRAMID, "PYRAMID" }
    };
};

}
