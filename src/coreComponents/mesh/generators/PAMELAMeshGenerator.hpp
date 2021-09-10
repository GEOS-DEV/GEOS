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
 * @file PAMELAMeshGenerator.hpp
 */

#ifndef GEOSX_MESH_GENERATORS_PAMELAMESHGENERATOR_HPP
#define GEOSX_MESH_GENERATORS_PAMELAMESHGENERATOR_HPP

#include "mesh/generators/MeshGeneratorBase.hpp"

namespace PAMELA
{
/// Forward declare Mesh class for unique_ptr member
class Mesh;
}

namespace geosx
{

/**
 *  @class PAMELAMeshGenerator
 *  @brief The PAMELAMeshGenerator class provides a class implementation of PAMELA generated meshes.
 */
class PAMELAMeshGenerator final : public MeshGeneratorBase
{
public:
/**
 * @brief Main constructor for MeshGenerator base class.
 * @param[in] name of the PAMELAMeshGenerator object
 * @param[in] parent the parent Group pointer for the MeshGenerator object
 */
  PAMELAMeshGenerator( const string & name,
                       Group * const parent );

/**
 * @brief Return the name of the PAMELAMeshGenerator in object Catalog.
 * @return string that contains the key name to PAMELAMeshGenerator in the Catalog
 */
  static string catalogName() { return "PAMELAMeshGenerator"; }

///@cond DO_NOT_DOCUMENT
  struct viewKeyStruct
  {
    constexpr static char const * filePathString() { return "file"; }
    constexpr static char const * scaleString() { return "scale"; }
    constexpr static char const * fieldsToImportString() { return "fieldsToImport"; }
    constexpr static char const * fieldNamesInGEOSXString() { return "fieldNamesInGEOSX"; }
    constexpr static char const * reverseZString() { return "reverseZ"; }
  };
/// @endcond

  /**
   * @brief Create a new geometric object (box, plane, etc) as a child of this group.
   * @param childKey the catalog key of the new geometric object to create
   * @param childName the name of the new geometric object in the repository
   * @return the group child
   */
  virtual Group * createChild( string const & childKey, string const & childName ) override;

  virtual void generateMesh( DomainPartition & domain ) override;

  virtual void importFields( DomainPartition & domain ) const override;

  virtual void freeResources() override;

protected:

  /**
   * @brief This function provides capability to post process input values prior to
   * any other initialization operations.
   */
  void postProcessInput() final;


private:

  /// Unique Pointer to the Mesh in the data structure of PAMELA.
  std::unique_ptr< PAMELA::Mesh > m_pamelaMesh;

  /// Names of the fields to be copied from PAMELA to GEOSX data structure
  string_array m_fieldsToImport;

  /// Path to the mesh file
  Path m_filePath;

  /// Scale factor that will be applied to the point coordinates
  real64 m_scale;

  /// String array of the GEOSX user declared fields
  string_array m_fieldNamesInGEOSX;

  /// z pointing direction flag, 0 (default) is upward, 1 is downward
  integer m_isZReverse;

  /// Map of cell block (subregion) names to PAMELA region names
  std::unordered_map< string, string > m_cellBlockRegions;
};

}

#endif /* GEOSX_MESH_GENERATORS_PAMELAMESHGENERATOR_HPP */
