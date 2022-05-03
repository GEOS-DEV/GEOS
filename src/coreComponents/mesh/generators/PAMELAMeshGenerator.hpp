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

#include "mesh/generators/ExternalMeshGeneratorBase.hpp"

#include "Elements/Polyhedron.hpp"

namespace PAMELA
{
/// Forward declare PAMELA types for use in member declarations
class Mesh;
template< typename > struct Part;
template< typename > struct SubPart;
}

namespace geosx
{

class ElementRegionBase;
class ElementSubRegionBase;

/**
 *  @class PAMELAMeshGenerator
 *  @brief The PAMELAMeshGenerator class provides a class implementation of PAMELA generated meshes.
 */
class PAMELAMeshGenerator final : public ExternalMeshGeneratorBase
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
  static string catalogName() { return "PAMELAMesh"; }

///@cond DO_NOT_DOCUMENT
  struct viewKeyStruct
  {
    constexpr static char const * filePathString() { return "file"; }
    constexpr static char const * scaleString() { return "scale"; }
    constexpr static char const * translateString() { return "translate"; }
    constexpr static char const * fieldsToImportString() { return "fieldsToImport"; }
    constexpr static char const * fieldNamesInGEOSXString() { return "fieldNamesInGEOSX"; }
  };
/// @endcond

  virtual void generateMesh( DomainPartition & domain ) override;

  virtual void importFields( DomainPartition & domain ) const override;

  virtual void freeResources() override;

private:

  void importFieldsOnSubRegion( PAMELA::Part< PAMELA::Polyhedron * > & srcRegion,
                                PAMELA::SubPart< PAMELA::Polyhedron * > & srcSubRegion,
                                ElementRegionBase const & dstRegion,
                                ElementSubRegionBase & dstSubRegion ) const;

  /// Unique Pointer to the Mesh in the data structure of PAMELA.
  std::unique_ptr< PAMELA::Mesh > m_pamelaMesh;

  /// Map of cell block (subregion) names to PAMELA region names
  std::unordered_map< string, string > m_cellBlockRegions;
};

}

#endif /* GEOSX_MESH_GENERATORS_PAMELAMESHGENERATOR_HPP */
