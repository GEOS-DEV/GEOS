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
 * @file ExternalMeshGeneratorBase.hpp
 */

#ifndef GEOSX_MESH_GENERATORS_EXTERNALMESHGENERATORBASE_HPP
#define GEOSX_MESH_GENERATORS_EXTERNALMESHGENERATORBASE_HPP

#include "mesh/generators/MeshGeneratorBase.hpp"

namespace geosx
{

/**
 * @brief Base class for external mesh generators (importers).
 */
class ExternalMeshGeneratorBase : public MeshGeneratorBase
{
public:

  /**
   * @brief Constructor.
   * @param[in] name name of the object
   * @param[in] parent the parent Group pointer
   */
  ExternalMeshGeneratorBase( const string & name,
                             Group * const parent );

protected:

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

  void postProcessInput() override;

  /// Path to the mesh file
  Path m_filePath;

  /// Translation vector that will be applied to the point coordinates (prior to scaling)
  R1Tensor m_translate;

  /// Scale factor that will be applied to the point coordinates (after translation)
  R1Tensor m_scale;

  /// Names of the fields to be copied from PAMELA to GEOSX data structure
  array1d< string > m_fieldsToImport;

  /// String array of the GEOSX user declared fields
  array1d< string > m_fieldNamesInGEOSX;

};

} // namespace geosx

#endif //GEOSX_MESH_GENERATORS_EXTERNALMESHGENERATORBASE_HPP
