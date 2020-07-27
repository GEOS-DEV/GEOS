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
 * @file GMSHMeshGenerator.hpp
 */

#ifndef GEOSX_MESHUTILITIES_GMSHMESHGENERATOR_HPP
#define GEOSX_MESHUTILITIES_GMSHMESHGENERATOR_HPP

#include "dataRepository/Group.hpp"
#include "codingUtilities/Utilities.hpp"
#include "codingUtilities/StringUtilities.hpp"


#include "MeshGeneratorBase.hpp"

namespace geosx
{


/**
 *  @class GMSHMeshGenerator
 *  @brief The GMSHMeshGenerator class provides a class implementation of GMSH generated meshes.
 */
class GMSHMeshGenerator : public MeshGeneratorBase
{
public:
/**
 * @brief Main constructor for MeshGenerator base class.
 * @param[in] name of the GMSHMeshGenerator object
 * @param[in] parent the parent Group pointer for the MeshGenerator object
 */
  GMSHMeshGenerator( const std::string & name,
                       Group * const parent );

  virtual ~GMSHMeshGenerator() override;

/**
 * @brief Return the name of the GMSHMeshGenerator in object Catalog.
 * @return string that contains the key name to GMSHMeshGenerator in the Catalog
 */
  static string CatalogName() { return "GMSHMeshGenerator"; }

///@cond DO_NOT_DOCUMENT
  struct viewKeyStruct
  {
    constexpr static auto filePathString = "file";
    constexpr static auto scaleString = "scale";
    constexpr static auto fieldsToImportString = "fieldsToImport";
    constexpr static auto fieldNamesInGEOSXString = "fieldNamesInGEOSX";
    constexpr static auto reverseZString = "reverseZ";
    constexpr static auto initNbOfProcString = "initNbOfProc";
  };
/// @endcond

  virtual void GenerateElementRegions( DomainPartition & domain ) override;

  /**
   * @brief Create a new geometric object (box, plane, etc) as a child of this group.
   * @param childKey the catalog key of the new geometric object to create
   * @param childName the name of the new geometric object in the repository
   * @return the group child
   */
  virtual Group * CreateChild( string const & childKey, string const & childName ) override;

  virtual void GenerateMesh( DomainPartition * const domain ) override;

  virtual void GetElemToNodesRelationInBox ( const std::string & elementType,
                                             const int index[],
                                             const int & iEle,
                                             int nodeIDInBox[],
                                             const int size ) override;

  virtual void RemapMesh ( dataRepository::Group * const domain ) override;

protected:

  /**
   * @brief This function provides capability to post process input values prior to
   * any other initialization operations.
   */
  void PostProcessInput() override final;

private:
  void GetLine(std::ifstream & fileStream, std::string& line);
private:
  /// Names of the fields to be copied from GMSH to GEOSX data structure
  string_array m_fieldsToImport;

  /// Path to the mesh file
  Path m_filePath;

  /// Scale factor that will be applied to the point coordinates
  real64 m_scale;

  /// String array of the GEOSX user decalred fields
  string_array m_fieldNamesInGEOSX;

  /// z pointing direction flag, 0 (default) is upward, 1 is downward
  int m_isZReverse;

  /// Number of processor on which the mesh will be loaded first before
  ///being distributed
  int m_initNbOfProc;

  /// To follow line number while reading the file
  globalIndex m_lineNumber;
  const std::unordered_map< localIndex, localIndex > m_gmshCellTypeToNbNodes
    =
    {
      { 4,4 },
      { 5,8 },
      { 6,6 },
      { 7,5 },
    };
};

}

#endif /* GEOSX_MESHUTILITIES_GMSHMESHGENERATOR_HPP */
