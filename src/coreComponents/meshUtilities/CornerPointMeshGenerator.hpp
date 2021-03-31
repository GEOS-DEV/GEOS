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
 * @file CornerPointMeshGenerator.hpp
 */

#ifndef GEOSX_MESHUTILITIES_CORNERPOINTGENERATOR_HPP
#define GEOSX_MESHUTILITIES_CORNERPOINTGENERATOR_HPP

#include "dataRepository/Group.hpp"
#include "codingUtilities/Utilities.hpp"
#include "codingUtilities/StringUtilities.hpp"
#include "meshUtilities/MeshGeneratorBase.hpp"
#include "meshUtilities/CPMesh/CPMeshBuilder.hpp"

namespace geosx
{

/**
 *  @class CornerPointMeshGenerator
 *  @brief Class that reads and processes a corner point grid
 */
class CornerPointMeshGenerator : public MeshGeneratorBase
{
public:
/**
 * @brief Main constructor for MeshGenerator base class.
 * @param[in] name of the CornerPointMeshGenerator object
 * @param[in] parent the parent Group pointer for the MeshGenerator object
 */
  CornerPointMeshGenerator( const string & name,
                            Group * const parent );

  virtual ~CornerPointMeshGenerator() override;

/**
 * @brief Return the name of the CornerPointMeshGenerator in object Catalog.
 * @return string that contains the key name to CornerPointMeshGenerator in the Catalog
 */
  static string catalogName() { return "CornerPointMeshGenerator"; }

///@cond DO_NOT_DOCUMENT
  struct viewKeyStruct
  {
    constexpr static char const * filePathString() { return "file"; }
  };
/// @endcond

  virtual void generateElementRegions( DomainPartition & domain ) override;

  /**
   * @brief Create a new geometric object (box, plane, etc) as a child of this group.
   * @param childKey the catalog key of the new geometric object to create
   * @param childName the name of the new geometric object in the repository
   * @return the group child
   */
  virtual Group * createChild( string const & childKey,
                               string const & childName ) override;

  virtual void generateMesh( DomainPartition & domain ) override;

  virtual void getElemToNodesRelationInBox ( const string & elementType,
                                             const int index[],
                                             const int & iEle,
                                             int nodeIDInBox[],
                                             const int size ) override;

  virtual void remapMesh ( dataRepository::Group & domain ) override;

protected:

  /**
   * @brief This function provides capability to post process input values prior to
   * any other initialization operations.
   */
  void postProcessInput() override final;


private:

  /// Driver class for the construction of the conformal corner-point mesh
  std::unique_ptr< CPMesh::CPMeshBuilder > m_cPMeshBuilder;

  /// Path to the mesh file
  Path m_filePath;

};

}

#endif /* GEOSX_MESHUTILITIES_CORNERPOINTGENERATOR_HPP */
