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
 * @file MeshGeneratorBase.h
 */

#ifndef MESHGENERATORBASE_H_
#define MESHGENERATORBASE_H_

#include "dataRepository/Group.hpp"
#include "codingUtilities/Utilities.hpp"
#include "common/DataTypes.hpp"

namespace geosx
{

namespace dataRepository
{}

class NodeManager;
class DomainPartition;

class MeshGeneratorBase : public dataRepository::Group
{
public:
  explicit MeshGeneratorBase( std::string const & name,
                              Group * const parent );

  virtual ~MeshGeneratorBase();

  static string CatalogName() { return "MeshGeneratorBase"; }

  virtual void GenerateElementRegions( DomainPartition & domain ) = 0;

  virtual void GenerateMesh( DomainPartition * const domain ) = 0;

  // virtual void GenerateNodesets( xmlWrapper::xmlNode const & targetNode,
  //                                NodeManager * nodeManager ) = 0;

  virtual void GetElemToNodesRelationInBox ( const std::string & elementType,
                                             const int index[],
                                             const int & iEle,
                                             int nodeIDInBox[],
                                             const int size ) = 0;

  virtual void RemapMesh ( dataRepository::Group * const domain ) = 0;

  int m_delayMeshDeformation = 0;

  using CatalogInterface = dataRepository::CatalogInterface< MeshGeneratorBase, std::string const &, Group * const >;
  static CatalogInterface::CatalogType & GetCatalog();

};
}

#endif /* MESHGENERATORBASE_H_ */
