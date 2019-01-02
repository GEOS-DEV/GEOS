/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/*
 * PAMELAMeshGenerator.cpp
 *
 *  Created on: Oct 08, 2018
 *      Author: Antoine Mazuyer
 */

#pragma once

#include "dataRepository/ManagedGroup.hpp"
#include "codingUtilities/Utilities.hpp"

//This is an include of PAMELA
#include "Mesh/Mesh.hpp"

#include "MeshGeneratorBase.hpp"

#ifdef USE_ATK
#include <slic/slic.hpp>
#endif

namespace geosx
{

namespace dataRepository
{
namespace keys
{
string const filePath = "file";
}
}


class PAMELAMeshGenerator : public MeshGeneratorBase
{
public:
  PAMELAMeshGenerator( const std::string& name,
                       ManagedGroup * const parent );

  virtual ~PAMELAMeshGenerator() override;

  static string CatalogName() { return "PAMELAMeshGenerator"; }


  virtual void GenerateElementRegions( DomainPartition& domain ) override;

  virtual ManagedGroup * CreateChild( string const & childKey, string const & childName ) override;

  virtual void GenerateMesh( dataRepository::ManagedGroup * const domain ) override;

  virtual void GetElemToNodesRelationInBox ( const std::string& elementType,
                                             const int index[],
                                             const int& iEle,
                                             int nodeIDInBox[],
                                             const int size ) override;

  virtual void RemapMesh ( dataRepository::ManagedGroup * const domain ) override;

  void ProcessInputFile_PostProcess() override final;

private:
  /// Mesh in the data structure of PAMELA.
  std::unique_ptr< PAMELA::Mesh >  m_pamelaMesh;

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
