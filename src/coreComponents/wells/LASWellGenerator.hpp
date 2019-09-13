/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
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

#ifndef GEOSX_CORECOMPONENTS_WELLS_LASWELLGENERATOR_HPP_
#define GEOSX_CORECOMPONENTS_WELLS_LASWELLGENERATOR_HPP_

#include "fileIO/las/LASFile.hpp"

#include "dataRepository/Group.hpp"

#include "meshUtilities/MeshGeneratorBase.hpp"

namespace geosx
{

class LASWellGenerator : public MeshGeneratorBase
{
  public:
  LASWellGenerator( const std::string& name,
                    Group * const parent );

  /**
   * @brief default destructor
   */
  virtual ~LASWellGenerator() override;

  /**
   * @return the name of this type in the catalog
   */  
  static string CatalogName() { return "LASWell"; }

  /// not implemented
  virtual void GenerateElementRegions( DomainPartition& domain ) override {}

  virtual Group * CreateChild( string const & childKey, 
                               string const & childName ) override;


  /**
   * @brief main function of this class: processes the well input and creates the globla well topology
   * @param domain the physical domain object
   */  
  virtual void GenerateMesh( DomainPartition * const domain ) override;

  /// not implemented 
  virtual void GetElemToNodesRelationInBox ( const std::string& elementType,
                                             const int index[],
                                             const int& iEle,
                                             int nodeIDInBox[],
                                             const int size) override {}

  /// not implemented
  virtual void RemapMesh ( dataRepository::Group * const domain ) override {}

  protected:
  void PostProcessInput() override final;

};

} // namespace

#endif
