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

#include "WellGeneratorBase.hpp"

namespace geosx
{

namespace dataRepository
{
namespace keys
{
string const fileName                     = "fileName";
string const geometryLogIndexInFile       = "geometryLogIndexInFile";
}
}

class LASWellGenerator : public WellGeneratorBase
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
  virtual void GenerateElementRegions( DomainPartition& GEOSX_UNUSED_ARG( domain ) ) override {}

  virtual void GeneratePolyLine() override final;

  /// not implemented 
  virtual void GetElemToNodesRelationInBox ( std::string const & GEOSX_UNUSED_ARG( elementType ),
                                             int const * GEOSX_UNUSED_ARG( index ),
                                             int const & GEOSX_UNUSED_ARG( iEle ),
                                             int * GEOSX_UNUSED_ARG( nodeIDInBox ),
                                             int const GEOSX_UNUSED_ARG( size )) override {}

  /// not implemented
  virtual void RemapMesh ( dataRepository::Group * const GEOSX_UNUSED_ARG( domain ) ) override {}

  protected:
  void PostProcessInput() override final;

  private:
  void GeneratePolyLineFromXYZ( LASFile const & lasFile );

  void GeneratePolyLineFromDepth( LASFile const & lasFile );

  /*!
   * @brief Get the factor to conver input unit into meters
   */
  double GetFactor( LASLine const & lasLine );
  private:

  /// Path to the LAS file
  string m_fileName;

  /// Index of the log to take to write the well geometry if
  /// the LAS file has several time the same sections
  integer m_logIndexToTakeForGeometry;

};

} // namespace

#endif
