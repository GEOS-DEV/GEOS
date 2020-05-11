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

#ifndef GEOSX_CORECOMPONENTS_WELLS_LASWELLGENERATOR_HPP_
#define GEOSX_CORECOMPONENTS_WELLS_LASWELLGENERATOR_HPP_


#include "dataRepository/Group.hpp"

#include "WellGeneratorBase.hpp"

namespace geosx
{

namespace dataRepository
{
namespace keys
{
string const fileName                = "fileName";
string const geometryLogIndexInFile  = "geometryLogIndexInFile";
}
}

class LASFile;
class LASLine;

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
  virtual void GenerateElementRegions( DomainPartition& GEOSX_UNUSED_PARAM( domain ) ) override {}

  virtual void GeneratePolyLine() override final;

  /// not implemented 
  virtual void GetElemToNodesRelationInBox ( std::string const & GEOSX_UNUSED_PARAM( elementType ),
                                             int const * GEOSX_UNUSED_PARAM( index ),
                                             int const & GEOSX_UNUSED_PARAM( iEle ),
                                             int * GEOSX_UNUSED_PARAM( nodeIDInBox ),
                                             int const GEOSX_UNUSED_PARAM( size )) override {}

  /// not implemented
  virtual void RemapMesh ( dataRepository::Group * const GEOSX_UNUSED_PARAM( domain ) ) override {}

  protected:
  void PostProcessInput() override final;

  private:
  void GeneratePolyLineFromXYZ( LASFile const & lasFile );

  void GeneratePolyLineFromDepth( LASFile const & lasFile );

  /*!
   * @brief Get the factor to conver input unit into meters
   */
  real64 GetFactor( LASLine const & lasLine );
  private:

  /// Path to the LAS file
  string m_fileName;

  /// Index of the log to take to write the well geometry if
  /// the LAS file has several time the same sections
  integer m_logIndexToTakeForGeometry;
};

} // namespace

#endif
