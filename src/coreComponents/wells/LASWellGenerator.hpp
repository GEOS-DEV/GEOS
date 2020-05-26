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

/*!
 * @brief Class handling the creation of the well into the GEOSX
 * data stucture from a LAS file import
 */
class LASWellGenerator : public WellGeneratorBase
{
public:
  /**
   * @brief Constructor.
   * @param name name of the object in the data hierarchy.
   * @param parent pointer to the parent group in the data hierarchy.
   */
  LASWellGenerator( const std::string & name,
                    Group * const parent );


  /**
   * @brief default destructor
   */
  virtual ~LASWellGenerator() override;

  /**
   * @return the name of this type in the catalog
   */
  static string CatalogName() { return "LASWell"; }

  /*!
   * @brief Generate the polyline from XML data
   */
  virtual void GeneratePolyLine() override final;

protected:

  /*!
   * @brief This function provides capability to post process input values prior to
   * any other initialization operations.
   */
  void PostProcessInput() override final;

private:
  /*!
   * @brief Generate the polyline if the LAS file contains XYZ information
   * @param[in] lasFile the LASFile reader
   */
  void GeneratePolyLineFromXYZ( LASFile const & lasFile );

  /*!
   * @brief Generate the polyline if the LAS file contains Depth information
   * @details usually for vertical wells
   * @param[in] lasFile the LASFile reader
   */
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
