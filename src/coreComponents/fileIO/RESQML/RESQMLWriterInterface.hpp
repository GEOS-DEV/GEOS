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

#ifndef GEOSX_FILEIO_VTK_RESQMLWRITERINTERFACE_HPP_
#define GEOSX_FILEIO_VTK_RESQMLWRITERINTERFACE_HPP_

#include "fileIO/vtk/VTKPolyDataWriterInterface.hpp"
#include "dataRepository/Group.hpp"

#include "fesapi/common/EpcDocument.h"
#include "fesapi/common/DataObjectRepository.h"

#include "fesapi/resqml2_0_1/ContinuousProperty.h"
#include "fesapi/resqml2_0_1/UnstructuredGridRepresentation.h"

#include <map>

namespace geosx
{

// using namespace dataRepository;

class RESQMLWriterInterface : private vtk::VTKPolyDataWriterInterface
{
public:
  /**
   * @brief Constructor
   * @param[in] outputName folder name in which all the epc and h5 files will be written
   */
  explicit RESQMLWriterInterface( string outputName );


  /**
   * @brief Sets the plot level
   * @details All fields have an associated plot level. If it is <= to \p plotLevel,
   * the field will be output.
   * @param[in] plotLevel the limit plotlevel
   */
  void setPlotLevel( integer plotLevel )
  {
    m_plotLevel = dataRepository::toPlotLevel( plotLevel );
  }

  /**
   * @brief Set the output directory name
   * @param[in] outputDir global output directory location
   * @param[in] outputName name of the VTK output subdirectory and corresponding PVD file
   */
  void setOutputLocation( string outputDir, string outputName )
  {
    m_outputDir = std::move( outputDir );
    m_outputName = std::move( outputName );
    //m_pvd.setFileName( joinPath( m_outputDir, m_outputName ) + ".epc" );
  }

  /**
   * @brief Set the flag to decide whether we only plot the fields specified by fieldNames, or if we also plot fields based on plotLevel
   * @param[in] onlyPlotSpecifiedFieldNames the flag
   */
  void setOnlyPlotSpecifiedFieldNamesFlag( integer const onlyPlotSpecifiedFieldNames )
  {
    m_onlyPlotSpecifiedFieldNames = onlyPlotSpecifiedFieldNames;
  }

  /**
   * @brief Set the names of the fields to output
   * @param[in] fieldNames the fields to output
   */
  void setFieldNames( arrayView1d< string const > const & fieldNames )
  {
    this->m_fieldNames.insert( fieldNames.begin(), fieldNames.end() );
  }

  void setParentRepresentation(const std::tuple<string,string>& parent)
  {
    std::cout << "set parent: " << std::get<0>(parent) << " " << std::get<1>(parent) << std::endl;
    m_parent = m_outputRepository->createPartial<RESQML2_0_1_NS::UnstructuredGridRepresentation>(std::get<0>(parent), std::get<1>(parent));
  }

  /**
   * @param[in] time the time step to be written
   * @param[in] cycle the current cycle of event
   * @param[in] domain the computation domain of this rank
   */
  void write( real64 time, integer cycle, DomainPartition const & domain );

  void generateOutput();


  private:
  

  void writeCellElementRegions( real64 const time,
                                ElementRegionManager const & elemManager,
                                NodeManager const & nodeManager,
                                time_t timestamp ) ;

  void writeElementFields( ElementRegionBase const & region, real64 const time, time_t timestamp ) ;

  void writeElementField( dataRepository::Group const & subRegions,
                   string const & field,
                   COMMON_NS::DataObjectRepository* outputRepository,
                   RESQML2_0_1_NS::UnstructuredGridRepresentation* parent,
                   EML2_NS::TimeSeries * timeSeries,
                   time_t timestamp) ;

  COMMON_NS::DataObjectRepository* m_outputRepository;

  RESQML2_0_1_NS::PropertyKind* m_propertyKind;

  EML2_NS::TimeSeries * m_timeSeries;

  RESQML2_0_1_NS::UnstructuredGridRepresentation* m_parent;

  std::map<string, RESQML2_0_1_NS::ContinuousProperty*> m_property_uuid;

};

}

#endif