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

#ifndef GEOSX_FILEIO_RESQML_RESQMLWRITERINTERFACE_HPP_
#define GEOSX_FILEIO_RESQML_RESQMLWRITERINTERFACE_HPP_

#include "fileIO/vtk/VTKPolyDataWriterInterface.hpp"
#include "dataRepository/Group.hpp"

#include "fesapi/common/EpcDocument.h"
#include "fesapi/common/DataObjectRepository.h"

#include "fesapi/resqml2_0_1/ContinuousProperty.h"
#include "fesapi/resqml2_0_1/UnstructuredGridRepresentation.h"
#include "fesapi/resqml2/SubRepresentation.h"

#include <vtkDataArray.h>

#include <map>
#include <unordered_set>
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

  /**
   * @brief Set the parameters to create the Parent Representation object
   * @param[in] parent A tuple of strings (UUID, Name)
   */
  void setParentRepresentation( const std::tuple< string, string > & parent )
  {
    m_parent = m_outputRepository->createPartial< RESQML2_0_1_NS::UnstructuredGridRepresentation >( std::get< 0 >( parent ), std::get< 1 >( parent ));
  }

  /**
   * @brief Main method of this class. Write all the files for one time step.
   * @param[in] time the time step to be written
   * @param[in] cycle the current cycle of event
   * @param[in] domain the computation domain of this rank
   */
  void write( real64 time, integer cycle, DomainPartition const & domain );

  /**
   * @brief Generates the output .epc and .hdf5 files from the data in the output repository
   */
  void generateOutput() const;

  /**
   * @brief Initialize HDFProxy to write numerical data
   */
  void initializeOutput();

  /**
   * @brief Generate a subrepresentation to map a property to the initial grid
   * @param[in] domain the computation domain of this rank
   */
  void generateSubRepresentations( DomainPartition const & domain );

private:

  /**
   * @brief Generate a subRepresentation for a field
   * @param[in] elemManager ElementRegion being written
   * @param[in] field field associated to the elements
   */
  void generateSubRepresentation( ElementRegionManager const & elemManager,
                                  string const & field );

private:

  /// Output repository of RESQML data objects
  COMMON_NS::DataObjectRepository * m_outputRepository;

  /// Property kind of output properties
  RESQML2_0_1_NS::PropertyKind * m_propertyKind;

  /// Time Series object to link to the output properties
  EML2_NS::TimeSeries * m_timeSeries;

  /// Parent representation of output properties
  RESQML2_0_1_NS::UnstructuredGridRepresentation * m_parent;

  /// Index the properties to reuse them accross the multiple regions subgroups
  std::map< string, RESQML2_0_1_NS::ContinuousProperty * > m_property_uuid;

  /// Permutation vector for each field
  std::map< string, array1d< int > > m_countPerProp;

  /// Regular fields to output
  std::unordered_set< string > m_regularFields;

  /// Keep track of subrepresentation throw the time: field name -> Subrep
  std::map< string, RESQML2_NS::SubRepresentation * > m_subrepresentations;

};

}

#endif
