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

/**
 * @file VTKOutput.cpp
 */

#include "VTKOutput.hpp"

namespace geosx
{

using namespace dataRepository;

VTKOutput::VTKOutput( string const & name,
                      Group * const parent ):
  OutputBase( name, parent ),
  m_plotFileRoot( name ),
  m_writeFaceMesh(),
  m_plotLevel(),
  m_fieldNames(),
  m_writer( getOutputDirectory() + '/' + m_plotFileRoot )
{
  registerWrapper( viewKeysStruct::plotFileRoot, &m_plotFileRoot ).
    setDefaultValue( m_plotFileRoot ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Name of the root file for this output." );

  registerWrapper( viewKeysStruct::writeFEMFaces, &m_writeFaceMesh ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "" );

  registerWrapper( viewKeysStruct::plotLevel, &m_plotLevel ).
    setApplyDefaultValue( 1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Level detail plot. Only fields with lower of equal plot level will be output." );

  registerWrapper( viewKeysStruct::fieldNames, &m_fieldNames ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Names of the fields to output. If this attribute is specified, GEOSX outputs all (and only) the fields specified by the user, regardless of their plotLevel" );

  registerWrapper( viewKeysStruct::binaryString, &m_writeBinaryData ).
    setApplyDefaultValue( m_writeBinaryData ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Output data format.  Valid options: ``" + EnumStrings< vtk::VTKOutputMode >::concat( "``, ``" ) + "``" );

  registerWrapper( viewKeysStruct::outputRegionTypeString, &m_outputRegionType ).
    setApplyDefaultValue( m_outputRegionType ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Output region types.  Valid options: ``" + EnumStrings< vtk::VTKRegionTypes >::concat( "``, ``" ) + "``" );
}

VTKOutput::~VTKOutput()
{}

void VTKOutput::postProcessInput()
{
  m_writer.setOutputLocation( getOutputDirectory(), m_plotFileRoot );
  m_writer.setFieldNames( m_fieldNames );

  GEOSX_LOG_RANK_0_IF( !m_fieldNames.empty(),
                       GEOSX_FMT(
                         "{} `{}`: found {} fields to plot. These fields will be output regardless of the plotLevel specified by the user. No other field will be output. Remove keyword `{}` from the XML file to output fields based on their plotLevel.",
                         catalogName(), getName(), std::to_string( m_fieldNames.size() ), viewKeysStruct::fieldNames ) );

}

bool VTKOutput::execute( real64 const time_n,
                         real64 const GEOSX_UNUSED_PARAM( dt ),
                         integer const cycleNumber,
                         integer const GEOSX_UNUSED_PARAM( eventCounter ),
                         real64 const GEOSX_UNUSED_PARAM ( eventProgress ),
                         DomainPartition & domain )
{
  m_writer.setOutputMode( m_writeBinaryData );
  m_writer.setOutputRegionType( m_outputRegionType );
  m_writer.setPlotLevel( m_plotLevel );
  m_writer.write( time_n, cycleNumber, domain );

  return false;
}


REGISTER_CATALOG_ENTRY( OutputBase, VTKOutput, string const &, Group * const )
} /* namespace geosx */
