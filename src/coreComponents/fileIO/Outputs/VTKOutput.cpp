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


#if defined(GEOSX_USE_PYGEOSX)
#include "fileIO/python/PyVTKOutputType.hpp"
#endif

namespace geosx
{

using namespace dataRepository;

VTKOutput::VTKOutput( string const & name,
                      Group * const parent ):
  OutputBase( name, parent ),
  m_plotFileRoot( name ),
  m_writeFaceMesh(),
  m_plotLevel(),
  m_onlyPlotSpecifiedFieldNames(),
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

  registerWrapper( viewKeysStruct::onlyPlotSpecifiedFieldNames, &m_onlyPlotSpecifiedFieldNames ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription(
    "If this flag is equal to 1, then we only plot the fields listed in `fieldNames`. Otherwise, we plot all the fields with the required `plotLevel`, plus the fields listed in `fieldNames`" );

  registerWrapper( viewKeysStruct::fieldNames, &m_fieldNames ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Names of the fields to output. If this attribute is specified, GEOSX outputs all the fields specified by the user, regardless of their `plotLevel`" );

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
  m_writer.setFieldNames( m_fieldNames.toViewConst() );
  m_writer.setOnlyPlotSpecifiedFieldNamesFlag( m_onlyPlotSpecifiedFieldNames );

  string const fieldNamesString = viewKeysStruct::fieldNames;
  string const onlyPlotSpecifiedFieldNamesString = viewKeysStruct::onlyPlotSpecifiedFieldNames;

  GEOSX_THROW_IF( ( m_onlyPlotSpecifiedFieldNames != 0 ) && m_fieldNames.empty(),
                  GEOSX_FMT( "{} `{}`: the flag `{}` is different from zero, but `{}` is empty, which is inconsistent",
                             catalogName(), getName(), onlyPlotSpecifiedFieldNamesString, fieldNamesString ),
                  InputError );

  GEOSX_LOG_RANK_0_IF( !m_fieldNames.empty() && ( m_onlyPlotSpecifiedFieldNames != 0 ),
                       GEOSX_FMT(
                         "{} `{}`: found {} fields to plot in `{}`. These fields will be output regardless of the `plotLevel` specified by the user. No other field will be output.",
                         catalogName(), getName(), std::to_string( m_fieldNames.size() ), fieldNamesString ) );

  GEOSX_LOG_RANK_0_IF( !m_fieldNames.empty() && ( m_onlyPlotSpecifiedFieldNames == 0 ),
                       GEOSX_FMT(
                         "{} `{}`: found {} fields to plot in `{}`, in addition to all fields with `plotLevel` smaller or equal to {}.",
                         catalogName(), getName(), std::to_string( m_fieldNames.size() ), fieldNamesString, m_plotLevel ) );
}


void VTKOutput::setPlotFileRoot( string const & root )
{
  m_plotFileRoot = root;
}


void VTKOutput::reinit()
{
  m_writer.clearData();

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

#if defined(GEOSX_USE_PYGEOSX)
PyTypeObject * VTKOutput::getPythonType() const
{
  return python::getPyVTKOutputType();
}
#endif

REGISTER_CATALOG_ENTRY( OutputBase, VTKOutput, string const &, Group * const )
} /* namespace geosx */
