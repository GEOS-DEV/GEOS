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

  registerWrapper( viewKeysStruct::binaryString, &m_writeBinaryData ).
    setApplyDefaultValue( m_writeBinaryData ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Output the data in binary format.  Valid options: ``" + EnumStrings< vtk::VTKOutputMode >::concat( "``, ``" ) + "``");

  registerWrapper( viewKeysStruct::outputRegionTypeString, &m_outputRegionType ).
    setApplyDefaultValue( m_outputRegionType ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Output region types.  Valid options: ``" + EnumStrings< vtk::VTKRegionTypes >::concat( "``, ``" ) + "``");
}

VTKOutput::~VTKOutput()
{}

void VTKOutput::postProcessInput()
{
  m_writer.setOutputLocation( getOutputDirectory(), m_plotFileRoot );
}

bool VTKOutput::execute( real64 const time_n,
                         real64 const GEOSX_UNUSED_PARAM( dt ),
                         integer const cycleNumber,
                         integer const GEOSX_UNUSED_PARAM( eventCounter ),
                         real64 const GEOSX_UNUSED_PARAM ( eventProgress ),
                         DomainPartition & domain )
{
  m_writer.setOutputMode(m_writeBinaryData);
  m_writer.setOutputRegionType(m_outputRegionType);
  m_writer.setPlotLevel( m_plotLevel );
  m_writer.write( time_n, cycleNumber, domain );

  return false;
}


REGISTER_CATALOG_ENTRY( OutputBase, VTKOutput, string const &, Group * const )
} /* namespace geosx */
