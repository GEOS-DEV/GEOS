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

/**
 * @file VTKOutput.cpp
 */

#include "VTKOutput.hpp"
#include "fileIO/vtk/VTKFile.hpp"
#include "managers/DomainPartition.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace cxx_utilities;

VTKOutput::VTKOutput( std::string const & name,
                        Group * const parent ):
  OutputBase( name, parent),
  m_plotFileRoot(),
  m_writeFaceMesh(),
  m_plotLevel(),
  m_vtkFile(name)
{
  registerWrapper(viewKeysStruct::plotFileRoot, &m_plotFileRoot, false )->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("");

  registerWrapper(viewKeysStruct::writeFEMFaces, &m_writeFaceMesh, false )->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("");

  registerWrapper(viewKeysStruct::plotLevel, &m_plotLevel, false )->
    setApplyDefaultValue(1)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("");

  registerWrapper(viewKeysStruct::binaryString, &m_writeBinaryData, false )->
    setApplyDefaultValue(1)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Output the data in binary format");

  m_vtkFile.SetPlotLevel( m_plotLevel );
  m_vtkFile.SetBinaryMode( m_writeBinaryData) ;

}

VTKOutput::~VTKOutput()
{}



void VTKOutput::Execute(real64 const time_n,
                         real64 const GEOSX_UNUSED_PARAM( dt ),
                         integer const GEOSX_UNUSED_PARAM( cycleNumber ),
                         integer const GEOSX_UNUSED_PARAM( eventCounter ),
                         real64 const GEOSX_UNUSED_PARAM( eventProgress ),
                         Group * domain)
{
  DomainPartition* domainPartition = Group::group_cast<DomainPartition*>(domain);
  m_vtkFile.Write( time_n, *domainPartition);
}


REGISTER_CATALOG_ENTRY( OutputBase, VTKOutput, std::string const &, Group * const )
} /* namespace geosx */
