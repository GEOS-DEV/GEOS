/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file ExternalMeshGeneratorBase.cpp
 */

#include "ExternalMeshGeneratorBase.hpp"

namespace geosx
{

using namespace dataRepository;

ExternalMeshGeneratorBase::ExternalMeshGeneratorBase( string const & name,
                                                      dataRepository::Group * const parent )
  : MeshGeneratorBase( name, parent )
{
  enableLogLevelInput();

  registerWrapper( viewKeyStruct::filePathString(), &m_filePath ).
    setInputFlag( InputFlags::REQUIRED ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Path to the mesh file" );

  registerWrapper( viewKeyStruct::translateString(), &m_translate ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( { 0.0, 0.0, 0.0 } ).
    setDescription( "Translate the coordinates of the vertices by a given vector (prior to scaling)" );

  registerWrapper( viewKeyStruct::scaleString(), &m_scale ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( { 1.0, 1.0, 1.0 } ).
    setDescription( "Scale the coordinates of the vertices by given scale factors (after translation)" );

  registerWrapper( viewKeyStruct::fieldsToImportString(), &m_fieldsToImport ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Fields to be imported from the external mesh file" );

  registerWrapper( viewKeyStruct::fieldNamesInGEOSXString(), &m_fieldNamesInGEOSX ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Names of fields in GEOSX to import into" );
}

void ExternalMeshGeneratorBase::postProcessInput()
{
  GEOSX_THROW_IF_NE_MSG( m_fieldsToImport.size(), m_fieldNamesInGEOSX.size(),
                         GEOSX_FMT( "Mesh '{}': attributes '{}' and '{}' must contain the same number of values",
                                    getName(),
                                    viewKeyStruct::fieldsToImportString(),
                                    viewKeyStruct::fieldNamesInGEOSXString() ),
                         InputError );
}

} // namespace geosx
