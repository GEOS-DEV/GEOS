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

  registerWrapper( viewKeyStruct::volumicFieldsToImportString(), &m_volumicFieldsToImport ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Volumic fields to be imported from the external mesh file" );

  registerWrapper( viewKeyStruct::volumicFieldsInGEOSXString(), &m_volumicFieldsInGEOSX ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Names of the volumic fields in GEOSX to import into" );

  registerWrapper( viewKeyStruct::surfacicFieldsToImportString(), &m_surfacicFieldsToImport ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Surfacic fields to be imported from the external mesh file" );

  registerWrapper( viewKeyStruct::surfacicFieldsInGEOSXString(), &m_surfacicFieldsInGEOSX ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Names of the surfacic fields in GEOSX to import into" );
}

void ExternalMeshGeneratorBase::postProcessInput()
{
  auto const checkSizes = [meshName = getName()]( arrayView1d< string const > from,
                                                  arrayView1d< string const > to,
                                                  string const & fromKey,
                                                  string const & toKey )
  {
    GEOSX_THROW_IF_NE_MSG( from.size(), to.size(),
                           "Mesh '" << meshName << "': attributes '" << fromKey << "' and '" << toKey << "' must contain the same number of values.",
                           InputError );
  };
  checkSizes( m_volumicFieldsToImport, m_volumicFieldsInGEOSX, viewKeyStruct::volumicFieldsToImportString(), viewKeyStruct::volumicFieldsInGEOSXString() );
  checkSizes( m_surfacicFieldsToImport, m_surfacicFieldsInGEOSX, viewKeyStruct::surfacicFieldsToImportString(), viewKeyStruct::surfacicFieldsInGEOSXString() );

  auto const checkDuplicates = [meshName = getName()]( arrayView1d< string const > v )
  {
    std::set< string > const tmp{ v.begin(), v.end() };
    bool const hasDuplicates = tmp.size() != LvArray::integerConversion< std::size_t >( v.size() );

    GEOSX_THROW_IF( hasDuplicates,
                    "Mesh '" << meshName << "': '" << stringutilities::join( v, ", " ) << "' already present in list of fields to import.",
                    InputError );
  };
  checkDuplicates( m_volumicFieldsInGEOSX );
  checkDuplicates( m_surfacicFieldsInGEOSX );

  // Building the fields mapping from the two separated input/output vectors.
  auto const buildMapping = [&]( arrayView1d< string const > from,
                                 arrayView1d< string const > to ) -> std::map< string, string >
  {
    std::map< string, string > mapping;
    for( int i = 0; i < from.size(); i++ )
    {
      mapping[from[i]] = to[i];
    }
    return mapping;
  };

  MeshGeneratorBase::m_volumicFields = buildMapping( m_volumicFieldsToImport.toViewConst(), m_volumicFieldsInGEOSX.toViewConst() );
  MeshGeneratorBase::m_surfacicFields = buildMapping( m_surfacicFieldsToImport.toViewConst(), m_surfacicFieldsInGEOSX.toViewConst() );
}

} // namespace geosx
