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
 * @file DeformationUpdateMPMEvent.cpp
 */

#include "DeformationUpdateMPMEvent.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geos
{

using namespace dataRepository;

DeformationUpdateMPMEvent::DeformationUpdateMPMEvent( const string & name,
                                                      Group * const parent ):
  MPMEventBase( name, parent ),
  m_relativeDeformation( 0 ),
  m_prescribedFTable( 0 ),
  m_prescribedBoundaryFTable( 0 ),
  m_stressControl(),
  m_fTableInterpType( mpm::InterpolationOption::Cosine ),
  m_stressTableInterpType( mpm::InterpolationOption::Cosine ),
  m_fTable(),
  m_stressTable()
{
  registerWrapper( viewKeyStruct::relativeDeformationString(), &m_relativeDeformation ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_relativeDeformation ).
    setDescription( "Determines where the deformation is absolute or relative to current configuration");

  registerWrapper( viewKeyStruct::prescribedFTableString(), &m_prescribedFTable ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setDescription( "Flag to enable prescribed F table" );

  registerWrapper( viewKeyStruct::prescribedBoundaryFTableString(), &m_prescribedBoundaryFTable ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setDescription( "Flag to enable prescribed boundary F table" );

  registerWrapper( viewKeyStruct::stressControlString(), &m_stressControl ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flags to enable stress control in x, y and z directions" );

  registerWrapper( viewKeyStruct::fTableInterpTypeString(), &m_fTableInterpType ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_fTableInterpType ).
    setDescription( "Interpolation type for F table" );

  registerWrapper( viewKeyStruct::stressTableInterpTypeString(), &m_stressTableInterpType ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_stressTableInterpType ).
    setDescription( "Interpolation type for stress table" );

  registerWrapper( viewKeyStruct::fTableString(), &m_fTable ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "F table" );
  
  registerWrapper( viewKeyStruct::stressTableString(), &m_stressTable ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Stress table" );
}

DeformationUpdateMPMEvent::~DeformationUpdateMPMEvent()
{}

void DeformationUpdateMPMEvent::postInputInitialization()
{
  MPMEventBase::postInputInitialization();

  if( m_stressControl.size() == 0 )
  {
    m_stressControl.resize( 3 );
    LvArray::tensorOps::fill< 3 >( m_stressControl, 0 );
  }
  else
  {
    GEOS_ERROR_IF( m_stressControl.size() != 3, "stressControl must be of length 3. " );
    
    if( m_stressControl[0] > 0 || m_stressControl[1] > 0 || m_stressControl[2] > 0 )
    {
      GEOS_ERROR_IF( m_stressTable.size( 0 ) == 0, "Stress table cannot be empty if stress control is enabled" );
      GEOS_ERROR_IF( m_stressTable.size( 1 ) == 0, "Stress table must have 4 columns" );

      for( int i = 0; i < m_stressTable.size( 0 ); ++i )
      {
        GEOS_ERROR_IF( m_stressTable[i][0] < 0.0, "Stress table times must be positive" );
      }
    }
  }

  if( m_prescribedBoundaryFTable == 1 || m_prescribedFTable == 1 )
  {
    // Reads the FTable directly from the xml
    int numRows = m_fTable.size( 0 );
    GEOS_ERROR_IF( numRows == 0, "Prescribed boundary deformation is enabled but no F table was specified." );
    for( int i = 0; i < numRows; ++i )
    {
      GEOS_ERROR_IF( m_fTable[i].size() != 4, "F table row " << i+1 << " must have 4 elements." );

      if( i == 0 )
      {
        GEOS_ERROR_IF( m_fTable[0][0] > 0.0, "F table times must be positive." );
      }
      else
      {
        GEOS_ERROR_IF( ( m_fTable[i][0] - m_fTable[i-1][0] ) < 0.0, "F table time entries must be monotonically increasing." );
      }

      // TODO: Add check that FTable times are monotonic

      for( int k=0; k<3; ++k )
      { // Stress-control = 1 overrides FTable control, so if we aren't doing stress control in a direction,
        // this checks that the table is well defined:
        if( m_stressControl[k] != 1 )
        {
          GEOS_ERROR_IF( m_fTable[i][k+1] < 0, "F table entries must be positive." );
        }
      }
    }
  }
}

REGISTER_CATALOG_ENTRY( MPMEventBase, DeformationUpdateMPMEvent, string const &, Group * const )

} /* namespace geos */
