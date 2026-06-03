/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include "CohesiveZoneRegionBase.hpp"

#include "common/TimingMacros.hpp"


namespace geos
{
using namespace dataRepository;


CohesiveZoneRegionBase::CohesiveZoneRegionBase( string const & name, Group * const parent ):
  ObjectManagerBase( name, parent ),
  m_enabled( 0 ),
  m_initialized( 0 ),
  m_czVolumeNormalization( 1 ),
  m_computeParticleSurfaceNormalsAndPositions( 0 ),
  m_normalsAndPositionsMethod( 0 ),
  m_czSurfaceDisplacementUpdate( 2 ),
  m_tag( 0 ),
  m_fieldA( 0 ),
  m_fieldB( 1 ),
  m_globalID(),
  m_referencePosition(),
  m_referencePartitioningSurfaceNormal(),
  m_referenceSurfaceNormal(),
  m_referenceArea()
{
//   setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );

  registerWrapper( "initialized", &m_initialized ).
    setInputFlag( InputFlags::FALSE ).
    setApplyDefaultValue( m_initialized ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Flag marking whether cohesive zone has already be initialized.");

  registerWrapper( "enabled", &m_enabled ).
    setInputFlag( InputFlags::FALSE ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Flag marking if cohesive zone is currently enabled" );

  registerWrapper( "constitutiveModel", &m_constitutiveModelName ).
    setInputFlag( InputFlags::REQUIRED).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Name of cohesive region constitutive model");

  registerWrapper( "czVolumeNormalization", &m_czVolumeNormalization ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_czVolumeNormalization ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Flag to perform volume normalization of cohesive zone area");

  registerWrapper( "computeParticleSurfaceNormalsAndPositions", &m_computeParticleSurfaceNormalsAndPositions ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_computeParticleSurfaceNormalsAndPositions ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Flag to compute particle surface normals and positions prior to initialization of cohesive zone");

  registerWrapper( "m_normalsAndPositionsMethod", &m_normalsAndPositionsMethod ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_normalsAndPositionsMethod ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Method for computing particle surface normals and positions");

  registerWrapper( viewKeyStruct::czSurfaceDisplacementUpdateString(), &m_czSurfaceDisplacementUpdate ).
    setInputFlag( InputFlags::FALSE ).
    setApplyDefaultValue( m_czSurfaceDisplacementUpdate ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Cohesive surface displacement update method: 0=TypeA, 1=TypeB, 2=Nodal" );

  registerWrapper( "tag", &m_tag ).
    setInputFlag( InputFlags::REQUIRED ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Tag ID");

  registerWrapper( "fieldA", &m_fieldA ).
    setInputFlag( InputFlags::FALSE ).
    setApplyDefaultValue( m_fieldA ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Index of field A");

  registerWrapper( "fieldB", &m_fieldB ).
    setInputFlag( InputFlags::FALSE ).
    setApplyDefaultValue( m_fieldB ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Index of field B");

 registerWrapper( viewKeyStruct::globalIDString(), &m_globalID ).
    setInputFlag( InputFlags::FALSE ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Array of the global indices for cohesive grid nodes" );

  registerWrapper( viewKeyStruct::referencePositionString(), &m_referencePosition ).
    setInputFlag( InputFlags::FALSE ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Reference cohesive grid node positions" ).
    reference().resizeDimension< 1 >( 3 );

  registerWrapper( viewKeyStruct::referencePartitioningSurfaceNormalString(), &m_referencePartitioningSurfaceNormal ).
    setInputFlag( InputFlags::FALSE ).
    setPlotLevel( PlotLevel::NOPLOT ).  
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Reference partitioning surface normal for cohesive grid nodes" ).
    reference().resizeDimension< 1 >( 3 );

  registerWrapper( viewKeyStruct::referenceSurfaceNormalString(), &m_referenceSurfaceNormal ).
    setInputFlag( InputFlags::FALSE ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Reference cohesive grid node surface normals" ).
    reference().resizeDimension< 1, 2 >( 2, 3 );

  registerWrapper( viewKeyStruct::referenceAreaString(), &m_referenceArea ).
    setInputFlag( InputFlags::FALSE ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Reference cohesive grid node areas" ).
    reference().resizeDimension< 1 >( 2 );
}


CohesiveZoneRegionBase::~CohesiveZoneRegionBase()
{}

}
