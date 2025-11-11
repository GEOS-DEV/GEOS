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

/**
 * @file CohesiveZoneRegion.cpp
 */

#include "CohesiveZoneRegion.hpp"
#include "common/TimingMacros.hpp"
#include "LvArray/src/SparsityPattern.hpp"
#include "common/MpiWrapper.hpp"

namespace geos
{
using namespace dataRepository;

CohesiveZoneRegion::CohesiveZoneRegion( string const & name, Group * const parent ):
  CohesiveZoneRegionBase( name, parent ),
  m_globalID(),
  m_referencePartitioningSurfaceNormal(),
  m_referenceSurfaceNormal(),
  m_referenceArea(),
  m_referencePosition(),
  m_maxNormalDisplacement(),
  m_maxTangentialDisplacement(),
  m_damage()
{  
  registerWrapper( viewKeyStruct::globalIDString(), &m_globalID ).
    setInputFlag( InputFlags::FALSE ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Array of the global indices for cohesive grid nodes" );

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
    reference().resizeDimension< 1 >( 2 );;

  registerWrapper( viewKeyStruct::referencePositionString(), &m_referencePosition ).
    setInputFlag( InputFlags::FALSE ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Reference cohesive grid node positions" ).
    reference().resizeDimension< 1 >( 3 );

  registerWrapper( viewKeyStruct::maxNormalDisplacementString(), &m_maxNormalDisplacement ).
    setInputFlag( InputFlags::FALSE ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Max cohesive normal displacement for each cohesive grid node" );

  registerWrapper( viewKeyStruct::maxTangentialDisplacementString(), &m_maxTangentialDisplacement ).
    setInputFlag( InputFlags::FALSE ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Max cohesive tangential displacement for each cohesive grid node" );

  registerWrapper( viewKeyStruct::damageString(), &m_damage ).
    setInputFlag( InputFlags::FALSE ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Cohesive grid node damages" );
}

CohesiveZoneRegion::~CohesiveZoneRegion()
{}

REGISTER_CATALOG_ENTRY( ObjectManagerBase, CohesiveZoneRegion, string const &, Group * const )

} /* namespace geos */
