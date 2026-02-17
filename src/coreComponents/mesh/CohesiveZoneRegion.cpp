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
  m_maxNormalDisplacement(),
  m_maxTangentialDisplacement(),
  m_damage(),
  m_temperature()
{
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
  
  registerWrapper( viewKeyStruct::temperatureString(), &m_temperature ).
    setInputFlag( InputFlags::FALSE ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Cohesive grid node temperatures" );
}

CohesiveZoneRegion::~CohesiveZoneRegion()
{}

REGISTER_CATALOG_ENTRY( ObjectManagerBase, CohesiveZoneRegion, string const &, Group * const )

} /* namespace geos */
