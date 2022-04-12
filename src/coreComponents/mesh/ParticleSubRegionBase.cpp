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
 * @file ParticleSubRegionBase.cpp
 */

#include "ParticleSubRegionBase.hpp"

namespace geosx
{
using namespace dataRepository;

ParticleSubRegionBase::ParticleSubRegionBase( string const & name, Group * const parent ): // @suppress("Class members should be properly initialized")
  ObjectManagerBase( name, parent ),
  m_constitutiveModels( groupKeyStruct::constitutiveModelsString(), this ),
  m_ghostRank(),
  m_particleID(),
  m_particleCenter(),
  m_particleVelocity(),
  m_particleVolume(),
  m_particleVolume0(),
  m_particleMass(),
  m_particleDeformationGradient(),
  m_particleRVectors(),
  m_particleRVectors0()
{
  registerGroup( groupKeyStruct::constitutiveModelsString(), &m_constitutiveModels ).
    setSizedFromParent( 1 );

  registerWrapper( viewKeyStruct::ghostRankString(), &m_ghostRank ).
    setPlotLevel( PlotLevel::LEVEL_1 );

  registerWrapper( viewKeyStruct::particleIDString(), &m_particleID ).
    setPlotLevel( PlotLevel::LEVEL_1 );

  registerWrapper( viewKeyStruct::particleCenterString(), &m_particleCenter ).
    setPlotLevel( PlotLevel::LEVEL_1 ).
    reference().resizeDimension< 1 >( 3 );

  registerWrapper( viewKeyStruct::particleVelocityString(), &m_particleVelocity ).
    setPlotLevel( PlotLevel::LEVEL_1 ).
    reference().resizeDimension< 1 >( 3 );

  registerWrapper( viewKeyStruct::particleVolumeString(), &m_particleVolume ).
    setPlotLevel( PlotLevel::LEVEL_1 );

  registerWrapper( viewKeyStruct::particleVolume0String(), &m_particleVolume0 ).
    setPlotLevel( PlotLevel::LEVEL_1 );

  registerWrapper( viewKeyStruct::particleMassString(), &m_particleMass ).
    setPlotLevel( PlotLevel::LEVEL_1 );

  // The only things that should be registered here are those that are read in from the input files. So e.g. particle mass shouldn't be here since it's not specified in any input file, whereas particle volume is.
  // A solver(?) should then on the first cycle register particle mass and calculate it from volume and density. Same idea for other similarly post-processed and/or solver-specific fields.

  // idk what I'm doing
//  registerWrapper( viewKeyStruct::particleDeformationGradientString(), &m_particleDeformationGradient ).
//    setPlotLevel( PlotLevel::LEVEL_1 ).
//    reference().resizeDimension< 1 >( 3 ).resizeDimension< 2 >( 3 );
}

ParticleSubRegionBase::~ParticleSubRegionBase()
{}

} /* namespace geosx */
