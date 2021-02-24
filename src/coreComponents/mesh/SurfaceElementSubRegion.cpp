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
 * @file SurfaceElementSubRegion.cpp
 */


#include "SurfaceElementSubRegion.hpp"

namespace geosx
{

using namespace dataRepository;

SurfaceElementSubRegion::SurfaceElementSubRegion( string const & name,
                                                  dataRepository::Group * const parent ):
  ElementSubRegionBase( name, parent ),
  m_surfaceElementsToCells(),
  m_toNodesRelation(),
  m_toEdgesRelation(),
  m_elementAperture(),
  m_elementArea()
{
  registerWrapper( viewKeyStruct::nodeListString(), &m_toNodesRelation ).
    setDescription( "Map to the nodes attached to each SurfaceElement." );

  registerWrapper( viewKeyStruct::edgeListString(), &m_toEdgesRelation ).
    setDescription( "Map to the edges attached to each SurfaceElement." );


  registerWrapper( viewKeyStruct::surfaceElementsToCellRegionsString(), &m_surfaceElementsToCells.m_toElementRegion ).
    setApplyDefaultValue( -1 ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "A map of face element local indices to the cell local indices" );

  registerWrapper( viewKeyStruct::surfaceElementsToCellSubRegionsString(), &m_surfaceElementsToCells.m_toElementSubRegion ).
    setApplyDefaultValue( -1 ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "A map of face element local indices to the cell local indices" );

  registerWrapper( viewKeyStruct::surfaceElementsToCellIndexString(), &m_surfaceElementsToCells.m_toElementIndex ).
    setApplyDefaultValue( -1 ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "A map of face element local indices to the cell local indices" );

  registerWrapper( viewKeyStruct::elementApertureString(), &m_elementAperture ).
    setApplyDefaultValue( 1.0e-5 ).
    setPlotLevel( dataRepository::PlotLevel::LEVEL_0 ).
    setDescription( "The aperture of each SurfaceElement." );

  registerWrapper( viewKeyStruct::elementAreaString(), &m_elementArea ).
    setApplyDefaultValue( -1.0 ).
    setPlotLevel( dataRepository::PlotLevel::LEVEL_2 ).
    setDescription( "The area of each SurfaceElement." );

  registerWrapper< real64_array >( viewKeyStruct::creationMassString() ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( dataRepository::PlotLevel::LEVEL_1 ).
    setDescription( "The amount of remaining mass that was introduced when the SurfaceElement was created." );

}

SurfaceElementSubRegion::~SurfaceElementSubRegion()
{}


} /* namespace geosx */
