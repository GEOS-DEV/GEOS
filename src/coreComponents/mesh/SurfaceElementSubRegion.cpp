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
 * @file SurfaceElementSubRegion.cpp
 */


#include "SurfaceElementSubRegion.hpp"
#include "ElementRegionManager.hpp"

namespace geos
{

using namespace dataRepository;

SurfaceElementSubRegion::SurfaceElementSubRegion( string const & name,
                                                  dataRepository::Group * const parent ):
  ElementSubRegionBase( name, parent ),
  m_2dElemToElems(),
  m_unmappedGlobalIndicesInToNodes(),
  m_toNodesRelation(),
  m_toEdgesRelation(),
  m_elementAperture(),
  m_elementArea(),
  m_normalVector(),
  m_tangentVector1(),
  m_tangentVector2()
{
  registerWrapper( viewKeyStruct::nodeListString(), &m_toNodesRelation ).
    setDescription( "Map to the nodes attached to each SurfaceElement." );

  registerWrapper( viewKeyStruct::edgeListString(), &m_toEdgesRelation ).
    setDescription( "Map to the edges attached to each SurfaceElement." );

  registerWrapper( viewKeyStruct::surfaceElementsToCellRegionsString(), &m_2dElemToElems.m_toElementRegion ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "A map of face element local indices to the cell local indices" );

  registerWrapper( viewKeyStruct::surfaceElementsToCellSubRegionsString(), &m_2dElemToElems.m_toElementSubRegion ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "A map of face element local indices to the cell local indices" );

  registerWrapper( viewKeyStruct::surfaceElementsToCellIndexString(), &m_2dElemToElems.m_toElementIndex ).
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

  registerWrapper( viewKeyStruct::normalVectorString(), &m_normalVector ).
    setDescription( "Unit normal vector to the embedded surface." );

  registerWrapper( viewKeyStruct::t1VectorString(), &m_tangentVector1 ).
    setDescription( "Unit vector in the first tangent direction to the embedded surface." );

  registerWrapper( viewKeyStruct::t2VectorString(), &m_tangentVector2 ).
    setDescription( "Unit vector in the second tangent direction to the embedded surface." );  

  excludeWrappersFromPacking( { viewKeyStruct::nodeListString(),
                                viewKeyStruct::edgeListString(),
                                viewKeyStruct::surfaceElementsToCellRegionsString(),
                                viewKeyStruct::surfaceElementsToCellSubRegionsString(),
                                viewKeyStruct::surfaceElementsToCellIndexString() } );

  // TODO there has to be a cleaner way than this.
  m_2dElemToElems.setElementRegionManager( dynamicCast< ElementRegionManager & >( getParent().getParent().getParent().getParent() ) );

}

SurfaceElementSubRegion::~SurfaceElementSubRegion()
{}

} /* namespace geos */
