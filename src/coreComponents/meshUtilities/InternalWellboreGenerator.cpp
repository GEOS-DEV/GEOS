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
 * @file InternalWellboreGenerator.cpp
 */

#include "InternalWellboreGenerator.hpp"
#include "managers/DomainPartition.hpp"
#include "mpiCommunications/PartitionBase.hpp"

namespace geosx
{
using namespace dataRepository;

InternalWellboreGenerator::InternalWellboreGenerator( string const & name, Group * const parent ):
  InternalMeshGenerator( name, parent )
{
  registerWrapper( viewKeyStruct::radiusString(), &m_radius ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Wellbore radius" );

  registerWrapper( viewKeyStruct::thetaString(), &m_theta ).
    setApplyDefaultValue( 360.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Tangent angle defining geometry size: 90 for quarter, 180 for half and 360 for full wellbore geometry" );

  registerWrapper( viewKeyStruct::rOutString(), &m_rOut ).
    setApplyDefaultValue( 1.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Farfield distance from wellbore center" );

  registerWrapper( viewKeyStruct::rElemsString(), &m_rElems ).
    setApplyDefaultValue( 10 ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Number of elements in the radial direction" );

  registerWrapper( viewKeyStruct::tElemsString(), &m_tElems ).
    setApplyDefaultValue( 40 ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Number of elements in the tangent direction" );

  registerWrapper( viewKeyStruct::rBiasString(), &m_rBias ).
    setApplyDefaultValue( -0.8 ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Bias of element sizes in the radial direction" );

  registerWrapper( viewKeyStruct::trajectoryString(), &m_trajectory ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Coordinates defining the wellbore trajectory" );
}

void InternalWellboreGenerator::postProcessInput()
{
  m_vertices[0].resize( 2 );
  m_vertices[1].resize( 2 );
  m_vertices[2].resize( 2 );
  m_nElems[0].resize( 1 );
  m_nElems[1].resize( 1 );
  m_nElemBias[0].resize( 1 );

  m_vertices[0][0]  = m_radius;
  m_vertices[0][1]  = m_rOut;
  m_vertices[1][0]  = 0;
  m_vertices[1][1]  = m_theta;
  m_vertices[2][0]  = m_trajectory[0][2];
  m_vertices[2][1]  = m_trajectory[1][2];
  m_nElems[0][0]    = m_rElems;
  m_nElems[1][0]    = m_tElems;
  m_nElemBias[0][0] = m_rBias;

  InternalMeshGenerator::postProcessInput();
}

void InternalWellboreGenerator::generateMesh( DomainPartition & domain )
{
  InternalMeshGenerator::generateMesh( domain );

  Group & meshBodies = domain.getGroup( string( "MeshBodies" ));
  MeshBody & meshBody = meshBodies.registerGroup< MeshBody >( this->getName() );
  MeshLevel & meshLevel0 = meshBody.registerGroup< MeshLevel >( string( "Level0" ));
  NodeManager & nodeManager = meshLevel0.getNodeManager();
  Group & nodeSets = nodeManager.sets();

  // Wellbore nodesets
  SortedArray< localIndex > & rnegNodes = nodeSets.registerWrapper< SortedArray< localIndex > >( string( "rneg" ) ).reference();
  SortedArray< localIndex > & rposNodes = nodeSets.registerWrapper< SortedArray< localIndex > >( string( "rpos" ) ).reference();
  SortedArray< localIndex > & tnegNodes = nodeSets.registerWrapper< SortedArray< localIndex > >( string( "tneg" ) ).reference();
  SortedArray< localIndex > & tposNodes = nodeSets.registerWrapper< SortedArray< localIndex > >( string( "tpos" ) ).reference();

  arrayView2d< real64, nodes::REFERENCE_POSITION_USD > const & X = nodeManager.referencePosition();
  for( int localNodeIndex=0; localNodeIndex<nodeManager.size(); ++localNodeIndex )
  {
    // Get Cartesian coordinates of a reference centered vertical wellbore
    real64 & xCoord = X( localNodeIndex, 0 );
    real64 & yCoord = X( localNodeIndex, 1 );
    real64 const & zCoord = X( localNodeIndex, 2 );

    // Compute cylindrical coordinates of a reference centered vertical wellbore
    real64 rCoord = sqrt( xCoord * xCoord + yCoord * yCoord );
    real64 tCoord;

    if( rCoord < 1e-10 )
    {
      tCoord = 0.0;
    }
    else if( yCoord >= 0.0 )
    {
      tCoord = acos( xCoord/rCoord );
    }
    else
    {
      tCoord = 2 * M_PI - acos( xCoord/rCoord );
    }

    tCoord *= 180.0 / M_PI;

    // Wellbore nodesets
    if( isEqual( rCoord, m_min[0], m_positionTolerance ) )
    {
      rnegNodes.insert( localNodeIndex );
    }

    if( isEqual( rCoord, m_max[0], m_positionTolerance ) )
    {
      rposNodes.insert( localNodeIndex );
    }

    if( isEqual( tCoord, m_min[1], m_positionTolerance ) )
    {
      tnegNodes.insert( localNodeIndex );
    }
    if( isEqual( tCoord, m_max[1], m_positionTolerance ) )
    {
      tposNodes.insert( localNodeIndex );
    }

    // Radial distance of the outer square boundary of a reference centered vertical wellbore
    real64 meshTheta = tCoord * M_PI / 180.0;
    int meshAxis = static_cast< int >(round( meshTheta * 2.0 / M_PI ));
    real64 meshPhi = fabs( meshTheta - meshAxis * M_PI / 2.0 );
    real64 meshRout = m_max[0] / cos( meshPhi );

    // Wellbore trajectory
    real64 xTopCenter = m_trajectory[0][0];
    real64 yTopCenter = m_trajectory[0][1];
    real64 zTop = m_min[2];

    real64 xBottomCenter = m_trajectory[1][0];
    real64 yBottomCenter = m_trajectory[1][1];
    real64 zBottom = m_max[2];

    real64 dx = xBottomCenter - xTopCenter;
    real64 dy = yBottomCenter - yTopCenter;
    real64 dz = zBottom - zTop;
    real64 dr = sqrt( dx*dx + dy*dy );
    real64 dl = sqrt( dr*dr + dz*dz );

    // Azimuth from x-axis
    real64 theta0;

    if( dr < 1e-10 )
    {
      theta0 = 0.0;
    }
    else if( dy>=0.0 )
    {
      theta0 = acos( dx/dr );
    }
    else
    {
      theta0 = 2.0 * M_PI - acos( dx/dr );
    }

    // The horizontal section of an inclined wellbore is an ellipse
    // The principle directions of this ellipse is defined by dTheta = 0, and PI/2
    real64 dTheta = meshTheta - theta0;
    real64 tanDTheta = tan( dTheta );

    // Transform radial coordinate regarding the elliptical shape of the wellbore horizontal section
    // This transformation ensures that the ourter square boundary is unchanged
    real64 transformCoeff = sqrt ( ( 1.0 + tanDTheta * tanDTheta )/( dz*dz/dl/dl + tanDTheta * tanDTheta ) );
    real64 rCoordTransform = rCoord * ( ( meshRout - rCoord ) / ( meshRout - m_min[0] ) * ( transformCoeff - 1.0 ) + 1.0 );

    // Compute transformed cartesian coordinates
    xCoord = rCoordTransform * cos( meshTheta );
    yCoord = rCoordTransform * sin( meshTheta );

    // Moving the coordinate in the horizontal plane with respect to the center of the wellbore section
    // This transformation ensures that the ourter square boundary is unchanged
    real64 zRatio = ( zCoord - zTop ) / ( zBottom -zTop );
    real64 xCenter = xTopCenter + (xBottomCenter - xTopCenter) * zRatio;
    real64 yCenter = yTopCenter + (yBottomCenter - yTopCenter) * zRatio;

    xCoord += xCenter * ( meshRout - rCoord ) / ( meshRout - m_min[0] * transformCoeff );
    yCoord += yCenter * ( meshRout - rCoord ) / ( meshRout - m_min[0] * transformCoeff );
  }
}

REGISTER_CATALOG_ENTRY( MeshGeneratorBase, InternalWellboreGenerator, string const &, Group * const )
} /* namespace geosx */
