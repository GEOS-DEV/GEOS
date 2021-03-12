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
  InternalMeshGenerator( name, parent ),
  m_trajectory(),
  m_cartesianOuterBoundary(),
  m_isFullAnnulus(false),
  m_autoSpaceRadialElems(),
  m_radialCoords( this->m_setCoords[0] )
{

  getWrapper< array1d< real64 > >( viewKeyStruct::xCoordsString() ).
    setInputFlag( InputFlags::FALSE );

  getWrapper< array1d< real64 > >( viewKeyStruct::yCoordsString() ).
    setInputFlag( InputFlags::FALSE );

  getWrapper< array1d< integer > >( viewKeyStruct::xElemsString() ).
    setInputFlag( InputFlags::FALSE );

  getWrapper< array1d< integer > >( viewKeyStruct::yElemsString() ).
    setInputFlag( InputFlags::FALSE );


  registerWrapper( viewKeyStruct::radiusString(), &(m_vertices[0]) ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Wellbore radius" );

  registerWrapper( viewKeyStruct::thetaString(), &(m_vertices[1]) ).
    setApplyDefaultValue( 360.0 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Tangent angle defining geometry size: 90 for quarter, 180 for half and 360 for full wellbore geometry" );

  registerWrapper( viewKeyStruct::rElemsString(), &(m_nElems[0]) ).
    setApplyDefaultValue( 10 ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Number of elements in the radial direction" );

  registerWrapper( viewKeyStruct::tElemsString(), &(m_nElems[1]) ).
    setApplyDefaultValue( 40 ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Number of elements in the tangent direction" );

  registerWrapper( viewKeyStruct::rBiasString(), &(m_nElemBias[0]) ).
    setApplyDefaultValue( -0.8 ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Bias of element sizes in the radial direction" );

  registerWrapper( viewKeyStruct::trajectoryString(), &m_trajectory ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Coordinates defining the wellbore trajectory" );


  registerWrapper( viewKeyStruct::cartesianOuterBoundaryString(), &m_cartesianOuterBoundary ).
    setApplyDefaultValue( 1000000 ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Enforce a Cartesian aligned outer boundary on the outer "
                    "block starting with the radial block specified in this value" );

  registerWrapper( viewKeyStruct::autoSpaceRadialElemsString(), &m_autoSpaceRadialElems ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Automatically set number and spacing of elements in the radial direction. "
                    "This overrides the values of nr!!" );

}

void InternalWellboreGenerator::postProcessInput()
{
  m_vertices[0].resize( 2 );
  m_vertices[1].resize( 2 );
  m_vertices[2].resize( 2 );

  m_nElems[0].resize( 1 );
  m_nElems[1].resize( 1 );
  m_nElems[2].resize( 1 );  //TODO extend to a trajectory defined by a polyline
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



  arrayView1d<real64 const> const theta = m_vertices[1];
  real64 const dTheta = theta.back() - theta[0];

  // enable full annulus corrections if the mesh is 360 degrees
  if( dTheta >= 360.0 - 1e-99 )
  {
    m_isFullAnnulus = true;
  }

  // automatically set radial coordinates
  if( m_autoSpaceRadialElems.size()>0 )
  {
    localIndex const numRadialBlocks = m_vertices[0].size()-1;

    for( localIndex iBlock=0; iBlock<numRadialBlocks; ++iBlock )
    {
      real64 const rInner = m_vertices[0][iBlock];
      real64 const rOuter = m_vertices[0][iBlock+1];

      if( iBlock==0 )
      {
        m_radialCoords.emplace_back( rInner );
      }

      if( m_autoSpaceRadialElems[iBlock] == 1 )
      {

        localIndex const startingIndex = m_radialCoords.size()-1;
        real64 const tElemSizeInner = rInner * ( 2 * M_PI * dTheta / 360 ) / m_nElems[1][0];
        real64 const tElemSizeOuter = rOuter * ( 2 * M_PI * dTheta / 360 ) / m_nElems[1][0];


        localIndex actualNumberOfRadialElements = 0;
        real64 scalingFactor = 0;

        printf( "  i   r_i      t_i     rip1_0    tip1_0    r_ip1\n" );
        for( localIndex i=0; i<m_nElems[0][iBlock]*10; ++i )
        {
          real64 const t_i =  m_radialCoords.back() * ( 2 * M_PI * dTheta / 360 ) / m_nElems[1][0];
          real64 const r_ip1_0 = ( m_radialCoords.back() +  t_i ) ;
          real64 const tElemSize_ip1_0 =  r_ip1_0 * ( 2 * M_PI * dTheta / 360 ) / m_nElems[1][0];


          constexpr real64 c = 0.5;
          real64 const r_ip1 = m_radialCoords.back() + ( 1.0 - c ) * t_i + c * tElemSize_ip1_0;

          printf( "%5ld %8.4f %8.4f %8.4f %8.4f %8.4f \n", i, m_radialCoords.back(), t_i, r_ip1_0, tElemSize_ip1_0, r_ip1 );


          // if the radius of the next layer is bigger than rOuter, we figure
          // out where to cut off the layer.
          if( r_ip1 > rOuter )
          {
            real64 const overshoot = r_ip1 - rOuter;
            real64 const undershoot = rOuter - m_radialCoords.back();

            if( overshoot < undershoot )
            {
              // use and append the overshot value
              m_radialCoords.emplace_back(r_ip1);
              actualNumberOfRadialElements = i+1;
              scalingFactor =  ( rOuter - rInner ) / ( r_ip1 - rInner ) ;
            }
            else
            {
              // use the undershot value...throw away the overshot value.
              actualNumberOfRadialElements = i;
              scalingFactor = ( rOuter - rInner ) / ( m_radialCoords.back() - rInner );
            }
            break;
          }
          else
          {
            m_radialCoords.emplace_back(r_ip1);
          }
        }
        std::cout<<actualNumberOfRadialElements<<", "<<scalingFactor<<std::endl;
        std::cout<<m_radialCoords.size()<<std::endl;

        m_nElems[0][iBlock] = actualNumberOfRadialElements;
        for( localIndex i=startingIndex; i<m_radialCoords.size(); ++i )
        {
          m_radialCoords[i] = ( m_radialCoords[i] - rInner ) * scalingFactor + rInner;
        }

      }
      else
      {
        real64 min = m_vertices[0][iBlock];
        real64 max = m_vertices[0][iBlock+1];
        real64 const h = (max-min) / m_nElems[0][iBlock];
        localIndex const numNodes = m_nElems[0][iBlock] + 1;
        for( localIndex i=1; i<numNodes; ++i )
        {
          m_radialCoords.emplace_back( min +  i * h );
        }
      }
    }
    std::cout<<m_radialCoords<<std::endl;
  }


  InternalMeshGenerator::postProcessInput();

}

void InternalWellboreGenerator::generateMesh( DomainPartition & domain )
{
  InternalMeshGenerator::generateMesh( domain );

//  Group & meshBodies = domain.getGroup( string( "MeshBodies" ));
//  MeshBody & meshBody = meshBodies.registerGroup< MeshBody >( this->getName() );
//  MeshLevel & meshLevel0 = meshBody.registerGroup< MeshLevel >( string( "Level0" ));
//  NodeManager & nodeManager = meshLevel0.getNodeManager();
//  Group & nodeSets = nodeManager.sets();
//
//  // Wellbore nodesets
//  // rneg, rpos, tneg and tpos are the named used by the end-used in the input files. Consider modifying them with care.
//  SortedArray< localIndex > & rnegNodes = nodeSets.registerWrapper< SortedArray< localIndex > >( string( "rneg" ) ).reference();
//  SortedArray< localIndex > & rposNodes = nodeSets.registerWrapper< SortedArray< localIndex > >( string( "rpos" ) ).reference();
//  SortedArray< localIndex > & tnegNodes = nodeSets.registerWrapper< SortedArray< localIndex > >( string( "tneg" ) ).reference();
//  SortedArray< localIndex > & tposNodes = nodeSets.registerWrapper< SortedArray< localIndex > >( string( "tpos" ) ).reference();
//
//  arrayView2d< real64, nodes::REFERENCE_POSITION_USD > const & X = nodeManager.referencePosition();
//  for( int localNodeIndex=0; localNodeIndex<nodeManager.size(); ++localNodeIndex )
//  {
//    // Get Cartesian coordinates of a reference centered vertical wellbore
//    real64 & xCoord = X( localNodeIndex, 0 );
//    real64 & yCoord = X( localNodeIndex, 1 );
//    real64 const & zCoord = X( localNodeIndex, 2 );
//
//    // Compute cylindrical coordinates of a reference centered vertical wellbore
//    real64 rCoord = sqrt( xCoord * xCoord + yCoord * yCoord );
//    real64 tCoord;
//
//    if( rCoord < m_coordinatePrecision )
//    {
//      tCoord = 0.0;
//    }
//    else if( yCoord >= 0.0 )
//    {
//      tCoord = acos( xCoord/rCoord );
//    }
//    else
//    {
//      tCoord = 2.0 * M_PI - acos( xCoord/rCoord );
//    }
//
//    tCoord *= 180.0 / M_PI;
//
//    // Wellbore nodesets
//    if( isEqual( rCoord, m_min[0], m_coordinatePrecision ) )
//    {
//      rnegNodes.insert( localNodeIndex );
//    }
//
//    if( isEqual( rCoord, m_max[0], m_coordinatePrecision ) )
//    {
//      rposNodes.insert( localNodeIndex );
//    }
//
//    if( isEqual( tCoord, m_min[1], m_coordinatePrecision ) )
//    {
//      tnegNodes.insert( localNodeIndex );
//    }
//    if( isEqual( tCoord, m_max[1], m_coordinatePrecision ) )
//    {
//      tposNodes.insert( localNodeIndex );
//    }
//
//    // Radial distance of the outer square boundary of a reference centered vertical wellbore
//    real64 meshTheta = tCoord * M_PI / 180.0;
//    int meshAxis = static_cast< int >(round( meshTheta * 2.0 / M_PI ));
//    real64 meshPhi = fabs( meshTheta - meshAxis * M_PI / 2.0 );
//    real64 meshRout = m_max[0] / cos( meshPhi );
//
//    // Wellbore trajectory
//    real64 xTopCenter = m_trajectory[0][0];
//    real64 yTopCenter = m_trajectory[0][1];
//    real64 zTop = m_min[2];
//
//    real64 xBottomCenter = m_trajectory[1][0];
//    real64 yBottomCenter = m_trajectory[1][1];
//    real64 zBottom = m_max[2];
//
//    real64 dx = xBottomCenter - xTopCenter;
//    real64 dy = yBottomCenter - yTopCenter;
//    real64 dz = zBottom - zTop;
//    real64 dr = sqrt( dx*dx + dy*dy );
//    real64 dl = sqrt( dr*dr + dz*dz );
//
//    // Azimuth of the wellbore from x-axis
//    real64 theta0;
//
//    if( dr < m_coordinatePrecision )
//    {
//      theta0 = 0.0;
//    }
//    else if( dy>=0.0 )
//    {
//      theta0 = acos( dx/dr );
//    }
//    else
//    {
//      theta0 = 2.0 * M_PI - acos( dx/dr );
//    }
//
//    // The horizontal section of an inclined wellbore is an ellipse
//    // The principle directions of this ellipse are defined by dTheta = 0, and PI/2
//    real64 dTheta = meshTheta - theta0;
//    real64 tanDTheta = tan( dTheta );
//
//    // Transform radial coordinate regarding the elliptical shape of the wellbore section in the horizontal plane
//    // This transformation ensures that the outer square boundary is unchanged
//    // TODO create a function in ComputationalGeometry class for this pure geometrical transformation
//    real64 transformCoeff = sqrt ( ( 1.0 + tanDTheta * tanDTheta )/( dz*dz/dl/dl + tanDTheta * tanDTheta ) );
//    real64 rCoordTransform = rCoord * ( ( meshRout - rCoord ) / ( meshRout - m_min[0] ) * ( transformCoeff - 1.0 ) + 1.0 );
//
//    // Compute transformed cartesian coordinates
//    xCoord = rCoordTransform * cos( meshTheta );
//    yCoord = rCoordTransform * sin( meshTheta );
//
//    // Moving the coordinate in the horizontal plane with respect to the center of the wellbore section
//    // This transformation ensures that the ourter square boundary is unchanged
//    real64 zRatio = ( zCoord - zTop ) / ( zBottom -zTop );
//    real64 xCenter = xTopCenter + (xBottomCenter - xTopCenter) * zRatio;
//    real64 yCenter = yTopCenter + (yBottomCenter - yTopCenter) * zRatio;
//
//    xCoord += xCenter * ( meshRout - rCoord ) / ( meshRout - m_min[0] * transformCoeff );
//    yCoord += yCenter * ( meshRout - rCoord ) / ( meshRout - m_min[0] * transformCoeff );
//  }
}

void InternalWellboreGenerator::reduceNumNodesForPeriodicBoundary( integer (&numNodesInDir)[3] )
{
  GEOSX_UNUSED_VAR(numNodesInDir);
  if( m_isFullAnnulus )
  {
    numNodesInDir[1] -= 1;
  }
}

void InternalWellboreGenerator::setNodeGlobalIndicesOnPeriodicBoundary( int (&index)[3],
                                             real64 (&minExtent)[3],
                                             real64 (&maxExtent)[3],
                                             arraySlice1d<real64 const> const & X,
                                             real64 const tol )
{
  GEOSX_UNUSED_VAR( minExtent );
  if( m_isFullAnnulus )
  {
    if( isEqual( X[ 1 ], maxExtent[1], tol ) )
    {
      index[1] = 0;
    }
  }
}

void InternalWellboreGenerator::setConnectivityForPeriodicBoundaries( integer const i,
                                                                      integer const j,
                                                                      integer const k,
                                                                      integer const iBlock,
                                                                      integer const jBlock,
                                                                      integer const kBlock,
                                                                      int (&globalIJK)[3],
                                                                      int const (&numElemsInDirForBlock)[3],
                                                                      integer const (&numNodesInDir)[3],
                                                                      int const (&firstElemIndexInPartition)[3],
                                                                      localIndex (&nodeOfBox)[8] )
{
  GEOSX_UNUSED_VAR( i, k, iBlock, kBlock );
  if( m_isFullAnnulus )
  {
    setConnectivityForPeriodicBoundary( 1,
                                        j,
                                        jBlock,
                                        globalIJK,
                                        numElemsInDirForBlock,
                                        numNodesInDir,
                                        firstElemIndexInPartition,
                                        nodeOfBox );
  }
}

void InternalWellboreGenerator::coordinateTransformation( NodeManager & nodeManager )
{
  arrayView2d< real64, nodes::REFERENCE_POSITION_USD > const & X = nodeManager.referencePosition();


  real64 const cartesianMappingInnerRadius = m_cartesianOuterBoundary<m_vertices[0].size() ?
                                             m_vertices[0][m_cartesianOuterBoundary] :
                                             1e99;

  // Map to radial mesh
  for( localIndex iN = 0; iN<nodeManager.size(); ++iN )
  {
    real64 meshTheta = X[iN][1] * M_PI / 180.0;
    int meshAxis = static_cast< int >(round( meshTheta * 2.0 / M_PI ));
    real64 meshPhi = fabs( meshTheta - meshAxis * M_PI / 2.0 );
    real64 meshRout = m_max[0] / cos( meshPhi );
    real64 meshRact;

    if( X[iN][0] > cartesianMappingInnerRadius )
    {
      meshRact = ( ( meshRout - m_min[0] ) / ( m_max[0] - m_min[0] ) ) * ( X[iN][0] - m_min[0] ) + m_min[0];
    }
    else
    {
      meshRact = X[iN][0];
    }

    X[iN][0] = meshRact * cos( meshTheta );
    X[iN][1] = meshRact * sin( meshTheta );

    // Add mapped values to nodesets
//    if( m_meshType == MeshType::CylindricalSquareBoundary )
//    {
//      if( isEqual( X[iN][0], -1 * m_max[0], m_coordinatePrecision ) )
//      {
//        xnegNodes.insert( iN );
//      }
//      if( isEqual( X[iN][0], m_max[0], m_coordinatePrecision ) )
//      {
//        xposNodes.insert( iN );
//      }
//      if( isEqual( X[iN][1], -1 * m_max[0], m_coordinatePrecision ) )
//      {
//        ynegNodes.insert( iN );
//      }
//      if( isEqual( X[iN][1], m_max[0], m_coordinatePrecision ) )
//      {
//        yposNodes.insert( iN );
//      }
//    }
  }
}


REGISTER_CATALOG_ENTRY( MeshGeneratorBase, InternalWellboreGenerator, string const &, Group * const )
} /* namespace geosx */
