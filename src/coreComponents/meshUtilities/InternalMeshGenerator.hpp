/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/*
 * InternalMeshGenerator.h
 *
 *  Created on: Nov 19, 2012
 *      Author: settgast1
 */

#ifndef INTERNALMESHGENERATOR_H_
#define INTERNALMESHGENERATOR_H_

#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

#include "dataRepository/ManagedGroup.hpp"
#include "codingUtilities/Utilities.hpp"
#include "MeshGeneratorBase.hpp"

namespace geosx
{

namespace dataRepository
{
namespace keys
{
string const xCoords = "xCoords";
string const yCoords = "yCoords";
string const zCoords = "zCoords";
string const xElems  = "nx";
string const yElems  = "ny";
string const zElems  = "nz";
string const xBias  = "xbias";
string const yBias  = "ybias";
string const zBias  = "zbias";
string const cellBlockNames = "cellBlockNames";
string const elementTypes = "elementTypes";
string const trianglePattern = "trianglePattern";
}
}

class NodeManager;
class DomainPartition;

class InternalMeshGenerator : public MeshGeneratorBase
{
public:
  InternalMeshGenerator( const std::string & name,
                         ManagedGroup * const parent );

  virtual ~InternalMeshGenerator() override;

  static string CatalogName() { return "InternalMesh"; }

  //  void ProcessInputFile( xmlWrapper::xmlNode const & targetNode ) override;
  //
  //

  virtual void GenerateElementRegions( DomainPartition & domain ) override;

  virtual ManagedGroup * CreateChild( string const & childKey, string const & childName ) override;

  virtual void GenerateMesh( DomainPartition * const domain ) override;

  // virtual void GenerateNodesets( xmlWrapper::xmlNode const & targetNode,
  //                                NodeManager * nodeManager ) override;

  virtual void GetElemToNodesRelationInBox( const std::string & elementType,
                                            const int index[],
                                            const int & iEle,
                                            int nodeIDInBox[],
                                            const int size ) override;

  virtual void RemapMesh( dataRepository::ManagedGroup * const domain ) override;

  //  int m_delayMeshDeformation;

protected:
  void PostProcessInput() override final;

private:

  int m_dim;
  array1d<real64> m_vertices[3];
  integer_array m_nElems[3];
  array1d<real64> m_nElemScaling[3];

  bool m_useBias = false;
  array1d<real64> m_nElemBias[3];

  string_array m_regionNames;

  realT m_min[3]; // Minimum extent of mesh dimensions
  realT m_max[3]; // Maximum extent of mesh dimensions

  //int m_numElems[3];
  integer_array m_firstElemIndexForBlock[3];
  integer_array m_lastElemIndexForBlock[3];



  //  realT m_wExtensionMin[3];
  //  realT m_wExtensionMax[3];
  //  int m_nExtensionLayersMin[3];
  //  int m_nExtensionLayersMax[3];
  //  realT m_extendedMin[3];
  //  realT m_extendedMax[3]; // This is the domain size after we apply n layers
  // of elements which are of the same size as the core elements.  We will move
  // these nodes to where they should be later when we finish the meshing.
  int m_numElemsTotal[3];
  //  realT m_commonRatioMin[3];
  //  realT m_commonRatioMax[3];

  string_array m_elementType;

  array1d<integer> m_numElePerBox;

  int m_trianglePattern;   // In pattern 0, half nodes have 4 edges and the
  // other half have 8; for Pattern 1, every node has
  // 6.

  realT m_fPerturb=0.0;
  int m_randSeed;

  int m_mapToRadial = 0;
  int m_meshAxis;
  realT m_meshTheta;
  realT m_meshPhi;
  realT m_meshRout;
  realT m_meshRact;

  realT m_skewAngle = 0;
  R1Tensor m_skewCenter = {0, 0, 0};

  std::string m_meshDx, m_meshDy, m_meshDz;

  inline globalIndex NodeGlobalIndex( const int index[3] )
  {
    globalIndex rval = 0;

    rval = index[0]*( m_numElemsTotal[1]+1 )*( m_numElemsTotal[2]+1 ) + index[1]*( m_numElemsTotal[2]+1 ) + index[2];
    return rval;
  }

  inline globalIndex ElemGlobalIndex( const int index[3] )
  {
    globalIndex rval = 0;

    rval = index[0]*m_numElemsTotal[1]*m_numElemsTotal[2] + index[1]*m_numElemsTotal[2] + index[2];
    return rval;
  }

  inline R1Tensor NodePosition( const int a[3], int trianglePattern )
  {
    R1Tensor X;
    realT xInterval( 0 );

    int xPosIndex = 0;
    if( trianglePattern == 1 )
    {
      int startingIndex = 0;
      int endingIndex = 0;
      int block = 0;
      for( block=0 ; block<m_nElems[0].size() ; ++block )
      {
        startingIndex = endingIndex;
        endingIndex = startingIndex + m_nElems[0][block];
      }
      xPosIndex = endingIndex;
    }

    for( int i=0 ; i<3 ; ++i )
    {

      int startingIndex = 0;
      int endingIndex = 0;
      int block = 0;
      for( block=0 ; block<m_nElems[i].size() ; ++block )
      {
        startingIndex = endingIndex;
        endingIndex = startingIndex + m_nElems[i][block];
        if( a[i]>=startingIndex && a[i]<=endingIndex )
        {
          break;
        }
      }
      realT min = m_vertices[i][block];
      realT max = m_vertices[i][block+1];


      X[i] = min + ( max-min ) * ( double( a[i] - startingIndex ) / m_nElems[i][block] );

      if( m_useBias && ( !isZero( m_nElemBias[i][block] ) ) & ( m_nElems[i][block]>1 ) )
      {
        GEOS_ERROR_IF( fabs( m_nElemBias[i][block] ) >= 1, "Mesh bias must between -1 and 1!" );

        realT len = max -  min;
        realT xmean = len / m_nElems[i][block];
        realT x0 = xmean * double( a[i] - startingIndex );
        realT chi = m_nElemBias[i][block]/( xmean/len - 1.0 );
        realT dx = -x0*chi + x0*x0*chi/len;
        X[i] += dx;
      }

      // This is for creating regular triangle pattern
      if( i==0 ) { xInterval = ( max-min ) / m_nElems[i][block]; }
      if( trianglePattern == 1 && i == 1 && a[1] % 2 == 1 && a[0] != 0 && a[0] != xPosIndex )
      {
        X[0] -= xInterval * 0.5;
      }
    }

    return X;
  }

  inline R1Tensor ElemCenterPosition( const int k[3] )
  {
    R1Tensor X;

    for( int i=0 ; i<3 ; ++i )
    {
      X[i] = m_min[i] + ( m_max[i]-m_min[i] ) * ( ( k[i] + 0.5 ) / m_numElemsTotal[i] );
    }

    return X;
  }

public:
  inline bool isRadial()
  {
    bool rval = ( m_mapToRadial > 0 );
    return rval;
  }

};
}

#endif /* INTERNALMESHGENERATOR_H_ */
