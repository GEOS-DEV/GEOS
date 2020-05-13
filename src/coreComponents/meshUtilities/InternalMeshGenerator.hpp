/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file InternalMeshGenerator.h
 */

#ifndef INTERNALMESHGENERATOR_H_
#define INTERNALMESHGENERATOR_H_

#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

#include "dataRepository/Group.hpp"
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
string const xBias  = "xBias";
string const yBias  = "yBias";
string const zBias  = "zBias";
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
                         Group * const parent );

  virtual ~InternalMeshGenerator() override;

  static string CatalogName() { return "InternalMesh"; }

//  void ProcessInputFile( xmlWrapper::xmlNode const & targetNode ) override;
//
//

  virtual void GenerateElementRegions( DomainPartition & domain ) override;

  virtual Group * CreateChild( string const & childKey, string const & childName ) override;

  virtual void GenerateMesh( DomainPartition * const domain ) override;

  // virtual void GenerateNodesets( xmlWrapper::xmlNode const & targetNode,
  //                                NodeManager * nodeManager ) override;

  virtual void GetElemToNodesRelationInBox ( const std::string & elementType,
                                             const int index[],
                                             const int & iEle,
                                             int nodeIDInBox[],
                                             const int size ) override;

  virtual void RemapMesh ( dataRepository::Group * const domain ) override;

//  int m_delayMeshDeformation;

protected:
  void PostProcessInput() override final;

private:

  int m_dim;
  array1d< real64 > m_vertices[3];
  integer_array m_nElems[3];
  array1d< real64 > m_nElemScaling[3];

  //bool m_useBias = false;
  array1d< real64 > m_nElemBias[3];

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

  array1d< integer > m_numElePerBox;

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

    rval = index[0]*(m_numElemsTotal[1]+1)*(m_numElemsTotal[2]+1) + index[1]*(m_numElemsTotal[2]+1) + index[2];
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
      for( block=0; block<m_nElems[0].size(); ++block )
      {
        startingIndex = endingIndex;
        endingIndex = startingIndex + m_nElems[0][block];
      }
      xPosIndex = endingIndex;
    }

    for( int i=0; i<3; ++i )
    {

      int startingIndex = 0;
      int endingIndex = 0;
      int block = 0;
      for( block=0; block<m_nElems[i].size(); ++block )
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


      X[i] = min + (max-min) * ( double( a[i] - startingIndex ) / m_nElems[i][block] );

      // First check if m_nElemBias contains values
      // Otherwise the next test will cause a segfault when looking for "block"
      if( m_nElemBias[i].size()>0 )
      {
        // Verify that the bias is non-zero and applied to more than one block:
        if( ( !isZero( m_nElemBias[i][block] ) ) && (m_nElems[i][block]>1))
        {
          GEOSX_ERROR_IF( fabs( m_nElemBias[i][block] ) >= 1, "Mesh bias must between -1 and 1!" );

          realT len = max -  min;
          realT xmean = len / m_nElems[i][block];
          realT x0 = xmean * double( a[i] - startingIndex );
          realT chi = m_nElemBias[i][block]/(xmean/len - 1.0);
          realT dx = -x0*chi + x0*x0*chi/len;
          X[i] += dx;
        }
      }

      // This is for creating regular triangle pattern
      if( i==0 ) xInterval = (max-min) / m_nElems[i][block];
      if( trianglePattern == 1 && i == 1 && a[1] % 2 == 1 && a[0] != 0 && a[0] != xPosIndex )
        X[0] -= xInterval * 0.5;
    }

    return X;
  }

  inline R1Tensor ElemCenterPosition( const int k[3] )
  {
    R1Tensor X;

    for( int i=0; i<3; ++i )
    {
      X[i] = m_min[i] + (m_max[i]-m_min[i]) * ( ( k[i] + 0.5 ) / m_numElemsTotal[i] );
    }

    return X;
  }

public:
  inline bool isRadial()
  {
    bool rval = (m_mapToRadial > 0);
    return rval;
  }

};
}

#endif /* INTERNALMESHGENERATOR_H_ */
