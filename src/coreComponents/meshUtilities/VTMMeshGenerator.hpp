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
 * VTMMeshGenerator.h
 *
 *  Created on: Aug 16, 2018
 *      Author: Antoine Mazuyer
 */

#ifndef VTMMESHGENERATOR_H_
#define VTMMESHGENERATOR_H_

#include "dataRepository/ManagedGroup.hpp"
#include "codingUtilities/Utilities.hpp"
#include "MeshGeneratorBase.hpp"

#include "fileIO/vtm/VtmFile.hpp"

#ifdef USE_ATK
#include <slic/slic.hpp>
#endif

namespace geosx
{

namespace dataRepository
{
namespace keys
{
string const filePath = "file";
}
}

class NodeManager;
class DomainPartition;

class VTMMeshGenerator : public MeshGeneratorBase
{
public:
  VTMMeshGenerator( const std::string& name,
                         ManagedGroup * const parent );

  virtual ~VTMMeshGenerator() override;

  static string CatalogName() { return "MeshFile"; }

  virtual void GenerateElementRegions( DomainPartition& domain ) override;

  virtual ManagedGroup * CreateChild( string const & childKey, string const & childName ) override;

  virtual void GenerateMesh( DomainPartition * const domain ) override;

  // virtual void GenerateNodesets( xmlWrapper::xmlNode const & targetNode,
  //                                NodeManager * nodeManager ) override;

  virtual void GetElemToNodesRelationInBox ( const std::string& elementType,
                                             const int index[],
                                             const int& iEle,
                                             int nodeIDInBox[],
                                             const int size) override;

  virtual void RemapMesh ( dataRepository::ManagedGroup * const domain ) override;

//  int m_delayMeshDeformation;

protected:
  void PostProcessInput() override final;

private:

  /// Contains the path to the VTM file
  string m_fileName;

  VtmFile m_vtmFile;

  inline globalIndex NodeGlobalIndex( const int index[3] )
  {
    globalIndex rval = 0;
      /*

    rval = index[0]*(m_numElemsTotal[1]+1)*(m_numElemsTotal[2]+1) + index[1]*(m_numElemsTotal[2]+1) + index[2];
    */
    return rval;
  }

  inline globalIndex ElemGlobalIndex( const int index[3] )
  {
    globalIndex rval = 0;
      /*

    rval = index[0]*m_numElemsTotal[1]*m_numElemsTotal[2] + index[1]*m_numElemsTotal[2] + index[2];
    */
    return rval;
  }

  inline R1Tensor NodePosition( const int a[3], int trianglePattern )
  {
    R1Tensor X;
      /*
    realT xInterval(0);

    int xPosIndex = 0;
    if (trianglePattern == 1)
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


      X[i] = min + (max-min) * ( double( a[i] - startingIndex ) / m_nElems[i][block] );

      if (( !isZero(m_nElemBias[i][block]) ) & (m_nElems[i][block]>1))
      {
        if (fabs(m_nElemBias[i][block]) >= 1)
        {
#ifdef USE_ATK
          SLIC_ERROR("Mesh bias must between -1 and 1!");
#endif
        }

        realT len = max -  min;
        realT xmean = len / m_nElems[i][block];
        realT x0 = xmean * double( a[i] - startingIndex );
        realT chi = m_nElemBias[i][block]/(xmean/len - 1.0);
        realT dx = -x0*chi + x0*x0*chi/len;
        X[i] += dx;
      }

      // This is for creating regular triangle pattern
      if (i==0) xInterval = (max-min) / m_nElems[i][block];
      if (trianglePattern == 1 && i == 1 && a[1] % 2 == 1 && a[0] != 0 && a[0] != xPosIndex)
        X[0] -= xInterval * 0.5;
    }

    */
    return X;
  }

  inline R1Tensor ElemCenterPosition( const int k[3] )
  {
    R1Tensor X;

    /*
    for( int i=0 ; i<3 ; ++i )
    {
      X[i] = m_min[i] + (m_max[i]-m_min[i]) * ( ( k[i] + 0.5 ) / m_numElemsTotal[i] );
    }
    */

    return X;
  }

public:
  inline bool isRadial()
  {
      /*
    bool rval = (m_mapToRadial > 0);
    return rval;
    */
      return false;
  }

};
}

#endif /* INTERNALMESHGENERATOR_H_ */
