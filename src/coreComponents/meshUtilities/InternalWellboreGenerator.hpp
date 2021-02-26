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
 * @file InternalWellboreGenerator.hpp
 */

#ifndef GEOSX_MESHUTILITIES_INTERNALWELLBOREGENERATOR_HPP_
#define GEOSX_MESHUTILITIES_INTERNALWELLBOREGENERATOR_HPP_

#include "dataRepository/Group.hpp"
#include "codingUtilities/Utilities.hpp"
#include "MeshGeneratorBase.hpp"

#include <stdlib.h>     
#include <time.h> 

namespace geosx
{

class NodeManager;
class DomainPartition;
/**
 * @class InternalWellboreGenerator
 * @brief The InternalWellboreGenerator class is a class handling GEOSX generated meshes.
 */
class InternalWellboreGenerator : public MeshGeneratorBase
{
public:

  /**
   * @brief Main constructor for InternalWellboreGenerator.
   * @param[in] name of the InternalWellboreGenerator
   * @param[in] parent point to the parent Group of the InternalWellboreGenerator
   */
  InternalWellboreGenerator( const string & name,
                             Group * const parent );

  virtual ~InternalWellboreGenerator() override = default;

  /**
   * @brief Return the name of the InternalWellboreGenerator in object Catalog.
   * @return string that contains the key name to InternalWellboreGenerator in the Catalog
   */
  static string catalogName() { return "InternalWellbore"; }

  struct viewKeyStruct
  {
    constexpr static char const * radialCoordsString() { return "radialCoords"; }
    constexpr static char const * tangentCoordsString() { return "tangentCoords"; }
    constexpr static char const * axialCoordsString() { return "axialCoords"; }
    constexpr static char const * radialElemsString() { return "nr"; }
    constexpr static char const * tangentElemsString() { return "nt"; }
    constexpr static char const * axialElemsString() { return "na"; }
    constexpr static char const * radialBiasString() { return "rBias"; }
    constexpr static char const * axialBiasString() { return "aBias"; }
    constexpr static char const * cellBlockNamesString() { return "cellBlockNames"; }
    constexpr static char const * elementTypesString() { return "elementTypes"; }
    constexpr static char const * trianglePatternString() { return "trianglePattern"; }
  };

  virtual void generateElementRegions( DomainPartition & domain ) override;

  /**
   * @brief Create a new geometric object (box, plane, etc) as a child of this group.
   * @param childKey the catalog key of the new geometric object to create
   * @param childName the name of the new geometric object in the repository
   * @return the group child
   */
  virtual Group * createChild( string const & childKey, string const & childName ) override;

  virtual void generateMesh( DomainPartition & domain ) override;

  // virtual void GenerateNodesets( xmlWrapper::xmlNode const & targetNode,
  //                                NodeManager * nodeManager ) override;

  virtual void getElemToNodesRelationInBox ( const string & elementType,
                                             const int index[],
                                             const int & iEle,
                                             int nodeIDInBox[],
                                             const int size ) override;

  virtual void remapMesh ( dataRepository::Group & domain ) override;

//  int m_delayMeshDeformation;

protected:

  void postProcessInput() override final;

private:

  /// Mesh number of dimension
  int m_dim;

  /// Array of vertex coordinates
  array1d< real64 > m_vertices[3];

  /// Ndim x nElem spatialized for element indexes
  integer_array m_nElems[3];

  /// Ndim x nElem spatialized array of element scaling factors
  array1d< real64 > m_nElemScaling[3];

  //bool m_useBias = false;
  /// Ndim x nElem spatialized array of element bias
  array1d< real64 > m_nElemBias[3];

  /// String array of region names
  string_array m_regionNames;

  /// Minimum extent of mesh dimensions
  real64 m_min[3];

  /// Maximum extent of mesh dimensions
  real64 m_max[3];

  //int m_numElems[3];
  /// Ndim x nBlock spatialized array of first elemnt index in the cellBlock
  integer_array m_firstElemIndexForBlock[3];

  /// Ndim x nBlock spatialized array of last elemnt index in the cellBlock
  integer_array m_lastElemIndexForBlock[3];

  /// Array of number of elements per direction
  int m_numElemsTotal[3];

  // String array listing the element type present
  string_array m_elementType;

  /// Array of number of element per box
  array1d< integer > m_numElePerBox;

  /**
   * @brief Member variable for triangle pattern seletion.
   * @note In pattern 0, half nodes have 4 edges and the other half have 8; for Pattern 1, every node has 6.
   */
  int m_trianglePattern;

  /// Node perturbation amplitude value
  real64 m_fPerturb=0.0;

  /// Random seed for generation of the node perturbation field
  int m_randSeed;

  /**
   * @brief Knob to map onto a radial mesh.
   * @note if 0 mesh is not radial, if positive mesh is, it larger than 1
   */
  int m_mapToRadial = 0;
  ///@cond DO_NOT_DOCUMENT
  /// axis index for cartesian to radial coordinates mapping
  // internal temp var
  int m_meshAxis;
  real64 m_meshTheta;
  real64 m_meshPhi;
  real64 m_meshRout;
  real64 m_meshRact;
  /// @endcond

  /// skew angle in radians for skewed mesh generation
  real64 m_skewAngle = 0;
  /// skew center for skew mesh generation
  real64 m_skewCenter[3] = { 0, 0, 0 };


  ///@cond DO_NOT_DOCUMENT
  //unused
  string m_meshDx, m_meshDy, m_meshDz;
  ///@endcond

/**
 * @brief Convert ndim node spatialized index to node global index.
 * @param[in] node ndim spatialized array index
 */
  inline globalIndex nodeGlobalIndex( const int index[3] )
  {
    globalIndex rval = 0;

    rval = index[0]*(m_numElemsTotal[1]+1)*(m_numElemsTotal[2]+1) + index[1]*(m_numElemsTotal[2]+1) + index[2];
    return rval;
  }

/**
 * @brief Convert ndim element spatialized index to element global index.
 * @param[in] element ndim spatialized array index
 */
  inline globalIndex elemGlobalIndex( const int index[3] )
  {
    globalIndex rval = 0;

    rval = index[0]*m_numElemsTotal[1]*m_numElemsTotal[2] + index[1]*m_numElemsTotal[2] + index[2];
    return rval;
  }

  /**
   * @brief Construct the node position for a spatially indexed node.
   * @tparam OUT_VECTOR type of output vector X
   * @param[in] a ndim spatial index for the considered node
   * @param[in] trianglePattern triangle pattern identifier
   * @param[out] X the node coordinates
   *
   * @note In pattern 0, half nodes have 4 edges and the other half have 8; for Pattern 1, every node has 6.
   */
  template< typename OUT_VECTOR >
  inline void getNodePosition( int const * a, int trianglePattern, OUT_VECTOR && X )
  {
    real64 xInterval( 0 );

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
      real64 min = m_vertices[i][block];
      real64 max = m_vertices[i][block+1];


      X[i] = min + (max-min) * ( double( a[i] - startingIndex ) / m_nElems[i][block] );

      // First check if m_nElemBias contains values
      // Otherwise the next test will cause a segfault when looking for "block"
      if( m_nElemBias[i].size()>0 )
      {
        // Verify that the bias is non-zero and applied to more than one block:
        if( ( !isZero( m_nElemBias[i][block] ) ) && (m_nElems[i][block]>1))
        {
          GEOSX_ERROR_IF( fabs( m_nElemBias[i][block] ) >= 1, "Mesh bias must between -1 and 1!" );

          real64 len = max -  min;
          real64 xmean = len / m_nElems[i][block];
          real64 x0 = xmean * double( a[i] - startingIndex );
          real64 chi = m_nElemBias[i][block]/(xmean/len - 1.0);
          real64 dx = -x0*chi + x0*x0*chi/len;
          X[i] += dx;
        }
      }

      // This is for creating regular triangle pattern
      if( i==0 ) xInterval = (max-min) / m_nElems[i][block];
      if( trianglePattern == 1 && i == 1 && a[1] % 2 == 1 && a[0] != 0 && a[0] != xPosIndex )
        X[0] -= xInterval * 0.5;
    }
  }

  /**
   * @brief
   * @tparam OUT_VECTOR type of output vector X
   * @param[in] k the ijk-index of the element
   * @param[out] X the element center coordinates
   */
  template< typename OUT_VECTOR >
  inline void getElemCenterPosition( const int k[3], OUT_VECTOR && X )
  {
    for( int i=0; i<3; ++i )
    {
      X[i] = m_min[i] + (m_max[i]-m_min[i]) * ( ( k[i] + 0.5 ) / m_numElemsTotal[i] );
    }
  }

public:
  /**
   * @brief Check if the mesh is a radial mesh.
   * @return true if the Internal mesh is radial, false else
   */
  inline bool isRadial()
  {
    bool rval = (m_mapToRadial > 0);
    return rval;
  }

};

} /* namespace geosx */

#endif /* GEOSX_MESHUTILITIES_INTERNALWELLBOREGENERATOR_HPP_ */

