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


#ifndef GEOSX_MESH_PARTICLEELEMENTSUBREGION_HPP_
#define GEOSX_MESH_PARTICLEELEMENTSUBREGION_HPP_

#include "mesh/generators/ParticleBlockABC.hpp"
#include "mesh/utilities/ComputationalGeometry.hpp"
#include "ParticleSubRegionBase.hpp"


namespace geosx
{

class MeshLevel;

/**
 * @class ParticleSubRegion
 * Class specializing the particle subregion for a cell particle.
 * This is the class used in the physics solvers to represent a collection of mesh cell particles
 */
class ParticleSubRegion : public ParticleSubRegionBase
{
public:

  /**
   * @brief Const getter for the catalog name.
   * @return the name of this type in the catalog
   */
  static const string catalogName()
  { return "ParticleSubRegion"; }

  /**
   * @copydoc catalogName()
   */
  virtual const string getCatalogName() const override final
  { return ParticleSubRegion::catalogName(); }

  /**
   * @name Constructor / Destructor
   */
  ///@{

  /**
   * @brief Constructor for this class.
   * @param[in] name the name of this object manager
   * @param[in] parent the parent Group
   */
  ParticleSubRegion( string const & name, Group * const parent );

  /**
   * @brief Destructor.
   */
  virtual ~ParticleSubRegion() override;

  ///@}

  /**
   * @name Helpers for ParticleSubRegion construction
   */
  ///@{

  /**
   * @brief Fill the ParticleSubRegion by copying those of the source ParticleBlock
   * @param particleBlock the ParticleBlock which properties (connectivity info) will be copied.
   */
  void copyFromParticleBlock( ParticleBlockABC & particleBlock );

  ///@}

  void updateRVectors(int const p,
                      LvArray::ArraySlice<double, 2, 1, long> const & p_F)
  {
    if(m_hasRVectors)
    {
//      real64 rTemp[3][3]; // Might be able to use this for an incremental update? Using I+L*dt?
//      for(int i=0; i<3; i++)
//      {
//        for(int j=0; j<3; j++)
//        {
//          rTemp[i][j] = m_particleRVectors[p][i][j];
//        }
//      }
      for(int i=0; i<3; i++)
      {
        for(int j=0; j<3; j++) // This seems wrong to me, but I checked it against a grid-imposed rigid body rotation velocity field and everything looked good...
        {
          m_particleRVectors[p][j][i] = p_F[0][i]*m_particleRVectors0[p][j][0] + p_F[1][i]*m_particleRVectors0[p][j][1] + p_F[2][i]*m_particleRVectors0[p][j][2];
        }
      }
    }
  }

  void getAllWeights(int const p,
                     LvArray::ArraySlice<double, 1, 0, long> const & p_x,
                     std::vector<real64> const & xMin,
                     std::vector<real64> const & hx,
                     std::vector<std::vector<std::vector<int>>> const & ijkMap,
                     arrayView2d< real64, nodes::REFERENCE_POSITION_USD > const & g_X,
                     std::vector<int> & nodeIDs,
                     std::vector<real64> & weights,
                     std::vector< std::vector<real64> > & gradWeights)
  {
    switch( m_particleType)
    {
      case ParticleType::SinglePoint:
      {
        // get cell IDs
        std::vector<int> cellID;
        cellID.resize(3);
        for(int i=0; i<3; i++)
        {
          cellID[i] = std::floor((p_x[i] - xMin[i])/hx[i]);
        }

        // get node IDs
        for(int i=0; i<2; i++)
        {
          for(int j=0; j<2; j++)
          {
            for(int k=0; k<2; k++)
            {
              nodeIDs.push_back(ijkMap[cellID[0]+i][cellID[1]+j][cellID[2]+k]);
            }
          }
        }

        // get weights and grad weights
        int corner = ijkMap[cellID[0]][cellID[1]][cellID[2]];
        auto corner_x = g_X[corner];

        real64 xRel = (p_x[0] - corner_x[0])/hx[0];
        real64 yRel = (p_x[1] - corner_x[1])/hx[1];
        real64 zRel = (p_x[2] - corner_x[2])/hx[2];

        for(int i=0; i<2; i++)
        {
          real64 xWeight = i*xRel + (1-i)*(1.0-xRel);
          real64 dxWeight = i/hx[0] - (1-i)/hx[0];
          for(int j=0; j<2; j++)
          {
            real64 yWeight = j*yRel + (1-j)*(1.0-yRel);
            real64 dyWeight = j/hx[1] - (1-j)/hx[1];
            for(int k=0; k<2; k++)
            {
              real64 zWeight = k*zRel + (1-k)*(1.0-zRel);
              real64 dzWeight = k/hx[2] - (1-k)/hx[2];
              weights.push_back(xWeight*yWeight*zWeight);
              gradWeights[0].push_back(dxWeight*yWeight*zWeight);
              gradWeights[1].push_back(xWeight*dyWeight*zWeight);
              gradWeights[2].push_back(xWeight*yWeight*dzWeight);
            }
          }
        }

        break;
      }

      case ParticleType::CPDI:
      {
        // Precalculated things
        int signs[8][3] = { { 1,  1,  1},
                            { 1,  1, -1},
                            { 1, -1,  1},
                            { 1, -1, -1},
                            {-1,  1,  1},
                            {-1,  1, -1},
                            {-1, -1,  1},
                            {-1, -1, -1} };

        real64 alpha[8][8];
        real64 one_over_V = 1.0 / m_particleVolume[p];
        real64 p_r1[3], p_r2[3], p_r3[3]; // allowing 1-indexed r-vectors to persist to torture future postdocs >:)

        for(int i=0; i<3; i++)
        {
          p_r1[i] = 0.5*m_particleRVectors[p][0][i];
          p_r2[i] = 0.5*m_particleRVectors[p][1][i];
          p_r3[i] = 0.5*m_particleRVectors[p][2][i];
        }

        alpha[0][0] = one_over_V * ( p_r1[2] * p_r2[1] - p_r1[1] * p_r2[2] - p_r1[2] * p_r3[1] + p_r2[2] * p_r3[1] + p_r1[1] * p_r3[2] - p_r2[1] * p_r3[2] );
        alpha[0][1] = one_over_V * ( -( p_r1[2] * p_r2[0] ) + p_r1[0] * p_r2[2] + p_r1[2] * p_r3[0] - p_r2[2] * p_r3[0] - p_r1[0] * p_r3[2] + p_r2[0] * p_r3[2] );
        alpha[0][2] = one_over_V * ( p_r1[1] * p_r2[0] - p_r1[0] * p_r2[1] - p_r1[1] * p_r3[0] + p_r2[1] * p_r3[0] + p_r1[0] * p_r3[1] - p_r2[0] * p_r3[1] );
        alpha[1][0] = one_over_V * ( p_r1[2] * p_r2[1] - p_r1[1] * p_r2[2] - p_r1[2] * p_r3[1] - p_r2[2] * p_r3[1] + p_r1[1] * p_r3[2] + p_r2[1] * p_r3[2] );
        alpha[1][1] = one_over_V * ( -( p_r1[2] * p_r2[0] ) + p_r1[0] * p_r2[2] + p_r1[2] * p_r3[0] + p_r2[2] * p_r3[0] - p_r1[0] * p_r3[2] - p_r2[0] * p_r3[2] );
        alpha[1][2] = one_over_V * ( p_r1[1] * p_r2[0] - p_r1[0] * p_r2[1] - p_r1[1] * p_r3[0] - p_r2[1] * p_r3[0] + p_r1[0] * p_r3[1] + p_r2[0] * p_r3[1] );
        alpha[2][0] = one_over_V * ( p_r1[2] * p_r2[1] - p_r1[1] * p_r2[2] + p_r1[2] * p_r3[1] - p_r2[2] * p_r3[1] - p_r1[1] * p_r3[2] + p_r2[1] * p_r3[2] );
        alpha[2][1] = one_over_V * ( -( p_r1[2] * p_r2[0] ) + p_r1[0] * p_r2[2] - p_r1[2] * p_r3[0] + p_r2[2] * p_r3[0] + p_r1[0] * p_r3[2] - p_r2[0] * p_r3[2] );
        alpha[2][2] = one_over_V * ( p_r1[1] * p_r2[0] - p_r1[0] * p_r2[1] + p_r1[1] * p_r3[0] - p_r2[1] * p_r3[0] - p_r1[0] * p_r3[1] + p_r2[0] * p_r3[1] );
        alpha[3][0] = one_over_V * ( p_r1[2] * p_r2[1] - p_r1[1] * p_r2[2] + p_r1[2] * p_r3[1] + p_r2[2] * p_r3[1] - p_r1[1] * p_r3[2] - p_r2[1] * p_r3[2] );
        alpha[3][1] = one_over_V * ( -( p_r1[2] * p_r2[0] ) + p_r1[0] * p_r2[2] - p_r1[2] * p_r3[0] - p_r2[2] * p_r3[0] + p_r1[0] * p_r3[2] + p_r2[0] * p_r3[2] );
        alpha[3][2] = one_over_V * ( p_r1[1] * p_r2[0] - p_r1[0] * p_r2[1] + p_r1[1] * p_r3[0] + p_r2[1] * p_r3[0] - p_r1[0] * p_r3[1] - p_r2[0] * p_r3[1] );
        alpha[4][0] = -alpha[2][0];
        alpha[4][1] = -alpha[2][1];
        alpha[4][2] = -alpha[2][2];
        alpha[5][0] = -alpha[3][0];
        alpha[5][1] = -alpha[3][1];
        alpha[5][2] = -alpha[3][2];
        alpha[6][0] = -alpha[0][0];
        alpha[6][1] = -alpha[0][1];
        alpha[6][2] = -alpha[0][2];
        alpha[7][0] = -alpha[1][0];
        alpha[7][1] = -alpha[1][1];
        alpha[7][2] = -alpha[1][2];

        // GEOS-to-GEOSX corner mapping, because I'm lazy
        int cornerMap[8] = {0, 4, 3, 7, 1, 5, 2, 6};

        // get cell IDs
        std::vector<std::vector<int>> cellID;
        cellID.resize(8); // CPDI can map to up to 8 cells
        for(int cell=0; cell<8; cell++)
        {
          cellID[cell].resize(3);
          for(int i=0; i<3; i++)
          {
            real64 CPDIcorner = p_x[i] + 0.5*(signs[cell][0]*m_particleRVectors[p][0][i] + signs[cell][1]*m_particleRVectors[p][1][i] + signs[cell][2]*m_particleRVectors[p][2][i]); // using full r-vectors, hence the 0.5
            cellID[cell][i] = std::floor((CPDIcorner - xMin[i])/hx[i]);
          }
        }

        // get node IDs
        for(size_t cell=0; cell<cellID.size(); cell++)
        {
          for(int i=0; i<2; i++)
          {
            for(int j=0; j<2; j++)
            {
              for(int k=0; k<2; k++)
              {
                nodeIDs.push_back(ijkMap[cellID[cell][0]+i][cellID[cell][1]+j][cellID[cell][2]+k]);
              }
            }
          }
        }

        // get weights and grad weights
        for(size_t cell=0; cell<cellID.size(); cell++)
        {
          int corner = ijkMap[cellID[cell][0]][cellID[cell][1]][cellID[cell][2]];
          auto corner_x = g_X[corner];

          real64 x, y, z;
          x = p_x[0] + 0.5*(signs[cell][0]*m_particleRVectors[p][0][0] + signs[cell][1]*m_particleRVectors[p][1][0] + signs[cell][2]*m_particleRVectors[p][2][0]);
          y = p_x[1] + 0.5*(signs[cell][0]*m_particleRVectors[p][0][1] + signs[cell][1]*m_particleRVectors[p][1][1] + signs[cell][2]*m_particleRVectors[p][2][1]);
          z = p_x[2] + 0.5*(signs[cell][0]*m_particleRVectors[p][0][2] + signs[cell][1]*m_particleRVectors[p][1][2] + signs[cell][2]*m_particleRVectors[p][2][2]);

          real64 xRel = (x - corner_x[0])/hx[0];
          real64 yRel = (y - corner_x[1])/hx[1];
          real64 zRel = (z - corner_x[2])/hx[2];

          for(int i=0; i<2; i++)
          {
            real64 xWeight = i*xRel + (1-i)*(1.0-xRel);
            for(int j=0; j<2; j++)
            {
              real64 yWeight = j*yRel + (1-j)*(1.0-yRel);
              for(int k=0; k<2; k++)
              {
                real64 zWeight = k*zRel + (1-k)*(1.0-zRel);
                real64 weight = xWeight*yWeight*zWeight;
                weights.push_back(0.125*weight); // note the built-in factor of 1/8 so we don't need it in the solver
                gradWeights[0].push_back(alpha[cornerMap[cell]][0]*weight);
                gradWeights[1].push_back(alpha[cornerMap[cell]][1]*weight);
                gradWeights[2].push_back(alpha[cornerMap[cell]][2]*weight);
              }
            }
          }
        }

        break;
      }

      default:
      {
        GEOSX_ERROR( "Invalid particle type: " << m_particleType );
      }
    }

  }

  /**
   * @name Overriding packing / Unpacking functions
   */
  ///@{

  virtual void viewPackingExclusionList( SortedArray< localIndex > & exclusionList ) const override;

  virtual localIndex packUpDownMapsSize( arrayView1d< localIndex const > const & packList ) const override;

  virtual localIndex packUpDownMaps( buffer_unit_type * & buffer,
                                     arrayView1d< localIndex const > const & packList ) const override;

  virtual localIndex unpackUpDownMaps( buffer_unit_type const * & buffer,
                                       array1d< localIndex > & packList,
                                       bool const overwriteUpMaps,
                                       bool const overwriteDownMaps ) override;

  ///@}

  /**
   * @name Miscellaneous
   */
  ///@{

  /**
   * @brief Helper function to apply a lambda function over all constructive groups
   * @tparam LAMBDA the type of the lambda function
   * @param lambda the lambda function
   */
  template< typename LAMBDA >
  void forMaterials( LAMBDA lambda )
  {

    for( auto & constitutiveGroup : m_constitutiveGrouping )
    {
      lambda( constitutiveGroup );
    }
  }

  ///@}

  /**
   * @brief struct to serve as a container for variable strings and keys
   * @struct viewKeyStruct
   */
  struct viewKeyStruct : public ParticleSubRegionBase::viewKeyStruct
  {
    /// @return String key for the constitutive point volume fraction
    static constexpr char const * constitutivePointVolumeFractionString() { return "ConstitutivePointVolumeFraction"; }
    /// @return String key for the derivatives of the shape functions with respect to the reference configuration
    static constexpr char const * dNdXString() { return "dNdX"; }
    /// @return String key for the derivative of the jacobian.
    static constexpr char const * detJString() { return "detJ"; }
    /// @return String key for the constitutive grouping
    static constexpr char const * constitutiveGroupingString() { return "ConstitutiveGrouping"; }
    /// @return String key for the constitutive map
    static constexpr char const * constitutiveMapString() { return "ConstitutiveMap"; }

    /// ViewKey for the constitutive grouping
    dataRepository::ViewKey constitutiveGrouping  = { constitutiveGroupingString() };
    /// ViewKey for the constitutive map
    dataRepository::ViewKey constitutiveMap       = { constitutiveMapString() };
  }
  /// viewKey struct for the ParticleSubRegion class
  m_ParticleBlockSubRegionViewKeys;

  virtual viewKeyStruct & viewKeys() override { return m_ParticleBlockSubRegionViewKeys; }
  virtual viewKeyStruct const & viewKeys() const override { return m_ParticleBlockSubRegionViewKeys; }

  /**
   * @brief Get the local indices of the nodes in a face of the particle.
   * @param[in] particleIndex The local index of the target particle.
   * @param[in] localFaceIndex The local index of the target face in the particle (this will be [0, numFacesInParticle[)
   * @param[out] nodeIndices A reference to the array of node indices of the face. Gets resized at the proper size.
   * @deprecated This method will be removed soon.
   */
  void getFaceNodes( localIndex const particleIndex,
                     localIndex const localFaceIndex,
                     array1d< localIndex > & nodeIndices ) const;





   /**
   * @brief @return The array of shape function derivatives.
   */
  array4d< real64 > & dNdX()
  { return m_dNdX; }

  /**
   * @brief @return The array of shape function derivatives.
   */
  arrayView4d< real64 const > dNdX() const
  { return m_dNdX; }

  /**
   * @brief @return The array of jacobian determinantes.
   */
  array2d< real64 > & detJ()
  { return m_detJ; }

  /**
   * @brief @return The array of jacobian determinantes.
   */
  arrayView2d< real64 const > detJ() const
  { return m_detJ; }

private:

  /// Map used for constitutive grouping
  map< string, localIndex_array > m_constitutiveGrouping;

  /// Array of constitutive point volume fraction
  array3d< real64 > m_constitutivePointVolumeFraction;

  /// Name of the properties registered from an external mesh
  string_array m_externalPropertyNames;

  /// The array of shape function derivaties.
  array4d< real64 > m_dNdX;

  /// The array of jacobian determinantes.
  array2d< real64 > m_detJ;

  /**
   * @brief Pack particle-to-node and particle-to-face maps
   * @tparam the flag for the bufferOps::Pack function
   * @param buffer the buffer used in the bufferOps::Pack function
   * @param packList the packList used in the bufferOps::Pack function
   * @return the pack size
   */
  template< bool DOPACK >
  localIndex packUpDownMapsPrivate( buffer_unit_type * & buffer,
                                    arrayView1d< localIndex const > const & packList ) const;


};

} /* namespace geosx */

#endif /* GEOSX_MESH_CELLELEMENTSUBREGION_HPP_ */
