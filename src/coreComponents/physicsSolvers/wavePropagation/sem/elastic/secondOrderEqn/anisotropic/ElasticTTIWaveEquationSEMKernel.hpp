/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file ElasticTTIWaveEquationSEMKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ElasticTTIWaveEquationSEMKERNEL_HPP_
#define GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ElasticTTIWaveEquationSEMKERNEL_HPP_

#include "finiteElement/kernelInterface/KernelBase.hpp"
#include "physicsSolvers/wavePropagation/shared/WaveSolverUtils.hpp"
#include "physicsSolvers/wavePropagation/sem/elastic/shared/ElasticFields.hpp"
#include "ElasticVTIFields.hpp"
#include "ElasticTTIFields.hpp"

namespace geos
{

/// Namespace to contain the elastic wave kernels.
namespace elasticTTIWaveEquationSEMKernels
{



/**
 * @brief Implements kernels for solving the elastic wave equations
 *   explicit central FD method and SEM in the Tilted Transverse Isotropic (TTI) case
 * @copydoc geos::finiteElement::KernelBase
 * @tparam SUBREGION_TYPE The type of subregion that the kernel will act on.
 *
 * ### ElasticTTIWaveEquationSEMKernel Description
 * Implements the KernelBase interface functions required for solving
 * the acoustic wave equations using the
 * "finite element kernel application" functions such as
 * geos::finiteElement::RegionBasedKernelApplication.
 *
 * The number of degrees of freedom per support point for both
 * the test and trial spaces are specified as `1`.
 */


template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class ExplicitElasticTTISEM : public finiteElement::KernelBase< SUBREGION_TYPE,
                                                                CONSTITUTIVE_TYPE,
                                                                FE_TYPE,
                                                                1,
                                                                1 >
{
public:

  /// Alias for the base class;
  using Base = finiteElement::KernelBase< SUBREGION_TYPE,
                                          CONSTITUTIVE_TYPE,
                                          FE_TYPE,
                                          1,
                                          1 >;

  /// Number of nodes per element...which is equal to the
  /// numTestSupportPointPerElem and numTrialSupportPointPerElem by definition.
  static constexpr int numNodesPerElem = Base::maxNumTestSupportPointsPerElem;

  using Base::numDofPerTestSupportPoint;
  using Base::numDofPerTrialSupportPoint;
  using Base::m_elemsToNodes;
  using Base::m_elemGhostRank;
  using Base::m_constitutiveUpdate;
  using Base::m_finiteElementSpace;

//*****************************************************************************
  /**
   * @brief Constructor
   * @copydoc geos::finiteElement::KernelBase::KernelBase
   * @param nodeManager Reference to the NodeManager object.
   * @param edgeManager Reference to the EdgeManager object.
   * @param faceManager Reference to the FaceManager object.
   * @param targetRegionIndex Index of the region the subregion belongs to.
   * @param dt The time interval for the step.
   */
  ExplicitElasticTTISEM( NodeManager & nodeManager,
                         EdgeManager const & edgeManager,
                         FaceManager const & faceManager,
                         localIndex const targetRegionIndex,
                         SUBREGION_TYPE const & elementSubRegion,
                         FE_TYPE const & finiteElementSpace,
                         CONSTITUTIVE_TYPE & inputConstitutiveType,
                         real64 const dt ):
    Base( elementSubRegion,
          finiteElementSpace,
          inputConstitutiveType ),
    m_nodeCoords( nodeManager.getField< fields::referencePosition32 >() ),
    m_ux_n( nodeManager.getField< fields::elasticfields::Displacementx_n >() ),
    m_uy_n( nodeManager.getField< fields::elasticfields::Displacementy_n >() ),
    m_uz_n( nodeManager.getField< fields::elasticfields::Displacementz_n >() ),
    m_stiffnessVectorx( nodeManager.getField< fields::elasticfields::StiffnessVectorx >() ),
    m_stiffnessVectory( nodeManager.getField< fields::elasticfields::StiffnessVectory >() ),
    m_stiffnessVectorz( nodeManager.getField< fields::elasticfields::StiffnessVectorz >() ),
    m_density( elementSubRegion.template getField< fields::elasticfields::ElasticDensity >() ),
    m_velocityVp( elementSubRegion.template getField< fields::elasticfields::ElasticVelocityVp >() ),
    m_velocityVs( elementSubRegion.template getField< fields::elasticfields::ElasticVelocityVs >() ),
    m_gamma( elementSubRegion.template getField< fields::elasticvtifields::Gamma >()),
    m_epsilon( elementSubRegion.template getField< fields::elasticvtifields::Epsilon >()),
    m_delta( elementSubRegion.template getField< fields::elasticvtifields::Delta >()),
    m_theta( elementSubRegion.template getField< fields::elasticttifields::Theta >()),
    m_phi( elementSubRegion.template getField< fields::elasticttifields::Phi >()),
    m_dt( dt )
  {
    GEOS_UNUSED_VAR( edgeManager );
    GEOS_UNUSED_VAR( faceManager );
    GEOS_UNUSED_VAR( targetRegionIndex );
  }



  //*****************************************************************************
  /**
   * @copydoc geos::finiteElement::KernelBase::StackVariables
   *
   * ### ExplicitElasticTTISEM Description
   * Adds a stack arrays for the nodal force, primary displacement variable, etc.
   */
  struct StackVariables : Base::StackVariables
  {
public:
    GEOS_HOST_DEVICE
    StackVariables():
      xLocal(),
      stiffnessVectorxLocal(),
      stiffnessVectoryLocal(),
      stiffnessVectorzLocal()
    {}
    /// C-array stack storage for element local the nodal positions.
    real64 xLocal[ 8 ][ 3 ]{};
    real32 stiffnessVectorxLocal[ numNodesPerElem ];
    real32 stiffnessVectoryLocal[ numNodesPerElem ];
    real32 stiffnessVectorzLocal[ numNodesPerElem ];
    real32 CVTI[6][6]{};
    real32 CTTI[6][6]{};

  };
  //***************************************************************************


  /**
   * @copydoc geos::finiteElement::KernelBase::setup
   *
   * Copies the primary variable, and position into the local stack array.
   */
  GEOS_HOST_DEVICE
  inline
  void setup( localIndex const k,
              StackVariables & stack ) const
  {
    for( localIndex a=0; a< 8; a++ )
    {
      localIndex const nodeIndex = m_elemsToNodes( k, FE_TYPE::meshIndexToLinearIndex3D( a ) );
      for( int i=0; i< 3; ++i )
      {
        stack.xLocal[ a ][ i ] = m_nodeCoords[ nodeIndex ][ i ];
      }
    }

    //Initialize VTI elasticity tensor
    for( localIndex i = 0; i < 6; i++ )
    {
      for( localIndex j = 0; j < 6; j++ )
      {
        stack.CVTI[i][j] = 0.0;
      }
    }


    stack.CVTI[0][0] = m_density[k] * pow( m_velocityVp[k], 2 ) * (1.0 + 2.0*m_epsilon[k]);
    stack.CVTI[1][1]= stack.CVTI[0][0];


    stack.CVTI[2][0] = m_density[k] *
                       sqrt( pow((pow( m_velocityVp[k],
                                       2 ) - pow( m_velocityVs[k], 2 )),
                                 2 ) + 2.0 * pow( m_velocityVp[k], 2 ) * m_delta[k] * (pow( m_velocityVp[k], 2 ) - pow( m_velocityVs[k], 2 )) ) - m_density[k] * pow(
      m_velocityVs[k], 2 );

    stack.CVTI[0][2] = stack.CVTI[2][0];

    stack.CVTI[1][2] = stack.CVTI[0][2];

    stack.CVTI[2][1] = stack.CVTI[1][2];

    stack.CVTI[2][2] = m_density[k] * pow( m_velocityVp[k], 2 );

    stack.CVTI[3][3] = m_density[k] * pow( m_velocityVs[k], 2 );

    stack.CVTI[4][4] = stack.CVTI[3][3];

    stack.CVTI[5][5] = m_density[k] * pow( m_velocityVs[k], 2 )*(1.0 + 2.0 * m_gamma[k]);

    stack.CVTI[1][0] = stack.CVTI[0][0] - 2.0 * stack.CVTI[5][5];
    stack.CVTI[0][1] = stack.CVTI[1][0];

    // Voigt array useful fot the CTTI tensor
    int Voigt[3][3];

    Voigt[0][0] = 0;
    Voigt[1][1] = 1;
    Voigt[2][2] = 2;
    Voigt[2][1] = 3;
    Voigt[1][2] = 3;
    Voigt[0][2] = 4;
    Voigt[2][0] = 4;
    Voigt[0][1] = 5;
    Voigt[1][0] = 5;

    //Rotation matrix
    real64 rotationMatrix[3][3];

    rotationMatrix[0][0] = LvArray::math::cos( m_theta[k] ) * LvArray::math::cos( m_phi[k] );
    rotationMatrix[0][1] = LvArray::math::cos( m_theta[k] )*LvArray::math::sin( m_phi[k] );
    rotationMatrix[0][2] = -LvArray::math::sin( m_theta[k] );
    rotationMatrix[1][0] = -LvArray::math::sin( m_phi[k] );
    rotationMatrix[1][1] = LvArray::math::cos( m_phi[k] );
    rotationMatrix[1][2] = 0.0;
    rotationMatrix[2][0] =  LvArray::math::sin( m_theta[k] ) * LvArray::math::cos( m_phi[k] );
    rotationMatrix[2][1] =  LvArray::math::sin( m_theta[k] ) * LvArray::math::sin( m_phi[k] );
    rotationMatrix[2][2] = LvArray::math::cos( m_theta[k] );


    // Build CTTI tensor
    real64 coeff;
    for( localIndex l = 0; l < 3; ++l )
    {
      for( localIndex m = 0; m < 3; ++m )
      {
        for( localIndex j = 0; j < 3; ++j )
        {
          for( localIndex i = 0; i < 3; ++i )
          {
            coeff=0.0;
            for( localIndex a = 0; a < 3; ++a )
            {
              for( localIndex b = 0; b < 3; ++b )
              {
                for( localIndex c = 0; c < 3; ++c )
                {
                  for( localIndex d = 0; d < 3; ++d )
                  {
                    coeff+= rotationMatrix[ d ][ i ] * rotationMatrix[ c ][ j ] * rotationMatrix[ b ][ m ] * rotationMatrix[ a ][ l ] * stack.CVTI[ Voigt[ d ][ c ] ][ Voigt[ b ][ a ] ];
                  }
                }
              }
            }
            stack.CTTI[ Voigt[ i ][ j ] ][ Voigt[ m ][ l ] ] = coeff;
          }
        }
      }
    }

  }

  /**
   * @copydoc geos::finiteElement::KernelBase::complete
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  real64 complete( localIndex const k,
                   StackVariables & stack ) const
  {
    for( int i=0; i<numNodesPerElem; i++ )
    {
      const localIndex nodeIndex = m_elemsToNodes( k, i );
      RAJA::atomicAdd< parallelDeviceAtomic >( &m_stiffnessVectorx[ nodeIndex ], stack.stiffnessVectorxLocal[ i ] );
      RAJA::atomicAdd< parallelDeviceAtomic >( &m_stiffnessVectory[ nodeIndex ], stack.stiffnessVectoryLocal[ i ] );
      RAJA::atomicAdd< parallelDeviceAtomic >( &m_stiffnessVectorz[ nodeIndex ], stack.stiffnessVectorzLocal[ i ] );
    }
    return 0;
  }

  /**
   * @copydoc geos::finiteElement::KernelBase::quadraturePointKernel
   *
   * ### ExplicitElasticTTISEM Description
   * Calculates stiffness vector
   *
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack ) const
  {

    m_finiteElementSpace.template computeFirstOrderStiffnessTerm<>( q, stack.xLocal, [&] ( int const i, int const j, real64 const val, real64 const (&J)[3][3], int const p, int const r )
    {

      real32 const Rxx_ij = val*(stack.CTTI[0][0]*J[p][0]*J[r][0]+stack.CTTI[5][0]*J[p][0]*J[r][1]+stack.CTTI[4][0]*J[p][0]*J[r][2]+
                                 stack.CTTI[0][5]*J[p][1]*J[r][0]+stack.CTTI[5][5]*J[p][1]*J[r][1]+stack.CTTI[4][5]*J[p][1]*J[r][2]+
                                 stack.CTTI[0][4]*J[p][2]*J[r][0]+stack.CTTI[5][4]*J[p][2]*J[r][1]+stack.CTTI[4][4]*J[p][2]*J[r][2]);

      real32 const Ryy_ij = val*(stack.CTTI[5][5]*J[p][0]*J[r][0]+stack.CTTI[1][5]*J[p][0]*J[r][1]+stack.CTTI[3][5]*J[p][0]*J[r][2]+
                                 stack.CTTI[5][2]*J[p][1]*J[r][0]+stack.CTTI[2][2]*J[p][1]*J[r][1]+stack.CTTI[3][2]*J[p][1]*J[r][2]+
                                 stack.CTTI[5][3]*J[p][2]*J[r][0]+stack.CTTI[2][3]*J[p][2]*J[r][1]+stack.CTTI[3][3]*J[p][2]*J[r][2]);

      real32 const Rzz_ij = val*(stack.CTTI[4][4]*J[p][0]*J[r][0]+stack.CTTI[3][4]*J[p][0]*J[r][1]+stack.CTTI[2][4]*J[p][0]*J[r][2]+
                                 stack.CTTI[4][3]*J[p][1]*J[r][0]+stack.CTTI[3][3]*J[p][1]*J[r][1]+stack.CTTI[2][3]*J[p][1]*J[r][2]+
                                 stack.CTTI[4][2]*J[p][2]*J[r][0]+stack.CTTI[3][2]*J[p][2]*J[r][1]+stack.CTTI[2][2]*J[p][2]*J[r][2]);

      real32 const Ryx_ij = val*(stack.CTTI[0][5]*J[p][0]*J[r][0]+stack.CTTI[5][5]*J[p][0]*J[r][1]+stack.CTTI[4][5]*J[p][0]*J[r][2]+
                                 stack.CTTI[0][1]*J[p][1]*J[r][0]+stack.CTTI[4][1]*J[p][1]*J[r][1]+stack.CTTI[4][1]*J[p][1]*J[r][2]+
                                 stack.CTTI[0][3]*J[p][2]*J[r][0]+stack.CTTI[5][3]*J[p][2]*J[r][1]+stack.CTTI[4][3]*J[p][2]*J[r][2]);

      real32 const Rxy_ij = val*(stack.CTTI[5][0]*J[p][0]*J[r][0]+stack.CTTI[1][0]*J[p][0]*J[r][1]+stack.CTTI[3][0]*J[p][0]*J[r][2]+
                                 stack.CTTI[5][5]*J[p][1]*J[r][0]+stack.CTTI[1][5]*J[p][1]*J[r][1]+stack.CTTI[3][5]*J[p][1]*J[r][2]+
                                 stack.CTTI[5][4]*J[p][2]*J[r][0]+stack.CTTI[1][4]*J[p][2]*J[r][1]+stack.CTTI[3][4]*J[p][2]*J[r][2]);

      real32 const Rzx_ij = val*(stack.CTTI[0][4]*J[p][0]*J[r][0]+stack.CTTI[5][4]*J[p][0]*J[r][1]+stack.CTTI[4][4]*J[p][0]*J[r][2]+
                                 stack.CTTI[0][3]*J[p][1]*J[r][0]+stack.CTTI[5][3]*J[p][1]*J[r][1]+stack.CTTI[4][3]*J[p][1]*J[r][2]+
                                 stack.CTTI[0][2]*J[p][2]*J[r][0]+stack.CTTI[4][2]*J[p][2]*J[r][1]+stack.CTTI[4][2]*J[p][2]*J[r][2]);

      real32 const Rxz_ij = val*(stack.CTTI[4][0]*J[p][0]*J[r][0]+stack.CTTI[3][0]*J[p][0]*J[r][1]+stack.CTTI[2][0]*J[p][0]*J[r][2]+
                                 stack.CTTI[4][5]*J[p][1]*J[r][0]+stack.CTTI[3][5]*J[p][1]*J[r][1]+stack.CTTI[2][5]*J[p][1]*J[r][2]+
                                 stack.CTTI[4][4]*J[p][2]*J[r][0]+stack.CTTI[3][4]*J[p][2]*J[r][1]+stack.CTTI[2][4]*J[p][2]*J[r][2]);

      real32 const Rzy_ij = val*(stack.CTTI[5][4]*J[p][0]*J[r][0]+stack.CTTI[1][4]*J[p][0]*J[r][1]+stack.CTTI[3][4]*J[p][0]*J[r][2]+
                                 stack.CTTI[5][3]*J[p][1]*J[r][0]+stack.CTTI[1][3]*J[p][1]*J[r][1]+stack.CTTI[3][3]*J[p][1]*J[r][2]+
                                 stack.CTTI[5][2]*J[p][2]*J[r][0]+stack.CTTI[1][2]*J[p][2]*J[r][1]+stack.CTTI[3][2]*J[p][2]*J[r][2]);

      real32 const Ryz_ij = val*(stack.CTTI[4][5]*J[p][0]*J[r][0]+stack.CTTI[3][5]*J[p][0]*J[r][1]+stack.CTTI[2][5]*J[p][0]*J[r][2]+
                                 stack.CTTI[4][1]*J[p][1]*J[r][0]+stack.CTTI[3][1]*J[p][1]*J[r][1]+stack.CTTI[2][1]*J[p][1]*J[r][2]+
                                 stack.CTTI[5][3]*J[p][2]*J[r][0]+stack.CTTI[3][3]*J[p][2]*J[r][1]+stack.CTTI[2][3]*J[p][2]*J[r][2]);

      real32 const localIncrementx = (Rxx_ij * m_ux_n[m_elemsToNodes( k, j )] + Rxy_ij*m_uy_n[m_elemsToNodes( k, j )] + Rxz_ij*m_uz_n[m_elemsToNodes( k, j )]);
      real32 const localIncrementy = (Ryx_ij * m_ux_n[m_elemsToNodes( k, j )] + Ryy_ij*m_uy_n[m_elemsToNodes( k, j )] + Ryz_ij*m_uz_n[m_elemsToNodes( k, j )]);
      real32 const localIncrementz = (Rzx_ij * m_ux_n[m_elemsToNodes( k, j )] + Rzy_ij*m_uy_n[m_elemsToNodes( k, j )] + Rzz_ij*m_uz_n[m_elemsToNodes( k, j )]);

      stack.stiffnessVectorxLocal[ i ] += localIncrementx;
      stack.stiffnessVectoryLocal[ i ] += localIncrementy;
      stack.stiffnessVectorzLocal[ i ] += localIncrementz;
    } );
  }


protected:
  /// The array containing the nodal position array.
  arrayView2d< WaveSolverBase::wsCoordType const, nodes::REFERENCE_POSITION_USD > const m_nodeCoords;

  /// The array containing the nodal displacement array in x direction.
  arrayView1d< real32 > const m_ux_n;

  /// The array containing the nodal displacement array in y direction.
  arrayView1d< real32 > const m_uy_n;

  /// The array containing the nodal displacement array in z direction.
  arrayView1d< real32 > const m_uz_n;

  /// The array containing the product of the stiffness matrix and the nodal displacement.
  arrayView1d< real32 > const m_stiffnessVectorx;

  /// The array containing the product of the stiffness matrix and the nodal displacement.
  arrayView1d< real32 > const m_stiffnessVectory;

  /// The array containing the product of the stiffness matrix and the nodal displacement.
  arrayView1d< real32 > const m_stiffnessVectorz;

  /// The array containing the density of the medium
  arrayView1d< real32 const > const m_density;

  /// The array containing the P-wavespeed
  arrayView1d< real32 const > const m_velocityVp;

  /// The array containing the S-wavespeed
  arrayView1d< real32 const > const m_velocityVs;

  ///The array containing the Thomsen constant gamma
  arrayView1d< real32 const > const m_gamma;

  ///The array containing the Thomsen constant epsilon
  arrayView1d< real32 const > const m_epsilon;

  ///The array containing the Thomsen constant delta
  arrayView1d< real32 const > const m_delta;

  ///The array containing tilt angle
  arrayView1d< real32 const > const m_theta;

  ///The array containing the azimut angle
  arrayView1d< real32 const > const m_phi;

  /// The time increment for this time integration step.
  real64 const m_dt;


};


/// The factory used to construct a ExplicitAcousticWaveEquation kernel.
using ExplicitElasticTTISEMFactory = finiteElement::KernelFactory< ExplicitElasticTTISEM,
                                                                   real64 >;

} // namespace elasticTTIWaveEquationSEMKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ElasticTTIWaveEquationSEMKERNEL_HPP_
