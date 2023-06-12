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
 * @file WaveSolverUtils.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_WAVESOLVERUTILS_HPP_
#define GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_WAVESOLVERUTILS_HPP_

#include "WaveSolverBase.hpp"

namespace geos
{

struct WaveSolverUtils
{

  GEOS_HOST_DEVICE
  static real32 evaluateRicker( real64 const & time_n, real32 const & f0, localIndex order )
  {
    real32 const o_tpeak = 1.0/f0;
    real32 pulse = 0.0;
    if((time_n <= -0.9*o_tpeak) || (time_n >= 2.9*o_tpeak))
    {
      return pulse;
    }

    constexpr real32 pi = M_PI;
    real32 const lam = (f0*pi)*(f0*pi);

    switch( order )
    {
      case 2:
      {
        pulse = 2.0*lam*(2.0*lam*(time_n-o_tpeak)*(time_n-o_tpeak)-1.0)*exp( -lam*(time_n-o_tpeak)*(time_n-o_tpeak));
      }
      break;
      case 1:
      {
        pulse = -2.0*lam*(time_n-o_tpeak)*exp( -lam*(time_n-o_tpeak)*(time_n-o_tpeak));
      }
      break;
      case 0:
      {
        pulse = -(time_n-o_tpeak)*exp( -2*lam*(time_n-o_tpeak)*(time_n-o_tpeak) );
      }
      break;
      default:
        GEOS_ERROR( "This option is not supported yet, rickerOrder must be 0, 1 or 2" );
    }

    return pulse;
  }

  static void computeSeismoTrace( real64 const time_n,
                                  real64 const dt,
                                  real64 const timeSeismo,
                                  localIndex iSeismo,
                                  arrayView2d< localIndex const > const receiverNodeIds,
                                  arrayView2d< real64 const > const receiverConstants,
                                  arrayView1d< localIndex const > const receiverIsLocal,
                                  localIndex const nsamplesSeismoTrace,
                                  localIndex const outputSeismoTrace,
                                  arrayView1d< real32 const > const var_np1,
                                  arrayView1d< real32 const > const var_n,
                                  arrayView2d< real32 > varAtReceivers )
  {
    real64 const time_np1 = time_n + dt;

    real32 const a1 = (LvArray::math::abs( dt ) < WaveSolverBase::epsilonLoc ) ? 1.0 : (time_np1 - timeSeismo)/dt;
    real32 const a2 = 1.0 - a1;

    if( nsamplesSeismoTrace > 0 )
    {
      forAll< WaveSolverBase::EXEC_POLICY >( receiverConstants.size( 0 ), [=] GEOS_HOST_DEVICE ( localIndex const ircv )
      {
        if( receiverIsLocal[ircv] == 1 )
        {
          varAtReceivers[iSeismo][ircv] = 0.0;
          real32 vtmp_np1 = 0.0;
          real32 vtmp_n = 0.0;
          for( localIndex inode = 0; inode < receiverConstants.size( 1 ); ++inode )
          {
            vtmp_np1 += var_np1[receiverNodeIds[ircv][inode]] * receiverConstants[ircv][inode];
            vtmp_n += var_n[receiverNodeIds[ircv][inode]] * receiverConstants[ircv][inode];
          }
          // linear interpolation between the pressure value at time_n and time_(n+1)
          varAtReceivers[iSeismo][ircv] = a1*vtmp_n + a2*vtmp_np1;
        }
      } );
    }

    // TODO DEBUG: the following output is only temporary until our wave propagation kernels are finalized.
    // Output will then only be done via the previous code.
    if( iSeismo == nsamplesSeismoTrace - 1 )
    {
      forAll< serialPolicy >( receiverConstants.size( 0 ), [=] ( localIndex const ircv )
      {
        if( outputSeismoTrace == 1 )
        {
          if( receiverIsLocal[ircv] == 1 )
          {
            // Note: this "manual" output to file is temporary
            //       It should be removed as soon as we can use TimeHistory to output data not registered on the mesh
            // TODO: remove saveSeismo and replace with TimeHistory
            std::ofstream f( GEOS_FMT( "seismoTraceReceiver{:03}.txt", ircv ), std::ios::app );
            for( localIndex iSample = 0; iSample < nsamplesSeismoTrace; ++iSample )
            {
              f << iSample << " " << varAtReceivers[iSample][ircv] << std::endl;
            }
            f.close();
          }
        }
      } );
    }
  }

  static void computeSeismoTracePOD( real64 const time_n,
                                     real64 const dt,
                                     real64 const timeSeismo,
                                     localIndex iSeismo,
				     arrayView2d< real64 const > const receiverConstants,
                                     arrayView1d< localIndex const > const receiverIsLocal,
                                     localIndex const nsamplesSeismoTrace,
                                     localIndex const outputSeismoTrace,
                                     arrayView1d< real32 const > const var_np1,
                                     arrayView1d< real32 const > const var_n,
                                     arrayView2d< real32 > varAtReceivers )
  {
    real64 const time_np1 = time_n + dt;

    real32 const a1 = (LvArray::math::abs( dt ) < WaveSolverBase::epsilonLoc ) ? 1.0 : (time_np1 - timeSeismo)/dt;
    real32 const a2 = 1.0 - a1;
    if( nsamplesSeismoTrace > 0 )
    {
      forAll< WaveSolverBase::EXEC_POLICY >( receiverConstants.size( 0 ), [=] GEOS_HOST_DEVICE ( localIndex const ircv )
	{
	  if( receiverIsLocal[ircv] == 1 )
	  {
	    varAtReceivers[iSeismo][ircv] = 0.0;
	    real32 vtmp_np1 = 0.0;
	    real32 vtmp_n = 0.0;
	    
	    for( localIndex m = 0; m<var_np1.size( 0 ); ++m )
	    {
	      vtmp_np1 += var_np1[m] * receiverConstants[ircv][m];
	      vtmp_n += var_n[m] * receiverConstants[ircv][m];
	    }
            // linear interpolation between the pressure value at time_n and time_(n+1)
	    varAtReceivers[iSeismo][ircv] = a1*vtmp_n + a2*vtmp_np1;
	  }
	} );
    }

    // TODO DEBUG: the following output is only temporary until our wave propagation kernels are finalized.
    // Output will then only be done via the previous code.
    if( iSeismo == nsamplesSeismoTrace - 1 )
    {
      forAll< serialPolicy >( receiverConstants.size( 0 ), [=] ( localIndex const ircv )
        {
	  if( outputSeismoTrace == 1 )
	  {
	    if( receiverIsLocal[ircv] == 1 )
	    {
	      // Note: this "manual" output to file is temporary
	      //       It should be removed as soon as we can use TimeHistory to output data not registered on the mesh
	      // TODO: remove saveSeismo and replace with TimeHistory
	      std::ofstream f( GEOS_FMT( "seismoTraceReceiver{:03}.txt", ircv ), std::ios::app );
	      for( localIndex iSample = 0; iSample < nsamplesSeismoTrace; ++iSample )
	      {
		f << iSample << " " << varAtReceivers[iSample][ircv] << std::endl;
	      }
	      f.close();
	    }
	  }
	} );
    }
  }

  static void compute2dVariableSeismoTrace( real64 const time_n,
                                            real64 const dt,
                                            real64 const timeSeismo,
                                            localIndex iSeismo,
                                            arrayView1d< localIndex const > const rcvElem,
                                            arrayView2d< real64 const > const receiverConstants,
                                            arrayView1d< localIndex const > const receiverIsLocal,
                                            localIndex const nsamplesSeismoTrace,
                                            localIndex const outputSeismoTrace,
                                            arrayView2d< real32 const > const var_np1,
                                            arrayView2d< real32 const > const var_n,
                                            arrayView2d< real32 > varAtReceivers )
  {
    real64 const time_np1 = time_n+dt;

    real32 const a1 = (dt < WaveSolverBase::epsilonLoc) ? 1.0 : (time_np1 - timeSeismo)/dt;
    real32 const a2 = 1.0 - a1;

    if( nsamplesSeismoTrace > 0 )
    {
      forAll< WaveSolverBase::EXEC_POLICY >( receiverConstants.size( 0 ), [=] GEOS_HOST_DEVICE ( localIndex const ircv )
      {
        if( receiverIsLocal[ircv] == 1 )
        {
          varAtReceivers[iSeismo][ircv] = 0.0;
          real32 vtmp_np1 = 0.0;
          real32 vtmp_n = 0.0;
          for( localIndex inode = 0; inode < receiverConstants.size( 1 ); ++inode )
          {
            vtmp_np1 += var_np1[rcvElem[ircv]][inode] * receiverConstants[ircv][inode];
            vtmp_n += var_n[rcvElem[ircv]][inode] * receiverConstants[ircv][inode];
          }
          // linear interpolation between the pressure value at time_n and time_(n+1)
          varAtReceivers[iSeismo][ircv] = a1*vtmp_n + a2*vtmp_np1;
        }
      } );
    }

    // TODO DEBUG: the following output is only temporary until our wave propagation kernels are finalized.
    // Output will then only be done via the previous code.
    if( iSeismo == nsamplesSeismoTrace - 1 )
    {
      forAll< serialPolicy >( receiverConstants.size( 0 ), [=] ( localIndex const ircv )
      {
        if( outputSeismoTrace == 1 )
        {
          if( receiverIsLocal[ircv] == 1 )
          {
            // Note: this "manual" output to file is temporary
            //       It should be removed as soon as we can use TimeHistory to output data not registered on the mesh
            // TODO: remove saveSeismo and replace with TimeHistory
            std::ofstream f( GEOS_FMT( "seismoTraceReceiver{:03}.txt", ircv ), std::ios::app );
            for( localIndex iSample = 0; iSample < nsamplesSeismoTrace; ++iSample )
            {
              f << iSample << " " << varAtReceivers[iSample][ircv] << std::endl;
            }
            f.close();
          }
        }
      } );
    }

  }


  /**
   * @brief Check if the source point is inside an element or not
   */

  /**
   * @brief Check if the source point is inside an element or not
   * @param numFacesPerElem number of face on an element
   * @param elemCenter array containing the center of the elements
   * @param faceNormal array containing the normal of all faces
   * @param faceCenter array containing the center of all faces
   * @param elemsToFaces map to get the global faces from element index and local face index
   * @param coords coordinate of the point
   * @return true if coords is inside the element
   */

  GEOS_HOST_DEVICE
  static bool
  locateSourceElement( real64 const numFacesPerElem,
                       real64 const (&elemCenter)[3],
                       arrayView2d< real64 const > const faceNormal,
                       arrayView2d< real64 const > const faceCenter,
                       arraySlice1d< localIndex const > const elemsToFaces,
                       real64 const (&coords)[3] )
  {
    //Loop over the element faces
    real64 tmpVector[3]{};
    for( localIndex kfe = 0; kfe < numFacesPerElem; ++kfe )
    {

      localIndex const iface = elemsToFaces[kfe];
      real64 faceCenterOnFace[3] = {faceCenter[iface][0],
                                    faceCenter[iface][1],
                                    faceCenter[iface][2]};
      real64 faceNormalOnFace[3] = {faceNormal[iface][0],
                                    faceNormal[iface][1],
                                    faceNormal[iface][2]};

      //Test to make sure if the normal is outwardly directed
      LvArray::tensorOps::copy< 3 >( tmpVector, faceCenterOnFace );
      LvArray::tensorOps::subtract< 3 >( tmpVector, elemCenter );
      if( LvArray::tensorOps::AiBi< 3 >( tmpVector, faceNormalOnFace ) < 0.0 )
      {
        LvArray::tensorOps::scale< 3 >( faceNormalOnFace, -1 );
      }

      // compute the vector face center to query point
      LvArray::tensorOps::subtract< 3 >( faceCenterOnFace, coords );
      localIndex const s = computationalGeometry::sign( LvArray::tensorOps::AiBi< 3 >( faceNormalOnFace, faceCenterOnFace ));

      // all dot products should be non-negative (we enforce outward normals)
      if( s < 0 )
      {
        return false;
      }

    }
    return true;
  }

/**
 * @brief Convert a mesh element point coordinate into a coordinate on the reference element
 * @tparam FE_TYPE finite element type
 * @param[in] coords coordinate of the point
 * @param[in] elemsToNodes map to obtaint global nodes from element index
 * @param[in] X array of mesh nodes coordinates
 * @param[out] coordsOnRefElem to contain the coordinate computed in the reference element
 */
  template< typename FE_TYPE >
  GEOS_HOST_DEVICE
  static void
  computeCoordinatesOnReferenceElement( real64 const (&coords)[3],
                                        arraySlice1d< localIndex const, cells::NODE_MAP_USD - 1 > const elemsToNodes,
                                        arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X,
                                        real64 (& coordsOnRefElem)[3] )
  {
    real64 xLocal[FE_TYPE::numNodes][3]{};
    for( localIndex a = 0; a < FE_TYPE::numNodes; ++a )
    {
      LvArray::tensorOps::copy< 3 >( xLocal[a], X[ elemsToNodes[a] ] );
    }
    // coordsOnRefElem = invJ*(coords-coordsNode_0)
    real64 invJ[3][3]{};
    FE_TYPE::invJacobianTransformation( 0, xLocal, invJ );
    for( localIndex i = 0; i < 3; ++i )
    {
      // init at (-1,-1,-1) as the origin of the referential elem
      coordsOnRefElem[i] = -1.0;
      for( localIndex j = 0; j < 3; ++j )
      {
        coordsOnRefElem[i] += invJ[i][j] * (coords[j] - xLocal[0][j]);
      }
    }
  }



};

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_WAVESOLVERUTILS_HPP_ */
