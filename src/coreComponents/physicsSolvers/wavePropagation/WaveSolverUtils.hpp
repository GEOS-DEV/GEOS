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

#ifndef GEOSX_PHYSICSSOLVERS_WAVEPROPAGATION_WAVESOLVERUTILS_HPP_
#define GEOSX_PHYSICSSOLVERS_WAVEPROPAGATION_WAVESOLVERUTILS_HPP_

namespace geosx
{

struct WaveSolverUtils
{

  GEOSX_HOST_DEVICE
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
        GEOSX_ERROR( "This option is not supported yet, rickerOrder must be 0, 1 or 2" );
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
    
    real32 const a1 = (std::abs( dt ) < 1e-8) ? 1.0 : (time_np1 - timeSeismo)/dt;
    real32 const a2 = 1.0 - a1;
  
    if( nsamplesSeismoTrace > 0 )
    {
      forAll< parallelDevicePolicy< 32 > >( receiverConstants.size( 0 ), [=] GEOSX_HOST_DEVICE ( localIndex const ircv )
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
            std::ofstream f( GEOSX_FMT( "seismoTraceReceiver{:03}.txt", ircv ), std::ios::app );
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


};

} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_WAVEPROPAGATION_WAVESOLVERUTILS_HPP_ */
