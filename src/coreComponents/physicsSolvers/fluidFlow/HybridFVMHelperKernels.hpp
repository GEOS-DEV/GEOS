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
 * @file CompositionalMultiphaseHybridFVMHelperKernels.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_HYBRIDFVMUPWINDINGHELPERKERNELS_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_HYBRIDFVMUPWINDINGHELPERKERNELS_HPP

#include "common/DataTypes.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "mesh/MeshLevel.hpp"

namespace geos
{

namespace hybridFVMKernels
{

/******************************** Kernel switches ********************************/

template< typename T, typename LAMBDA >
void kernelLaunchSelectorFaceSwitch( T value, LAMBDA && lambda )
{
  static_assert( std::is_integral< T >::value, "kernelLaunchSelectorFaceSwitch: type should be integral" );

  switch( value )
  {
    case 4:
    { lambda( std::integral_constant< T, 4 >() ); return;}
    case 5:
    { lambda( std::integral_constant< T, 5 >() ); return;}
    case 6:
    { lambda( std::integral_constant< T, 6 >() ); return;}
    case 7:
    { lambda( std::integral_constant< T, 7 >() ); return;}
    case 8:
    { lambda( std::integral_constant< T, 8 >() ); return;}
    case 9:
    { lambda( std::integral_constant< T, 9 >() ); return;}
    case 10:
    { lambda( std::integral_constant< T, 10 >() ); return;}
    case 11:
    { lambda( std::integral_constant< T, 11 >() ); return;}
    /*
       This is disabled for now to pass the CI tests for the CUDA builds
       case 12:
       { lambda( std::integral_constant< T, 12 >() ); return;}
       case 13:
       { lambda( std::integral_constant< T, 13 >() ); return;}
     */
    default: GEOS_ERROR( "Unknown numFacesInElem value: " << value );
  }
}

template< typename KERNELWRAPPER, typename INNER_PRODUCT, typename ... ARGS >
void simpleKernelLaunchSelector( localIndex numFacesInElem, ARGS && ... args )
{
  kernelLaunchSelectorFaceSwitch( numFacesInElem, [&] ( auto NUM_FACES )
  {
    KERNELWRAPPER::template launch< INNER_PRODUCT, NUM_FACES() >( std::forward< ARGS >( args )... );
  } );
}

/******************************** CellConnectivity ********************************/

struct CellConnectivity
{

  GEOS_HOST_DEVICE
  static bool
  isNeighborFound( localIndex const (&localIds)[3],
                   localIndex const ifaceLoc,
                   arrayView2d< localIndex const > const & elemRegionList,
                   arrayView2d< localIndex const > const & elemSubRegionList,
                   arrayView2d< localIndex const > const & elemList,
                   SortedArrayView< localIndex const > const & regionFilter,
                   arraySlice1d< localIndex const > const & elemToFaces,
                   localIndex ( & neighborIds )[3] )
  {
    neighborIds[0] = localIds[0];
    neighborIds[1] = localIds[1];
    neighborIds[2] = localIds[2];

    // the face has at most two adjacent elements
    // one of these two elements is the current element indexed by er, esr, ei
    // but here we are interested in the indices of the other element
    // this other element is "the neighbor" for this one-sided face
    for( localIndex k=0; k<elemRegionList.size( 1 ); ++k )
    {

      localIndex const erNeighbor  = elemRegionList[elemToFaces[ifaceLoc]][k];
      localIndex const esrNeighbor = elemSubRegionList[elemToFaces[ifaceLoc]][k];
      localIndex const eiNeighbor  = elemList[elemToFaces[ifaceLoc]][k];

      // this element is not the current element
      // we have found the neighbor or we are at the boundary
      if( erNeighbor != localIds[0] || esrNeighbor != localIds[1] || eiNeighbor != localIds[2] )
      {
        bool const onBoundary       = (erNeighbor == -1 || esrNeighbor == -1 || eiNeighbor == -1);
        bool const neighborInTarget = regionFilter.contains( erNeighbor );

        // if not on boundary, save the element indices
        if( !onBoundary && neighborInTarget )
        {
          neighborIds[0] = erNeighbor;
          neighborIds[1] = esrNeighbor;
          neighborIds[2] = eiNeighbor;
        }
        // if the face is on the boundary, use the properties of the local elem
      }
    }
    return !( localIds[0] == neighborIds[0] &&
              localIds[1] == neighborIds[1] &&
              localIds[2] == neighborIds[2] );
  }

};


} // namespace hybridFVMUpwindingKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEHYBRIDFVMUPWINDINGHELPERKERNELS_HPP
