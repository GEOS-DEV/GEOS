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

#ifndef GEOSX_MESH_MESHFORLOOPINTERFACE_HPP
#define GEOSX_MESH_MESHFORLOOPINTERFACE_HPP

#include "finiteElement/FiniteElementDiscretization.hpp"
#include "finiteElement/FiniteElementDiscretizationManager.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"

#include "common/DataTypes.hpp"
#include "mesh/MeshLevel.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "finiteElement/ElementLibrary/FiniteElementBase.h"

namespace geosx
{

template< class POLICY=serialPolicy, typename LAMBDA=void >
void forAllElemsInMesh( MeshLevel const * const mesh, LAMBDA && lambdaBody )
{

  ElementRegionManager const * const elemManager = mesh->getElemManager();

  for( localIndex er=0; er<elemManager->numRegions(); ++er )
  {
    ElementRegionBase const * const elemRegion = elemManager->GetRegion( er );
    for( localIndex esr=0; esr<elemRegion->numSubRegions(); ++esr )
    {
      ElementSubRegionBase const * const elementSubRegion = elemRegion->GetSubRegion( esr );

      forall_in_range< POLICY >( 0, elementSubRegion->size(),
                                 [=]( localIndex index ) mutable -> void
      {
        lambdaBody( er, esr, index );
      } );
    }
  }
}

template< typename NUMBER=real64, class POLICY=serialPolicy, typename LAMBDA=void >
std::pair< NUMBER, localIndex >
minloc_in_range( localIndex const begin, const localIndex end, LAMBDA && body )
{
  NUMBER minDist = std::numeric_limits< NUMBER >::max();
  localIndex minDistId = -1;

  // TODO: make this a RAJA loop
  for( localIndex ei = begin; ei < end; ++ei )
  {
    NUMBER dist = body( ei );

    if( dist < minDist )
    {
      minDist   = dist;
      minDistId = ei;
    }
  }

  return std::make_pair( minDist, minDistId );
}


template< typename NUMBER=real64, class POLICY=serialPolicy, typename LAMBDA=void >
std::pair< NUMBER, std::tuple< localIndex, localIndex, localIndex > >
minLocOverElemsInMesh( MeshLevel const * const mesh, LAMBDA && lambdaBody )
{
  NUMBER minVal = std::numeric_limits< NUMBER >::max();
  localIndex minReg = -1, minSubreg = -1, minIndex = -1;

  ElementRegionManager const * const elemManager = mesh->getElemManager();

  for( localIndex er=0; er<elemManager->numRegions(); ++er )
  {
    ElementRegionBase const * const elemRegion = elemManager->GetRegion( er );

    elemRegion->forElementSubRegionsIndex< CellElementSubRegion >( [&]( localIndex const esr, CellElementSubRegion const & subRegion )
    {
      auto ebody = [=]( localIndex index )
      {
        return lambdaBody( er, esr, index );
      };

      auto ret = minloc_in_range< NUMBER, POLICY >( 0, subRegion.size(), ebody );

      if( ret.first < minVal )
      {
        minVal    = ret.first;
        minReg    = er;
        minSubreg = esr;
        minIndex  = ret.second;
      }
    } );
  }

  return std::make_pair( minVal, std::make_tuple( minReg, minSubreg, minIndex ));
}

} // namespace geosx

#endif // GEOSX_MESH_MESHFORLOOPINTERFACE_HPP
