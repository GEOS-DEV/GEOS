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
 * @file AcousticMatricesSEMKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICVTISEMKERNEL_HPP_
#define GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICVTISEMKERNEL_HPP_

#include "mainInterface/ProblemManager.hpp"
#include "mesh/ElementType.hpp"
#include "AcousticVTIFields.hpp"

namespace geos
{
  struct AcousticVTISEM
  {

    static void precomputeSurfaceFieldIndicator( DomainPartition & domain )
    {
      real64 const time = 0.0;
      FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();
      FunctionManager const & functionManager = FunctionManager::getInstance();

      FaceManager & faceManager = domain.getMeshBody( 0 ).getMeshLevel( m_discretizationName ).getFaceManager();
      NodeManager & nodeManager = domain.getMeshBody( 0 ).getMeshLevel( m_discretizationName ).getNodeManager();

      ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();

      /// array of indicators: 1 if a face is on on lateral surface; 0 otherwise
      arrayView1d< localIndex > const lateralSurfaceFaceIndicator = faceManager.getField< fields::acousticvtifields::AcousticLateralSurfaceFaceIndicator >();
      /// array of indicators: 1 if a node is on on lateral surface; 0 otherwise
      arrayView1d< localIndex > const lateralSurfaceNodeIndicator = nodeManager.getField< fields::acousticvtifields::AcousticLateralSurfaceNodeIndicator >();

      /// array of indicators: 1 if a face is on on bottom surface; 0 otherwise
      arrayView1d< localIndex > const bottomSurfaceFaceIndicator = faceManager.getField< fields::acousticvtifields::AcousticBottomSurfaceFaceIndicator >();
      /// array of indicators: 1 if a node is on on bottom surface; 0 otherwise
      arrayView1d< localIndex > const bottomSurfaceNodeIndicator = nodeManager.getField< fields::acousticvtifields::AcousticBottomSurfaceNodeIndicator >();

      // Lateral surfaces
      fsManager.apply< FaceManager >( time,
                                      domain.getMeshBody( 0 ).getMeshLevel( m_discretizationName ),
                                      viewKeyStruct::lateralSurfaceString(),
                                      [&]( FieldSpecificationBase const & bc,
                                          string const &,
                                          SortedArrayView< localIndex const > const & targetSet,
                                          FaceManager &,
                                          string const & )
      {
        string const & functionName = bc.getFunctionName();

        if( functionName.empty() || functionManager.getGroup< FunctionBase >( functionName ).isFunctionOfTime() == 2 )
        {
          for( localIndex i = 0; i < targetSet.size(); ++i )
          {
            localIndex const kf = targetSet[ i ];
            lateralSurfaceFaceIndicator[kf] = 1;

            localIndex const numNodes = faceToNodeMap.sizeOfArray( kf );
            for( localIndex a=0; a < numNodes; ++a )
            {
              localIndex const dof = faceToNodeMap( kf, a );
              lateralSurfaceNodeIndicator[dof] = 1;
            }
          }
        }
        else
        {
          GEOS_ERROR( "This option is not supported yet" );
        }
      } );
    }
  };
}

#endif //GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICVTISEMKERNEL_HPP_