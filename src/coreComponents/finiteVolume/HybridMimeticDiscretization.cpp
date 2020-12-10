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
 * @file HybridMimeticDiscretization.cpp
 *
 */

#include "HybridMimeticDiscretization.hpp"

#include "finiteVolume/MimeticInnerProductDispatch.hpp"
#include "finiteVolume/mimeticInnerProducts/TPFAInnerProduct.hpp"
#include "finiteVolume/mimeticInnerProducts/QuasiTPFAInnerProduct.hpp"
#include "finiteVolume/mimeticInnerProducts/QuasiRTInnerProduct.hpp"
#include "finiteVolume/mimeticInnerProducts/SimpleInnerProduct.hpp"
#include "finiteVolume/mimeticInnerProducts/BdVLMInnerProduct.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace mimeticInnerProduct;

HybridMimeticDiscretization::HybridMimeticDiscretization( std::string const & name,
                                                          Group * const parent )
  : FluxApproximationBase( name, parent )
{
  registerWrapper( viewKeyStruct::innerProductTypeString, &m_innerProductType )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Type of inner product used in the hybrid FVM solver" );
}

void HybridMimeticDiscretization::InitializePostInitialConditions_PreSubGroups( Group * const rootGroup )
{
  FluxApproximationBase::InitializePostInitialConditions_PreSubGroups( rootGroup );

  std::unique_ptr< MimeticInnerProductBase > newMimeticIP = factory( m_innerProductType );

  registerWrapper< MimeticInnerProductBase >( viewKeyStruct::innerProductString, std::move( newMimeticIP ) )->
    setRestartFlags( dataRepository::RestartFlags::NO_WRITE );
}

std::unique_ptr< MimeticInnerProductBase >
HybridMimeticDiscretization::factory( string const & mimeticInnerProductType ) const
{
  std::unique_ptr< MimeticInnerProductBase > rval;
  if( mimeticInnerProductType == MimeticInnerProductTypeStrings::TPFA )
  {
    rval = std::make_unique< TPFAInnerProduct >();
  }
  else if( mimeticInnerProductType == MimeticInnerProductTypeStrings::QuasiTPFA )
  {
    rval = std::make_unique< QuasiTPFAInnerProduct >();
  }
  else if( mimeticInnerProductType == MimeticInnerProductTypeStrings::QuasiRT )
  {
    rval = std::make_unique< QuasiRTInnerProduct >();
  }
  else if( mimeticInnerProductType == MimeticInnerProductTypeStrings::Simple )
  {
    rval = std::make_unique< SimpleInnerProduct >();
  }
  else if( mimeticInnerProductType == MimeticInnerProductTypeStrings::BdVLM )
  {
    rval = std::make_unique< BdVLMInnerProduct >();
  }
  else
  {
    GEOSX_ERROR( "Key value of "<< mimeticInnerProductType <<" does not have an associated mimetic inner product." );
  }
  return rval;
}

void HybridMimeticDiscretization::addToFractureStencil( MeshLevel & mesh,
                                                        string const & faceElementRegionName,
                                                        bool const initFlag ) const
{
  GEOSX_UNUSED_VAR( mesh );
  GEOSX_UNUSED_VAR( faceElementRegionName );
  GEOSX_UNUSED_VAR( initFlag );
  GEOSX_ERROR( "addToFractureStencil: This function is not implemented for the hybrid mimetic discretization" );
}

void HybridMimeticDiscretization::addEDFracToFractureStencil( MeshLevel & mesh,
                                                              string const & embeddedSurfaceRegionName ) const
{
  GEOSX_UNUSED_VAR( mesh );
  GEOSX_UNUSED_VAR( embeddedSurfaceRegionName );
  GEOSX_ERROR( "addEDFracToFractureStencil: This function is not implemented for the hybrid mimetic discretization" );
}


void HybridMimeticDiscretization::registerCellStencil( Group & stencilGroup ) const
{
  GEOSX_UNUSED_VAR( stencilGroup );
}

void HybridMimeticDiscretization::computeCellStencil( MeshLevel & mesh ) const
{
  GEOSX_UNUSED_VAR( mesh );
}

void HybridMimeticDiscretization::registerFractureStencil( Group & stencilGroup ) const
{
  GEOSX_UNUSED_VAR( stencilGroup );
}


void HybridMimeticDiscretization::registerBoundaryStencil( Group & stencilGroup,
                                                           string const & setName ) const
{
  GEOSX_UNUSED_VAR( stencilGroup );
  GEOSX_UNUSED_VAR( setName );
}

void HybridMimeticDiscretization::computeBoundaryStencil( MeshLevel & mesh,
                                                          string const & setName,
                                                          SortedArrayView< localIndex const > const & faceSet ) const
{
  GEOSX_UNUSED_VAR( mesh );
  GEOSX_UNUSED_VAR( setName );
  GEOSX_UNUSED_VAR( faceSet );
}

REGISTER_CATALOG_ENTRY( FluxApproximationBase, HybridMimeticDiscretization, std::string const &, Group * const )
}
