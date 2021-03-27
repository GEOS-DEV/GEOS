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
#include "mesh/DomainPartition.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace mimeticInnerProduct;

HybridMimeticDiscretization::HybridMimeticDiscretization( string const & name,
                                                          Group * const parent )
  : Group( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );

  // will need to add a fieldName when hybrid FVM can properly enforce (non-zero) face boundary conditions

  // needed for the multiplier
  registerWrapper( viewKeyStruct::coeffNameString(), &m_coeffName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of coefficient field" );

  registerWrapper( viewKeyStruct::innerProductTypeString(), &m_innerProductType ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Type of inner product used in the hybrid FVM solver" );
}

void HybridMimeticDiscretization::initializePostInitialConditionsPreSubGroups()
{
  Group::initializePostInitialConditionsPreSubGroups();

  std::unique_ptr< MimeticInnerProductBase > newMimeticIP = factory( m_innerProductType );

  registerWrapper< MimeticInnerProductBase >( viewKeyStruct::innerProductString(), std::move( newMimeticIP ) ).
    setRestartFlags( dataRepository::RestartFlags::NO_WRITE );
}

void HybridMimeticDiscretization::registerDataOnMesh( Group & meshBodies )
{
  meshBodies.forSubGroups< MeshBody >( [&]( MeshBody & meshBody )
  {
    meshBody.forSubGroups< MeshLevel >( [&]( MeshLevel & mesh )
    {
      FaceManager & faceManager = mesh.getFaceManager();
      faceManager.registerWrapper< array1d< real64 > >( m_coeffName + viewKeyStruct::transMultiplierString() ).
        setApplyDefaultValue( 1.0 ).
        setPlotLevel( PlotLevel::LEVEL_0 ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds the transmissibility multipliers" );

    } );
  } );
}

HybridMimeticDiscretization::CatalogInterface::CatalogType &
HybridMimeticDiscretization::getCatalog()
{
  static HybridMimeticDiscretization::CatalogInterface::CatalogType catalog;
  return catalog;
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

REGISTER_CATALOG_ENTRY( HybridMimeticDiscretization, HybridMimeticDiscretization, string const &, Group * const )

}
